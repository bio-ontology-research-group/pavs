import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';

const API = import.meta.env.VITE_API_URL ?? '';

const PREFIXES = `PREFIX pavs:  <http://pavs.kaust.edu.sa/ontology/>
PREFIX pav:   <http://pavs.kaust.edu.sa/data/>
PREFIX hgnc:  <http://identifiers.org/hgnc.symbol/>
PREFIX hp:    <http://purl.obolibrary.org/obo/HP_>
PREFIX dc:    <http://purl.org/dc/terms/>
PREFIX rdfs:  <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd:   <http://www.w3.org/2001/XMLSchema#>

`;

// Each entry: i18nKey maps to t('sparql.queries.<key>')
const EXAMPLE_QUERIES: { key: string; query: string }[] = [
  {
    key: 'commonPhenotypes',
    query: `${PREFIXES}SELECT ?label (COUNT(DISTINCT ?case) AS ?n) WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          pavs:hasPhenotype ?hpo .
  }
  GRAPH <http://pavs.kaust.edu.sa/graph/hpo-ic> {
    ?hpo rdfs:label ?label .
  }
}
GROUP BY ?label
ORDER BY DESC(?n)
LIMIT 20`,
  },
  {
    key: 'sumf1Cases',
    query: `${PREFIXES}SELECT ?id ?hgvsC ?acmg ?disease WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:hasVariant ?v .
    ?v pavs:affectsGene hgnc:SUMF1 .
    OPTIONAL { ?v pavs:hgvsC ?hgvsC }
    OPTIONAL { ?v pavs:acmgClass ?acmg }
    OPTIONAL { ?case pavs:diseaseLabel ?disease }
  }
}
ORDER BY ?id`,
  },
  {
    key: 'pathogenicVariants',
    query: `${PREFIXES}SELECT ?id ?gene ?hgvsC ?acmg ?disease WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:diseaseLabel ?disease ;
          pavs:hasVariant ?v .
    ?v pavs:acmgClass ?acmg ;
       pavs:affectsGene ?gUri .
    FILTER(CONTAINS(LCASE(?acmg), "pathogenic"))
    OPTIONAL { ?v pavs:hgvsC ?hgvsC }
    BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene)
  }
}
ORDER BY ?id
LIMIT 50`,
  },
  {
    key: 'constrainedGenes',
    query: `${PREFIXES}SELECT ?gene ?loeuf WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/genes> {
    ?geneUri pavs:loeuf ?loeuf .
    BIND(STRAFTER(STR(?geneUri), "hgnc.symbol/") AS ?gene)
  }
  FILTER(xsd:decimal(?loeuf) < 0.15)
}
ORDER BY xsd:decimal(?loeuf)
LIMIT 30`,
  },
  {
    key: 'diseaseHpo',
    query: `${PREFIXES}SELECT ?hpo ?label WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/hpoa> {
    ?assoc pavs:disease <https://omim.org/entry/272200> ;
           pavs:hpoTerm ?hpoTerm .
    BIND(REPLACE(STR(?hpoTerm), ".*HP_", "HP:") AS ?hpo)
  }
  OPTIONAL {
    GRAPH <http://pavs.kaust.edu.sa/graph/hpo-ic> {
      ?hpoTerm rdfs:label ?label .
    }
  }
}`,
  },
  {
    key: 'consanguineous',
    query: `${PREFIXES}SELECT ?id ?gene ?disease WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:consanguinity "Offspring of consanguineous parents" .
    OPTIONAL { ?case pavs:diseaseLabel ?disease }
    OPTIONAL {
      ?case pavs:hasVariant ?v .
      ?v pavs:affectsGene ?gUri .
      BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene)
    }
  }
}
LIMIT 50`,
  },
];

const SparqlExplorer: React.FC = () => {
  const { t } = useTranslation();
  const [query, setQuery] = useState(EXAMPLE_QUERIES[0].query);
  const [vars, setVars] = useState<string[]>([]);
  const [rows, setRows] = useState<Record<string, string>[]>([]);
  const [count, setCount] = useState<number | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const runQuery = async () => {
    setLoading(true);
    setError('');
    setRows([]);
    setVars([]);
    setCount(null);
    try {
      const res = await axios.post(`${API}/api/sparql`, { query });
      setVars(res.data.vars);
      setRows(res.data.rows);
      setCount(res.data.count);
    } catch (e: any) {
      setError(e.response?.data?.detail || e.message || 'Query failed');
    } finally {
      setLoading(false);
    }
  };

  const handleKeyDown = (e: React.KeyboardEvent<HTMLTextAreaElement>) => {
    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
      e.preventDefault();
      runQuery();
    }
  };

  const downloadResults = (format: 'json' | 'tsv') => {
    if (format === 'json') {
      const blob = new Blob([JSON.stringify({ vars, rows }, null, 2)], { type: 'application/json' });
      triggerDownload(blob, 'sparql_results.json');
    } else {
      const lines = [vars.join('\t'), ...rows.map(r => vars.map(v => r[v] ?? '').join('\t'))];
      const blob = new Blob([lines.join('\n')], { type: 'text/tab-separated-values' });
      triggerDownload(blob, 'sparql_results.tsv');
    }
  };

  // Resolve the public base URL at runtime so it reflects the actual host
  // (works for both localhost dev and production deployments).
  const publicBase = API || (typeof window !== 'undefined' ? window.location.origin : '');
  const apiEndpoint = `${publicBase}/api/sparql`;

  return (
    <div className="sparql-explorer">
      <h2>{t('sparql.title')}</h2>
      <p className="sparql-hint">{t('sparql.hint')}</p>

      {/* Endpoint info box */}
      <div className="sparql-endpoint-box">
        <span className="sparql-endpoint-label">{t('sparql.endpoint')}:</span>
        <code className="sparql-endpoint-url">{apiEndpoint}</code>
        <span className="sparql-endpoint-hint">{t('sparql.endpointHint')}</span>
        <code className="sparql-curl-example">{`curl '${apiEndpoint}?format=tsv&query=SELECT+...' > results.tsv`}</code>
      </div>

      <div className="sparql-examples">
        <strong>{t('sparql.examples')}:</strong>
        {EXAMPLE_QUERIES.map(ex => (
          <button
            key={ex.key}
            className="sparql-example-btn"
            onClick={() => { setQuery(ex.query); setRows([]); setVars([]); setCount(null); setError(''); }}
          >
            {t(`sparql.queries.${ex.key}`)}
          </button>
        ))}
      </div>

      <div className="sparql-editor-wrap">
        <textarea
          className="sparql-editor"
          value={query}
          onChange={e => setQuery(e.target.value)}
          onKeyDown={handleKeyDown}
          spellCheck={false}
          rows={14}
        />
        <button
          className="btn-primary sparql-run-btn"
          onClick={runQuery}
          disabled={loading}
        >
          {loading ? t('sparql.running') : t('sparql.run')}
        </button>
      </div>

      {error && <div className="error-msg">{error}</div>}

      {count !== null && (
        <div className="sparql-results">
          <div className="sparql-results-header">
            <p className="result-count">
              {t('sparql.rows', { count })}
              {count === 2000 ? t('sparql.limitReached') : ''}
            </p>
            {rows.length > 0 && (
              <div className="sparql-download-btns">
                <button className="sparql-dl-btn" onClick={() => downloadResults('json')}>
                  ⬇ {t('sparql.downloadJson')}
                </button>
                <button className="sparql-dl-btn" onClick={() => downloadResults('tsv')}>
                  ⬇ {t('sparql.downloadTsv')}
                </button>
              </div>
            )}
          </div>
          {rows.length > 0 && (
            <div className="results-table-wrapper">
              <table className="results-table sparql-results-table">
                <thead>
                  <tr>
                    {vars.map(v => <th key={v}>{v}</th>)}
                  </tr>
                </thead>
                <tbody>
                  {rows.map((row, idx) => (
                    <tr key={idx}>
                      {vars.map(v => (
                        <td key={v}>
                          {isUri(row[v])
                            ? <a href={row[v]} target="_blank" rel="noopener noreferrer"
                                 className="sparql-uri">{shortenUri(row[v])}</a>
                            : (row[v] ?? '')}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

function triggerDownload(blob: Blob, filename: string) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function isUri(val: string | undefined): boolean {
  return !!val && (val.startsWith('http://') || val.startsWith('https://'));
}

function shortenUri(uri: string): string {
  const known: [string, string][] = [
    ['http://pavs.kaust.edu.sa/data/', 'pav:'],
    ['http://pavs.kaust.edu.sa/ontology/', 'pavs:'],
    ['http://identifiers.org/hgnc.symbol/', 'hgnc:'],
    ['http://purl.obolibrary.org/obo/HP_', 'HP:'],
    ['https://omim.org/entry/', 'OMIM:'],
    ['http://purl.org/dc/terms/', 'dc:'],
    ['http://www.w3.org/2000/01/rdf-schema#', 'rdfs:'],
  ];
  for (const [prefix, short] of known) {
    if (uri.startsWith(prefix)) return short + uri.slice(prefix.length);
  }
  const last = uri.split(/[/#]/).filter(Boolean).pop() || uri;
  return last.length < 60 ? last : uri.slice(0, 60) + '…';
}

export default SparqlExplorer;
