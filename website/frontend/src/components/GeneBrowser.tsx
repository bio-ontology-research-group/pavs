import React, { useState, useEffect, useRef } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL ?? '';

interface GeneListItem {
  symbol: string;
  saudi_count: number;
}

interface GeneDetail {
  symbol: string;
  ncbiGeneId?: string;
  pli?: string;
  loeuf?: string;
  goProcess?: string;
  expressedIn?: string;
  relatedDisease?: string;
  hgnc_url?: string;
  [key: string]: string | undefined;
}

interface GeneCaseResult {
  id?: string;
  case?: string;
  gene?: string;
  disease?: string;
  source?: string;
  isSaudi?: string;
  hgvsC?: string;
  acmg?: string;
  saudiAF?: string;
  rsId?: string;
  clinvarId?: string;
  clinvarSig?: string;
  togovar_url?: string;
}

interface HpoTerm {
  id: string;
  label: string;
}

interface DiseaseInfo {
  disease_id: string;
  label: string;
  omim_url: string;
  clinvar_url: string;
  hpo_terms: HpoTerm[];
}

interface Props {
  initialGene?: string;
}

// ─── Expression bar chart ────────────────────────────────────────────────────

interface Tissue { name: string; val: number }

const ExpressionChart: React.FC<{ expressedIn: string }> = ({ expressedIn }) => {
  const [tooltip, setTooltip] = useState<{ x: number; y: number; tissue: string; val: number } | null>(null);

  const tissues: Tissue[] = expressedIn
    .split('|')
    .map(t => { const [name, v] = t.split(':'); return { name: name.replace(/_/g, ' '), val: Number(v) }; })
    .filter(t => !isNaN(t.val));

  if (!tissues.length) return null;

  const maxVal = Math.max(...tissues.map(t => t.val));
  const barW = 28;
  const gap = 3;
  const chartH = 140;
  const labelH = 160;
  const svgW = tissues.length * (barW + gap);
  const svgH = chartH + labelH;

  return (
    <div className="expression-chart-wrap">
      <svg width={svgW} height={svgH} className="expression-svg" style={{ overflow: 'visible' }}>
        {tissues.map((t, i) => {
          const x = i * (barW + gap);
          const barH = maxVal > 0 ? (t.val / maxVal) * chartH : 0;
          const y = chartH - barH;
          const cx = x + barW / 2;
          return (
            <g
              key={t.name + i}
              onMouseEnter={e => setTooltip({ x: e.clientX, y: e.clientY, tissue: t.name, val: t.val })}
              onMouseLeave={() => setTooltip(null)}
              style={{ cursor: 'default' }}
            >
              <rect x={x} y={y} width={barW} height={barH}
                fill="var(--color-primary)" opacity="0.75" rx="2" />
              <text
                x={cx} y={chartH + 6}
                fontSize="10" textAnchor="end" fill="var(--color-text-light)"
                transform={`rotate(-45, ${cx}, ${chartH + 6})`}
              >
                {t.name}
              </text>
            </g>
          );
        })}
        {/* Y-axis label */}
        <text x={-4} y={8} fontSize="9" textAnchor="end" fill="var(--color-text-light)">100%</text>
        <text x={-4} y={chartH / 2 + 4} fontSize="9" textAnchor="end" fill="var(--color-text-light)">50%</text>
        <line x1={0} y1={0} x2={0} y2={chartH} stroke="var(--color-border)" strokeWidth="1" />
        <line x1={0} y1={chartH / 2} x2={svgW} y2={chartH / 2} stroke="var(--color-border)" strokeWidth="0.5" strokeDasharray="3,3" />
      </svg>
      {tooltip && (
        <div
          className="expression-tooltip"
          style={{ position: 'fixed', left: tooltip.x + 12, top: tooltip.y - 36, pointerEvents: 'none' }}
        >
          <strong>{tooltip.tissue}</strong>: {tooltip.val} (normalized {maxVal > 0 ? Math.round((tooltip.val / maxVal) * 100) : 0}%)
        </div>
      )}
    </div>
  );
};

// ─── Main component ──────────────────────────────────────────────────────────

/** Localize ACMG class labels and normalize for CSS classes. */
function localizeAcmg(acmg: string, t: any): { label: string; className: string } {
  const raw = acmg.toLowerCase().replace(/[^a-z]+/g, '_');
  const key = raw === 'uncertain_significance' ? 'uncertain_significance' :
              raw === 'pathogenic' ? 'pathogenic' :
              raw === 'likely_pathogenic' ? 'likely_pathogenic' :
              raw === 'likely_benign' ? 'likely_benign' :
              raw === 'benign' ? 'benign' : raw;

  const label = t(`acmg.${key}`, { defaultValue: acmg });
  return { label, className: `acmg-${key.replace(/_/g, '-')}` };
}

const GeneBrowser: React.FC<Props> = ({ initialGene = '' }) => {
  const { t } = useTranslation();

  const [geneList, setGeneList] = useState<GeneListItem[]>([]);
  const [geneListLoading, setGeneListLoading] = useState(true);
  const [filter, setFilter] = useState('');

  const [selectedGene, setSelectedGene] = useState<string>(initialGene);
  const [geneDetail, setGeneDetail] = useState<GeneDetail | null>(null);
  const [geneCases, setGeneCases] = useState<GeneCaseResult[]>([]);
  const [geneDiseases, setGeneDiseases] = useState<DiseaseInfo[]>([]);
  const [detailLoading, setDetailLoading] = useState(false);
  const [detailError, setDetailError] = useState('');

  const [showDDD, setShowDDD] = useState(false);
  const [showLiterature, setShowLiterature] = useState(false);
  const [expandedDiseases, setExpandedDiseases] = useState<Set<string>>(new Set());

  const selectedRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    axios.get(`${API}/api/genes`)
      .then(res => setGeneList(res.data))
      .catch(err => console.error('Failed to load gene list:', err))
      .finally(() => setGeneListLoading(false));
  }, []);

  useEffect(() => {
    if (initialGene) {
      setSelectedGene(initialGene);
      setFilter('');
    }
  }, [initialGene]);

  useEffect(() => {
    if (!selectedGene) return;
    setDetailLoading(true);
    setDetailError('');
    setGeneDetail(null);
    setGeneCases([]);
    setGeneDiseases([]);
    setExpandedDiseases(new Set());

    const params = { include_ddd: showDDD, include_literature: showLiterature };

    Promise.all([
      axios.get(`${API}/api/gene/${encodeURIComponent(selectedGene)}`),
      axios.get(`${API}/api/gene/${encodeURIComponent(selectedGene)}/cases`, { params }),
      axios.get(`${API}/api/gene/${encodeURIComponent(selectedGene)}/diseases`),
    ])
      .then(([detailRes, casesRes, diseasesRes]) => {
        setGeneDetail(detailRes.data);
        setGeneCases(casesRes.data);
        setGeneDiseases(diseasesRes.data);
      })
      .catch(err => setDetailError(err.response?.data?.detail || err.message || 'Failed to load gene'))
      .finally(() => setDetailLoading(false));
  }, [selectedGene]);

  // Re-fetch cases only when toggles change (not diseases/detail)
  useEffect(() => {
    if (!selectedGene) return;
    const params = { include_ddd: showDDD, include_literature: showLiterature };
    axios.get(`${API}/api/gene/${encodeURIComponent(selectedGene)}/cases`, { params })
      .then(res => setGeneCases(res.data))
      .catch(() => {});
  }, [showDDD, showLiterature]);

  useEffect(() => {
    if (selectedRef.current) {
      selectedRef.current.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    }
  }, [selectedGene, geneList]);

  const filteredGenes = filter
    ? geneList.filter(g => g.symbol.toLowerCase().includes(filter.toLowerCase()))
    : geneList;

  // Client-side cohort filter (DDD cases with isSaudi=false can appear in graph/cases)
  const displayedCases = geneCases.filter(c => {
    const isSaudi = c.isSaudi === '1' || c.isSaudi === 'true';
    const isDDD = c.source?.toLowerCase().includes('ddd');
    if (isSaudi) return true;
    if (isDDD) return showDDD;
    return showLiterature;
  });

  const toggleDisease = (id: string) => {
    setExpandedDiseases(prev => {
      const next = new Set(prev);
      next.has(id) ? next.delete(id) : next.add(id);
      return next;
    });
  };

  return (
    <div className="gene-browser-layout">
      {/* Left panel */}
      <div className="gene-list-panel">
        <div className="gene-list-filter">
          <input
            type="text"
            value={filter}
            onChange={e => setFilter(e.target.value)}
            placeholder={t('gene.filterPlaceholder')}
            className="text-input"
          />
        </div>
        <div className="gene-list-meta">
          {geneListLoading
            ? <span className="loading" style={{ padding: '8px 14px', display: 'block' }}>{t('search.loading')}</span>
            : <span>{filteredGenes.length} {t('gene.genesLabel')}</span>
          }
        </div>
        <div className="gene-list-scroll">
          {filteredGenes.map(g => (
            <div
              key={g.symbol}
              ref={selectedGene === g.symbol ? selectedRef : null}
              className={`gene-list-item ${selectedGene === g.symbol ? 'active' : ''}`}
              onClick={() => setSelectedGene(g.symbol)}
            >
              <span className="gene-symbol">{g.symbol}</span>
              <span className="gene-count">{g.saudi_count}</span>
            </div>
          ))}
        </div>
      </div>

      {/* Right panel */}
      <div className="gene-detail-panel">
        {!selectedGene && <div className="gene-select-prompt">{t('gene.selectPrompt')}</div>}
        {selectedGene && detailLoading && <div className="loading">{t('search.loading')}</div>}
        {selectedGene && detailError && <div className="error-msg">{detailError}</div>}

        {selectedGene && geneDetail && !detailLoading && (
          <>
            {/* Gene detail card */}
            <div className="gene-detail-card">
              <h3>{geneDetail.symbol}</h3>
              <dl className="prop-list">
                {geneDetail.ncbiGeneId && (
                  <><dt>NCBI Gene</dt>
                  <dd><a href={`https://www.ncbi.nlm.nih.gov/gene/${geneDetail.ncbiGeneId.replace('http://identifiers.org/ncbigene/', '')}`}
                         target="_blank" rel="noopener noreferrer">
                    {geneDetail.ncbiGeneId.replace('http://identifiers.org/ncbigene/', '')}
                  </a></dd></>
                )}
                {geneDetail.pli && (
                  <><dt>
                    <span title="Probability of being Loss-of-function Intolerant. Values &gt; 0.9 indicate the gene is highly intolerant to heterozygous loss-of-function variants (gnomAD).">
                      gnomAD pLI ⓘ
                    </span>
                  </dt>
                  <dd>{parseFloat(geneDetail.pli).toExponential(2)}</dd></>
                )}
                {geneDetail.loeuf && (
                  <><dt>
                    <span title="Loss-of-function Observed/Expected Upper bound Fraction. Lower values (< 0.35) indicate stronger selective constraint against loss-of-function variants (gnomAD).">
                      gnomAD LOEUF ⓘ
                    </span>
                  </dt>
                  <dd>{parseFloat(geneDetail.loeuf).toFixed(3)}</dd></>
                )}
                {geneDetail.goProcess && (
                  <><dt>GO Biological Process</dt>
                  <dd>
                    <ul className="go-list">
                      {geneDetail.goProcess.split('|').map(term => (
                        <li key={term}>
                          <a
                            href={`https://www.ebi.ac.uk/QuickGO/search/${encodeURIComponent(term)}`}
                            target="_blank" rel="noopener noreferrer"
                          >
                            {term}
                          </a>
                        </li>
                      ))}
                    </ul>
                  </dd></>
                )}
                {geneDetail.expressedIn && (
                  <><dt>GTEx expression</dt>
                  <dd style={{ overflowX: 'auto', overflowY: 'visible', paddingBottom: '8px' }}>
                    <div style={{ minWidth: 'max-content', overflowY: 'visible' }}>
                      <ExpressionChart expressedIn={geneDetail.expressedIn} />
                    </div>
                  </dd></>
                )}
              </dl>
              <div className="gene-links">
                {geneDetail.hgnc_url && (
                  <a href={geneDetail.hgnc_url} target="_blank" rel="noopener noreferrer" className="ext-link">{t('links.hgnc')}</a>
                )}
                <a href={`https://www.uniprot.org/uniprotkb?query=gene%3A${geneDetail.symbol}&reviewed=true&taxonomy_id=9606`}
                   target="_blank" rel="noopener noreferrer" className="ext-link">{t('links.uniprot')}</a>
                <a href={`https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${geneDetail.symbol}`}
                   target="_blank" rel="noopener noreferrer" className="ext-link">{t('links.ensembl')}</a>
                <a href={`https://decipher.sanger.ac.uk/search?q=${geneDetail.symbol}`}
                   target="_blank" rel="noopener noreferrer" className="ext-link">{t('links.decipher')}</a>
                <a href={`https://omim.org/search?search=${geneDetail.symbol}`}
                   target="_blank" rel="noopener noreferrer" className="ext-link">{t('links.omim')}</a>
                <a href={`https://www.ncbi.nlm.nih.gov/clinvar/?term=${geneDetail.symbol}[gene]`}
                   target="_blank" rel="noopener noreferrer" className="ext-link">ClinVar</a>
              </div>
            </div>

            {/* Associated diseases with HPO phenotypes */}
            {geneDiseases.length > 0 && (
              <div className="gene-diseases-section">
                <h4>Associated Diseases</h4>
                {geneDiseases.map(d => (
                  <div key={d.disease_id} className="gene-disease-item">
                    <div className="gene-disease-header">
                      <button
                        className="gene-disease-toggle"
                        onClick={() => toggleDisease(d.disease_id)}
                        aria-expanded={expandedDiseases.has(d.disease_id)}
                      >
                        {expandedDiseases.has(d.disease_id) ? '▼' : '▶'}
                      </button>
                      <span className="gene-disease-label">{d.label}</span>
                      <span className="gene-disease-id">{d.disease_id}</span>
                      <a href={d.omim_url} target="_blank" rel="noopener noreferrer" className="ext-link">OMIM</a>
                      <a href={d.clinvar_url} target="_blank" rel="noopener noreferrer" className="ext-link">ClinVar</a>
                    </div>
                    {expandedDiseases.has(d.disease_id) && d.hpo_terms.length > 0 && (
                      <div className="gene-disease-hpo">
                        {d.hpo_terms.map(hp => (
                          <a
                            key={hp.id}
                            href={`https://hpo.jax.org/browse/term/${hp.id}`}
                            target="_blank" rel="noopener noreferrer"
                            className="hpo-tag-chip"
                            style={{ display: 'inline-flex', margin: '2px' }}
                          >
                            {hp.label || hp.id} ({hp.id})
                          </a>
                        ))}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            )}

            {/* Source toggles */}
            <div className="gene-source-toggles">
              <strong>{t('gene.showCasesFrom')}:</strong>
              <label className="toggle-label">
                <input type="checkbox" checked disabled readOnly />
                {' '}{t('badges.saudi')}
              </label>
              <label className="toggle-label">
                <input type="checkbox" checked={showDDD} onChange={e => setShowDDD(e.target.checked)} />
                {' '}{t('badges.ddd')}
              </label>
              <label className="toggle-label">
                <input type="checkbox" checked={showLiterature} onChange={e => setShowLiterature(e.target.checked)} />
                {' '}{t('badges.literature')}
              </label>
            </div>

            {/* Cases table */}
            {displayedCases.length > 0 ? (
              <div className="results-table-wrapper">
                <p className="result-count">{displayedCases.length} {t('gene.casesLabel')}</p>
                <table className="results-table">
                  <thead>
                    <tr>
                      <th>{t('results.caseId')}</th>
                      <th>{t('results.gene')}</th>
                      <th>{t('results.hgvsC')}</th>
                      <th>{t('results.acmg')}</th>
                      <th>{t('results.saudiAF')}</th>
                      <th>{t('results.disease')}</th>
                      <th>{t('results.source')}</th>
                      <th>{t('results.links')}</th>
                    </tr>
                  </thead>
                  <tbody>
                    {displayedCases.map((r, idx) => (
                      <tr key={idx}>
                        <td>
                          <a href={`/case/${encodeURIComponent(r.id || r.case || '')}`} className="case-link">
                            {r.id || r.case || '—'}
                          </a>
                        </td>
                        <td>
                          {r.gene
                            ? <a href={`/?gene=${encodeURIComponent(r.gene)}`} className="gene-link bold">{r.gene}</a>
                            : '—'}
                        </td>
                        <td className="mono">{r.hgvsC || '—'}</td>
                        <td>
                          {r.acmg ? (() => {
                            const { label, className } = localizeAcmg(r.acmg, t);
                            return (
                              <span className={`acmg-badge ${className}`}>
                                {label}
                              </span>
                            );
                          })() : '—'}
                        </td>
                        <td>
                          {r.saudiAF === '0' || r.saudiAF === '0.0'
                            ? <span className="novel-badge">{t('results.novel')}</span>
                            : (r.saudiAF || '—')}
                        </td>
                        <td className="disease-cell">{r.disease || '—'}</td>
                        <td><SourceBadge source={r.source || ''} /></td>
                        <td className="links-cell">
                          {r.togovar_url && (
                            <a href={r.togovar_url} target="_blank" rel="noopener noreferrer"
                               className="ext-link togovar-link" title={t('links.togovar')}>TGV</a>
                          )}
                          {r.rsId && (
                            <a href={`https://www.ncbi.nlm.nih.gov/snp/${r.rsId}`}
                               target="_blank" rel="noopener noreferrer" className="ext-link">dbSNP</a>
                          )}
                          {r.clinvarId && (
                            <a href={`https://www.ncbi.nlm.nih.gov/clinvar/variation/${r.clinvarId}`}
                               target="_blank" rel="noopener noreferrer" className="ext-link">ClinVar</a>
                          )}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            ) : (
              <p className="result-count" style={{ marginTop: '12px' }}>{t('search.noResults')}</p>
            )}
          </>
        )}
      </div>
    </div>
  );
};

export default GeneBrowser;
