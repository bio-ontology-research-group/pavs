import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';

const API = import.meta.env.VITE_API_URL || 'http://localhost:8000';

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

const GeneBrowser: React.FC = () => {
  const { t } = useTranslation();
  const [query, setQuery] = useState('');
  const [result, setResult] = useState<GeneDetail | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const handleSearch = async () => {
    if (!query.trim()) return;
    setLoading(true);
    setError('');
    try {
      const res = await axios.get(`${API}/api/gene/${query.trim()}`);
      setResult(res.data);
    } catch (e: any) {
      setError(e.response?.status === 404 ? `Gene "${query}" not found` : e.message);
      setResult(null);
    } finally {
      setLoading(false);
    }
  };

  const topTissues = (expressed: string) =>
    expressed.split('|').slice(0, 5).join(' · ');

  return (
    <div className="search-panel">
      <h2>{t('nav.gene')}</h2>
      <div className="form-group inline">
        <input
          type="text" value={query}
          onChange={e => setQuery(e.target.value)}
          onKeyDown={e => e.key === 'Enter' && handleSearch()}
          placeholder="Gene symbol, e.g. SUMF1"
          className="text-input"
        />
        <button className="btn-primary" onClick={handleSearch} disabled={loading}>
          {loading ? t('search.loading') : t('search.search')}
        </button>
      </div>

      {error && <div className="error-msg">{error}</div>}

      {result && (
        <div className="gene-detail-card">
          <h3>{result.symbol}</h3>
          <dl className="prop-list">
            {result.ncbiGeneId && (
              <><dt>NCBI Gene</dt>
              <dd><a href={`https://www.ncbi.nlm.nih.gov/gene/${result.ncbiGeneId.replace('http://identifiers.org/ncbigene/', '')}`}
                     target="_blank" rel="noopener noreferrer">
                {result.ncbiGeneId.replace('http://identifiers.org/ncbigene/', '')}
              </a></dd></>
            )}
            {result.pli && <><dt>gnomAD pLI</dt><dd>{parseFloat(result.pli).toExponential(2)}</dd></>}
            {result.loeuf && <><dt>gnomAD LOEUF</dt><dd>{parseFloat(result.loeuf).toFixed(3)}</dd></>}
            {result.goProcess && (
              <><dt>GO Biological Process</dt><dd className="go-cell">{result.goProcess}</dd></>
            )}
            {result.expressedIn && (
              <><dt>Top expression</dt><dd>{topTissues(result.expressedIn)}</dd></>
            )}
          </dl>
          <div className="gene-links">
            {result.hgnc_url && (
              <a href={result.hgnc_url} target="_blank" rel="noopener noreferrer" className="ext-link">
                {t('links.hgnc')}
              </a>
            )}
            <a href={`https://www.uniprot.org/uniprotkb?query=gene%3A${result.symbol}&reviewed=true&taxonomy_id=9606`}
               target="_blank" rel="noopener noreferrer" className="ext-link">
              {t('links.uniprot')}
            </a>
            <a href={`https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${result.symbol}`}
               target="_blank" rel="noopener noreferrer" className="ext-link">
              {t('links.ensembl')}
            </a>
            <a href={`https://decipher.sanger.ac.uk/search?q=${result.symbol}`}
               target="_blank" rel="noopener noreferrer" className="ext-link">
              {t('links.decipher')}
            </a>
            <a href={`https://omim.org/search?search=${result.symbol}`}
               target="_blank" rel="noopener noreferrer" className="ext-link">
              {t('links.omim')}
            </a>
          </div>
        </div>
      )}
    </div>
  );
};

export default GeneBrowser;
