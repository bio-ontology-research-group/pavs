import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL || 'http://localhost:8000';

type SearchTab = 'gene' | 'rsid' | 'hgvs' | 'acmg';

interface VariantResult {
  id?: string;
  case?: string;
  gene?: string;
  disease?: string;
  source?: string;
  isSaudi?: string;
  hgvsC?: string;
  hgvsG?: string;
  acmg?: string;
  saudiAF?: string;
  rsId?: string;
  togovar_url?: string;
}

const VariantLookup: React.FC = () => {
  const { t } = useTranslation();
  const [activeTab, setActiveTab] = useState<SearchTab>('gene');
  const [gene, setGene] = useState('');
  const [rsid, setRsid] = useState('');
  const [hgvs, setHgvs] = useState('');
  const [acmgClass, setAcmgClass] = useState('Pathogenic');
  const [results, setResults] = useState<VariantResult[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const handleSearch = async () => {
    setLoading(true);
    setError('');
    setResults([]);
    try {
      const params: Record<string, string> = {};
      if (activeTab === 'gene') params.gene = gene;
      else if (activeTab === 'rsid') params.rsid = rsid;
      else if (activeTab === 'hgvs') params.hgvs = hgvs;
      else if (activeTab === 'acmg') params.acmg = acmgClass;

      const res = await axios.get(`${API}/api/search/variant`, { params });
      setResults(res.data);
    } catch (e: any) {
      setError(e.response?.data?.detail || e.message || 'Search failed');
    } finally {
      setLoading(false);
    }
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') handleSearch();
  };

  return (
    <div className="search-panel">
      <h2>{t('nav.variant')}</h2>

      <div className="tab-bar">
        {(['gene', 'rsid', 'hgvs', 'acmg'] as SearchTab[]).map(tab => (
          <button key={tab} className={`tab-btn ${activeTab === tab ? 'active' : ''}`}
            onClick={() => setActiveTab(tab)}>
            {t(`variant.by${tab.charAt(0).toUpperCase() + tab.slice(1)}`)}
          </button>
        ))}
      </div>

      <div className="form-group">
        {activeTab === 'gene' && (
          <input type="text" value={gene} onChange={e => setGene(e.target.value)}
            onKeyDown={handleKeyDown}
            placeholder={t('variant.genePlaceholder')} className="text-input" />
        )}
        {activeTab === 'rsid' && (
          <input type="text" value={rsid} onChange={e => setRsid(e.target.value)}
            onKeyDown={handleKeyDown}
            placeholder={t('variant.rsidPlaceholder')} className="text-input" />
        )}
        {activeTab === 'hgvs' && (
          <input type="text" value={hgvs} onChange={e => setHgvs(e.target.value)}
            onKeyDown={handleKeyDown}
            placeholder={t('variant.hgvsPlaceholder')} className="text-input" />
        )}
        {activeTab === 'acmg' && (
          <select value={acmgClass} onChange={e => setAcmgClass(e.target.value)} className="text-input">
            {['Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign'].map(c => (
              <option key={c} value={c}>{c}</option>
            ))}
          </select>
        )}
      </div>

      <button className="btn-primary" onClick={handleSearch} disabled={loading}>
        {loading ? t('search.loading') : t('search.search')}
      </button>

      {error && <div className="error-msg">{error}</div>}

      {results.length > 0 && (
        <div className="results-table-wrapper">
          <p className="result-count">{results.length} results</p>
          <table className="results-table">
            <thead>
              <tr>
                <th>Case ID</th>
                <th>{t('results.gene')}</th>
                <th>{t('results.hgvsC')}</th>
                <th>{t('results.acmg')}</th>
                <th>{t('results.saudiAF')}</th>
                <th>{t('results.source')}</th>
                <th>Links</th>
              </tr>
            </thead>
            <tbody>
              {results.map((r, idx) => (
                <tr key={idx}>
                  <td>
                    <a href={`/case/${encodeURIComponent(r.id || r.case || '')}`} className="case-link">
                      {r.id || r.case || '—'}
                    </a>
                  </td>
                  <td>
                    {r.gene ? (
                      <a href={`/gene/${r.gene}`} className="gene-link">{r.gene}</a>
                    ) : '—'}
                  </td>
                  <td className="mono">{r.hgvsC || '—'}</td>
                  <td>
                    <span className={`acmg-badge acmg-${(r.acmg || '').toLowerCase().replace(/\s+/g, '-')}`}>
                      {r.acmg || '—'}
                    </span>
                  </td>
                  <td>
                    {r.saudiAF === '0' || r.saudiAF === '0.0'
                      ? <span className="novel-badge">{t('results.novel')}</span>
                      : (r.saudiAF || '—')}
                  </td>
                  <td><SourceBadge source={r.source || ''} /></td>
                  <td className="links-cell">
                    {r.togovar_url && (
                      <a href={r.togovar_url} target="_blank" rel="noopener noreferrer"
                         className="ext-link togovar-link" title={t('links.togovar')}>
                        TGV
                      </a>
                    )}
                    {r.rsId && (
                      <a href={`https://www.ncbi.nlm.nih.gov/snp/${r.rsId}`}
                         target="_blank" rel="noopener noreferrer" className="ext-link" title={t('links.dbsnp')}>
                        dbSNP
                      </a>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
};

export default VariantLookup;
