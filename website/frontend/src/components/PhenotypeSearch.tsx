import React, { useState, useCallback } from 'react';
import AsyncSelect from 'react-select/async';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL || 'http://localhost:8000';

interface HpoOption { value: string; label: string }
interface CaseResult {
  id: string;
  gene: string;
  disease: string;
  source: string;
  score: number;
  is_saudi: boolean;
}

const PhenotypeSearch: React.FC = () => {
  const { t } = useTranslation();
  const [selectedHpos, setSelectedHpos] = useState<HpoOption[]>([]);
  const [method, setMethod] = useState<'lin' | 'resnik'>('lin');
  const [includeDisease, setIncludeDisease] = useState(false);
  const [includeNonSaudi, setIncludeNonSaudi] = useState(false);
  const [limit, setLimit] = useState(20);
  const [results, setResults] = useState<CaseResult[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const loadHpoOptions = useCallback(async (inputValue: string): Promise<HpoOption[]> => {
    if (inputValue.length < 2) return [];
    try {
      const res = await axios.get(`${API}/api/search/hpo`, { params: { q: inputValue } });
      return res.data.map((d: any) => ({ value: d.id, label: `${d.id} — ${d.label}` }));
    } catch {
      return [];
    }
  }, []);

  const handleSearch = async () => {
    if (!selectedHpos.length) return;
    setLoading(true);
    setError('');
    try {
      const res = await axios.post(`${API}/api/search/phenotype`, {
        hpo_ids: selectedHpos.map(h => h.value),
        method,
        limit,
        include_disease_phenotypes: includeDisease,
        include_non_saudi: includeNonSaudi,
      });
      setResults(res.data);
    } catch (e: any) {
      setError(e.message || 'Search failed');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="search-panel">
      <h2>{t('nav.phenotype')}</h2>

      <div className="form-group">
        <label>{t('search.placeholder')}</label>
        <AsyncSelect
          isMulti
          cacheOptions
          loadOptions={loadHpoOptions}
          onChange={(opts) => setSelectedHpos(opts as HpoOption[])}
          value={selectedHpos}
          placeholder={t('search.placeholder')}
          noOptionsMessage={() => t('search.noResults')}
          loadingMessage={() => t('search.loading')}
          className="hpo-select"
        />
      </div>

      <div className="search-options">
        <div className="form-group">
          <label>{t('search.method')}</label>
          <select value={method} onChange={e => setMethod(e.target.value as any)}>
            <option value="lin">{t('search.lin')}</option>
            <option value="resnik">{t('search.resnik')}</option>
          </select>
        </div>

        <div className="form-group">
          <label>
            <input type="checkbox" checked={includeDisease}
              onChange={e => setIncludeDisease(e.target.checked)} />
            {' '}{t('search.includeDisease')}
          </label>
        </div>

        <div className="form-group">
          <label>
            <input type="checkbox" checked={includeNonSaudi}
              onChange={e => setIncludeNonSaudi(e.target.checked)} />
            {' '}{t('search.includeNonSaudi')}
          </label>
        </div>

        <div className="form-group">
          <label>{t('search.limit')}</label>
          <select value={limit} onChange={e => setLimit(Number(e.target.value))}>
            {[10, 20, 50, 100].map(n => <option key={n} value={n}>{n}</option>)}
          </select>
        </div>
      </div>

      <div className="search-actions">
        <button className="btn-primary" onClick={handleSearch} disabled={loading || !selectedHpos.length}>
          {loading ? t('search.loading') : t('search.search')}
        </button>
        <button className="btn-secondary" onClick={() => { setSelectedHpos([]); setResults([]); }}>
          {t('search.clear')}
        </button>
      </div>

      {error && <div className="error-msg">{error}</div>}

      {results.length > 0 && (
        <div className="results-table-wrapper">
          <table className="results-table">
            <thead>
              <tr>
                <th>{t('results.score')}</th>
                <th>ID</th>
                <th>{t('results.gene')}</th>
                <th>{t('results.disease')}</th>
                <th>{t('results.source')}</th>
                <th></th>
              </tr>
            </thead>
            <tbody>
              {results.map(r => (
                <tr key={r.id}>
                  <td><span className="score-badge">{r.score.toFixed(4)}</span></td>
                  <td>
                    <a href={`/case/${encodeURIComponent(r.id)}`} className="case-link">
                      {r.id}
                    </a>
                  </td>
                  <td>{r.gene || '—'}</td>
                  <td className="disease-cell">{r.disease || '—'}</td>
                  <td><SourceBadge source={r.source} /></td>
                  <td>
                    <a href={`${API}/api/phenopacket/${encodeURIComponent(r.id)}/download`}
                       className="download-link" title={t('results.download')} download>
                      ↓
                    </a>
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

export default PhenotypeSearch;
