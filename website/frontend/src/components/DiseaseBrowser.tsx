import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL ?? '';

interface DiseaseResult {
  case?: string;
  id?: string;
  disease?: string;
  source?: string;
  gene?: string;
  isSaudi?: string;
}

const DiseaseBrowser: React.FC = () => {
  const { t } = useTranslation();
  const [query, setQuery] = useState('');
  const [results, setResults] = useState<DiseaseResult[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const handleSearch = async () => {
    if (!query.trim()) return;
    setLoading(true);
    setError('');
    try {
      const res = await axios.get(`${API}/api/search/disease`, { params: { q: query.trim() } });
      setResults(res.data);
    } catch (e: any) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="search-panel">
      <h2>{t('nav.disease')}</h2>
      <div className="form-group inline">
        <input
          type="text" value={query}
          onChange={e => setQuery(e.target.value)}
          onKeyDown={e => e.key === 'Enter' && handleSearch()}
          placeholder="Disease name or OMIM ID…"
          className="text-input"
        />
        <button className="btn-primary" onClick={handleSearch} disabled={loading}>
          {loading ? t('search.loading') : t('search.search')}
        </button>
      </div>

      {error && <div className="error-msg">{error}</div>}

      {results.length > 0 && (
        <div className="results-table-wrapper">
          <p className="result-count">{results.length} cases</p>
          <table className="results-table">
            <thead>
              <tr>
                <th>Case ID</th>
                <th>{t('results.disease')}</th>
                <th>{t('results.gene')}</th>
                <th>{t('results.source')}</th>
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
                  <td>{r.disease || '—'}</td>
                  <td>{r.gene || '—'}</td>
                  <td><SourceBadge source={r.source || ''} /></td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
};

export default DiseaseBrowser;
