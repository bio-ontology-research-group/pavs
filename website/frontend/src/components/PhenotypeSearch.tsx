import React, { useState, useCallback, useRef, useEffect } from 'react';
import AsyncSelect from 'react-select/async';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL ?? '';
const PAGE_SIZE = 20;

interface HpoOption { value: string; label: string }
interface CaseResult {
  id: string;
  gene: string;
  disease: string;
  suggested_disease: string;
  source: string;
  score: number;
  is_saudi: boolean;
}

const PhenotypeSearch: React.FC = () => {
  const { t } = useTranslation();
  const [selectedHpos, setSelectedHpos] = useState<HpoOption[]>([]);
  const [method, setMethod] = useState<'lin' | 'resnik'>('lin');
  const [includeDisease, setIncludeDisease] = useState(false);
  const [includeSaudi, setIncludeSaudi] = useState(true);
  const [includeDDD, setIncludeDDD] = useState(false);
  const [includeLiterature, setIncludeLiterature] = useState(false);
  const [includeClinVar, setIncludeClinVar] = useState(false);
  const [fetchLimit, setFetchLimit] = useState(200);
  const [results, setResults] = useState<CaseResult[]>([]);
  const [displayCount, setDisplayCount] = useState(PAGE_SIZE);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const sentinelRef = useRef<HTMLDivElement>(null);

  // Infinite scroll: when sentinel enters viewport, show more rows
  useEffect(() => {
    if (!sentinelRef.current) return;
    const observer = new IntersectionObserver(
      (entries) => {
        if (entries[0].isIntersecting) {
          setDisplayCount(prev => Math.min(prev + PAGE_SIZE, results.length));
        }
      },
      { threshold: 0.1 }
    );
    observer.observe(sentinelRef.current);
    return () => observer.disconnect();
  }, [results.length]);

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
    setDisplayCount(PAGE_SIZE);
    try {
      const res = await axios.post(`${API}/api/search/phenotype`, {
        hpo_ids: selectedHpos.map(h => h.value),
        method,
        limit: fetchLimit,
        include_disease_phenotypes: includeDisease,
        include_saudi: includeSaudi,
        include_ddd: includeDDD,
        include_literature: includeLiterature,
        include_clinvar: includeClinVar,
      });
      setResults(res.data);
    } catch (e: any) {
      setError(e.message || 'Search failed');
    } finally {
      setLoading(false);
    }
  };

  const isClinVarCase = (id: string) => id.startsWith('ClinVar:');
  const visible = results.slice(0, displayCount);
  const hasMore = displayCount < results.length;

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
          <label style={{ fontWeight: 600 }}>{t('search.cohorts')}</label>
          <div style={{ display: 'flex', gap: '1.2rem', flexWrap: 'wrap', marginTop: '0.3rem' }}>
            <label>
              <input type="checkbox" checked={includeSaudi}
                onChange={e => setIncludeSaudi(e.target.checked)} />
              {' '}{t('search.includeSaudi')}
            </label>
            <label>
              <input type="checkbox" checked={includeDDD}
                onChange={e => setIncludeDDD(e.target.checked)} />
              {' '}{t('search.includeDDD')}
            </label>
            <label>
              <input type="checkbox" checked={includeLiterature}
                onChange={e => setIncludeLiterature(e.target.checked)} />
              {' '}{t('search.includeLiterature')}
            </label>
            <label>
              <input type="checkbox" checked={includeClinVar}
                onChange={e => setIncludeClinVar(e.target.checked)} />
              {' '}{t('search.includeClinVar')}
            </label>
          </div>
        </div>

        <div className="form-group">
          <label>{t('search.limit')}</label>
          <select value={fetchLimit} onChange={e => setFetchLimit(Number(e.target.value))}>
            {[50, 100, 200, 500].map(n => <option key={n} value={n}>{n}</option>)}
          </select>
        </div>
      </div>

      <div className="search-actions">
        <button className="btn-primary" onClick={handleSearch} disabled={loading || !selectedHpos.length}>
          {loading ? t('search.loading') : t('search.search')}
        </button>
        <button className="btn-secondary" onClick={() => { setSelectedHpos([]); setResults([]); setDisplayCount(PAGE_SIZE); }}>
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
              {visible.map(r => (
                <tr key={r.id}>
                  <td><span className="score-badge">{r.score.toFixed(4)}</span></td>
                  <td>
                    {isClinVarCase(r.id) ? (
                      <span className="case-id">{r.id}</span>
                    ) : (
                      <a href={`/case/${encodeURIComponent(r.id)}`} className="case-link">
                        {r.id}
                      </a>
                    )}
                  </td>
                  <td>{r.gene || '—'}</td>
                  <td className="disease-cell">
                    {r.disease
                      ? r.disease
                      : r.suggested_disease
                        ? <span title={t('case.suggestedDiseaseNote')} style={{ color: '#b07000', fontStyle: 'italic' }}>⚠ {r.suggested_disease}</span>
                        : '—'}
                  </td>
                  <td><SourceBadge source={r.source} /></td>
                  <td>
                    {!isClinVarCase(r.id) && (
                      <a href={`${API}/api/phenopacket/${encodeURIComponent(r.id)}/download`}
                         className="download-link" title={t('results.download')} download>
                        ↓
                      </a>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>

          <div style={{ padding: '0.5rem 0', color: '#888', fontSize: '0.85rem', textAlign: 'center' }}>
            {t('search.showing', { count: visible.length, total: results.length })}
          </div>

          {/* Sentinel for infinite scroll */}
          <div ref={sentinelRef} style={{ height: 1 }} />

          {hasMore && (
            <div style={{ textAlign: 'center', padding: '0.5rem' }}>
              <button className="btn-secondary"
                onClick={() => setDisplayCount(c => Math.min(c + PAGE_SIZE, results.length))}>
                {t('search.showMore')}
              </button>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default PhenotypeSearch;
