import React, { useState, useRef } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL ?? '';

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
  // Pathogenicity predictions
  consequence?: string;
  vepImpact?: string;
  sift?: string;
  polyphen?: string;
  gnomadAF?: string;
  zygosity?: string;
}

// Small inline info tooltip
function InfoTip({ text }: { text: string }) {
  const [open, setOpen] = useState(false);
  const ref = useRef<HTMLSpanElement>(null);
  return (
    <span
      ref={ref}
      style={{ position: 'relative', display: 'inline-block', marginLeft: '0.25rem', verticalAlign: 'middle' }}
      onMouseEnter={() => setOpen(true)}
      onMouseLeave={() => setOpen(false)}
      onFocus={() => setOpen(true)}
      onBlur={() => setOpen(false)}
    >
      <span
        tabIndex={0}
        role="button"
        aria-label="More information"
        style={{
          display: 'inline-flex', alignItems: 'center', justifyContent: 'center',
          width: '1.1em', height: '1.1em', borderRadius: '50%',
          border: '1.5px solid #7aa', color: '#4a9', background: '#f0faf8',
          fontSize: '0.68em', fontWeight: 700, fontStyle: 'normal',
          cursor: 'help', userSelect: 'none', lineHeight: 1,
        }}
      >i</span>
      {open && (
        <span style={{
          position: 'absolute', bottom: 'calc(100% + 6px)', left: '50%',
          transform: 'translateX(-50%)',
          background: '#1a2a2a', color: '#e8f4f0', borderRadius: '5px',
          padding: '0.45rem 0.65rem', fontSize: '0.75rem', lineHeight: 1.45,
          width: '240px', whiteSpace: 'normal', zIndex: 200,
          boxShadow: '0 3px 10px rgba(0,0,0,0.3)',
          pointerEvents: 'none',
        }}>
          {text}
        </span>
      )}
    </span>
  );
}

// Colour-coded consequence label
function ConsequenceBadge({ val }: { val: string }) {
  const high = ['stop_gained', 'frameshift_variant', 'splice_acceptor_variant',
                 'splice_donor_variant', 'start_lost', 'stop_lost'];
  const moderate = ['missense_variant', 'inframe_insertion', 'inframe_deletion',
                    'protein_altering_variant'];
  const cls = high.some(h => val.includes(h))
    ? 'badge-high'
    : moderate.some(m => val.includes(m))
    ? 'badge-moderate'
    : 'badge-low';
  const label = val.replace(/_variant$/, '').replace(/_/g, ' ');
  return <span className={`consequence-badge ${cls}`}>{label}</span>;
}

// SIFT / PolyPhen formatted label
function PredBadge({ label, val }: { label: string; val: string }) {
  const lower = val.toLowerCase();
  const bad = ['deleterious', 'probably_damaging', 'probably damaging'].some(k => lower.includes(k));
  const mid = ['tolerated', 'possibly_damaging', 'possibly damaging'].some(k => lower.includes(k));
  const cls = bad ? 'badge-high' : mid ? 'badge-moderate' : 'badge-low';
  return (
    <span>
      <strong>{label}:</strong>{' '}
      <span className={`consequence-badge ${cls}`}>{val.replace(/_/g, ' ')}</span>
    </span>
  );
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
  const [expandedRows, setExpandedRows] = useState<Set<number>>(new Set());

  const handleSearch = async () => {
    setLoading(true);
    setError('');
    setResults([]);
    setExpandedRows(new Set());
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

  const toggleRow = (idx: number) => {
    setExpandedRows(prev => {
      const next = new Set(prev);
      next.has(idx) ? next.delete(idx) : next.add(idx);
      return next;
    });
  };

  const hasPathogenicity = (r: VariantResult) =>
    !!(r.consequence || r.vepImpact || r.sift || r.polyphen || r.gnomadAF || r.zygosity);

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
                <th style={{ width: '24px' }}></th>
                <th>{t('results.caseId')}</th>
                <th>{t('results.gene')}</th>
                <th>{t('results.hgvsC')}</th>
                <th>{t('results.acmg')}</th>
                <th>{t('results.saudiAF')}</th>
                <th>{t('results.source')}</th>
                <th>{t('results.links')}</th>
              </tr>
            </thead>
            <tbody>
              {results.map((r, idx) => (
                <React.Fragment key={idx}>
                  <tr className={expandedRows.has(idx) ? 'variant-row-expanded' : ''}>
                    <td>
                      {hasPathogenicity(r) && (
                        <button
                          className="row-expand-btn"
                          onClick={() => toggleRow(idx)}
                          title="Show pathogenicity predictions"
                          aria-expanded={expandedRows.has(idx)}
                        >
                          {expandedRows.has(idx) ? '▼' : '▶'}
                        </button>
                      )}
                    </td>
                    <td>
                      <a href={`/case/${encodeURIComponent(r.id || r.case || '')}`} className="case-link">
                        {r.id || r.case || '—'}
                      </a>
                    </td>
                    <td>
                      {r.gene ? (
                        <a href={`/?gene=${encodeURIComponent(r.gene)}`} className="gene-link">{r.gene}</a>
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
                  {expandedRows.has(idx) && hasPathogenicity(r) && (
                    <tr className="variant-detail-row">
                      <td></td>
                      <td colSpan={7}>
                        <div className="variant-detail-panel">
                          {r.consequence && (
                            <div className="variant-detail-item">
                              <strong>VEP Consequence:</strong>{' '}
                              <ConsequenceBadge val={r.consequence} />
                            </div>
                          )}
                          {r.vepImpact && (
                            <div className="variant-detail-item">
                              <strong>Impact:
                                <InfoTip text="VEP impact rating predicts the severity of a variant's effect. HIGH = likely loss-of-function (e.g. stop-gain, frameshift). MODERATE = protein-altering (e.g. missense). LOW = likely tolerated (e.g. synonymous). MODIFIER = non-coding or uncertain effect." />
                              </strong>{' '}
                              <span className={`consequence-badge ${
                                r.vepImpact === 'HIGH' ? 'badge-high'
                                : r.vepImpact === 'MODERATE' ? 'badge-moderate'
                                : r.vepImpact === 'LOW' ? 'badge-low'
                                : ''
                              }`}>{r.vepImpact}</span>
                            </div>
                          )}
                          {r.sift && (
                            <div className="variant-detail-item">
                              <PredBadge label="SIFT" val={r.sift} />
                            </div>
                          )}
                          {r.polyphen && (
                            <div className="variant-detail-item">
                              <PredBadge label="PolyPhen-2" val={r.polyphen} />
                            </div>
                          )}
                          {r.gnomadAF && (
                            <div className="variant-detail-item">
                              <strong>gnomAD AF:
                                <InfoTip text="Allele frequency in gnomAD (Genome Aggregation Database) — a large reference population of ~800 000 exomes and genomes. Very low or absent values indicate a rare or novel variant. Displayed as a decimal (e.g. 0.0001 = 0.01%)." />
                              </strong>{' '}{r.gnomadAF}
                            </div>
                          )}
                          {r.zygosity && (
                            <div className="variant-detail-item">
                              <strong>Zygosity:</strong> {r.zygosity}
                            </div>
                          )}
                          {r.hgvsG && (
                            <div className="variant-detail-item mono" style={{ fontSize: '0.8em' }}>
                              <strong>g.</strong> {r.hgvsG}
                            </div>
                          )}
                        </div>
                      </td>
                    </tr>
                  )}
                </React.Fragment>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
};

export default VariantLookup;
