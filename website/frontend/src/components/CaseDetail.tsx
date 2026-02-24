import React, { useEffect, useState } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';
import HpoTag from './HpoTag';

const API = import.meta.env.VITE_API_URL || 'http://localhost:8000';

interface Variant {
  gene?: string;
  hgvsC?: string;
  hgvsG?: string;
  hgvsP?: string;
  rsId?: string;
  zygosity?: string;
  zygosity_label?: string;
  zygosity_ar?: string;
  acmg?: string;
  consequence?: string;
  consequence_label?: string;
  consequence_ar?: string;
  vepImpact?: string;
  sift?: string;
  sift_label?: string;
  sift_ar?: string;
  sift_score?: string;
  polyphen?: string;
  polyphen_label?: string;
  polyphen_ar?: string;
  polyphen_score?: string;
  gnomadAF?: string;
  gnomadAF_MID?: string;
  gnomadAF_SAS?: string;
  caddPhred?: string;
  revelScore?: string;
  alphaMissense?: string;
  spliceAI?: string;
  saudiAF?: string;
  saudiAC?: string;
  clinvarId?: string;
  clinvarSig?: string;
  togovar_url?: string;
}

interface HpoTerm {
  id: string;
  label: string;
}

interface SuggestedDisease {
  id: string;
  label: string;
  url: string;
  gene: string;
}

interface CaseData {
  id: string;
  properties: Record<string, string>;
  variants: Variant[];
  phenotypes?: HpoTerm[];
  excluded_phenotypes?: HpoTerm[];
  diseases?: { id: string; label?: string; url?: string }[];
  suggested_diseases?: SuggestedDisease[];
  phenopacket_download_url?: string;
  field_metadata?: Record<string, { en: string; ar: string }>;
}

interface Props {
  caseId: string;
}

/** Simple info icon with a title tooltip. */
const InfoIcon: React.FC<{ text?: string }> = ({ text }) => {
  if (!text) return null;
  return (
    <span className="info-icon" title={text} style={{ marginLeft: '4px', cursor: 'help', opacity: 0.6 }}>
      ⓘ
    </span>
  );
};

/** Derive a publication link from the source file name or the case's pmid property. */
function getPublication(source: string, pmid?: string): { label: string; url: string } | null {
  if (pmid) {
    return { label: `PMID:${pmid}`, url: `https://pubmed.ncbi.nlm.nih.gov/${pmid}` };
  }
  // e.g. "ahmed-pmid28454995" → PMID 28454995
  const pmidMatch = source.match(/pmid(\d+)/i);
  if (pmidMatch) {
    return { label: `PMID:${pmidMatch[1]}`, url: `https://pubmed.ncbi.nlm.nih.gov/${pmidMatch[1]}` };
  }
  // e.g. "PMC6562004" or "PMC7082194"
  const pmcMatch = source.match(/^(PMC\d+)$/i);
  if (pmcMatch) {
    const pmcId = pmcMatch[1].toUpperCase();
    return { label: pmcId, url: `https://www.ncbi.nlm.nih.gov/pmc/articles/${pmcId}/` };
  }
  if (source === 'ddd-diagnoses') {
    return { label: 'DECIPHER/DDD', url: 'https://www.deciphergenomics.org/' };
  }
  return null;
}

/** Localize ACMG class labels and normalize for CSS classes. */
function localizeAcmg(acmg: string, t: any): { label: string; className: string } {
  const raw = acmg.toLowerCase().replace(/[^a-z]+/g, '_');
  // Handle common variants like UNCERTAIN_SIGNIFICANCE vs uncertain_significance
  const key = raw === 'uncertain_significance' ? 'uncertain_significance' :
              raw === 'pathogenic' ? 'pathogenic' :
              raw === 'likely_pathogenic' ? 'likely_pathogenic' :
              raw === 'likely_benign' ? 'likely_benign' :
              raw === 'benign' ? 'benign' : raw;

  const label = t(`acmg.${key}`, { defaultValue: acmg });
  return { label, className: `acmg-${key.replace(/_/g, '-')}` };
}

const CaseDetail: React.FC<Props> = ({ caseId }) => {
  const { t, i18n } = useTranslation();
  const isAr = i18n.language === 'ar';
  const [data, setData] = useState<CaseData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');

  useEffect(() => {
    setLoading(true);
    axios.get(`${API}/api/case/${encodeURIComponent(caseId)}`)
      .then(res => { setData(res.data); setLoading(false); })
      .catch(e => { setError(e.message); setLoading(false); });
  }, [caseId]);

  if (loading) return <div className="loading">{t('search.loading')}</div>;
  if (error) return <div className="error-msg">{error}</div>;
  if (!data) return null;

  const props = data.properties;
  const source = props['source'] || 'unknown';
  const pubLink = getPublication(source, props['pmid']);
  const phenotypes = data.phenotypes || [];
  const excludedPhenotypes = data.excluded_phenotypes || [];
  const diseases = data.diseases || [];

  return (
    <div className="case-detail">
      <div className="case-header">
        <h2>{caseId}</h2>
        <SourceBadge source={source} />
        {(props['isSaudi'] === 'true' || props['isSaudi'] === '1') && (
          <span className="ancestry-badge" title="Middle Eastern / Saudi Arabia">
            🇸🇦 {t('case.middleEastern')}
          </span>
        )}
        <div className="case-actions">
          <a href={`${API}/api/phenopacket/${encodeURIComponent(caseId)}/download`}
             className="btn-primary" download>
            ↓ {t('results.download')}
          </a>
        </div>
      </div>

      <section className="detail-section">
        <h3>{t('case.demographics')}</h3>
        <dl className="prop-list">
          {props['sex'] && (
            <>
              <dt>{t('case.sex')}</dt>
              <dd>
                {props['sex'].includes('C16576') ? t('case.female') : 
                 props['sex'].includes('C20197') ? t('case.male') : 
                 props['sex'] === 'Male' ? t('case.male') :
                 props['sex'] === 'Female' ? t('case.female') :
                 props['sex']}
              </dd>
            </>
          )}
          {props['age'] && (
            <>
              <dt>{t('case.age')}</dt>
              <dd>
                {props['age']
                  .replace(/^P/, '')
                  .replace(/Y/, ` ${t('case.years')}`)
                  .replace(/M/, ` ${t('case.months')}`)
                  .replace(/D/, ` ${t('case.days')}`)}
              </dd>
            </>
          )}
          {props['consanguinity'] && (
            <>
              <dt>{t('case.consanguinity')}</dt>
              <dd>
                {isAr ? (props['consanguinity_ar'] || props['consanguinity_label'] || props['consanguinity']) : (props['consanguinity_label'] || props['consanguinity'])}
              </dd>
            </>
          )}
          {props['familyHistory'] && (
            <>
              <dt>{t('case.familyHistory')}</dt>
              <dd>{props['familyHistory']}</dd>
            </>
          )}
          {props['progressStatus'] && (
            <>
              <dt>{t('case.progressStatus')}</dt>
              <dd>
                {props['progressStatus'] === 'Solved' || props['progressStatus'] === 'SOLVED' ? t('case.solved') :
                 props['progressStatus'] === 'Unsolved' || props['progressStatus'] === 'UNSOLVED' ? t('case.unsolved') :
                 props['progressStatus'] === 'Unknown' || props['progressStatus'] === 'UNKNOWN_PROGRESS' || props['progressStatus'] === 'IN_PROGRESS' ? t('case.unknown') :
                 props['progressStatus']}
              </dd>
            </>
          )}
          {props['diseaseLabel'] && diseases.length === 0 && <><dt>{t('case.disease')}</dt><dd>{props['diseaseLabel']}</dd></>}
          {diseases.length === 0 && data.suggested_diseases && data.suggested_diseases.length > 0 && (
            <>
              <dt>{t('case.disease')}</dt>
              <dd>
                <div className="suggested-disease-block">
                  <span className="suggested-label">{t('case.suggestedDiseaseNote')}</span>
                  <ul className="suggested-disease-list">
                    {data.suggested_diseases.map(sd => (
                      <li key={sd.id}>
                        <a href={sd.url} target="_blank" rel="noopener noreferrer">
                          {sd.label}
                        </a>
                        <span className="suggested-id"> ({sd.id})</span>
                        <span className="suggested-gene"> — {t('case.basedOnGene')} {sd.gene}</span>
                      </li>
                    ))}
                  </ul>
                </div>
              </dd>
            </>
          )}
          {props['pmid'] && (
            <><dt>{t('case.pmid')}</dt>
            <dd><a href={`https://pubmed.ncbi.nlm.nih.gov/${props['pmid']}`} target="_blank" rel="noopener noreferrer">
              {props['pmid']}
            </a></dd></>
          )}
        </dl>
      </section>

      {diseases.length > 0 && (
        <section className="detail-section">
          <h3>{t('case.diseases')}</h3>
          <ul className="disease-list">
            {diseases.map((d, idx) => (
              <li key={idx} className="disease-item">
                <a href={d.url} target="_blank" rel="noopener noreferrer" className="disease-link">
                  <strong>{d.id}</strong> {d.label ? `— ${d.label}` : ''}
                </a>
              </li>
            ))}
          </ul>
        </section>
      )}

      {phenotypes.length > 0 && (
        <section className="detail-section">
          <h3>{t('case.phenotypes')} ({phenotypes.length})</h3>
          <div className="hpo-tags">
            {phenotypes.map(hp => (
              <HpoTag key={hp.id} id={hp.id} label={hp.label} />
            ))}
          </div>
        </section>
      )}

      {excludedPhenotypes.length > 0 && (
        <section className="detail-section">
          <h3>{t('case.excludedPhenotypes')} ({excludedPhenotypes.length})</h3>
          <div className="hpo-tags">
            {excludedPhenotypes.map(hp => (
              <HpoTag key={hp.id} id={hp.id} label={hp.label} excluded />
            ))}
          </div>
        </section>
      )}

      {data.variants.length > 0 && (
        <section className="detail-section">
          <h3>{t('case.variants')}</h3>
          {data.variants.map((v, idx) => (
            <div key={idx} className="variant-card">
              <div className="variant-header">
                {v.gene && (
                  <div style={{ display: 'flex', alignItems: 'baseline', gap: '8px' }}>
                    <a href={`/?gene=${encodeURIComponent(v.gene)}`} className="gene-link bold">{v.gene}</a>
                    <a href={`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${v.gene}`}
                       target="_blank" rel="noopener noreferrer" className="gene-ext-link" title="Open in GeneCards">
                      ↗
                    </a>
                  </div>
                )}
                {v.acmg && (() => {
                  const { label, className } = localizeAcmg(v.acmg, t);
                  return (
                    <span className={`acmg-badge ${className}`}>
                      {label}
                    </span>
                  );
                })()}
                <div style={{ marginLeft: 'auto', display: 'flex', gap: '6px', alignItems: 'center' }}>
                  <SourceBadge source={source} />
                  {pubLink && (
                    <a href={pubLink.url} target="_blank" rel="noopener noreferrer" className="pub-link">
                      {pubLink.label}
                    </a>
                  )}
                </div>
              </div>
              <dl className="prop-list">
                {v.hgvsC && <><dt>{t('results.hgvsC')}</dt><dd className="mono">{v.hgvsC}</dd></>}
                {v.hgvsP && <><dt>HGVS p.</dt><dd className="mono">{v.hgvsP}</dd></>}
                {v.hgvsG && <><dt>{t('results.hgvsG')}</dt><dd className="mono">{v.hgvsG}</dd></>}
                {(v.zygosity_label || v.zygosity) && (
                  <>
                    <dt>
                      {t('results.zygosity')}
                      <InfoIcon text={data.field_metadata?.zygosity?.[isAr ? 'ar' : 'en']} />
                    </dt>
                    <dd>
                      {isAr ? (v.zygosity_ar || v.zygosity_label || v.zygosity?.split('/').pop()) : (v.zygosity_label || v.zygosity?.split('/').pop())}
                    </dd>
                  </>
                )}
                {(v.consequence_label || v.consequence) && (
                  <>
                    <dt>
                      {t('results.consequence')}
                      <InfoIcon text={data.field_metadata?.consequence?.[isAr ? 'ar' : 'en']} />
                    </dt>
                    <dd>
                      {isAr ? (v.consequence_ar || v.consequence_label || v.consequence) : (v.consequence_label || v.consequence)}
                      {isAr && v.consequence_label && <span className="translation-sub"> ({v.consequence_label})</span>}
                    </dd>
                  </>
                )}
                {v.vepImpact && (
                  <>
                    <dt>{t('results.vepImpact')}</dt>
                    <dd>{v.vepImpact}</dd>
                  </>
                )}
                {(v.sift_label || v.sift) && (
                  <>
                    <dt>
                      {t('results.sift')}
                      <InfoIcon text={data.field_metadata?.sift?.[isAr ? 'ar' : 'en']} />
                    </dt>
                    <dd>
                      {isAr ? (v.sift_ar || v.sift_label || v.sift) : (v.sift_label || v.sift)}
                      {v.sift_score && <span className="score-value"> ({v.sift_score})</span>}
                    </dd>
                  </>
                )}
                {(v.polyphen_label || v.polyphen) && (
                  <>
                    <dt>
                      {t('results.polyphen')}
                      <InfoIcon text={data.field_metadata?.polyphen?.[isAr ? 'ar' : 'en']} />
                    </dt>
                    <dd>
                      {isAr ? (v.polyphen_ar || v.polyphen_label || v.polyphen) : (v.polyphen_label || v.polyphen)}
                      {v.polyphen_score && <span className="score-value"> ({v.polyphen_score})</span>}
                    </dd>
                  </>
                )}
                {v.gnomadAF !== undefined && (
                  <><dt>{t('results.gnomadAF')}</dt><dd>{v.gnomadAF || '0'}</dd></>
                )}
                {v.gnomadAF_MID && <><dt>gnomAD MID</dt><dd>{v.gnomadAF_MID}</dd></>}
                {v.gnomadAF_SAS && <><dt>gnomAD SAS</dt><dd>{v.gnomadAF_SAS}</dd></>}
                {v.caddPhred && <><dt>CADD Phred</dt><dd>{v.caddPhred}</dd></>}
                {v.revelScore && <><dt>REVEL</dt><dd>{v.revelScore}</dd></>}
                {v.alphaMissense && <><dt>AlphaMissense</dt><dd>{v.alphaMissense}</dd></>}
                {v.spliceAI && <><dt>SpliceAI Max</dt><dd>{v.spliceAI}</dd></>}
                {v.saudiAF !== undefined && (
                  <><dt>{t('results.saudiAF')}</dt>
                  <dd>
                    {v.saudiAF === '0' || v.saudiAF === '0.0'
                      ? <strong className="novel-badge">{t('results.novel')}</strong>
                      : v.saudiAF}
                    {v.saudiAC && <span className="ac-count"> (AC: {v.saudiAC})</span>}
                  </dd></>
                )}
              </dl>
              <div className="variant-links">
                {v.togovar_url && (
                  <a href={v.togovar_url} target="_blank" rel="noopener noreferrer" className="ext-link togovar-link">
                    {t('links.togovar')}
                  </a>
                )}
                {v.rsId && (
                  <a href={`https://www.ncbi.nlm.nih.gov/snp/${v.rsId}`}
                     target="_blank" rel="noopener noreferrer" className="ext-link">
                    {t('links.dbsnp')}: {v.rsId}
                  </a>
                )}
                {v.clinvarId && (
                  <a href={`https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.clinvarId}/`}
                     target="_blank" rel="noopener noreferrer" className="ext-link">
                    {t('links.clinvar')}: {v.clinvarSig || v.clinvarId}
                  </a>
                )}
                {v.gene && (
                  <a href={`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${v.gene}`}
                     target="_blank" rel="noopener noreferrer" className="ext-link">
                    GeneCards
                  </a>
                )}
              </div>
            </div>
          ))}
        </section>
      )}
    </div>
  );
};

export default CaseDetail;
