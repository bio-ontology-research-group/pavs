import React, { useEffect, useState } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL || 'http://localhost:8000';

interface Variant {
  gene?: string;
  hgvsC?: string;
  hgvsG?: string;
  hgvsP?: string;
  rsId?: string;
  zygosity?: string;
  acmg?: string;
  consequence?: string;
  sift?: string;
  polyphen?: string;
  gnomadAF?: string;
  saudiAF?: string;
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
  suggested_diseases?: SuggestedDisease[];
  phenopacket_download_url?: string;
}

interface Props {
  caseId: string;
}

const HP_LINK = (id: string) => `https://hpo.jax.org/browse/term/${id}`;

const CaseDetail: React.FC<Props> = ({ caseId }) => {
  const { t } = useTranslation();
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
  const phenotypes = data.phenotypes || [];
  const excludedPhenotypes = data.excluded_phenotypes || [];

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
          {props['sex'] && <><dt>{t('case.sex')}</dt><dd>{props['sex'].includes('C16576') ? 'Female' : props['sex'].includes('C20197') ? 'Male' : props['sex']}</dd></>}
          {props['consanguinity'] && <><dt>{t('case.consanguinity')}</dt><dd>{props['consanguinity']}</dd></>}
          {props['progressStatus'] && <><dt>{t('case.progressStatus')}</dt><dd>{props['progressStatus']}</dd></>}
          {props['diseaseLabel'] && <><dt>{t('case.disease')}</dt><dd>{props['diseaseLabel']}</dd></>}
          {!props['diseaseLabel'] && data.suggested_diseases && data.suggested_diseases.length > 0 && (
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

      {phenotypes.length > 0 && (
        <section className="detail-section">
          <h3>{t('case.phenotypes')} ({phenotypes.length})</h3>
          <div className="hpo-tags">
            {phenotypes.map(hp => (
              <a key={hp.id} href={HP_LINK(hp.id)} target="_blank" rel="noopener noreferrer" className="hpo-tag">
                {hp.label ? `${hp.label} (${hp.id})` : hp.id}
              </a>
            ))}
          </div>
        </section>
      )}

      {excludedPhenotypes.length > 0 && (
        <section className="detail-section">
          <h3>{t('case.excludedPhenotypes')} ({excludedPhenotypes.length})</h3>
          <div className="hpo-tags excluded">
            {excludedPhenotypes.map(hp => (
              <a key={hp.id} href={HP_LINK(hp.id)} target="_blank" rel="noopener noreferrer" className="hpo-tag excluded">
                {hp.label ? `${hp.label} (${hp.id})` : hp.id}
              </a>
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
                  <a href={`/gene/${v.gene}`} className="gene-link bold">{v.gene}</a>
                )}
                {v.acmg && (
                  <span className={`acmg-badge acmg-${(v.acmg).toLowerCase().replace(/\s+/g, '-')}`}>
                    {v.acmg}
                  </span>
                )}
              </div>
              <dl className="prop-list">
                {v.hgvsC && <><dt>{t('results.hgvsC')}</dt><dd className="mono">{v.hgvsC}</dd></>}
                {v.hgvsP && <><dt>HGVS p.</dt><dd className="mono">{v.hgvsP}</dd></>}
                {v.hgvsG && <><dt>{t('results.hgvsG')}</dt><dd className="mono">{v.hgvsG}</dd></>}
                {v.zygosity && <><dt>{t('results.zygosity')}</dt><dd>{v.zygosity.split('/').pop()}</dd></>}
                {v.consequence && <><dt>VEP consequence</dt><dd>{v.consequence}</dd></>}
                {v.sift && <><dt>SIFT</dt><dd>{v.sift}</dd></>}
                {v.polyphen && <><dt>PolyPhen-2</dt><dd>{v.polyphen}</dd></>}
                {v.gnomadAF !== undefined && (
                  <><dt>{t('results.gnomadAF')}</dt><dd>{v.gnomadAF || '0'}</dd></>
                )}
                {v.saudiAF !== undefined && (
                  <><dt>{t('results.saudiAF')}</dt>
                  <dd>
                    {v.saudiAF === '0' || v.saudiAF === '0.0'
                      ? <strong className="novel-badge">{t('results.novel')}</strong>
                      : v.saudiAF}
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
