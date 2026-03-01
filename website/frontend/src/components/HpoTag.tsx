import React, { useState, useEffect, useRef } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';

const API = import.meta.env.VITE_API_URL ?? '';

const HP_LINK = (id: string) => `https://hpo.jax.org/browse/term/${id}`;

interface HpoEnrichment {
  id: string;
  definition: string;
  layperson_synonyms: string[];
  synonyms: string[];
  arabic_label: string;
  arabic_layperson: string;
  arabic_definition: string;
}

interface Props {
  id: string;
  label: string;
  excluded?: boolean;
}

const HpoTag: React.FC<Props> = ({ id, label, excluded = false }) => {
  const { t, i18n } = useTranslation();
  const [expanded, setExpanded] = useState(false);
  const [enrichment, setEnrichment] = useState<HpoEnrichment | null>(null);
  const [loading, setLoading] = useState(false);

  const handleToggle = async (e: React.MouseEvent) => {
    e.preventDefault();
    e.stopPropagation();

    if (!expanded && !enrichment && !loading) {
      setLoading(true);
      try {
        const res = await axios.get(`${API}/api/hpo/${encodeURIComponent(id)}`);
        setEnrichment(res.data);
      } catch {
        // Show expand without enrichment data (network error etc.)
        setEnrichment({ id, definition: '', layperson_synonyms: [], synonyms: [], arabic_label: '', arabic_layperson: '', arabic_definition: '' });
      } finally {
        setLoading(false);
      }
    }
    setExpanded(prev => !prev);
  };

  const isAr = i18n.language === 'ar';

  // Auto-fetch enrichment when language switches to Arabic so the chip label updates immediately
  const fetchedRef = useRef(false);
  useEffect(() => {
    if (!isAr || fetchedRef.current) return;
    fetchedRef.current = true;
    if (enrichment || loading) return;
    setLoading(true);
    axios.get(`${API}/api/hpo/${encodeURIComponent(id)}`)
      .then(res => setEnrichment(res.data))
      .catch(() => setEnrichment({ id, definition: '', layperson_synonyms: [], synonyms: [], arabic_label: '', arabic_layperson: '', arabic_definition: '' }))
      .finally(() => setLoading(false));
  }, [isAr, id]); // eslint-disable-line react-hooks/exhaustive-deps

  // In Arabic mode, show Arabic label if available
  const chipLabel = isAr && enrichment?.arabic_label ? enrichment.arabic_label : label;

  // Layperson label for the current language
  const layLabel = enrichment && (
    isAr ? enrichment.arabic_layperson : enrichment.layperson_synonyms?.[0]
  );

  // Definition for the current language
  const definition = enrichment && (
    isAr ? enrichment.arabic_definition : enrichment.definition
  );

  // Synonyms (always English for now; shown only when expanded and not Arabic mode)
  const synonyms = enrichment?.synonyms ?? [];

  return (
    <span className={`hpo-tag-wrapper ${expanded ? 'hpo-expanded' : ''}`}>
      <span className={`hpo-tag-chip ${excluded ? 'excluded' : ''}`}>
        <a href={HP_LINK(id)} target="_blank" rel="noopener noreferrer">
          {chipLabel ? `${chipLabel} (${id})` : id}
        </a>
        <button
          className="hpo-info-btn"
          onClick={handleToggle}
          title={expanded ? t('hpo.collapse') : t('hpo.expand')}
          aria-expanded={expanded}
        >
          {loading ? '…' : expanded ? '▲' : 'ⓘ'}
        </button>
      </span>

      {expanded && enrichment && (
        <span className="hpo-tag-details">
          {layLabel && (
            <span className="hpo-detail-row">
              <span className="hpo-detail-key">{t('hpo.layTerm')}:</span>
              <span className="hpo-detail-val hpo-lay">{layLabel}</span>
            </span>
          )}
          {definition && (
            <span className="hpo-detail-row">
              <span className="hpo-detail-key">{t('hpo.definition')}:</span>
              <span className="hpo-detail-val">{definition}</span>
            </span>
          )}
          {!isAr && synonyms.length > 0 && (
            <span className="hpo-detail-row">
              <span className="hpo-detail-key">{t('hpo.alsoKnownAs')}:</span>
              <span className="hpo-detail-val hpo-synonyms">{synonyms.slice(0, 4).join(' · ')}</span>
            </span>
          )}
          {isAr && enrichment.arabic_label && !isAr && (
            /* Arabic label shown in chip when in AR mode; show EN label here when in AR mode */
            <span className="hpo-detail-row">
              <span className="hpo-detail-key">EN:</span>
              <span className="hpo-detail-val">{label}</span>
            </span>
          )}
        </span>
      )}
    </span>
  );
};

export default HpoTag;
