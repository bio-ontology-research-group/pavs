import React from 'react';
import { useTranslation } from 'react-i18next';

interface Props {
  source: string;
}

const SOURCE_CONFIG: Record<string, { label: string; className: string }> = {
  'ahmed-variants':       { label: 'Saudi cohort', className: 'badge-saudi' },
  'ahmed-pmid28454995':   { label: 'Saudi cohort', className: 'badge-saudi' },
  'fawzan-variants':      { label: 'Saudi cohort', className: 'badge-saudi' },
  'marwa-variants':       { label: 'Saudi cohort', className: 'badge-saudi' },
  'PMC6562004':           { label: 'Saudi cohort', className: 'badge-saudi' },
  'PMC7082194':           { label: 'Saudi cohort', className: 'badge-saudi' },
  'ddd-diagnoses':        { label: 'DDD', className: 'badge-ddd' },
  'Literature':           { label: 'Literature', className: 'badge-literature' },
  'ClinVar':              { label: 'ClinVar', className: 'badge-clinvar' },
};

const SourceBadge: React.FC<Props> = ({ source }) => {
  const { t } = useTranslation();
  const cfg = SOURCE_CONFIG[source] || { label: source, className: 'badge-unknown' };
  return (
    <span className={`source-badge ${cfg.className}`} title={source}>
      {cfg.label}
    </span>
  );
};

export default SourceBadge;
