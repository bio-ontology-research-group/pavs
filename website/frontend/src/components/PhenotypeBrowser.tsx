import React, { useState, useEffect, useCallback } from 'react';
import { useTranslation } from 'react-i18next';
import axios from 'axios';
import SourceBadge from './SourceBadge';

const API = import.meta.env.VITE_API_URL || 'http://localhost:8000';

// HP:0000118 = Phenotypic abnormality (root for browsable tree)
const ROOT_ID = 'HP:0000118';

interface HpoNode {
  id: string;
  label: string;
  definition?: string;
  arabic_label?: string;
  case_count: number;
  saudi_case_count: number;
  child_count: number;
}

interface HpoEnrichment {
  id: string;
  definition: string;
  layperson_synonyms: string[];
  synonyms: string[];
  arabic_label: string;
  arabic_layperson: string;
  arabic_definition: string;
}

interface CaseResult {
  id: string;
  gene?: string;
  disease?: string;
  source?: string;
  is_saudi?: boolean;
  hpo_ids?: string[];
}

// ─── Tree node component ──────────────────────────────────────────────────────

interface TreeNodeProps {
  node: HpoNode;
  depth: number;
  selected: HpoNode | null;
  onSelect: (node: HpoNode) => void;
  childrenMap: Record<string, HpoNode[]>;
  expanded: Set<string>;
  loadingNodes: Set<string>;
  onToggle: (id: string) => void;
  isAr: boolean;
}

const TreeNode: React.FC<TreeNodeProps> = ({
  node, depth, selected, onSelect,
  childrenMap, expanded, loadingNodes, onToggle, isAr,
}) => {
  const hasChildren = node.child_count > 0;
  const isExpanded = expanded.has(node.id);
  const isLoading = loadingNodes.has(node.id);
  const isSelected = selected?.id === node.id;
  const children = childrenMap[node.id] || [];
  const displayLabel = isAr && node.arabic_label ? node.arabic_label : node.label;

  return (
    <div className="hpo-tree-node">
      <div
        className={`hpo-tree-row ${isSelected ? 'hpo-tree-selected' : ''}`}
        style={{ paddingLeft: `${depth * 16 + 8}px` }}
      >
        <button
          className="hpo-tree-toggle"
          onClick={e => { e.stopPropagation(); if (hasChildren) onToggle(node.id); }}
          aria-label={isExpanded ? 'Collapse' : 'Expand'}
          style={{ visibility: hasChildren ? 'visible' : 'hidden' }}
        >
          {isLoading ? '⋯' : isExpanded ? '▼' : '▶'}
        </button>
        <span
          className="hpo-tree-label"
          onClick={() => onSelect(node)}
          title={node.definition}
        >
          {displayLabel}
          {node.saudi_case_count > 0 && (
            <span className="hpo-tree-count"> ({node.saudi_case_count})</span>
          )}
        </span>
        <span className="hpo-tree-id">{node.id}</span>
      </div>
      {isExpanded && children.map(child => (
        <TreeNode
          key={child.id}
          node={child}
          depth={depth + 1}
          selected={selected}
          onSelect={onSelect}
          childrenMap={childrenMap}
          expanded={expanded}
          loadingNodes={loadingNodes}
          onToggle={onToggle}
          isAr={isAr}
        />
      ))}
    </div>
  );
};

// ─── Case table ──────────────────────────────────────────────────────────────

const CaseTable: React.FC<{ cases: CaseResult[] }> = ({ cases }) => (
  <div className="gene-cases-table-wrap">
    <table className="gene-cases-table">
      <thead>
        <tr>
          <th>Case ID</th>
          <th>Gene</th>
          <th>Disease</th>
          <th>Source</th>
          <th>HPO</th>
        </tr>
      </thead>
      <tbody>
        {cases.map(c => (
          <tr key={c.id}>
            <td><a href={`/case/${encodeURIComponent(c.id)}`}>{c.id}</a></td>
            <td>
              {c.gene
                ? <a href={`/?gene=${encodeURIComponent(c.gene)}`} className="gene-link">{c.gene}</a>
                : '—'}
            </td>
            <td>{c.disease || '—'}</td>
            <td><SourceBadge source={c.source || ''} /></td>
            <td className="pb-hpo-count">{c.hpo_ids?.length ?? 0}</td>
          </tr>
        ))}
      </tbody>
    </table>
  </div>
);

// ─── Main component ──────────────────────────────────────────────────────────

const PhenotypeBrowser: React.FC = () => {
  const { t, i18n } = useTranslation();
  const isAr = i18n.language === 'ar';

  // Tree state
  const [childrenMap, setChildrenMap] = useState<Record<string, HpoNode[]>>({});
  const [expanded, setExpanded] = useState<Set<string>>(new Set());
  const [loadingNodes, setLoadingNodes] = useState<Set<string>>(new Set());

  // Selected term + cases
  const [selected, setSelected] = useState<HpoNode | null>(null);
  const [termEnrichment, setTermEnrichment] = useState<HpoEnrichment | null>(null);
  const [termCases, setTermCases] = useState<CaseResult[]>([]);
  const [casesLoading, setCasesLoading] = useState(false);

  // Toggles for case sources
  const [includeChildren, setIncludeChildren] = useState(true);
  const [includeDDD, setIncludeDDD] = useState(false);
  const [includeLiterature, setIncludeLiterature] = useState(false);

  const fetchChildren = useCallback(async (hp_id: string) => {
    if (childrenMap[hp_id] !== undefined) return; // already fetched
    setLoadingNodes(prev => new Set(prev).add(hp_id));
    try {
      const res = await axios.get(`${API}/api/hpo-children/${encodeURIComponent(hp_id)}`);
      setChildrenMap(prev => ({ ...prev, [hp_id]: res.data }));
    } catch {
      setChildrenMap(prev => ({ ...prev, [hp_id]: [] }));
    } finally {
      setLoadingNodes(prev => { const next = new Set(prev); next.delete(hp_id); return next; });
    }
  }, [childrenMap]);

  // Load root children on mount
  useEffect(() => {
    fetchChildren(ROOT_ID).then(() => {
      setExpanded(new Set([ROOT_ID]));
    });
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  const handleToggle = useCallback(async (id: string) => {
    if (!expanded.has(id)) {
      // Expand: fetch children first if needed
      await fetchChildren(id);
      setExpanded(prev => new Set(prev).add(id));
    } else {
      setExpanded(prev => { const next = new Set(prev); next.delete(id); return next; });
    }
  }, [expanded, fetchChildren]);

  const fetchCases = useCallback(async (
    node: HpoNode,
    children: boolean,
    ddd: boolean,
    lit: boolean,
  ) => {
    setCasesLoading(true);
    setTermCases([]);
    try {
      const res = await axios.get(`${API}/api/phenotype/${encodeURIComponent(node.id)}/cases`, {
        params: { include_children: children, include_ddd: ddd, include_literature: lit, limit: 1000 },
      });
      setTermCases(res.data);
    } catch {
      setTermCases([]);
    } finally {
      setCasesLoading(false);
    }
  }, []);

  const handleSelect = useCallback((node: HpoNode) => {
    setSelected(node);
    setTermEnrichment(null);
    fetchCases(node, includeChildren, includeDDD, includeLiterature);
    axios.get(`${API}/api/hpo/${encodeURIComponent(node.id)}`)
      .then(res => setTermEnrichment(res.data))
      .catch(() => {});
  }, [includeChildren, includeDDD, includeLiterature, fetchCases]);

  const handleToggleOption = (
    setter: React.Dispatch<React.SetStateAction<boolean>>,
    currentValue: boolean,
    key: 'children' | 'ddd' | 'lit',
  ) => {
    const newVal = !currentValue;
    setter(newVal);
    if (selected) {
      const c = key === 'children' ? newVal : includeChildren;
      const d = key === 'ddd' ? newVal : includeDDD;
      const l = key === 'lit' ? newVal : includeLiterature;
      fetchCases(selected, c, d, l);
    }
  };

  // Root nodes: children of HP:0000118
  const rootNodes = childrenMap[ROOT_ID] || [];
  const rootLoading = loadingNodes.has(ROOT_ID);

  const saudiCases = termCases.filter(r => r.is_saudi);
  const dddCases = termCases.filter(r => !r.is_saudi && r.source?.toLowerCase().includes('ddd'));
  const litCases = termCases.filter(r => !r.is_saudi && !r.source?.toLowerCase().includes('ddd'));

  return (
    <div className="hpo-browser-layout">
      {/* Left: HPO tree */}
      <div className="hpo-tree-panel">
        <div className="hpo-tree-header">
          <h3>{t('phenotypeBrowser.treeTitle')}</h3>
          <p className="hpo-tree-hint">{t('phenotypeBrowser.treeHint')}</p>
        </div>
        <div className="hpo-tree-scroll">
          {rootLoading && <div className="loading" style={{ padding: '12px' }}>{t('search.loading')}</div>}
          {rootNodes.map(node => (
            <TreeNode
              key={node.id}
              node={node}
              depth={0}
              selected={selected}
              onSelect={handleSelect}
              childrenMap={childrenMap}
              expanded={expanded}
              loadingNodes={loadingNodes}
              onToggle={handleToggle}
              isAr={isAr}
            />
          ))}
        </div>
      </div>

      {/* Right: selected term detail + cases */}
      <div className="hpo-detail-panel">
        {!selected ? (
          <p className="gene-select-prompt">{t('phenotypeBrowser.selectPrompt')}</p>
        ) : (
          <>
            <div className="hpo-term-card">
              <h3>
                {isAr && (termEnrichment?.arabic_label || selected.arabic_label)
                  ? (termEnrichment?.arabic_label || selected.arabic_label)
                  : selected.label}
                <span className="hpo-term-id"> {selected.id}</span>
              </h3>
              {isAr && termEnrichment?.arabic_layperson && (
                <p className="hpo-term-lay">
                  <span className="hpo-detail-key">{t('hpo.layTerm')}:</span>{' '}
                  {termEnrichment.arabic_layperson}
                </p>
              )}
              {(() => {
                const def = isAr
                  ? termEnrichment?.arabic_definition
                  : (termEnrichment?.definition || selected.definition);
                return def ? <p className="hpo-term-def">{def}</p> : null;
              })()}
              <div className="hpo-term-stats">
                <span>{t('phenotypeBrowser.casesCount')}: <strong>{termCases.length}</strong></span>
                {selected.child_count > 0 && (
                  <span> · {selected.child_count} {t('phenotypeBrowser.subcategories')}</span>
                )}
              </div>
              <a
                href={`https://hpo.jax.org/browse/term/${selected.id}`}
                target="_blank" rel="noopener noreferrer"
                className="ext-link"
              >
                HPO Browser ↗
              </a>
            </div>

            {/* Toggles */}
            <div className="gene-source-toggles" style={{ marginTop: '12px' }}>
              <label className="toggle-label">
                <input type="checkbox" checked={includeChildren}
                  onChange={() => handleToggleOption(setIncludeChildren, includeChildren, 'children')} />
                {' '}{t('phenotypeBrowser.includeChildren')}
              </label>
              <label className="toggle-label">
                <input type="checkbox" checked={includeDDD}
                  onChange={() => handleToggleOption(setIncludeDDD, includeDDD, 'ddd')} />
                {' '}DDD
              </label>
              <label className="toggle-label">
                <input type="checkbox" checked={includeLiterature}
                  onChange={() => handleToggleOption(setIncludeLiterature, includeLiterature, 'lit')} />
                {' '}{t('badges.literature')}
              </label>
            </div>

            {casesLoading && <div className="loading">{t('search.loading')}</div>}

            {!casesLoading && (
              <div className="pb-sections">
                {saudiCases.length > 0 && (
                  <section className="pb-section">
                    <h4>Saudi cohort ({saudiCases.length})</h4>
                    <CaseTable cases={saudiCases} />
                  </section>
                )}
                {includeDDD && dddCases.length > 0 && (
                  <section className="pb-section">
                    <h4>DDD ({dddCases.length})</h4>
                    <CaseTable cases={dddCases} />
                  </section>
                )}
                {includeLiterature && litCases.length > 0 && (
                  <section className="pb-section">
                    <h4>Literature ({litCases.length})</h4>
                    <CaseTable cases={litCases} />
                  </section>
                )}
                {termCases.length === 0 && (
                  <p className="gene-select-prompt">{t('phenotypeBrowser.noCases')}</p>
                )}
              </div>
            )}
          </>
        )}
      </div>
    </div>
  );
};

export default PhenotypeBrowser;
