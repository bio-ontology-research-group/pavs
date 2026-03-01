import React, { useState, useEffect } from 'react';
import { useTranslation } from 'react-i18next';
import './i18n';
import './App.css';
import PhenotypeSearch from './components/PhenotypeSearch';
import VariantLookup from './components/VariantLookup';
import DiseaseBrowser from './components/DiseaseBrowser';
import GeneBrowser from './components/GeneBrowser';
import PhenotypeBrowser from './components/PhenotypeBrowser';
import SparqlExplorer from './components/SparqlExplorer';
import AboutPage from './components/AboutPage';
import CaseDetail from './components/CaseDetail';

type TabId = 'phenotype' | 'variant' | 'disease' | 'gene' | 'phenotype-browser' | 'sparql' | 'about';

const API = import.meta.env.VITE_API_URL ?? '';

const App: React.FC = () => {
  const { t, i18n } = useTranslation();
  const [activeTab, setActiveTab] = useState<TabId>('phenotype');
  const [caseId, setCaseId] = useState<string | null>(null);
  const [initialGene, setInitialGene] = useState<string>('');

  // Simple client-side routing for /case/:id; gene navigation via ?gene= param
  useEffect(() => {
    const path = window.location.pathname;
    const caseMatch = path.match(/^\/case\/(.+)/);
    if (caseMatch) {
      setCaseId(decodeURIComponent(caseMatch[1]));
      return;
    }
    // Check for ?gene= URL search param (used by gene links in all views)
    const params = new URLSearchParams(window.location.search);
    const gene = params.get('gene');
    if (gene) {
      setActiveTab('gene');
      setInitialGene(decodeURIComponent(gene));
    }
  }, []);

  const switchLang = () => {
    const next = i18n.language === 'ar' ? 'en' : 'ar';
    i18n.changeLanguage(next);
    document.documentElement.dir = next === 'ar' ? 'rtl' : 'ltr';
    document.documentElement.lang = next;
    localStorage.setItem('i18nextLng', next);
  };

  // Apply RTL on init
  useEffect(() => {
    const lang = i18n.language;
    document.documentElement.dir = lang === 'ar' ? 'rtl' : 'ltr';
    document.documentElement.lang = lang;
  }, [i18n.language]);

  // Tab click: always navigate away from case view and update URL
  const handleTabClick = (tabId: TabId) => {
    setActiveTab(tabId);
    if (caseId) {
      setCaseId(null);
      window.history.pushState({}, '', '/');
    }
  };

  const tabs: { id: TabId; label: string }[] = [
    { id: 'phenotype',         label: t('nav.phenotype') },
    { id: 'variant',           label: t('nav.variant') },
    { id: 'disease',           label: t('nav.disease') },
    { id: 'gene',              label: t('nav.gene') },
    { id: 'phenotype-browser', label: t('nav.phenotypeBrowser') },
    { id: 'sparql',            label: t('nav.sparql') },
    { id: 'about',             label: t('nav.about') },
  ];

  return (
    <div className="app">
      <nav className="navbar">
        <a href="/" className="nav-brand" onClick={(e) => {
          e.preventDefault();
          setActiveTab('phenotype');
          setCaseId(null);
          window.history.pushState({}, '', '/');
        }}>
          <img src="/logo.svg" alt="PAVS Logo" height="36" />
        </a>
        <div className="nav-tabs">
          {tabs.map(tab => (
            <button
              key={tab.id}
              className={`nav-tab ${!caseId && activeTab === tab.id ? 'active' : ''}`}
              onClick={() => handleTabClick(tab.id)}
            >
              {tab.label}
            </button>
          ))}
        </div>
        <div className="nav-right">
          <a href={`${API}/api/phenopackets/download-all`}
             className="download-all-btn" title="Download all Saudi phenopackets" download>
            ↓ {t('results.downloadAll')}
          </a>
          <button className="lang-switcher" onClick={switchLang} title="Switch language / تغيير اللغة">
            {i18n.language === 'ar' ? 'EN' : 'ع'}
          </button>
        </div>
      </nav>

      <main className="main-content">
        {caseId ? (
          <>
            <button className="back-btn" onClick={() => { setCaseId(null); window.history.back(); }}>← Back</button>
            <CaseDetail caseId={caseId} />
          </>
        ) : (
          <>
            {activeTab === 'phenotype'         && <PhenotypeSearch />}
            {activeTab === 'variant'           && <VariantLookup />}
            {activeTab === 'disease'           && <DiseaseBrowser />}
            {activeTab === 'gene'              && <GeneBrowser initialGene={initialGene} />}
            {activeTab === 'phenotype-browser' && <PhenotypeBrowser />}
            {activeTab === 'sparql'            && <SparqlExplorer />}
            {activeTab === 'about'             && <AboutPage />}
          </>
        )}
      </main>

      <footer className="footer">
        <p>
          PAVS — Phenotype-Associated Variants in Saudi Arabia ·{' '}
          <a href="https://www.kaust.edu.sa" target="_blank" rel="noopener noreferrer">KAUST</a>
        </p>
      </footer>
    </div>
  );
};

export default App;
