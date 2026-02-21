import React, { useState, useEffect } from 'react';
import { useTranslation } from 'react-i18next';
import './i18n';
import './App.css';
import PhenotypeSearch from './components/PhenotypeSearch';
import VariantLookup from './components/VariantLookup';
import DiseaseBrowser from './components/DiseaseBrowser';
import GeneBrowser from './components/GeneBrowser';
import AboutPage from './components/AboutPage';
import CaseDetail from './components/CaseDetail';

type TabId = 'phenotype' | 'variant' | 'disease' | 'gene' | 'about';

const API = import.meta.env.VITE_API_URL || 'http://localhost:8000';

const App: React.FC = () => {
  const { t, i18n } = useTranslation();
  const [activeTab, setActiveTab] = useState<TabId>('phenotype');
  const [caseId, setCaseId] = useState<string | null>(null);

  // Simple client-side routing for /case/:id and /gene/:symbol
  useEffect(() => {
    const path = window.location.pathname;
    const caseMatch = path.match(/^\/case\/(.+)/);
    if (caseMatch) {
      setCaseId(decodeURIComponent(caseMatch[1]));
      return;
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

  if (caseId) {
    return (
      <div className="app">
        <nav className="navbar">
          <a href="/" className="nav-brand">
            <img src="/logo.svg" alt="PAVS" height="36" />
          </a>
          <button className="lang-switcher" onClick={switchLang} title="Switch language">
            {i18n.language === 'ar' ? 'EN' : 'ع'}
          </button>
        </nav>
        <main className="main-content">
          <button className="back-btn" onClick={() => window.history.back()}>← Back</button>
          <CaseDetail caseId={caseId} />
        </main>
      </div>
    );
  }

  const tabs: { id: TabId; label: string }[] = [
    { id: 'phenotype', label: t('nav.phenotype') },
    { id: 'variant',   label: t('nav.variant') },
    { id: 'disease',   label: t('nav.disease') },
    { id: 'gene',      label: t('nav.gene') },
    { id: 'about',     label: t('nav.about') },
  ];

  return (
    <div className="app">
      <nav className="navbar">
        <div className="nav-brand">
          <img src="/logo.svg" alt="PAVS" height="36" />
          <span className="nav-title">PAVS</span>
        </div>
        <div className="nav-tabs">
          {tabs.map(tab => (
            <button
              key={tab.id}
              className={`nav-tab ${activeTab === tab.id ? 'active' : ''}`}
              onClick={() => setActiveTab(tab.id)}
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
        {activeTab === 'phenotype' && <PhenotypeSearch />}
        {activeTab === 'variant'   && <VariantLookup />}
        {activeTab === 'disease'   && <DiseaseBrowser />}
        {activeTab === 'gene'      && <GeneBrowser />}
        {activeTab === 'about'     && <AboutPage />}
      </main>

      <footer className="footer">
        <p>
          PAVS — Phenotypic and Variant Standardization ·{' '}
          <a href="https://www.kaust.edu.sa" target="_blank" rel="noopener noreferrer">KAUST</a>
        </p>
      </footer>
    </div>
  );
};

export default App;
