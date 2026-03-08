import React, { useEffect, useState } from 'react';
import ReactMarkdown from 'react-markdown';
import { useTranslation } from 'react-i18next';

const AboutPage: React.FC = () => {
  const { i18n } = useTranslation();
  const [content, setContent] = useState('');
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const lang = i18n.language === 'ar' ? 'ar' : 'en';
    fetch(`/api/about?lang=${lang}`)
      .then(r => r.text())
      .then(text => { setContent(text); setLoading(false); })
      .catch(() => { setContent('# About\n\nContent not available.'); setLoading(false); });
  }, [i18n.language]);

  if (loading) return <div className="loading">Loading…</div>;

  return (
    <div className="about-page">
      <ReactMarkdown>{content}</ReactMarkdown>
    </div>
  );
};

export default AboutPage;
