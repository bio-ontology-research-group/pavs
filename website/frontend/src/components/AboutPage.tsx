import React, { useEffect, useState } from 'react';
import ReactMarkdown from 'react-markdown';

const AboutPage: React.FC = () => {
  const [content, setContent] = useState('');
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    fetch('/about.md')
      .then(r => r.text())
      .then(text => { setContent(text); setLoading(false); })
      .catch(() => { setContent('# About\n\nContent not available.'); setLoading(false); });
  }, []);

  if (loading) return <div className="loading">Loading…</div>;

  return (
    <div className="about-page">
      <ReactMarkdown>{content}</ReactMarkdown>
    </div>
  );
};

export default AboutPage;
