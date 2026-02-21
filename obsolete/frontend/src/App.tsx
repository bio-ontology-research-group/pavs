import React, { useState, useEffect } from 'react';
import { Container, Nav, Navbar, Tab, Tabs, Form, Button, Table, Badge, Card, Spinner, Modal } from 'react-bootstrap';
import AsyncSelect from 'react-select/async';
import axios from 'axios';
import 'bootstrap/dist/css/bootstrap.min.css';

const API_BASE = '/api';

interface HPOTerm {
  id: string;
  name: string;
}

interface PhenopacketSummary {
  id: number;
  external_id: string;
  gene_symbol: string | null;
  disease_id: string | null;
  disease_name: string | null;
  disease_is_suggested: number;
  source: string;
  score?: number;
}

interface PhenopacketDetail extends PhenopacketSummary {
  raw_json: any;
}

function App() {
  const [selectedHPO, setSelectedHPO] = useState<any[]>([]);
  const [results, setResults] = useState<PhenopacketSummary[]>([]);
  const [loading, setLoading] = useState(false);
  const [includeDiseasePhenotypes, setIncludeDiseasePhenotypes] = useState(false);
  const [geneSearch, setGeneSearch] = useState('');
  const [genes, setGenes] = useState<string[]>([]);
  const [diseases, setDiseases] = useState<any[]>([]);
  
  const [geneFilter, setGeneFilter] = useState('');
  const [diseaseFilter, setDiseaseFilter] = useState('');
  
  const [activeTab, setActiveTab] = useState('phenotype');
  const [showDetail, setShowDetail] = useState(false);
  const [selectedPacket, setSelectedPacket] = useState<PhenopacketDetail | null>(null);
  const [detailLoading, setDetailLoading] = useState(false);

  useEffect(() => {
    fetchMetadata();
  }, []);

  const fetchMetadata = async () => {
    try {
      const [genesRes, diseasesRes] = await Promise.all([
        axios.get(`${API_BASE}/genes`),
        axios.get(`${API_BASE}/diseases`)
      ]);
      setGenes(genesRes.data);
      setDiseases(diseasesRes.data);
    } catch (err) {
      console.error("Error fetching metadata", err);
    }
  };

  const loadHPOOptions = async (inputValue: string) => {
    if (inputValue.length < 2) return [];
    try {
      const res = await axios.get(`${API_BASE}/search/hpo?q=${inputValue}`);
      return res.data.map((t: HPOTerm) => ({ value: t.id, label: `${t.name} (${t.id})` }));
    } catch (err) {
      return [];
    }
  };

  const handlePhenotypeSearch = async () => {
    setLoading(true);
    try {
      const hpo_ids = selectedHPO.map(h => h.value);
      const res = await axios.post(`${API_BASE}/search/phenotype`, {
        hpo_ids,
        similarity_method: 'resnik',
        limit: 50,
        include_disease_phenotypes: includeDiseasePhenotypes
      });
      setResults(res.data);
      setActiveTab('results');
    } catch (err) {
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const handleGeneSearch = async (gene?: string) => {
    const term = gene || geneSearch;
    setLoading(true);
    try {
      const res = await axios.get(`${API_BASE}/search/genotype?gene=${term}`);
      setResults(res.data);
      setActiveTab('results');
    } catch (err) {
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const handleDiseaseSearch = async (diseaseId: string) => {
    setLoading(true);
    try {
      const res = await axios.get(`${API_BASE}/search/disease?disease_id=${diseaseId}`);
      setResults(res.data);
      setActiveTab('results');
    } catch (err) {
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const handleRowClick = async (id: number) => {
    setDetailLoading(true);
    setShowDetail(true);
    try {
      const res = await axios.get(`${API_BASE}/phenopacket/${id}`);
      setSelectedPacket(res.data);
    } catch (err) {
      console.error("Error fetching phenopacket detail", err);
    } finally {
      setDetailLoading(false);
    }
  };

  return (
    <div className="App">
      <Navbar bg="dark" variant="dark" expand="lg">
        <Container>
          <Navbar.Brand href="#home">Saudi Phenopackets Search</Navbar.Brand>
        </Container>
      </Navbar>

      <Container className="mt-4">
        <Tabs activeKey={activeTab} onSelect={(k: any) => setActiveTab(k)} id="search-tabs" className="mb-3">
          <Tab eventKey="phenotype" title="Search by Phenotype">
            <Card className="p-3">
              <Form.Group className="mb-3">
                <Form.Label>Enter HPO Terms</Form.Label>
                <AsyncSelect
                  isMulti
                  cacheOptions
                  loadOptions={loadHPOOptions}
                  onChange={(val: any) => setSelectedHPO(val)}
                  placeholder="Type to search phenotypes (e.g. 'Microcephaly')..."
                />
              </Form.Group>
              <Form.Group className="mb-3">
                <Form.Check 
                  type="checkbox"
                  label="Include phenotypes from disease clinical profile (HPOA)"
                  checked={includeDiseasePhenotypes}
                  onChange={(e) => setIncludeDiseasePhenotypes(e.target.checked)}
                />
                <Form.Text className="text-muted">
                  When enabled, search will also match against typical phenotypes associated with the patient's diagnosis.
                </Form.Text>
              </Form.Group>
              <Button onClick={handlePhenotypeSearch} disabled={loading}>
                {loading ? <Spinner size="sm" /> : 'Rank by Similarity'}
              </Button>
            </Card>
          </Tab>
          
          <Tab eventKey="genotype" title="Search by Genotype">
            <Card className="p-3">
              <Form.Group className="mb-3">
                <Form.Label>Gene Symbol</Form.Label>
                <Form.Control 
                  type="text" 
                  placeholder="e.g. ABCA4" 
                  value={geneSearch}
                  onChange={(e) => setGeneSearch(e.target.value)}
                />
              </Form.Group>
              <Button onClick={() => handleGeneSearch()} disabled={loading}>
                {loading ? <Spinner size="sm" /> : 'Search Gene'}
              </Button>
            </Card>
          </Tab>

          <Tab eventKey="genes" title={`Genes (${genes.length})`}>
             <Card className="p-3">
                <Form.Control 
                  type="text" 
                  placeholder="Filter genes..." 
                  className="mb-3"
                  value={geneFilter}
                  onChange={(e) => setGeneFilter(e.target.value)}
                />
                <div style={{maxHeight: '400px', overflowY: 'auto'}}>
                   <ListGroup 
                     genes={genes.filter(g => g.toLowerCase().includes(geneFilter.toLowerCase()))} 
                     onSelect={(g) => {setGeneSearch(g); handleGeneSearch(g);}} 
                   />
                </div>
             </Card>
          </Tab>

          <Tab eventKey="diseases" title={`Diseases (${diseases.length})`}>
             <Card className="p-3">
                <Form.Control 
                  type="text" 
                  placeholder="Filter diseases..." 
                  className="mb-3"
                  value={diseaseFilter}
                  onChange={(e) => setDiseaseFilter(e.target.value)}
                />
                <div style={{maxHeight: '400px', overflowY: 'auto'}}>
                   {diseases
                     .filter(d => 
                       d.name.toLowerCase().includes(diseaseFilter.toLowerCase()) || 
                       d.id.toLowerCase().includes(diseaseFilter.toLowerCase())
                     )
                     .map(d => (
                       <div key={d.id} className="small mb-1 text-primary" style={{cursor: 'pointer', textDecoration: 'underline'}} onClick={() => handleDiseaseSearch(d.id)}>
                          {d.name} ({d.id})
                       </div>
                     ))
                   }
                </div>
             </Card>
          </Tab>

          <Tab eventKey="results" title={`Results (${results.length})`}>
            <div className="mt-2">
              <Table striped bordered hover responsive>
                <thead>
                  <tr>
                    <th>ID</th>
                    <th>Gene</th>
                    <th>Disease</th>
                    <th>Source</th>
                    <th>Score</th>
                  </tr>
                </thead>
                <tbody>
                  {results.map((r) => (
                    <tr key={r.id} onClick={() => handleRowClick(r.id)} style={{cursor: 'pointer'}}>
                      <td>{r.external_id}</td>
                      <td>{r.gene_symbol || 'N/A'}</td>
                      <td>
                        {r.disease_is_suggested === 1 ? (
                          <span>
                            {String(r.disease_name || r.disease_id).split('|').map((d, idx) => (
                              <div key={idx}>
                                <em>{d}</em>
                                <Badge bg="info" className="ms-1" style={{fontSize: '0.6rem'}}>Suggested</Badge>
                              </div>
                            ))}
                          </span>
                        ) : (
                          r.disease_name || r.disease_id || 'N/A'
                        )}
                      </td>
                      <td>
                        <Badge bg={r.source.includes('Saudi') ? 'primary' : 'secondary'}>
                          {r.source}
                        </Badge>
                      </td>
                      <td>{r.score ? r.score.toFixed(4) : '-'}</td>
                    </tr>
                  ))}
                  {results.length === 0 && (
                    <tr>
                      <td colSpan={5} className="text-center text-muted">No results to display. Start a search above.</td>
                    </tr>
                  )}
                </tbody>
              </Table>
            </div>
          </Tab>
        </Tabs>

        <Modal show={showDetail} onHide={() => {setShowDetail(false); setSelectedPacket(null);}} size="lg">
          <Modal.Header closeButton>
            <Modal.Title>
              Phenopacket: {selectedPacket?.external_id || 'Loading...'}
            </Modal.Title>
          </Modal.Header>
          <Modal.Body>
            {detailLoading ? (
              <div className="text-center p-5">
                <Spinner animation="border" />
              </div>
            ) : selectedPacket ? (
              <Tabs defaultActiveKey="overview" id="detail-tabs">
                <Tab eventKey="overview" title="Overview">
                  <div className="mt-3">
                    <PhenopacketView data={selectedPacket.raw_json} />
                  </div>
                </Tab>
                <Tab eventKey="json" title="Raw JSON">
                  <div className="mt-3">
                    <pre style={{backgroundColor: '#f8f9fa', padding: '15px', borderRadius: '5px', maxHeight: '500px', overflow: 'auto'}}>
                      {JSON.stringify(selectedPacket.raw_json, null, 2)}
                    </pre>
                  </div>
                </Tab>
              </Tabs>
            ) : (
              <p>No data available.</p>
            )}
          </Modal.Body>
          <Modal.Footer>
            <Button variant="secondary" onClick={() => setShowDetail(false)}>
              Close
            </Button>
          </Modal.Footer>
        </Modal>
      </Container>
    </div>
  );
}

const ListGroup = ({genes, onSelect}: {genes: string[], onSelect: (g: string) => void}) => (
  <div>
    {genes.map(g => (
      <div key={g} className="small mb-1" style={{cursor: 'pointer'}} onClick={() => onSelect(g)}>
        {g}
      </div>
    ))}
  </div>
)

const PhenopacketView = ({ data }: { data: any }) => {
  const subject = data.subject || {};
  const phenotypes = data.phenotypicFeatures || [];
  const interpretations = data.interpretations || [];

  return (
    <div>
      <section className="mb-4">
        <h6 className="border-bottom pb-1 text-primary">Subject</h6>
        <div className="row small">
          <div className="col-md-4"><strong>ID:</strong> {subject.id}</div>
          <div className="col-md-4"><strong>Sex:</strong> {subject.sex}</div>
          <div className="col-md-4"><strong>Age:</strong> {subject.timeAtLastEncounter?.age?.iso8601duration || 'N/A'}</div>
        </div>
      </section>

      <section className="mb-4">
        <h6 className="border-bottom pb-1 text-primary">Phenotypic Features</h6>
        <div className="d-flex flex-wrap gap-2">
          {phenotypes.map((p: any, i: number) => (
            <Badge key={i} bg="light" text="dark" className="fw-normal border">
              {p.type?.label} ({p.type?.id})
            </Badge>
          ))}
          {phenotypes.length === 0 && <span className="text-muted small">None recorded</span>}
        </div>
      </section>

      <section className="mb-4">
        <h6 className="border-bottom pb-1 text-primary">Interpretations</h6>
        {interpretations.map((interp: any, i: number) => (
          <div key={i} className="card p-2 mb-2 bg-light small">
            <div className="d-flex justify-content-between">
              <strong>{interp.diagnosis?.disease?.label || 'Diagnosis'}</strong>
              <Badge bg="success">{interp.progressStatus}</Badge>
            </div>
            {interp.diagnosis?.genomicInterpretations?.map((gi: any, j: number) => {
              const vd = gi.variantInterpretation?.variationDescriptor || {};
              return (
                <div key={j} className="mt-2 ps-2 border-start border-3">
                  <div><strong>Gene:</strong> {vd.geneContext?.symbol}</div>
                  <div><strong>Variant:</strong> {vd.expressions?.[0]?.value}</div>
                  <div><strong>Significance:</strong> {gi.variantInterpretation?.clinicalSignificance}</div>
                </div>
              );
            })}
          </div>
        ))}
        {interpretations.length === 0 && <span className="text-muted small">No interpretations</span>}
      </section>
    </div>
  );
};

export default App;