-- PAVS RDF bulk load script for Virtuoso isql
-- Run with: isql virtuoso:1111 dba pavs_dba /load/load_ttl.sql

-- Clear any previous load list entries for our files
DELETE FROM DB.DBA.LOAD_LIST WHERE LL_FILE LIKE '/rdf_import/%';

-- Clear existing graphs so we start clean
SPARQL CLEAR SILENT GRAPH <http://pavs.kaust.edu.sa/graph/cases>;
SPARQL CLEAR SILENT GRAPH <http://pavs.kaust.edu.sa/graph/genes>;
SPARQL CLEAR SILENT GRAPH <http://pavs.kaust.edu.sa/graph/hpoa>;
SPARQL CLEAR SILENT GRAPH <http://pavs.kaust.edu.sa/graph/hpo-ic>;
SPARQL CLEAR SILENT GRAPH <http://pavs.kaust.edu.sa/graph/literature>;

-- Register TTL files for bulk loading (ld_dir queues them; rdf_loader_run processes)
ld_dir('/rdf_import', 'cases.ttl',      'http://pavs.kaust.edu.sa/graph/cases');
ld_dir('/rdf_import', 'genes.ttl',      'http://pavs.kaust.edu.sa/graph/genes');
ld_dir('/rdf_import', 'hpoa.ttl',       'http://pavs.kaust.edu.sa/graph/hpoa');
ld_dir('/rdf_import', 'hpo_ic.ttl',     'http://pavs.kaust.edu.sa/graph/hpo-ic');
ld_dir('/rdf_import', 'literature.ttl', 'http://pavs.kaust.edu.sa/graph/literature');

-- Show what is queued
SELECT LL_FILE, LL_GRAPH FROM DB.DBA.LOAD_LIST;

-- Run the loader (processes all queued files)
rdf_loader_run();

-- Persist to disk
checkpoint;

-- Show result counts per graph
SPARQL SELECT ?g (COUNT(*) AS ?n) WHERE { GRAPH ?g { ?s ?p ?o } } GROUP BY ?g ORDER BY ?g;

exit;
