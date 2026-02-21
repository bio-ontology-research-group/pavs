
@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
@Grab(group='org.slf4j', module='slf4j-simple', version='1.7.36')

import java.util.*
import java.util.concurrent.*
import java.util.stream.*
import org.openrdf.model.URI
import slib.sglib.io.loader.bio.obo.OboLoader
import slib.sglib.model.graph.G
import slib.sglib.model.impl.graph.memory.GraphMemory
import slib.sglib.model.impl.repo.URIFactoryMemory
import slib.sglib.model.repo.URIFactory
import slib.sml.sm.core.engine.SM_Engine
import slib.sml.sm.core.metrics.ic.utils.ICconf
import slib.sml.sm.core.utils.SMconf
import slib.sml.sm.core.metrics.utils.SMConstants

// Paths
def hpoOboPath = "ontology/hp.obo"
def genesToPhenotypePath = "data/genes_to_phenotype.txt"
def casesPath = "data/combined_annotated.tsv"
def outputPath = "analysis/phenotype_similarity_results.csv"

println "Loading HPO Ontology..."
URIFactory factory = URIFactoryMemory.getSingleton()
URI graphURI = factory.getURI("http://hp/")
G graph = new GraphMemory(graphURI)
OboLoader loader = new OboLoader()
loader.populateGraph(hpoOboPath, graph)
SM_Engine engine = new SM_Engine(graph)

println "Loading Gene to Phenotype mappings..."
def geneToHpos = [:]
new File(genesToPhenotypePath).eachLine { line ->
    if (line.startsWith("#") || line.startsWith("ncbi_gene_id")) return
    def parts = line.split("\t")
    if (parts.size() < 3) return
    def gene = parts[1]
    def hpoId = parts[2].replace(":", "_")
    def uri = factory.getURI("http://purl.obolibrary.org/obo/" + hpoId)
    if (graph.containsVertex(uri)) {
        if (!geneToHpos[gene]) geneToHpos[gene] = [] as Set
        geneToHpos[gene] << uri
    }
}

println "Loading cases..."
def cases = []
def headers = []
new File(casesPath).eachLine { line, idx ->
    def parts = line.split("\t")
    if (idx == 1) {
        headers = parts
        return
    }
    def row = [:]
    for (int i = 0; i < [headers.size(), parts.size()].min(); i++) {
        row[headers[i]] = parts[i]
    }
    if (row.source_file != "ddd-variants" && row.phenotypes_present_ids && row.gene_symbol) {
        if (geneToHpos.containsKey(row.gene_symbol)) {
            cases << row
        }
    }
}

// Pre-parse case phenotype sets
def caseToHpos = [:]
cases.each { row ->
    def hpoStr = row.phenotypes_present_ids
    def uris = hpoStr.split(/[|,]/).collect { it.trim() }
        .findAll { it.startsWith("HP:") }
        .collect { factory.getURI("http://purl.obolibrary.org/obo/" + it.replace(":", "_")) }
        .findAll { graph.containsVertex(it) } as Set
    if (uris) caseToHpos[row.pavs_id] = uris
}
cases = cases.findAll { caseToHpos.containsKey(it.pavs_id) }

println "Cases with valid HPOs: ${cases.size()}"

// Configure SML
ICconf icConf = new ICconf("Sanchez", SMConstants.IC_ID_SANCHEZ)
SMconf smConfLin = new SMconf("Lin", SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN)
smConfLin.setICconf(icConf)

// Pre-calculate unique terms
def allTerms = (caseToHpos.values().flatten() + geneToHpos.values().flatten()) as Set
println "Unique terms involved: ${allTerms.size()}"

println "Pre-calculating pairwise Lin similarity..."
def termIndices = [:]
def termList = allTerms.toList()
termList.eachWithIndex { term, i -> termIndices[term] = i }
float[] simMatrix = new float[allTerms.size() * allTerms.size()]

IntStream.range(0, termList.size()).parallel().forEach { i ->
    URI u1 = termList[i]
    for (int j = i; j < termList.size(); j++) {
        URI u2 = termList[j]
        try {
            float s = (float) engine.compare(smConfLin, u1, u2)
            simMatrix[i * termList.size() + j] = s
            simMatrix[j * termList.size() + i] = s
        } catch (Exception e) {}
    }
}

println "Starting BMA analysis..."
def results = Collections.synchronizedList([])
def allGenes = geneToHpos.keySet().toList()

// BMA Helper using simMatrix
def getBMA = { set1, set2 ->
    if (!set1 || !set2) return 0.0
    
    double sum1 = 0
    for (t1 in set1) {
        int idx1 = termIndices[t1]
        float maxSim = 0
        for (t2 in set2) {
            int idx2 = termIndices[t2]
            float s = simMatrix[idx1 * termList.size() + idx2]
            if (s > maxSim) maxSim = s
        }
        sum1 += maxSim
    }
    
    double sum2 = 0
    for (t2 in set2) {
        int idx2 = termIndices[t2]
        float maxSim = 0
        for (t1 in set1) {
            int idx1 = termIndices[t1]
            float s = simMatrix[idx1 * termList.size() + idx2]
            if (s > maxSim) maxSim = s
        }
        sum2 += maxSim
    }
    
    return (sum1 / set1.size() + sum2 / set2.size()) / 2.0
}

def start = System.currentTimeMillis()
int total = cases.size()
int count = 0

cases.parallelStream().forEach { row ->
    def set1 = caseToHpos[row.pavs_id]
    def correctGene = row.gene_symbol
    def similarities = []
    
    for (gene in allGenes) {
        def set2 = geneToHpos[gene]
        double sim = getBMA(set1, set2)
        similarities << [gene: gene, sim: sim]
    }
    
    similarities.sort { a, b -> b.sim <=> a.sim ?: a.gene <=> b.gene }
    
    def rank = -1
    double score = 0
    for (int i = 0; i < similarities.size(); i++) {
        if (similarities[i].gene == correctGene) {
            rank = i + 1
            score = similarities[i].sim
            break
        }
    }
    
    results << [
        pavs_id: row.pavs_id,
        gene: correct_gene,
        rank: rank,
        score: score,
        status: row.result_status,
        num_candidates: similarities.size()
    ]
    
    synchronized(this) {
        count++
        if (count % 100 == 0) {
            System.err.print "\rProcessed $count/$total cases... (" + (int)(count*100.0/total) + "%)"
        }
    }
}

println "\nSaving results to $outputPath"
new File("analysis").mkdirs()
new File(outputPath).withWriter { writer ->
    writer.writeLine("pavs_id,gene,rank,score,status,num_candidates")
    results.each { r ->
        writer.writeLine("${r.pavs_id},${r.gene},${r.rank},${r.score},${r.status},${r.num_candidates}")
    }
}

def end = System.currentTimeMillis()
println "Analysis complete in ${(end - start)/1000.0} seconds."
