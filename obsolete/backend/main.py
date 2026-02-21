from fastapi import FastAPI, Depends, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select
from sqlalchemy import or_, func
from typing import List, Optional
import os

from .database import engine, Base, get_db, AsyncSessionLocal
from .models import Phenopacket
from .schemas import PhenotypeSearchRequest, PhenopacketSummary, PhenopacketDetail
from .loader import load_data
from .similarity import calculate_similarity, search_hpo_terms, get_hpo_set

app = FastAPI(title="Saudi Phenopackets Search")

# In-memory cache for faster similarity search
packet_cache = []

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.on_event("startup")
async def startup():
    print("STARTING UP BACKEND...")
    async with engine.begin() as conn:
        # Drop and recreate tables to handle schema changes in prototype
        await conn.run_sync(Base.metadata.drop_all)
        await conn.run_sync(Base.metadata.create_all)
    
    async with AsyncSessionLocal() as db:
        await load_data(db)
        # Load all into memory for fast search
        result = await db.execute(select(Phenopacket))
        all_packets = result.scalars().all()
        for p in all_packets:
            packet_cache.append({
                "id": p.id,
                "external_id": p.external_id,
                "gene_symbol": p.gene_symbol,
                "disease_id": p.disease_id,
                "disease_name": p.disease_name,
                "disease_is_suggested": p.disease_is_suggested,
                "source": p.source,
                "hpo_set": get_hpo_set(p.phenotypes),
                "disease_hpo_set": get_hpo_set(p.disease_phenotypes)
            })
    print(f"Loaded {len(packet_cache)} phenopackets into memory cache.")

@app.get("/api/search/hpo")
async def hpo_autocomplete(q: str = Query(..., min_length=2)):
    return search_hpo_terms(q)

@app.post("/api/search/phenotype", response_model=List[PhenopacketSummary])
async def search_by_phenotype(request: PhenotypeSearchRequest, db: AsyncSession = Depends(get_db)):
    scored_results = []
    for packet in packet_cache:
        target_hpo_set = packet["hpo_set"]
        if request.include_disease_phenotypes and packet["disease_hpo_set"]:
            target_hpo_set = packet["hpo_set"] | packet["disease_hpo_set"]

        if not target_hpo_set:
            continue
            
        score = calculate_similarity(request.hpo_ids, list(target_hpo_set), method=request.similarity_method)
        if score > 0:
            scored_results.append({
                "id": packet["id"],
                "external_id": packet["external_id"],
                "gene_symbol": packet["gene_symbol"],
                "disease_id": packet["disease_id"],
                "disease_name": packet["disease_name"],
                "disease_is_suggested": packet.get("disease_is_suggested", 0),
                "source": packet["source"],
                "score": score
            })
    
    # Sort by score descending
    scored_results.sort(key=lambda x: x["score"], reverse=True)
    return scored_results[:request.limit]

@app.get("/api/search/genotype", response_model=List[PhenopacketSummary])
async def search_by_genotype(gene: Optional[str] = None, db: AsyncSession = Depends(get_db)):
    query = select(Phenopacket)
    if gene:
        query = query.where(Phenopacket.gene_symbol.ilike(f"%{gene}%"))
    
    result = await db.execute(query.limit(100))
    return result.scalars().all()

@app.get("/api/search/disease", response_model=List[PhenopacketSummary])
async def search_by_disease(disease_id: str, db: AsyncSession = Depends(get_db)):
    query = select(Phenopacket).where(Phenopacket.disease_id == disease_id)
    result = await db.execute(query.limit(100))
    return result.scalars().all()

@app.get("/api/genes", response_model=List[str])
async def list_genes(db: AsyncSession = Depends(get_db)):
    result = await db.execute(select(Phenopacket.gene_symbol).distinct().where(Phenopacket.gene_symbol != None))
    return sorted([r for r in result.scalars().all()])

@app.get("/api/diseases", response_model=List[dict])
async def list_diseases(db: AsyncSession = Depends(get_db)):
    result = await db.execute(select(Phenopacket.disease_id, Phenopacket.disease_name).distinct().where(Phenopacket.disease_id != None))
    return [{"id": r[0], "name": r[1]} for r in result.all()]

@app.get("/api/phenopacket/{id}", response_model=PhenopacketDetail)
async def get_phenopacket(id: int, db: AsyncSession = Depends(get_db)):
    result = await db.execute(select(Phenopacket).where(Phenopacket.id == id))
    packet = result.scalars().first()
    if not packet:
        raise HTTPException(status_code=404, detail="Phenopacket not found")
    return packet

# Serve React static files (optional integration)
from fastapi.staticfiles import StaticFiles
if os.path.exists("frontend/build"):
    app.mount("/", StaticFiles(directory="frontend/build", html=True), name="static")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
