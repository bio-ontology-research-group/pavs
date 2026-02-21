from pydantic import BaseModel
from typing import List, Optional, Any

class PhenotypeSearchRequest(BaseModel):
    hpo_ids: List[str]
    similarity_method: str = "resnik"  # resnik, lin, etc.
    limit: int = 20
    include_disease_phenotypes: bool = False

class PhenopacketSummary(BaseModel):
    id: int
    external_id: str
    gene_symbol: Optional[str]
    disease_id: Optional[str]
    disease_name: Optional[str]
    disease_is_suggested: int = 0
    source: str
    score: Optional[float] = None

    class Config:
        from_attributes = True

class PhenopacketDetail(PhenopacketSummary):
    raw_json: Any
    disease_phenotypes: Optional[str] = None
