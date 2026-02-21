from sqlalchemy import Column, Integer, String, Text, JSON
from .database import Base

class Phenopacket(Base):
    __tablename__ = "phenopackets"

    id = Column(Integer, primary_key=True, index=True)
    external_id = Column(String, index=True)
    gene_symbol = Column(String, index=True)
    disease_id = Column(String, index=True)
    disease_name = Column(String)
    disease_is_suggested = Column(Integer, default=0) # 0 = confirmed, 1 = suggested
    phenotypes = Column(Text)  # Comma separated HPO IDs for quick filtering
    disease_phenotypes = Column(Text) # Comma separated HPO IDs from disease associations
    source = Column(String, index=True)  # 'Saudi' or 'Global'
    file_path = Column(String)
    raw_json = Column(JSON)
