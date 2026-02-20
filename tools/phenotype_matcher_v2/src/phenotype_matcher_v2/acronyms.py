"""
Common medical acronyms and their expansions.
Copied from v1 and kept as-is.
"""

from typing import Dict, List

# Acronym → List of possible expansions (HPO labels or search terms)
MEDICAL_ACRONYMS: Dict[str, List[str]] = {
    "ASD": [
        "Atrial septal defect",
        "Autism spectrum disorder",
        "Autistic behavior",
    ],
    "AD": [
        "Autosomal dominant",
        "Alzheimer disease",
    ],
    "AR": [
        "Autosomal recessive",
        "Aortic regurgitation",
    ],
    "BPD": [
        "Bronchopulmonary dysplasia",
        "Bipolar disorder",
    ],
    "CHD": [
        "Congenital heart disease",
        "Congenital heart defect",
    ],
    "CKD": [
        "Chronic kidney disease",
    ],
    "CP": [
        "Cerebral palsy",
        "Cleft palate",
    ],
    "DD": [
        "Developmental delay",
        "Global developmental delay",
    ],
    "DM": [
        "Diabetes mellitus",
        "Myotonic dystrophy",
    ],
    "FTT": [
        "Failure to thrive",
    ],
    "GDD": [
        "Global developmental delay",
    ],
    "GER": [
        "Gastroesophageal reflux",
    ],
    "HD": [
        "Huntington disease",
        "Hip dysplasia",
    ],
    "HOB": [
        "Head of bed",
    ],
    "ID": [
        "Intellectual disability",
    ],
    "IQ": [
        "Intelligence quotient",
    ],
    "MD": [
        "Muscular dystrophy",
    ],
    "MR": [
        "Mental retardation",
        "Intellectual disability",
    ],
    "MS": [
        "Multiple sclerosis",
        "Mitral stenosis",
    ],
    "OCD": [
        "Obsessive-compulsive disorder",
    ],
    "PDA": [
        "Patent ductus arteriosus",
    ],
    "PICU": [
        "Pediatric intensive care unit",
    ],
    "SCD": [
        "Sickle cell disease",
        "Sudden cardiac death",
    ],
    "SD": [
        "Standard deviation",
        "Speech delay",
    ],
    "TOF": [
        "Tetralogy of Fallot",
    ],
    "VSD": [
        "Ventricular septal defect",
    ],
}

# Acronyms that are NOT phenotypes
NON_PHENOTYPE_ACRONYMS = {
    "PICU", "ICU", "NICU", "ER", "OR", "CT", "MRI",
    "EEG", "ECG", "EKG", "IQ", "SD", "HOB", "NPO",
    "PRN", "BID", "TID", "QID",
}


def expand_acronym(acronym: str) -> List[str]:
    """Return possible expansions for a medical acronym, or [] if non-phenotype."""
    acronym_upper = acronym.upper().strip()
    if acronym_upper in NON_PHENOTYPE_ACRONYMS:
        return []
    return MEDICAL_ACRONYMS.get(acronym_upper, [acronym])


def is_likely_acronym(text: str) -> bool:
    """Return True if *text* looks like a medical acronym."""
    text = text.strip()
    if len(text) <= 5 and text.isupper():
        return True
    if text.upper() in MEDICAL_ACRONYMS:
        return True
    return False


def should_expand_for_search(text: str) -> bool:
    """Return True if the acronym is ambiguous and all expansions should be searched."""
    text = text.strip()
    acronym_upper = text.upper()
    if acronym_upper in MEDICAL_ACRONYMS:
        return len(MEDICAL_ACRONYMS[acronym_upper]) > 1
    return False
