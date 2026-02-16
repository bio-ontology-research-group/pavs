"""
Common medical acronyms and their expansions for disambiguation.

These are hard-coded mappings for ambiguous acronyms that may not
be retrieved well by semantic search alone.
"""

from typing import Dict, List


# Acronym → List of possible expansions (HPO labels or search terms)
MEDICAL_ACRONYMS: Dict[str, List[str]] = {
    # A
    "ASD": [
        "Atrial septal defect",  # Cardiac
        "Autism spectrum disorder",  # Neurodevelopmental
        "Autistic behavior",  # Alternative neurodevelopmental term
    ],
    "AD": [
        "Autosomal dominant",
        "Alzheimer disease",
    ],
    "AR": [
        "Autosomal recessive",
        "Aortic regurgitation",
    ],
    # B
    "BPD": [
        "Bronchopulmonary dysplasia",  # Respiratory
        "Bipolar disorder",  # Psychiatric
    ],
    # C
    "CHD": [
        "Congenital heart disease",  # Cardiac
        "Congenital heart defect",  # Alternative
    ],
    "CKD": [
        "Chronic kidney disease",
    ],
    "CP": [
        "Cerebral palsy",
        "Cleft palate",
    ],
    # D
    "DD": [
        "Developmental delay",
        "Global developmental delay",
    ],
    "DM": [
        "Diabetes mellitus",
        "Myotonic dystrophy",
    ],
    # F
    "FTT": [
        "Failure to thrive",
    ],
    # G
    "GDD": [
        "Global developmental delay",
    ],
    "GER": [
        "Gastroesophageal reflux",
    ],
    # H
    "HD": [
        "Huntington disease",
        "Hip dysplasia",
    ],
    "HOB": [
        "Head of bed",  # Not a phenotype - should be filtered
    ],
    # I
    "ID": [
        "Intellectual disability",
    ],
    "IQ": [
        "Intelligence quotient",  # Not a phenotype
    ],
    # M
    "MD": [
        "Muscular dystrophy",
    ],
    "MR": [
        "Mental retardation",  # Legacy term
        "Intellectual disability",  # Modern term
    ],
    "MS": [
        "Multiple sclerosis",
        "Mitral stenosis",
    ],
    # O
    "OCD": [
        "Obsessive-compulsive disorder",
    ],
    # P
    "PDA": [
        "Patent ductus arteriosus",  # Cardiac
    ],
    "PICU": [
        "Pediatric intensive care unit",  # Location, not phenotype
    ],
    # S
    "SCD": [
        "Sickle cell disease",
        "Sudden cardiac death",
    ],
    "SD": [
        "Standard deviation",  # Not a phenotype
        "Speech delay",
    ],
    # T
    "TOF": [
        "Tetralogy of Fallot",
    ],
    # V
    "VSD": [
        "Ventricular septal defect",
    ],
}


# Acronyms that are NOT phenotypes (locations, measurements, etc.)
NON_PHENOTYPE_ACRONYMS = {
    "PICU",
    "ICU",
    "NICU",
    "ER",
    "OR",
    "CT",
    "MRI",
    "EEG",
    "ECG",
    "EKG",
    "IQ",
    "SD",
    "HOB",
    "NPO",
    "PRN",
    "BID",
    "TID",
    "QID",
}


def expand_acronym(acronym: str) -> List[str]:
    """
    Get possible expansions for a medical acronym.

    Args:
        acronym: Medical acronym (e.g., "ASD", "DD")

    Returns:
        List of possible expansions, or [acronym] if unknown
    """
    acronym_upper = acronym.upper().strip()

    # Check if it's a non-phenotype
    if acronym_upper in NON_PHENOTYPE_ACRONYMS:
        return []

    # Return expansions or original if not found
    return MEDICAL_ACRONYMS.get(acronym_upper, [acronym])


def is_likely_acronym(text: str) -> bool:
    """
    Check if text is likely a medical acronym.

    Args:
        text: Text to check

    Returns:
        True if likely an acronym
    """
    text = text.strip()

    # Short, all caps or mixed case
    if len(text) <= 5 and text.isupper():
        return True

    # Known acronym
    if text.upper() in MEDICAL_ACRONYMS:
        return True

    return False


def should_expand_for_search(text: str) -> bool:
    """
    Determine if we should expand acronym for RAG search.

    For ambiguous acronyms, we want to search for ALL expansions
    and then disambiguate based on context.

    Args:
        text: Phenotype term

    Returns:
        True if we should expand and search for all variants
    """
    text = text.strip()

    # Check if it's a known ambiguous acronym
    acronym_upper = text.upper()
    if acronym_upper in MEDICAL_ACRONYMS:
        expansions = MEDICAL_ACRONYMS[acronym_upper]
        # If multiple expansions, it's ambiguous
        return len(expansions) > 1

    return False
