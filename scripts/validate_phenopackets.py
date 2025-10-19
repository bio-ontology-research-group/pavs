#!/usr/bin/env python3
"""
Phenopackets validator following phenopacket-tools validation workflow
"""

import json
import sys
from typing import List, Dict, Any
from dataclasses import dataclass


@dataclass
class ValidationError:
    field: str
    message: str
    severity: str = "ERROR"


class PhenopacketValidator:
    """Validates Phenopackets according to v2 schema requirements"""
    
    def __init__(self):
        self.errors: List[ValidationError] = []
        self.warnings: List[ValidationError] = []
    
    def validate_phenopacket(self, pp: Dict[str, Any], pp_id: str = "") -> bool:
        """Validate a single phenopacket"""
        self.errors = []
        self.warnings = []
        
        # Required fields
        self._check_required_field(pp, "id", pp_id)
        self._check_required_field(pp, "metaData", pp_id)
        
        # Subject validation
        if "subject" in pp:
            self._validate_subject(pp["subject"], pp_id)
        else:
            self.warnings.append(
                ValidationError("subject", f"[{pp_id}] Subject is recommended but missing", "WARNING")
            )
        
        # Phenotypic features validation
        if "phenotypicFeatures" in pp:
            self._validate_phenotypic_features(pp["phenotypicFeatures"], pp_id)
        
        # Interpretations validation
        if "interpretations" in pp:
            self._validate_interpretations(pp["interpretations"], pp_id)
        
        # Diseases validation
        if "diseases" in pp:
            self._validate_diseases(pp["diseases"], pp_id)
        
        # MetaData validation
        if "metaData" in pp:
            self._validate_metadata(pp["metaData"], pp_id)
        
        return len(self.errors) == 0
    
    def _check_required_field(self, obj: Dict, field: str, context: str):
        """Check if required field exists"""
        if field not in obj or not obj[field]:
            self.errors.append(
                ValidationError(field, f"[{context}] Required field '{field}' is missing or empty")
            )
    
    def _validate_subject(self, subject: Dict, context: str):
        """Validate subject"""
        self._check_required_field(subject, "id", f"{context}.subject")
        
        # Check sex is valid enum
        if "sex" in subject:
            valid_sex = ["UNKNOWN_SEX", "FEMALE", "MALE", "OTHER_SEX"]
            if subject["sex"] not in valid_sex:
                self.errors.append(
                    ValidationError("subject.sex", 
                                  f"[{context}] Invalid sex value: {subject['sex']}. Must be one of {valid_sex}")
                )
    
    def _validate_phenotypic_features(self, features: List[Dict], context: str):
        """Validate phenotypic features"""
        for i, feature in enumerate(features):
            if "type" not in feature:
                self.errors.append(
                    ValidationError(f"phenotypicFeatures[{i}]", 
                                  f"[{context}] Feature missing required 'type' field")
                )
                continue
            
            feature_type = feature["type"]
            if "id" not in feature_type:
                self.errors.append(
                    ValidationError(f"phenotypicFeatures[{i}].type", 
                                  f"[{context}] Feature type missing required 'id' field")
                )
            
            # Check HPO ID format
            if "id" in feature_type and not feature_type["id"].startswith("HP:"):
                self.warnings.append(
                    ValidationError(f"phenotypicFeatures[{i}].type.id", 
                                  f"[{context}] ID {feature_type['id']} doesn't appear to be an HPO term",
                                  "WARNING")
                )
    
    def _validate_interpretations(self, interpretations: List[Dict], context: str):
        """Validate interpretations"""
        for i, interp in enumerate(interpretations):
            self._check_required_field(interp, "id", f"{context}.interpretations[{i}]")
            
            # progressStatus should be valid enum
            if "progressStatus" in interp:
                valid_status = ["UNKNOWN_PROGRESS", "IN_PROGRESS", "COMPLETED", "SOLVED", "UNSOLVED"]
                if interp["progressStatus"] not in valid_status:
                    self.errors.append(
                        ValidationError(f"interpretations[{i}].progressStatus",
                                      f"[{context}] Invalid progressStatus: {interp['progressStatus']}")
                    )
            
            # Validate diagnosis if present
            if "diagnosis" in interp:
                diagnosis = interp["diagnosis"]
                if "genomicInterpretations" in diagnosis:
                    self._validate_genomic_interpretations(
                        diagnosis["genomicInterpretations"], 
                        f"{context}.interpretations[{i}].diagnosis"
                    )
    
    def _validate_genomic_interpretations(self, genomic_interps: List[Dict], context: str):
        """Validate genomic interpretations"""
        for i, gi in enumerate(genomic_interps):
            self._check_required_field(gi, "subjectOrBiosampleId", f"{context}.genomicInterpretations[{i}]")
            self._check_required_field(gi, "interpretationStatus", f"{context}.genomicInterpretations[{i}]")
    
    def _validate_diseases(self, diseases: List[Dict], context: str):
        """Validate diseases"""
        for i, disease in enumerate(diseases):
            if "term" not in disease:
                self.errors.append(
                    ValidationError(f"diseases[{i}]", 
                                  f"[{context}] Disease missing required 'term' field")
                )
    
    def _validate_metadata(self, metadata: Dict, context: str):
        """Validate metadata"""
        self._check_required_field(metadata, "created", f"{context}.metaData")
        self._check_required_field(metadata, "phenopacketSchemaVersion", f"{context}.metaData")
        
        # Resources is recommended
        if "resources" not in metadata or not metadata["resources"]:
            self.warnings.append(
                ValidationError("metaData.resources", 
                              f"[{context}] Resources field is recommended", "WARNING")
            )


def validate_phenopackets_file(file_path: str) -> bool:
    """Validate phenopackets from JSON file"""
    print(f"Validating {file_path}...")
    print("=" * 60)
    
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"❌ Invalid JSON: {e}")
        return False
    
    # Handle both single phenopacket and array
    phenopackets = data if isinstance(data, list) else [data]
    
    validator = PhenopacketValidator()
    total_errors = 0
    total_warnings = 0
    valid_count = 0
    
    for i, pp in enumerate(phenopackets):
        pp_id = pp.get('id', f'phenopacket_{i+1}')
        is_valid = validator.validate_phenopacket(pp, pp_id)
        
        if is_valid:
            valid_count += 1
        
        total_errors += len(validator.errors)
        total_warnings += len(validator.warnings)
        
        # Print first few errors/warnings
        if i < 5:  # Only show details for first 5
            if validator.errors:
                print(f"\n❌ Errors in {pp_id}:")
                for error in validator.errors[:3]:  # Show first 3 errors
                    print(f"   - {error.field}: {error.message}")
            
            if validator.warnings:
                print(f"\n⚠️  Warnings in {pp_id}:")
                for warning in validator.warnings[:3]:
                    print(f"   - {warning.field}: {warning.message}")
    
    # Summary
    print("\n" + "=" * 60)
    print("Validation Summary")
    print("=" * 60)
    print(f"Total phenopackets: {len(phenopackets)}")
    print(f"Valid: {valid_count} ({valid_count/len(phenopackets)*100:.1f}%)")
    print(f"With errors: {len(phenopackets) - valid_count} ({(len(phenopackets)-valid_count)/len(phenopackets)*100:.1f}%)")
    print(f"Total errors: {total_errors}")
    print(f"Total warnings: {total_warnings}")
    
    if valid_count == len(phenopackets):
        print("\n✅ All phenopackets are valid!")
        return True
    else:
        print(f"\n⚠️  {len(phenopackets) - valid_count} phenopackets have validation errors")
        return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python validate_phenopackets_proper.py <phenopackets.json>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    is_valid = validate_phenopackets_file(file_path)
    
    sys.exit(0 if is_valid else 1)