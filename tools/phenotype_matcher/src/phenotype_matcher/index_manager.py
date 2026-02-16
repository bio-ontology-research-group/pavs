"""
Index management for the Phenotype Matcher tool.

Manages embedding indices for different models (HPO, MONDO, OMIM).
Allows building, rebuilding, and switching between indices.
"""

import os
import pickle
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
import numpy as np
from pyhpo import Ontology
import pronto
from sentence_transformers import SentenceTransformer


# Available embedding models
EMBEDDING_MODELS = {
    "fast": "all-MiniLM-L6-v2",
    "balanced": "all-mpnet-base-v2",
    "accurate": "BAAI/bge-large-en-v1.5",
    "medical": "pritamdeka/S-PubMedBert-MS-MARCO",
    "biobert": "dmis-lab/biobert-base-cased-v1.2",
}


class IndexManager:
    """
    Manages embedding indices for different ontologies and models.

    Indices are organized by:
    - Ontology type (HPO, MONDO, OMIM)
    - Embedding model

    Directory structure:
    cache_dir/
        hpo/
            embeddings_fast.pkl
            embeddings_balanced.pkl
            embeddings_accurate.pkl
        mondo/
            embeddings_fast.pkl
            embeddings_balanced.pkl
        omim/
            embeddings_fast.pkl
            embeddings_balanced.pkl
    """

    def __init__(self, cache_dir: str = "data/graph_rag_cache", device: str = "cpu"):
        """
        Initialize index manager.

        Args:
            cache_dir: Base directory for storing indices
            device: Compute device (cpu or cuda)
        """
        self.cache_dir = Path(cache_dir)
        self.device = device
        self.logger = logging.getLogger(__name__)

        # Create directory structure
        self.hpo_dir = self.cache_dir / "hpo"
        self.mondo_dir = self.cache_dir / "mondo"
        self.omim_dir = self.cache_dir / "omim"

        for directory in [self.hpo_dir, self.mondo_dir, self.omim_dir]:
            directory.mkdir(parents=True, exist_ok=True)

    def get_index_path(self, ontology: str, model_name: str) -> Path:
        """Get path to index file for ontology and model."""
        # Resolve model name from preset
        full_model_name = EMBEDDING_MODELS.get(model_name, model_name)
        safe_name = full_model_name.replace("/", "_")

        if ontology.lower() == "hpo":
            return self.hpo_dir / f"embeddings_{safe_name}.pkl"
        elif ontology.lower() == "mondo":
            return self.mondo_dir / f"embeddings_{safe_name}.pkl"
        elif ontology.lower() == "omim":
            return self.omim_dir / f"embeddings_{safe_name}.pkl"
        else:
            raise ValueError(f"Unknown ontology: {ontology}")

    def index_exists(self, ontology: str, model_name: str) -> bool:
        """Check if index exists for ontology and model."""
        index_path = self.get_index_path(ontology, model_name)
        return index_path.exists()

    def list_indices(self) -> Dict[str, List[str]]:
        """
        List all existing indices.

        Returns:
            Dictionary mapping ontology -> list of models with indices
        """
        indices = {
            "hpo": [],
            "mondo": [],
            "omim": [],
        }

        # Check each ontology directory
        for ontology, directory in [
            ("hpo", self.hpo_dir),
            ("mondo", self.mondo_dir),
            ("omim", self.omim_dir),
        ]:
            if directory.exists():
                for file in directory.glob("embeddings_*.pkl"):
                    # Extract model name from filename
                    model_part = file.stem.replace("embeddings_", "")
                    indices[ontology].append(model_part)

        return indices

    def load_ontology_terms(self, ontology: str) -> List[Dict[str, Any]]:
        """
        Load terms from specified ontology.

        Args:
            ontology: Ontology type (hpo, mondo, omim)

        Returns:
            List of term dictionaries
        """
        terms = []

        if ontology.lower() == "hpo":
            self.logger.info("Loading HPO ontology...")
            try:
                _ = Ontology()  # Initialize pyhpo
                for term in Ontology:
                    terms.append(
                        {
                            "id": term.id,
                            "name": term.name,
                            "definition": term.definition or "",
                            "synonyms": list(term.synonym)
                            if hasattr(term, "synonym")
                            else [],
                            "parents": [p.id for p in term.parents],
                            "children": [c.id for c in term.children],
                            "source": "HPO",
                        }
                    )
                self.logger.info(f"Loaded {len(terms)} HPO terms")
            except Exception as e:
                self.logger.error(f"Failed to load HPO: {e}")

        elif ontology.lower() == "mondo":
            self.logger.info("Loading MONDO ontology...")
            mondo_path = "ontology/mondo.obo"
            if os.path.exists(mondo_path):
                try:
                    ms = pronto.Ontology(mondo_path)
                    for term in ms.terms():
                        if term.id.startswith("MONDO:"):
                            terms.append(
                                {
                                    "id": term.id,
                                    "name": term.name,
                                    "definition": str(term.definition)
                                    if term.definition
                                    else "",
                                    "synonyms": [s.description for s in term.synonyms],
                                    "parents": [
                                        p.id
                                        for p in term.superclasses(distance=1).to_set()
                                    ],
                                    "children": [
                                        c.id
                                        for c in term.subclasses(distance=1).to_set()
                                    ],
                                    "source": "MONDO",
                                }
                            )
                    self.logger.info(f"Loaded {len(terms)} MONDO terms")
                except Exception as e:
                    self.logger.error(f"Failed to load MONDO: {e}")
            else:
                self.logger.warning(f"MONDO ontology not found at {mondo_path}")

        elif ontology.lower() == "omim":
            self.logger.info("Loading OMIM data...")
            # OMIM is loaded from phenotype.hpoa file
            omim_path = "data/phenotype.hpoa"
            if os.path.exists(omim_path):
                try:
                    # Parse phenotype.hpoa for OMIM disease terms
                    with open(omim_path, "r") as f:
                        for line in f:
                            if line.startswith("#"):
                                continue
                            parts = line.strip().split("\t")
                            if len(parts) >= 2:
                                omim_id = parts[0]
                                disease_name = parts[1] if len(parts) > 1 else ""

                                if omim_id.startswith("OMIM:"):
                                    terms.append(
                                        {
                                            "id": omim_id,
                                            "name": disease_name,
                                            "definition": "",
                                            "synonyms": [],
                                            "parents": [],
                                            "children": [],
                                            "source": "OMIM",
                                        }
                                    )

                    # Remove duplicates
                    seen = set()
                    unique_terms = []
                    for term in terms:
                        if term["id"] not in seen:
                            seen.add(term["id"])
                            unique_terms.append(term)
                    terms = unique_terms

                    self.logger.info(f"Loaded {len(terms)} OMIM terms")
                except Exception as e:
                    self.logger.error(f"Failed to load OMIM: {e}")
            else:
                self.logger.warning(f"OMIM data not found at {omim_path}")

        return terms

    def build_index(
        self,
        ontology: str,
        model_name: str,
        force_rebuild: bool = False,
    ) -> bool:
        """
        Build embedding index for specified ontology and model.

        Args:
            ontology: Ontology type (hpo, mondo, omim)
            model_name: Embedding model name or preset
            force_rebuild: Rebuild even if index exists

        Returns:
            True if successful, False otherwise
        """
        index_path = self.get_index_path(ontology, model_name)

        # Check if already exists
        if index_path.exists() and not force_rebuild:
            self.logger.info(f"Index already exists: {index_path}")
            return True

        # Load ontology terms
        self.logger.info(f"Building index for {ontology} with model {model_name}...")
        terms = self.load_ontology_terms(ontology)

        if not terms:
            self.logger.error(f"No terms loaded for {ontology}")
            return False

        # Resolve model name
        full_model_name = EMBEDDING_MODELS.get(model_name, model_name)

        # Load sentence transformer
        try:
            self.logger.info(
                f"Loading model: {full_model_name} on device: {self.device}"
            )
            model = SentenceTransformer(full_model_name, device=self.device)
            self.logger.info(f"Model loaded successfully")
        except Exception as e:
            self.logger.error(f"Failed to load model: {e}")
            if self.device == "cuda":
                self.logger.warning("CUDA failed, trying CPU...")
                try:
                    model = SentenceTransformer(full_model_name, device="cpu")
                    self.device = "cpu"
                    self.logger.info("Successfully loaded model on CPU")
                except Exception as e2:
                    self.logger.error(f"Failed on CPU too: {e2}")
                    return False
            else:
                return False

        # Prepare texts for embedding
        # CRITICAL: Include ALL synonyms for comprehensive matching
        # We include every synonym (no limit) because:
        # - Mean synonyms per term: 1.21
        # - Median: 1
        # - Max: 28 (only 6 terms have >20 synonyms)
        # This ensures the tool can match layperson terms, alternative phrasings,
        # and clinical synonyms without any loss of coverage.
        texts = []
        for t in terms:
            txt = f"Name: {t['name']}."
            if t["definition"]:
                txt += f" Definition: {t['definition']}."
            if t["synonyms"]:
                # Include ALL synonyms (no filtering - 100% coverage)
                synonyms_text = ", ".join(t["synonyms"])
                txt += f" Synonyms: {synonyms_text}."
            texts.append(txt)

        # Compute embeddings
        self.logger.info(f"Computing embeddings for {len(terms)} terms...")

        # Estimate time based on device and term count
        # Rough estimates: CPU ~100 terms/sec, CUDA ~500 terms/sec
        if self.device == "cuda":
            est_time_sec = len(terms) / 500
        else:
            est_time_sec = len(terms) / 100

        est_min = int(est_time_sec / 60)
        if est_min > 0:
            self.logger.info(f"Estimated time: ~{est_min} minutes")
        else:
            self.logger.info(f"Estimated time: <1 minute")

        self.logger.info(f"Model: {full_model_name}, Device: {self.device}")
        self.logger.info(f"Press Ctrl+C to cancel if needed")

        try:
            embeddings = model.encode(
                texts,
                show_progress_bar=True,
                convert_to_numpy=True,
                normalize_embeddings=True,
                batch_size=32,  # Process in smaller batches
            )
        except KeyboardInterrupt:
            self.logger.warning("Embedding computation interrupted by user")
            return False
        except Exception as e:
            self.logger.error(f"Failed to compute embeddings: {e}")
            import traceback

            self.logger.debug(traceback.format_exc())
            return False

        # Save to cache
        try:
            with open(index_path, "wb") as f:
                pickle.dump(
                    {
                        "embeddings": embeddings,
                        "terms": terms,
                        "model": full_model_name,
                        "ontology": ontology,
                    },
                    f,
                )
            self.logger.info(f"Index saved to: {index_path}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to save index: {e}")
            return False

    def build_all_indices(
        self,
        models: Optional[List[str]] = None,
        ontologies: Optional[List[str]] = None,
        force_rebuild: bool = False,
    ) -> Dict[str, Dict[str, bool]]:
        """
        Build indices for all combinations of ontologies and models.

        Args:
            models: List of models to build (default: all presets)
            ontologies: List of ontologies to build (default: all)
            force_rebuild: Rebuild even if indices exist

        Returns:
            Dictionary mapping ontology -> model -> success status
        """
        if models is None:
            models = list(EMBEDDING_MODELS.keys())

        if ontologies is None:
            ontologies = ["hpo", "mondo", "omim"]

        results = {}
        total = len(ontologies) * len(models)
        current = 0

        for ontology in ontologies:
            results[ontology] = {}
            for model in models:
                current += 1
                self.logger.info(f"[{current}/{total}] Building {ontology} / {model}")
                success = self.build_index(ontology, model, force_rebuild)
                results[ontology][model] = success

        return results

    def delete_index(self, ontology: str, model_name: str) -> bool:
        """Delete index for specified ontology and model."""
        index_path = self.get_index_path(ontology, model_name)

        if index_path.exists():
            try:
                index_path.unlink()
                self.logger.info(f"Deleted index: {index_path}")
                return True
            except Exception as e:
                self.logger.error(f"Failed to delete index: {e}")
                return False
        else:
            self.logger.warning(f"Index does not exist: {index_path}")
            return False

    def get_index_info(
        self, ontology: str, model_name: str
    ) -> Optional[Dict[str, Any]]:
        """Get information about an index."""
        index_path = self.get_index_path(ontology, model_name)

        if not index_path.exists():
            return None

        try:
            with open(index_path, "rb") as f:
                data = pickle.load(f)

            return {
                "path": str(index_path),
                "ontology": data.get("ontology"),
                "model": data.get("model"),
                "num_terms": len(data.get("terms", [])),
                "embedding_shape": data.get("embeddings").shape,
                "size_mb": index_path.stat().st_size / (1024 * 1024),
            }
        except Exception as e:
            self.logger.error(f"Failed to read index info: {e}")
            return None


def main():
    """CLI for index management."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Phenotype Matcher Index Manager",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "command",
        choices=["build", "rebuild", "list", "info", "delete"],
        help="Command to execute",
    )
    parser.add_argument(
        "--ontology",
        choices=["hpo", "mondo", "omim", "all"],
        default="all",
        help="Ontology to operate on (default: all)",
    )
    parser.add_argument(
        "--model",
        choices=list(EMBEDDING_MODELS.keys()) + ["all"],
        default="balanced",
        help="Embedding model (default: balanced)",
    )
    parser.add_argument(
        "--cache-dir",
        type=str,
        default="data/graph_rag_cache",
        help="Cache directory (default: data/graph_rag_cache)",
    )
    parser.add_argument(
        "--device",
        choices=["cpu", "cuda"],
        default="cpu",
        help="Compute device (default: cpu)",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    # Initialize manager
    manager = IndexManager(cache_dir=args.cache_dir, device=args.device)

    # Execute command
    if args.command == "build":
        if args.ontology == "all" or args.model == "all":
            ontologies = (
                ["hpo", "mondo", "omim"] if args.ontology == "all" else [args.ontology]
            )
            models = (
                list(EMBEDDING_MODELS.keys()) if args.model == "all" else [args.model]
            )
            results = manager.build_all_indices(
                models=models, ontologies=ontologies, force_rebuild=False
            )
            print("\nBuild Results:")
            for ont, model_results in results.items():
                print(f"\n{ont.upper()}:")
                for model, success in model_results.items():
                    status = "✓" if success else "✗"
                    print(f"  {status} {model}")
        else:
            success = manager.build_index(
                args.ontology, args.model, force_rebuild=False
            )
            if success:
                print(f"✓ Successfully built index for {args.ontology} / {args.model}")
            else:
                print(f"✗ Failed to build index for {args.ontology} / {args.model}")

    elif args.command == "rebuild":
        if args.ontology == "all" or args.model == "all":
            ontologies = (
                ["hpo", "mondo", "omim"] if args.ontology == "all" else [args.ontology]
            )
            models = (
                list(EMBEDDING_MODELS.keys()) if args.model == "all" else [args.model]
            )
            results = manager.build_all_indices(
                models=models, ontologies=ontologies, force_rebuild=True
            )
            print("\nRebuild Results:")
            for ont, model_results in results.items():
                print(f"\n{ont.upper()}:")
                for model, success in model_results.items():
                    status = "✓" if success else "✗"
                    print(f"  {status} {model}")
        else:
            success = manager.build_index(args.ontology, args.model, force_rebuild=True)
            if success:
                print(
                    f"✓ Successfully rebuilt index for {args.ontology} / {args.model}"
                )
            else:
                print(f"✗ Failed to rebuild index for {args.ontology} / {args.model}")

    elif args.command == "list":
        indices = manager.list_indices()
        print("\nExisting Indices:")
        for ontology, models in indices.items():
            print(f"\n{ontology.upper()}:")
            if models:
                for model in sorted(models):
                    info = manager.get_index_info(
                        ontology, model.replace("embeddings_", "").replace(".pkl", "")
                    )
                    if info:
                        print(
                            f"  ✓ {model} ({info['num_terms']} terms, {info['size_mb']:.1f} MB)"
                        )
                    else:
                        print(f"  ✓ {model}")
            else:
                print("  (none)")

    elif args.command == "info":
        if args.ontology == "all":
            print("Error: --ontology must be specified for 'info' command")
            return

        info = manager.get_index_info(args.ontology, args.model)
        if info:
            print(f"\nIndex Information:")
            print(f"  Ontology: {info['ontology']}")
            print(f"  Model: {info['model']}")
            print(f"  Terms: {info['num_terms']}")
            print(f"  Embedding shape: {info['embedding_shape']}")
            print(f"  Size: {info['size_mb']:.1f} MB")
            print(f"  Path: {info['path']}")
        else:
            print(f"Index not found for {args.ontology} / {args.model}")

    elif args.command == "delete":
        if args.ontology == "all":
            print("Error: --ontology must be specified for 'delete' command")
            return

        success = manager.delete_index(args.ontology, args.model)
        if success:
            print(f"✓ Deleted index for {args.ontology} / {args.model}")
        else:
            print(f"✗ Failed to delete index for {args.ontology} / {args.model}")


if __name__ == "__main__":
    main()
