"""
Command-line interface for Phenotype Matcher v2.

Subcommands:
    match   — match one or more text strings (stdin, file, or positional args)
    batch   — process a TSV/CSV file column by column
    info    — look up an ontology term by ID or label
    check   — verify configuration, data files, and Python packages
"""

import argparse
import json
import os
import sys
from typing import List, Optional


# ---------------------------------------------------------------------------
# Shared argument setup
# ---------------------------------------------------------------------------

def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add ontology/model arguments shared by all subcommands."""
    g = parser.add_argument_group("ontology paths")
    g.add_argument("--hpo", default="ontology/hp.obo", metavar="PATH",
                   help="Path to hp.obo  [default: ontology/hp.obo]")
    g.add_argument("--mondo", default="ontology/mondo.obo", metavar="PATH",
                   help="Path to mondo.obo  [default: ontology/mondo.obo]")
    g.add_argument("--hpoa", default="data/phenotype.hpoa", metavar="PATH",
                   help="Path to phenotype.hpoa  [default: data/phenotype.hpoa]")

    g2 = parser.add_argument_group("model selection")
    g2.add_argument("--model", default="sapbert",
                    choices=["sapbert", "fast", "balanced", "accurate", "medical"],
                    help="Embedding model preset  [default: sapbert]")
    g2.add_argument("--llm", default="deepseek",
                    choices=["deepseek", "gpt4oss", "accurate", "gemini"],
                    help="LLM model preset  [default: deepseek]")
    g2.add_argument("--device", default="cpu",
                    help="Torch device ('cpu' or 'cuda')  [default: cpu]")

    g3 = parser.add_argument_group("runtime")
    g3.add_argument("--cache-dir", default="data/phenotype_matcher_v2_cache",
                    metavar="DIR", help="Embedding cache directory")
    g3.add_argument("--api-key", default=None, metavar="KEY",
                    help="OpenRouter API key (overrides OPENROUTER_API_KEY env var)")
    g3.add_argument("--debug", action="store_true",
                    help="Print progress messages to stderr")
    g3.add_argument("--match-cache-size", type=int, default=10_000, metavar="N",
                    help="Match-level cache size (0 = disabled)  [default: 10000]")


def _build_matcher(args: argparse.Namespace):
    """Construct a PhenotypeMatcher from parsed CLI arguments."""
    from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig
    cfg = MatcherConfig(
        hpo_path=args.hpo,
        mondo_path=args.mondo,
        hpoa_path=args.hpoa,
        embedding_model=args.model,
        llm_model=args.llm,
        device=args.device,
        cache_dir=args.cache_dir,
        api_key=args.api_key or os.environ.get("OPENROUTER_API_KEY"),
        debug=args.debug,
        match_cache_size=getattr(args, "match_cache_size", 10_000),
    )
    return PhenotypeMatcher(cfg)


# ---------------------------------------------------------------------------
# Output formatters
# ---------------------------------------------------------------------------

def _result_to_json_obj(text: str, out) -> dict:
    d = out.to_dict()
    d["input"] = text
    return d


def _result_to_tsv_row(text: str, out) -> str:
    hpo_present  = ";".join(out.get_hpo_ids())
    hpo_excluded = ";".join(out.get_hpo_ids(excluded=True))
    severity     = ";".join(
        f"{p.hpo_id}={p.severity_id}"
        for p in out.phenotypes
        if p.severity_id
    )
    mondo_ids   = ";".join(out.get_mondo_ids())
    omim_ids    = ";".join(out.get_omim_ids())
    orphan_ids  = ";".join(out.get_orphanet_ids())
    return "\t".join([text, hpo_present, hpo_excluded, severity,
                      mondo_ids, omim_ids, orphan_ids])


_TSV_HEADER = "\t".join([
    "input", "hpo_present", "hpo_excluded", "severity",
    "mondo_ids", "omim_ids", "orphanet_ids",
])


def _result_to_pretty(text: str, out) -> str:
    lines = [f"Input : {text!r}", ""]

    present  = out.get_present_phenotypes()
    excluded = out.get_excluded_phenotypes()
    diseases = out.diseases

    if present:
        lines.append("Present phenotypes:")
        for p in present:
            sev = f"  [{p.severity_label}]" if p.severity_label else ""
            lines.append(f"  {p.hpo_id:<14}  {p.label}{sev}")

    if excluded:
        lines.append("Excluded phenotypes:")
        for p in excluded:
            lines.append(f"  {p.hpo_id:<14}  {p.label}  [EXCLUDED]")

    if diseases:
        lines.append("Diseases:")
        for d in diseases:
            parts = []
            if d.mondo_id:
                parts.append(f"{d.mondo_id} ({d.mondo_label})")
            for oid in d.omim_phenotype_ids:
                parts.append(oid)
            for rid in d.orphanet_ids:
                parts.append(rid)
            lines.append("  " + "  ".join(parts))

    if not present and not excluded and not diseases:
        lines.append("  (no matches)")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Subcommand: match
# ---------------------------------------------------------------------------

def cmd_match(args: argparse.Namespace) -> None:
    """Match one or more text strings and print results."""
    texts: List[str] = []

    if args.input:
        src = sys.stdin if args.input == "-" else open(args.input, encoding="utf-8")
        texts = [line.rstrip("\n") for line in src if line.strip()]
        if args.input != "-":
            src.close()

    texts.extend(args.text)

    if not texts:
        print("error: provide text arguments or --input FILE", file=sys.stderr)
        sys.exit(1)

    matcher = _build_matcher(args)

    if args.format == "tsv":
        print(_TSV_HEADER)

    json_results = []

    match_fn = matcher.match_narrative if getattr(args, "narrative", False) else matcher.match

    for i, text in enumerate(texts):
        out = match_fn(text)

        # Apply output filters
        if args.present_only:
            out.phenotypes = out.get_present_phenotypes()
        if args.hpo_only:
            out.diseases = []
        if args.diseases_only:
            out.phenotypes = []

        if args.format == "json":
            json_results.append(_result_to_json_obj(text, out))
        elif args.format == "tsv":
            print(_result_to_tsv_row(text, out))
        else:  # pretty
            print(_result_to_pretty(text, out))
            if i < len(texts) - 1:
                print()

    if args.format == "json":
        if len(json_results) == 1:
            print(json.dumps(json_results[0], indent=2))
        else:
            print(json.dumps(json_results, indent=2))


# ---------------------------------------------------------------------------
# Subcommand: batch
# ---------------------------------------------------------------------------

def _apply_result(row: dict, result) -> None:
    """Write matcher output fields into a row dict in-place."""
    if result is None:
        row.update({
            "hpo_present": "", "hpo_excluded": "", "severity": "",
            "mondo_ids": "", "omim_ids": "", "orphanet_ids": "",
        })
    else:
        row["hpo_present"]  = ";".join(result.get_hpo_ids())
        row["hpo_excluded"] = ";".join(result.get_hpo_ids(excluded=True))
        row["severity"]     = ";".join(
            f"{p.hpo_id}={p.severity_id}"
            for p in result.phenotypes if p.severity_id
        )
        row["mondo_ids"]    = ";".join(result.get_mondo_ids())
        row["omim_ids"]     = ";".join(result.get_omim_ids())
        row["orphanet_ids"] = ";".join(result.get_orphanet_ids())


def cmd_batch(args: argparse.Namespace) -> None:
    """Process a TSV/CSV file, matching text in the specified column."""
    import csv
    from concurrent.futures import ThreadPoolExecutor, as_completed

    delim = "\t" if args.sep == "tab" else ","
    matcher = _build_matcher(args)
    workers = args.workers

    with open(args.input_file, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if reader.fieldnames is None:
            print("error: input file has no header row", file=sys.stderr)
            sys.exit(1)

        col = args.column
        if col not in reader.fieldnames:
            available = ", ".join(reader.fieldnames)
            print(f"error: column {col!r} not found. Available: {available}",
                  file=sys.stderr)
            sys.exit(1)

        out_fields = list(reader.fieldnames) + [
            "hpo_present", "hpo_excluded", "severity",
            "mondo_ids", "omim_ids", "orphanet_ids",
        ]
        rows = list(reader)

    out_fh = (
        open(args.output, "w", encoding="utf-8", newline="")
        if args.output else sys.stdout
    )
    writer = csv.DictWriter(
        out_fh, fieldnames=out_fields, delimiter="\t",
        extrasaction="ignore", lineterminator="\n",
    )
    writer.writeheader()

    total = len(rows)
    match_fn = matcher.match_narrative if getattr(args, "narrative", False) else matcher.match

    if workers > 1:
        # Parallel processing: submit all non-empty rows to the thread pool,
        # collect results keyed by row index, then write in original order.
        results_map: dict = {}

        def _match_row(idx: int) -> tuple:
            row = rows[idx]
            text = row.get(col, "").strip()
            if not text:
                return idx, None
            if args.debug:
                print(f"[submitting {idx + 1}/{total}] {text[:70]}", file=sys.stderr)
            return idx, match_fn(text)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_to_idx = {pool.submit(_match_row, i): i for i in range(total)}
            done = 0
            for future in as_completed(future_to_idx):
                idx, result = future.result()
                results_map[idx] = result
                done += 1
                if args.debug or (done % 10 == 0):
                    print(f"[{done}/{total} done]", file=sys.stderr)

        for i, row in enumerate(rows):
            _apply_result(row, results_map[i])
            writer.writerow(row)

    else:
        # Sequential processing (original behaviour)
        for i, row in enumerate(rows, 1):
            text = row.get(col, "").strip()

            if not text:
                _apply_result(row, None)
                writer.writerow(row)
                continue

            if args.debug or (i % 10 == 0):
                print(f"[{i}/{total}] {text[:70]}", file=sys.stderr)

            result = match_fn(text)
            _apply_result(row, result)
            writer.writerow(row)

    if args.output:
        out_fh.close()

    print(f"\nDone. {total} rows processed.", file=sys.stderr)


# ---------------------------------------------------------------------------
# Subcommand: info
# ---------------------------------------------------------------------------

def cmd_info(args: argparse.Namespace) -> None:
    """Look up an ontology term by ID or label string."""
    from phenotype_matcher_v2.ontology_index import OntologyIndex

    idx = OntologyIndex(
        hpo_path=args.hpo,
        mondo_path=args.mondo,
        hpoa_path=args.hpoa,
        cache_dir=args.cache_dir,
        debug=args.debug,
    )

    query = args.query.strip()
    known_prefixes = ("HP:", "MONDO:", "OMIM:", "ORPHA:")

    if any(query.upper().startswith(p) for p in known_prefixes):
        # Ontology ID lookup
        term_id = query
        label = idx.primary_label.get(term_id)
        if label is None:
            print(f"Term not found: {term_id}", file=sys.stderr)
            sys.exit(1)

        print(f"ID     : {term_id}")
        print(f"Label  : {label}")

        synonyms = sorted(
            lbl for lbl, ids in idx.exact_map.items()
            if term_id in ids and lbl != label.lower()
        )
        if synonyms:
            print(f"Synonyms ({len(synonyms)}):")
            for s in synonyms:
                print(f"  {s}")

        if term_id.startswith("HP:"):
            if term_id in idx.phenotype_ids:
                ns = "HPO phenotype (descendant of HP:0000118)"
            elif term_id in idx.modifier_ids:
                ns = "HPO modifier / severity (descendant of HP:0012824)"
            else:
                ns = "HPO (other)"
            print(f"Namespace: {ns}")
        elif term_id.startswith("MONDO:"):
            print("Namespace: MONDO")
            orpha = idx.mondo_to_orphanet.get(term_id, [])
            if orpha:
                print(f"Orphanet xrefs: {', '.join(orpha)}")
        elif term_id.startswith("OMIM:"):
            print("Namespace: OMIM")
        elif term_id.startswith("ORPHA:"):
            print("Namespace: Orphanet")

    else:
        # Label / free-text lookup
        query_lower = query.lower()
        exact = idx.exact_map.get(query_lower, [])

        if exact:
            print(f"Exact match for {query!r}:")
            for term_id in exact:
                primary = idx.primary_label.get(term_id, term_id)
                print(f"  {term_id:<22}  {primary}")
        else:
            # Partial match
            partial = sorted(
                (lbl, ids)
                for lbl, ids in idx.exact_map.items()
                if query_lower in lbl
            )
            if not partial:
                print(f"No terms found for: {query!r}", file=sys.stderr)
                sys.exit(1)
            print(f"No exact match for {query!r}. Partial matches (up to 20):")
            for lbl, ids in partial[:20]:
                id_str = ", ".join(ids)
                print(f"  {lbl!r:<40}  {id_str}")


# ---------------------------------------------------------------------------
# Subcommand: check
# ---------------------------------------------------------------------------

def cmd_check(args: argparse.Namespace) -> None:
    """Verify configuration, data files, and installed packages."""
    issues = []

    def ok(msg: str) -> None:
        print(f"  [OK]  {msg}")

    def warn(msg: str) -> None:
        print(f"  [--]  {msg}")
        issues.append(msg)

    print("Ontology files:")
    for label, path in [
        ("hp.obo   (HPO)", args.hpo),
        ("mondo.obo (MONDO)", args.mondo),
        ("phenotype.hpoa (OMIM/Orphanet)", args.hpoa),
    ]:
        if os.path.isfile(path):
            size_kb = os.path.getsize(path) // 1024
            ok(f"{label}: {path}  ({size_kb:,} KB)")
        else:
            warn(f"{label}: NOT FOUND at {path}")

    print("\nAPI key:")
    key = args.api_key or os.environ.get("OPENROUTER_API_KEY", "")
    if key:
        ok(f"OPENROUTER_API_KEY set  ({key[:8]}...)")
    else:
        warn("OPENROUTER_API_KEY not set — LLM steps will be skipped")

    print("\nEmbedding cache:")
    from phenotype_matcher_v2.matcher import EMBEDDING_MODELS
    model_name = EMBEDDING_MODELS.get(args.model, args.model)
    cache_file = os.path.join(
        args.cache_dir,
        f"embeddings_{model_name.replace('/', '_')}.pkl",
    )
    if os.path.isfile(cache_file):
        size_mb = os.path.getsize(cache_file) // (1024 * 1024)
        ok(f"{args.model} ({model_name}): cached at {cache_file}  ({size_mb:,} MB)")
    else:
        warn(f"{args.model} ({model_name}): not yet built — will build on first 'match' call")

    print("\nPython packages:")
    packages = [
        ("pronto",               "pronto"),
        ("ahocorasick",          "pyahocorasick"),
        ("rapidfuzz",            "rapidfuzz"),
        ("spacy",                "spacy"),
        ("sentence_transformers","sentence-transformers"),
        ("numpy",                "numpy"),
        ("requests",             "requests"),
        ("torch",                "torch"),
    ]
    for mod, label in packages:
        try:
            pkg = __import__(mod)
            ver = getattr(pkg, "__version__", "?")
            ok(f"{label}  ({ver})")
        except ImportError:
            warn(f"{label}: NOT INSTALLED")

    print()
    if issues:
        print(f"Found {len(issues)} issue(s):")
        for issue in issues:
            print(f"  - {issue}")
        sys.exit(1)
    else:
        print("Configuration OK.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="phenotype-matcher",
        description=(
            "Map clinical phenotype text to HPO / MONDO / OMIM / Orphanet IDs.\n\n"
            "Subcommands:\n"
            "  match   Match one or more text strings\n"
            "  batch   Process a TSV/CSV file column by column\n"
            "  info    Look up an ontology term by ID or label\n"
            "  check   Verify configuration, data files, and packages"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ------------------------------------------------------------------ match
    p_match = subparsers.add_parser(
        "match",
        help="Match one or more text strings",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Match clinical phenotype text to ontology IDs.\n\n"
            "Examples:\n"
            "  phenotype-matcher match \"seizures, hypotonia\"\n"
            "  phenotype-matcher match \"no fever\" \"severe ataxia\" --format pretty\n"
            "  phenotype-matcher match --input terms.txt --format tsv > results.tsv\n"
            "  echo \"intellectual disability\" | phenotype-matcher match --input -"
        ),
    )
    p_match.add_argument("text", nargs="*", help="Text string(s) to match")
    p_match.add_argument("--input", "-i", metavar="FILE",
                         help="Read one text per line from FILE (use - for stdin)")
    p_match.add_argument("--format", "-f",
                         choices=["json", "tsv", "pretty"], default="json",
                         help="Output format  [default: json]")
    p_match.add_argument("--present-only", action="store_true",
                         help="Suppress excluded (negated) phenotypes from output")
    p_match.add_argument("--hpo-only", action="store_true",
                         help="Suppress disease matches from output")
    p_match.add_argument("--diseases-only", action="store_true",
                         help="Suppress phenotype matches from output")
    p_match.add_argument("--narrative", action="store_true",
                         help="Use NER-based narrative mode for long clinical text")
    _add_common_args(p_match)
    p_match.set_defaults(func=cmd_match)

    # ------------------------------------------------------------------ batch
    p_batch = subparsers.add_parser(
        "batch",
        help="Process a TSV/CSV file, one row per record",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Process a delimited file, matching text in a named column.\n"
            "All original columns are preserved; six new columns are appended:\n"
            "  hpo_present, hpo_excluded, severity, mondo_ids, omim_ids, orphanet_ids\n\n"
            "Severity entries are formatted as HPO_ID=SEVERITY_ID pairs.\n\n"
            "Examples:\n"
            "  phenotype-matcher batch patients.tsv --column phenotypes --output results.tsv\n"
            "  phenotype-matcher batch data.csv --sep comma --column description --output out.tsv"
        ),
    )
    p_batch.add_argument("input_file", metavar="INPUT",
                         help="Input TSV or CSV file")
    p_batch.add_argument("--column", "-c", required=True,
                         help="Name of the column containing phenotype text")
    p_batch.add_argument("--output", "-o", metavar="FILE",
                         help="Output TSV file  [default: stdout]")
    p_batch.add_argument("--sep", choices=["tab", "comma"], default="tab",
                         help="Input file delimiter  [default: tab]")
    p_batch.add_argument("--workers", type=int, default=1, metavar="N",
                         help="Parallel worker threads for LLM I/O  [default: 1]")
    p_batch.add_argument("--narrative", action="store_true",
                         help="Use NER-based narrative mode for long clinical text")
    _add_common_args(p_batch)
    p_batch.set_defaults(func=cmd_batch)

    # ------------------------------------------------------------------- info
    p_info = subparsers.add_parser(
        "info",
        help="Look up an ontology term by ID or label",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Look up a term by ontology ID or free-text label.\n"
            "Shows canonical label, all synonyms, namespace, and cross-references.\n\n"
            "Examples:\n"
            "  phenotype-matcher info HP:0001250\n"
            "  phenotype-matcher info MONDO:0100079\n"
            "  phenotype-matcher info OMIM:261600\n"
            "  phenotype-matcher info \"seizure\"\n"
            "  phenotype-matcher info \"dravet\"       # partial label search"
        ),
    )
    p_info.add_argument("query", help="Term ID (e.g. HP:0001250) or label string")
    _add_common_args(p_info)
    p_info.set_defaults(func=cmd_info)

    # ------------------------------------------------------------------ check
    p_check = subparsers.add_parser(
        "check",
        help="Verify configuration, data files, and Python packages",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Verify that all required files and packages are present.\n\n"
            "Examples:\n"
            "  phenotype-matcher check\n"
            "  phenotype-matcher check --hpo /data/hp.obo --model fast"
        ),
    )
    _add_common_args(p_check)
    p_check.set_defaults(func=cmd_check)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
