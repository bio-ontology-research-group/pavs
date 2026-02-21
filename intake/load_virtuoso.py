#!/usr/bin/env python3
"""
load_virtuoso.py — Bulk-load TTL files into a Virtuoso triple store.

Uses /sparql-auth with HTTP Digest authentication (dba user) to send
CLEAR and LOAD statements. Virtuoso requires Digest, not Basic auth.

Usage:
    uv run python scripts/load_virtuoso.py \
        --endpoint http://localhost:8890 \
        --input-dir rdf_output/ \
        [--virtuoso-import-dir /rdf_import]
"""

import argparse
import json
import time
import urllib.parse
import urllib.request
import urllib.error
from pathlib import Path

GRAPH_MAP = {
    "cases.ttl":      "http://pavs.kaust.edu.sa/graph/cases",
    "genes.ttl":      "http://pavs.kaust.edu.sa/graph/genes",
    "hpoa.ttl":       "http://pavs.kaust.edu.sa/graph/hpoa",
    "hpo_ic.ttl":     "http://pavs.kaust.edu.sa/graph/hpo-ic",
    "literature.ttl": "http://pavs.kaust.edu.sa/graph/literature",
    "ddd.ttl":        "http://pavs.kaust.edu.sa/graph/ddd",
    "clinvar.ttl":    "http://pavs.kaust.edu.sa/graph/clinvar",
}


def make_opener(user: str, password: str, endpoint: str) -> urllib.request.OpenerDirector:
    """Build an opener with HTTP Digest auth (required by Virtuoso /sparql-auth)."""
    passwd_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    passwd_mgr.add_password(None, endpoint, user, password)
    handler = urllib.request.HTTPDigestAuthHandler(passwd_mgr)
    return urllib.request.build_opener(handler)


def sparql_update(opener, auth_url: str, query: str) -> bool:
    data = urllib.parse.urlencode({"update": query}).encode("utf-8")
    req = urllib.request.Request(
        auth_url, data=data,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    try:
        with opener.open(req, timeout=600) as resp:
            resp.read()
        return True
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")[:300]
        print(f"  HTTP {e.code}: {body}")
        return False
    except Exception as e:
        print(f"  Error: {e}")
        return False


def sparql_count(sparql_url: str, graph: str) -> int:
    query = f"SELECT (COUNT(*) AS ?n) WHERE {{ GRAPH <{graph}> {{ ?s ?p ?o }} }}"
    url = sparql_url + "?" + urllib.parse.urlencode({
        "query": query, "format": "application/sparql-results+json"
    })
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            result = json.loads(resp.read())
            return int(result["results"]["bindings"][0]["n"]["value"])
    except Exception:
        return -1


def load_file(opener, auth_url: str, sparql_url: str,
              fname: str, graph: str, virtuoso_dir: str) -> bool:
    file_uri = f"file://{virtuoso_dir.rstrip('/')}/{fname}"
    print(f"  {fname} → <{graph}>")
    print(f"    CLEAR …", end=" ", flush=True)
    sparql_update(opener, auth_url, f"CLEAR SILENT GRAPH <{graph}>")
    print(f"LOAD …", end=" ", flush=True)
    ok = sparql_update(opener, auth_url, f"LOAD <{file_uri}> INTO GRAPH <{graph}>")
    if ok:
        n = sparql_count(sparql_url, graph)
        print(f"OK  ({n:,} triples)")
    else:
        print("FAILED")
    return ok


def wait_for_virtuoso(endpoint: str, max_seconds: int = 120):
    print(f"  Waiting for Virtuoso at {endpoint} …", flush=True)
    deadline = time.time() + max_seconds
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(f"{endpoint}/sparql", timeout=5):
                print("  Virtuoso is up!", flush=True)
                return
        except Exception:
            time.sleep(2)
    print(f"  ERROR: Virtuoso not available after {max_seconds}s")
    raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description="Load TTL files into Virtuoso")
    parser.add_argument("--endpoint", default="http://localhost:8890")
    parser.add_argument("--input-dir", default="rdf_output/")
    parser.add_argument("--virtuoso-import-dir", default="/rdf_import",
                        help="Path to TTL files as seen inside the Virtuoso container")
    parser.add_argument("--user", default="dba")
    parser.add_argument("--password", default="pavs_dba")
    parser.add_argument("--wait", action="store_true")
    args = parser.parse_args()

    if args.wait:
        wait_for_virtuoso(args.endpoint)

    sparql_url  = f"{args.endpoint}/sparql"
    auth_url    = f"{args.endpoint}/sparql-auth"
    opener      = make_opener(args.user, args.password, args.endpoint)
    input_dir   = Path(args.input_dir)

    print(f"\nLoading TTL from {input_dir}  (Virtuoso sees: {args.virtuoso_import_dir})")
    print(f"Endpoint: {auth_url}  (Digest auth as {args.user})")
    print("=" * 60)

    failed = []
    for fname, graph in GRAPH_MAP.items():
        if not (input_dir / fname).exists():
            print(f"  {fname} — not found, skipping")
            continue
        ok = load_file(opener, auth_url, sparql_url, fname, graph, args.virtuoso_import_dir)
        if not ok:
            failed.append(fname)

    print()
    if failed:
        print(f"WARNING: {len(failed)} file(s) failed: {failed}")
        raise SystemExit(1)
    else:
        print("All files loaded successfully.")


if __name__ == "__main__":
    main()
