#!/bin/bash

# PAVS Website Start Script
# Usage: ./start.sh [OPTIONS]

COMPOSE_FILE="docker-compose-sparql.yml"
INPUT_FILE="data/combined_annotated-gpt4oss.tsv"
RELOAD=false

show_help() {
    echo "PAVS Website Management Script"
    echo "Usage: ./start.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --reload             Force regeneration of RDF data from the default TSV."
    echo "  --input, -i FILE     Specify a TSV file to reload from (implies --reload)."
    echo "  --help, -h           Show this help message."
    echo ""
    echo "Examples:"
    echo "  ./start.sh                    # Start services with existing data"
    echo "  ./start.sh --reload           # Refresh data from default TSV and start"
    echo "  ./start.sh -i data/my.tsv     # Refresh data from specific file and start"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --reload)
      RELOAD=true
      shift
      ;;
    --input|-i)
      INPUT_FILE="$2"
      RELOAD=true
      shift 2
      ;;
    --help|-h)
      show_help
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      show_help
      exit 1
      ;;
  esac
done

if [ "$RELOAD" = true ]; then
  if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
  fi
  echo "=== Reloading data (TSV -> RDF) from: $INPUT_FILE ==="
  # Run the full preparation pipeline
  # We use --skip-load because the Docker 'loader' service handles Virtuoso ingestion
  uv run python intake/prepare_all.py --force --skip-load --input "$INPUT_FILE"
fi

echo "=== Starting PAVS services with Docker Compose ==="
# Start virtuoso first to ensure health checks pass
docker compose -f "$COMPOSE_FILE" up -d virtuoso

echo "=== Waiting for Virtuoso to be healthy ==="
# The loader, backend, and frontend services will start based on dependencies in the compose file
if [ "$RELOAD" = true ]; then
    # Force loader to run to pick up the new TTL files
    docker compose -f "$COMPOSE_FILE" up -d loader backend frontend
else
    docker compose -f "$COMPOSE_FILE" up -d backend frontend
fi

echo "=== PAVS is starting ==="
echo "Frontend:       http://localhost:3000"
echo "Backend API:    http://localhost:8000"
echo "Virtuoso:       http://localhost:8890/sparql"
echo "========================="
echo "Use ./stop.sh to shut down."
