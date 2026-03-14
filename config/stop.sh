#!/bin/bash

# PAVS Website Stop Script
# Usage: ./stop.sh

COMPOSE_FILE="docker-compose-sparql.yml"

echo "=== Stopping PAVS services with Docker Compose ==="
docker compose -f "$COMPOSE_FILE" down

echo "=== PAVS is stopped ==="
