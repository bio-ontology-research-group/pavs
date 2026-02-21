# Stage 1: Build React Frontend
FROM node:20-slim AS frontend-builder
WORKDIR /app/frontend
COPY frontend/package*.json ./
RUN npm install
COPY frontend/ ./
RUN npm run build

# Stage 2: Backend & Final Image
FROM python:3.10-slim
WORKDIR /app

# Install system dependencies
# We keep this extremely minimal to avoid filling up the small /var partition.
RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install uv for faster dependency management
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uv/bin/uv

# Install Python dependencies
COPY requirements.txt ./
# Surgical install to avoid building psycopg2 from source
RUN /uv/bin/uv pip install --system psycopg2-binary && \
    /uv/bin/uv pip install --system hgvs --no-deps && \
    /uv/bin/uv pip install --system -r requirements.txt

# Copy backend code
COPY backend/ ./backend/
COPY ontology/ ./ontology/
COPY data/ ./data/
COPY phenopackets/ ./phenopackets/

# Copy frontend build from Stage 1
COPY --from=frontend-builder /app/frontend/build ./frontend/build

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV HOST=0.0.0.0
ENV PORT=8000

# Expose port
EXPOSE 8000

# Start command
CMD ["python", "-m", "backend.main"]
