FROM python:3.12-slim

WORKDIR /app

RUN pip install --no-cache-dir fastapi uvicorn pydantic pandas

# Copy application code
COPY backend_sparql/ ./backend_sparql/
COPY scripts/ ./scripts/

EXPOSE 8000

CMD ["uvicorn", "backend_sparql.main:app", "--host", "0.0.0.0", "--port", "8000"]
