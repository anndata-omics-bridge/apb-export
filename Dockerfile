# apb-export-web, the upload-and-convert page, with its anndata_bridge siblings built from
# source. The build context is the folder holding the checkouts side by side, as the
# anndata_bridge workspace and the publish workflow lay them out:
#   docker build -f apb-export/Dockerfile -t apb-export ..
FROM ghcr.io/astral-sh/uv:0.9.8-python3.13-trixie-slim

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never \
    UV_PROJECT_ENVIRONMENT=/app/.venv \
    PATH=/app/.venv/bin:$PATH

COPY apb2 /src/apb2
COPY apb-catalog /src/apb-catalog
COPY apb-fasta /src/apb-fasta
COPY protein_fasta /src/protein_fasta
COPY prozor /src/prozor
COPY apb-export /src/apb-export
WORKDIR /src/apb-export
RUN --mount=type=cache,target=/root/.cache/uv uv sync --no-dev --extra web --no-editable

WORKDIR /app
# The store holds options, examples and job folders. Uploads stream through TMPDIR, so the
# deployment points it at disk, not at a memory-backed /tmp.
ENTRYPOINT ["apb-export-web", "--host", "0.0.0.0", "--port", "8000", "--store", "/data/store"]
EXPOSE 8000
