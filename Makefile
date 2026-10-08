VENV_BIN := .venv/bin
STUDIO ?= ../apb_studio
MODULES ?= ../apb-proteobench/src/apb_proteobench/data/modules
WEB_EXAMPLES ?= build/web-examples
IMAGE ?= ghcr.io/anndata-omics-bridge/apb-export
IMAGE_TAG ?= local

.DEFAULT_GOAL := help
.PHONY: help sync format format-check lint typecheck deps test docs build schema web web-examples image check clean

help:  ## Show developer commands
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

sync:  ## Synchronize the locked development environment
	uv sync --group dev --group docs

format:  ## Format and autofix source and tests
	$(VENV_BIN)/ruff format src tests
	$(VENV_BIN)/ruff check --fix src tests

format-check:  ## Check formatting without changing files
	$(VENV_BIN)/ruff format --check src tests

lint:  ## Run Ruff lint checks
	$(VENV_BIN)/ruff check src tests

typecheck:  ## Run standard Pyright in strict mode
	$(VENV_BIN)/pyright

deps:  ## Validate dependency declarations
	$(VENV_BIN)/deptry .

test:  ## Run tests with branch coverage
	$(VENV_BIN)/pytest --cov --cov-branch

docs:  ## Build the documentation with warnings as errors
	uv run --group docs zensical build --clean --strict

build:  ## Build and validate source and wheel distributions
	uv build
	$(VENV_BIN)/twine check dist/*

schema:  ## Rewrite the export-rule JSON Schemas from the models
	$(VENV_BIN)/python -m apb_export.export_rules.schema_artifact

web:  ## Test, type-check and rebuild the web page bundle into the package
	cd web && npm ci && npm test && npm run build

web-examples:  ## Collect the routine corpus, scripts/web_examples.csv and module SDRFs into one folder the server mounts
	$(VENV_BIN)/python scripts/web_examples.py $(STUDIO)/corpuses/routine.csv $(STUDIO)/test_data_download $(MODULES) $(WEB_EXAMPLES) --extra scripts/web_examples.csv --workspace .. --fastas $(STUDIO)/workflow_tables/workflow_proteobench.csv

image:  ## Build the web app's image from this checkout and its sibling packages
	docker build -f Dockerfile -t $(IMAGE):$(IMAGE_TAG) ..

check:  ## Run every merge-blocking quality gate
	$(MAKE) format-check lint typecheck deps test docs build

clean:  ## Remove generated build and quality artifacts
	$(VENV_BIN)/python -c "import shutil; [shutil.rmtree(path, ignore_errors=True) for path in ('build', 'dist', 'public', '.pytest_cache', '.ruff_cache')]"
