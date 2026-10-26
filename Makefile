VENV_BIN := .venv/bin

.DEFAULT_GOAL := help
.PHONY: help sync format format-check lint typecheck deps test docs build schema check check-full clean

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

check:  ## Run every merge-blocking quality gate
	$(MAKE) format-check lint typecheck deps test docs build

check-full: check  ## Run every pre-push quality gate

clean:  ## Remove generated build and quality artifacts
	$(VENV_BIN)/python -c "import shutil; [shutil.rmtree(path, ignore_errors=True) for path in ('build', 'dist', 'public', '.pytest_cache', '.ruff_cache')]"
