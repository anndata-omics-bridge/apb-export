# Deployment

The public instance is [apps-dev.bfabric.org/apb-export](https://apps-dev.bfabric.org/apb-export/), run by FGCZ behind its app portal. A release is an image tag; the examples are a data folder the container mounts.

## Release the image

1. Commit, then tag: `git tag vX.Y.Z && git push origin vX.Y.Z`
2. The `publish` workflow builds linux/arm64 and linux/amd64 from this repository and its sibling packages at their `main`, and pushes `ghcr.io/anndata-omics-bridge/apb-export:vX.Y.Z`; it takes about five minutes
3. For a local test build: `make image`, with the siblings checked out beside this repository

## Stage the examples

```bash
make web-examples WEB_EXAMPLES=build/web-examples-vX.Y.Z
```

The folder holds `examples.csv`, the datasets under `data/` and the module SDRFs under `sdrf/`, about 1.8 GB: APB Studio's routine corpus plus the rows of `scripts/web_examples.csv`. Stage into a new folder each time, since staging only adds files. Copy it to the host with `rsync -rlt --delete`.

## Run the container

| Item | Value |
| --- | --- |
| Entrypoint | `apb-export-web --host 0.0.0.0 --port 8000 --store /data/store` |
| Arguments | `--max-upload-mb 2048 --ttl-hours 24 --timeout-minutes 60 --examples /examples/examples.csv --examples-root /examples/data --examples-sdrf /examples/sdrf` |
| `/data/store` | writable: options, examples index, previews, job folders |
| `/data/tmp` | writable, set as `TMPDIR`; uploads stream through it, so it must be disk |
| `/examples` | the staged folder, read-only |
| Health | `GET /api/options.json` |

Behind a path prefix the proxy strips the prefix and redirects `/<prefix>` to `/<prefix>/`; the page builds only relative URLs, so the app needs no root path.

## The FGCZ instance

Its descriptor is `portal/apb-export` in FGCZ's web-apps repository, whose README holds the host recipe. A new release is a new `IMAGE_TAG` in that folder's `.env`, then `make deploy` on the host; new examples are a new staged folder copied over, then `make deploy`.
