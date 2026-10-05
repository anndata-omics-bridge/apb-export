# The Python tools that read apb-export's files: msmu, ProteoPy and AlphaPeptTools, pinned to the
# export rules' target_version_pattern. prolfquapp is checked in its own published image, not
# here; apb-export itself depends on none of them.
FROM ghcr.io/astral-sh/uv:0.9.8-python3.12-bookworm-slim

# msmu needs inmoose, which compiles from source where no wheel exists, as on linux/arm64.
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential \
 && rm -rf /var/lib/apt/lists/*

ENV VIRTUAL_ENV=/opt/consumers \
    PATH=/opt/consumers/bin:$PATH
RUN uv venv /opt/consumers \
 && uv pip install "msmu>=0.4,<0.5" "proteopy>=0.1,<0.2" "alphapepttools>=0.4,<0.5"

COPY check_msmu.py check_proteopy.py check_alphapepttools.py /checks/
