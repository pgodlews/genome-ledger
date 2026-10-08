# Genome Ledger in a container. Supported platform: Linux x86-64 only — several of the
# tools `setup` downloads ship x86-64 Linux binaries and nothing else.
#
# The image holds the system tools and the Python environment. Reference data, annotation
# caches and the downloaded tools (GATK, VEP, PharmCAT, …) are NOT in it: `setup` fetches
# them into the mounted refs directory, exactly as it does on a normal install.
FROM ghcr.io/astral-sh/uv:python3.12-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
      bcftools samtools tabix bwa fastp default-jre-headless perl git build-essential \
      unzip curl ca-certificates zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
      libssl-dev libdbi-perl libarchive-zip-perl libwww-perl libjson-perl liblist-moreutils-perl cpanminus \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
# The environment is built here, as root, and only read at run time: the container runs as
# your own user id (see compose.yaml), which could not write it and does not need to.
ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    UV_NO_SYNC=1 \
    PATH=/opt/venv/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    GENOMES_ROOT=/data \
    HOME=/home/app \
    UV_CACHE_DIR=/data/tools/.uv-cache
COPY pyproject.toml uv.lock ./
RUN UV_NO_SYNC=0 uv sync --frozen --no-dev
COPY . .
# A home directory any user id can use (git, java and uv all want one).
RUN mkdir -p /home/app && chmod 1777 /home/app

ENTRYPOINT ["python", "run.py"]
CMD ["--help"]
