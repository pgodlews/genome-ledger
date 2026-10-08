# Running with Docker Compose

A container image with the system tools and the Python environment, and a compose file
that keeps your data in three host directories. **Supported on Linux x86-64 only**: several
of the tools `setup` downloads ship x86-64 Linux binaries and nothing else. On a Mac, use
the native install in the README.

## Try the demo

```bash
docker compose up demo
```

Open http://127.0.0.1:8765/. This builds the synthetic family's reports inside the
container and serves them; nothing is written to your disk and nothing else needs
installing. Stop it with Ctrl-C.

## The three directories

| Variable | Default | Holds | Treat it as |
|---|---|---|---|
| `GENOMES_DIR` | `./data/genomes` | raw reads, CRAMs, calls, the database, snapshots, `family/family.tsv` | Irreplaceable. Back it up. Never expose it. |
| `REFS_DIR` | `./data/refs` | reference genome and index, annotation cache, downloaded tools | Re-downloadable. Tens of GB. No backup needed. |
| `REPORTS_DIR` | `./data/reports` | `md/`, `json/`, `html/` (see [report-formats.md](report-formats.md)) | The only part you read. Point other tools here. |

Inside the container they appear as one tree under `/data`, which is why every `run.py`
command works unchanged.

## First-time setup

Create the directories yourself, and tell Compose to run as you. Docker creates a missing
mount directory as root, and the pipeline — which runs as your user — then cannot write to
it, so do this before the first command:

```bash
mkdir -p data/genomes/{refs,vep_cache,tools,reports} \
         data/refs/{refs,vep_cache,tools} data/reports/html
printf 'PUID=%s\nPGID=%s\n' "$(id -u)" "$(id -g)" > .env
docker compose build pipeline
docker compose run --rm pipeline setup     # ~40 GB into REFS_DIR; 1-3 h depending on CPU and bandwidth; once
```

(The four empty folders inside `data/genomes` are only mount points for the other two
directories. `data/reports/html` is one too: `serve` mounts it on its own, so if `serve` is
started before any report exists, Docker creates it as root and every later report render
fails with "Permission denied".) To keep the data elsewhere, add `GENOMES_DIR=…`,
`REFS_DIR=…` and `REPORTS_DIR=…` to `.env` and create those paths instead, including
`$REPORTS_DIR/html`.

## Everyday use

Any `run.py` command runs through the `pipeline` service; paths are as the container sees
them, so put input files under `data/genomes/` and refer to them as `/data/…`:

```bash
docker compose run --rm pipeline ingest /data/incoming/maria.vcf.gz --sample Maria --relation mother --sex female
docker compose run --rm pipeline normalize Maria
docker compose run --rm pipeline scan
```

Read the reports:

```bash
docker compose up -d serve        # http://127.0.0.1:8765/
```

The `serve` container mounts `reports/html` only, read-only. It cannot see a genome, the
database or even the markdown reports — but the pages it serves are themselves sensitive:
names, family relationships and genetic results. The server has no authentication, so it
listens on this machine only. If others should reach it, put an authenticating reverse proxy
in front, and give access only to the people whose results it shows.

Check for new database releases once a day, and scan when there is one:

```bash
docker compose up -d scheduler
```

Settings go in `.env`: `NOTIFY_CMD` (a command run with the scan's summary), `PRS_S1_GATE`
(`0` to show incidental findings to everyone), `SERVE_PORT`.

## What does not work in a container

- **`sv` and `hla-type`.** Those stages start Docker containers of their own (Manta,
  arcasHLA). From inside a container that needs the host's Docker socket, which is
  root-equivalent on the host, so it is not wired up. Run those two stages from a native
  install against the same data directory.
- **`interpret`.** It sends jobs to a GPU box over SSH and needs your SSH key and host
  configuration in the container; not set up by default.
- **Apple Silicon and other non-x86-64 hosts.**

## How the image is built

`Dockerfile` installs the system packages listed in the README, then the locked Python
environment. The tools `setup` downloads (GATK, VEP, PharmCAT, ExpansionHunter, plink2, …)
are not baked in: they land in `REFS_DIR/tools`, exactly as on a native install, so the
image stays small and a tool upgrade does not need a rebuild. The trade-off is that the
image alone does not pin those versions; `pipeline/config.py` does.
