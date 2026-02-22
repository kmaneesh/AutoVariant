# AutoVariant

Search a VCF by gene symbol or location and get the 7-step variant message (WES Variant Details SOP). No AI; deterministic only.

## CLI

```bash
source .venv/bin/activate
python scripts/auto_variant.py FAH --vcf path/to/file.vcf
python scripts/auto_variant.py 80472479 --vcf path/to/file.vcf.gz
```

## Web app

1. Put VCF files (`.vcf` or `.vcf.gz`) in the **`data/`** folder. The app scans `data/` by default when run locally.

2. From the project root, run:

   ```bash
   source .venv/bin/activate
   uvicorn app.main:app --reload --port 9000
   ```

3. Open http://127.0.0.1:9000 — enter a gene (e.g. `FAH`) or location (e.g. `80472479`), choose a VCF from the dropdown, and click Search. The 7-step message is shown in the result box.

## Docker

Run the app on port 9000 with Docker Compose. The **`./data`** folder on the host is always mounted at **`/data`** in the container; the app is configured to use `/data` there.

```bash
# Put VCF files in ./data, then:
docker compose up --build
```

Open http://127.0.0.1:9000. No `.env` needed; the compose file sets `VCF_DIR=/data` and mounts `./data:/data`.

The app reads the selected VCF from the server path (no upload). `.vcf` and `.vcf.gz` are supported (via cyvcf2 when available).
