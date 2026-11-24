# OpenFE Zenodo reference freeze

This directory freezes complete API metadata for the 38 records returned by the
exact Zenodo query `"OpenFE 2024 Industry Benchmarking Project"` on 2026-09-22.
No raw result archives or trajectory payloads are downloaded by default. Several
records exceed multiple GB and must be acquired only after selecting the exact
public Ross category/system records needed for a table.

`manifest.tsv` records every DOI, file name, file size, and API checksum.
`records.json` preserves the complete searchable API response per record. The
raw result payload is intentionally not treated as part of Gate 0.
