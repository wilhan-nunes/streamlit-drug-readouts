# Standalone web app changelog

This log covers `drug_readouts.html`, not the original Streamlit application.

## Unreleased — under development

### Added

- Self-contained browser app with embedded drug ontology, local table uploads, and optional public GNPS2 task retrieval.
- Welcome guidance, collapsible horizontal inputs, editable annotations, result tables, and TSV downloads.
- Top pharmacological classes chart with percentages and a sample-count axis.
- UpSet-style co-occurrence plot with SVG download, plus the classification overview visualization.
- Credits and citations, including the 2026 MetaboApps paper.
- Offline-use guidance and a visible development-status notice.
- A deterministic build script that embeds the current drug and analog CSV data from `data/` into the standalone HTML, with a stale-output check and an editable HTML source.
- Optional checksum-verified downloads or explicit build inputs for the drug and analog metadata CSVs from the GNPS Drug Library Zenodo record; the welcome link now points to that record.
- An in-app **Fetch latest metadata** option that reads Zenodo's latest published dataset through its CORS-enabled API and caches successful refreshes in the browser for later offline use, with a control to clear the cache and return to embedded data.

### Notes

- Offline analysis requires previously downloaded FBMN quantification and annotation tables. Task-ID retrieval and external links require internet access.
- The standalone browser port remains under development; validate important results against source data and the original workflow.
