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
- A **Feature table + query MGF** input runs the library search step on a local feature table and query MGF, with no GNPS2 task ID. It shares the **Library search** panel, Zenodo reference retrieval and engine hook with Everything Bagel tasks; the engine receives `taskId: null`.
- Download progress for the EB feature MGF and the Zenodo reference MGFs: the status bar shows MB received, percent, speed and time left, with a progress bar beneath it. Cached reference files are reported as such.
- One tab per input source (**GNPS2 task**, **FBMN tables**, **Feature table + query MGF**) replaces the source drop-down. Each tab opens with a short how-to, and keyboard arrows move between tabs. **About the input files** is kept and lists each source separately.
- Drop boxes turn green with a ✓ once they hold a file.
- A **Files** box in the Library search panel shows whether the query spectra and each reference MGF are in memory, cached in this browser, read from disk, or still to be downloaded. **Check files** compares the cache with Zenodo's latest version and gives the download size.
- The Streamlit app's link parameters `task_id`, `threshold` and `blank_ids` prefill the inputs. As in Streamlit, nothing runs until **Run analysis** is clicked, and an out-of-range threshold is ignored with a warning.

### Changed

- Peak-area threshold and blank IDs moved from the GNPS2 box into a shared **Analysis settings** section, since they apply to every source. **Load example** moved into the GNPS2 task tab, and **Run analysis** moved below the settings.
- The citation now points to the published GNPS Drug Library paper (*Nature Communications* 16, 10600, 2025; doi:10.1038/s41467-025-65993-5) instead of the bioRxiv preprint. The Streamlit app's welcome page was updated to match.
- Less text on the page: each tab's explanation is now a **How to use this tab** hover/focus popup (Escape closes it; on phones it opens as a bottom sheet). The drug metadata card is gone: **Fetch latest metadata** sits beside **Run analysis**, **Use embedded metadata** (formerly **Use embedded & clear cache**) appears only while Zenodo metadata is in use, the active version is in the Fetch button's tooltip, and the card's explanation is an ⓘ popup.
- The welcome card and the development notice span the full page width, matching the other cards.
- EB-only names in the search-preparation UI and code became generic (`Prepare library search files`, `prepareSearch`, `runSearch`); `window.DrugReadoutsEBEngine` is unchanged. Prepare and Run are disabled while the other is working, so the reference files are never downloaded twice at once.
- The page carries `<meta name="status" content="beta">`, so app indexes such as apps.gnps2.org list it as in development.
- The builder writes the generated page with mode 0644 instead of the temporary file's 0600, and its stale-output error and the page's header comment no longer name repository-specific paths, so they read correctly in vendored copies.

### Notes

- Offline analysis requires previously downloaded FBMN quantification and annotation tables. Task-ID retrieval and external links require internet access.
- The standalone browser port remains under development; validate important results against source data and the original workflow.
