# Drug Readout standalone web app

**Status: under development.** This is a browser-only version of the Streamlit Drug Readouts app. It is intended for exploration and feedback; verify important interpretations against the original workflow and source data.

## Download and use offline

Download [drug_readouts.html](drug_readouts.html) and open it in a current Chrome, Edge, Firefox, or Safari browser. The HTML contains its styles, scripts, and compressed drug ontology, so no Python installation, web server, or internet connection is needed when using uploaded tables.

For offline analysis, obtain the two result tables from a GNPS2 Feature-Based Molecular Networking (FBMN) job while online, then:

1. Open `drug_readouts.html` locally and select the **FBMN tables** tab.
2. Provide `featuretable_reformated.csv` as the quantification table and `merged_results_with_gnps.tsv` as the annotation table.
3. Set the peak-area threshold and, if needed, a regular expression for blank sample IDs. Click **Run analysis**.
4. Review the summaries, top-class chart, UpSet plot, classification overview, and tables. Use the download buttons to save results.

The quantification table needs a `row ID` column and sample columns containing `.mzML` or `.mzXML`. The annotation table needs `#Scan#`, `SpectrumID`, `Compound_Name`, `MQScore`, and `SharedPeaks`. The two tables should come from the same FBMN task.

You can instead use the **GNPS2 task** tab: enter a task ID and run without downloading tables first. On every run, the page reads `submission_parameters.yaml` from that job and uses its `workflowname` to choose the FBMN or Everything Bagel loader. Tasks with missing metadata or unsupported workflows stop with an error instead of guessing result paths. That tab, its **Load example** shortcut, and external documentation links require internet access. Uploaded files and analysis results stay in your browser; the task-ID option requests public result files from the selected GNPS2 host.

## Link parameters

The page accepts the same query parameters as the Streamlit app:

| Parameter | Effect |
|---|---|
| `task_id` | GNPS2 task ID; selects the **GNPS2 task** tab |
| `threshold` | Peak area threshold, 100 to 5,000,000,000 (default 1000); an invalid value is ignored with a warning |
| `blank_ids` | Blank sample pattern, e.g. `BLANK\|QC` |

For example: `drug_readouts.html?task_id=4d99fc25d84143bdbbf2dd07bf044e5e&threshold=1000&blank_ids=blank`. As in the Streamlit app, these only prefill the inputs; the analysis starts when you click **Run analysis**.

## Library search preparation

Library search runs on query spectra from one of two sources:

- **Everything Bagel task.** In the **GNPS2 task** tab, enter an Everything Bagel task ID, and select **Run analysis**. Once the job metadata identifies it as Everything Bagel, the **Library search** panel appears.
- **Your own files, no task ID.** Use the **Feature table + query MGF** tab. The feature table needs a `row ID` column and sample columns containing `.mzML` or `.mzXML`, such as an EB `aligned_features.csv`. The query MGF's `SCANS` values must match those row IDs. The MGF is read from disk when needed and never uploaded.


For an Everything Bagel task, **Prepare library search files** downloads `nf_output/feature_finding/feature_finding_results/aligned_features.csv` and `nf_output/feature_finding/aligned_features_filled.mgf` from GNPS2. For both sources it checks that MGF `SCANS` values match feature table `row ID` values. It then retrieves the latest version of `GNPS_Drug_Library_Spectra_Drugs_and_Metabolites.mgf` from Zenodo, plus `GNPS_Drug_Library_Spectra_Drug_Analogs_Updated.mgf` when the analog option is selected. It checks file size and the MGF header, and caches complete reference files in IndexedDB when the browser has space. These files are large (hundreds of MB), so preparation starts only on request. While they download, the status bar shows each file's progress (MB received, percent, speed and time left) above a progress bar. Later preparations reuse the cached copies. The **Files** box in the Library search panel shows, for the query spectra and each reference MGF, whether it is already in memory, cached in this browser, read from a local file, or still to be downloaded (and how large). It reads only the browser cache. **Check files** also asks Zenodo for the latest version, to tell whether a cached copy is current. The task MGF stays in memory for the current page session.

Library search is pending the external engine. Until it is installed, **Run analysis** in either library-search mode reports that dependency rather than displaying fabricated matches. The engine integration point is `window.DrugReadoutsEBEngine.search({featureMgf, libraries, taskId, zenodoRecordId})`, where `featureMgf` and each `libraries[].blob` are `Blob` objects (the uploaded query MGF is a `File`), and `taskId` is `null` for a local library search. The engine must return an array of `{featureId, spectrumId, compoundName, score, matchedPeaks}` hits. Feature IDs must equal MGF `SCANS` and table `row ID`; spectrum IDs must match the drug/analog metadata identifiers. Once the engine is supplied, these hits flow through the existing readout and visualization code.

The page starts with its embedded drug metadata. **Fetch latest metadata**, next to **Run analysis**, checks the latest version of the [GNPS Drug Library Zenodo record](https://zenodo.org/records/17232042), downloads its two metadata CSVs through Zenodo's browser-readable API, and caches them in this browser using IndexedDB. A cached copy is reused on later visits, including offline visits. When Zenodo metadata is in use, a **Use embedded metadata** button appears beside it; it returns to the built-in dataset and clears the cache. Hover over the Fetch button to see which version is in use, or over the ⓘ beside it for this explanation. After changing metadata, rerun the analysis; existing results are not silently overwritten. Zenodo metadata may differ from the embedded version and change readouts. The browser cache belongs to this browser/site, not to the HTML file: moving the HTML to another device does not carry the refreshed dataset with it.

## Input tabs

The **Inputs** panel has one tab per source: **GNPS2 task**, **FBMN tables** and **Feature table + query MGF**. Each tab has a **How to use this tab** popup (hover over it or focus it with the keyboard; Escape closes it) explaining what it needs, and **About the input files** at the bottom repeats the file requirements for all three. **Run analysis** uses the selected tab. Files added in other tabs are kept, so you can switch back without re-adding them. Peak-area threshold and blank IDs are under **Analysis settings**, which applies to every tab. A drop box turns green with a ✓ once it holds a file. Use the arrow keys to move between tabs.

## Working with results

- Top pharmacological classes shows sample counts and percentages relative to samples with any detected class.
- Class summary includes an UpSet-style co-occurrence plot. Its top bars count exact class combinations, its left bars count samples containing each class, and connected dots identify combination membership. The plot can be saved as SVG.
- Edit blue classification cells in Feature annotations, then select **Apply edits & recalculate** to refresh every result. You can also remove low-confidence annotations and recalculate.
- Tables can be downloaded as TSV. Results are session-local; save downloads before closing or reloading the page.

This browser port is under development, and its outputs should be checked before research or clinical interpretation. For project background and the original Streamlit setup, see the [repository README](../README.md). Development notes are in [CHANGELOG.md](CHANGELOG.md).

## Rebuild after ontology updates

The editable page is [source/drug_readouts.html](source/drug_readouts.html). The distributed `drug_readouts.html` is generated; edit the source page, not the generated copy. The builder embeds the two CSVs from the repository's `data/` directory as deterministic gzip/base64 blocks, following the other standalone MetaboApps build pattern.

From the repository root, run:

```bash
python3 webapp/build_drug_readouts.py
python3 webapp/build_drug_readouts.py --check
```

Run the first command whenever either CSV or the editable HTML changes. `--check` exits with an error if the generated page is stale. The script uses only Python's standard library and makes no network requests. To build a copy for another app host without changing this repository's output, pass `--output /path/to/drug_readouts.html`.

The source CSVs can also be obtained from the [GNPS Drug Library Zenodo record](https://zenodo.org/records/17232042). To download and embed that record's two metadata CSVs in one step, run:

```bash
python3 webapp/build_drug_readouts.py --zenodo
```

The builder verifies the downloads against Zenodo's published checksums and does not replace the repository's `data/` files. After a `--zenodo` build, use `--check --zenodo` to check against the same inputs; a plain `--check` checks against the repository defaults.

Alternatively, download `GNPS_Drug_Library_Metadata_Drugs.csv` and `GNPS_Drug_Library_Metadata_Drug_Analogs_Updated_With_Source.csv` yourself (not the much larger MGF spectral files). The Zenodo analog CSV has a different filename from the repository's default analog CSV, so pass both paths explicitly when building from downloaded files:

```bash
python3 webapp/build_drug_readouts.py \
  --drugs /path/to/GNPS_Drug_Library_Metadata_Drugs.csv \
  --analogs /path/to/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated_With_Source.csv
```

The builder checks that required columns are present and records each input file's name and SHA-256 digest in the generated HTML. Review changes in classifications when switching dataset versions.
