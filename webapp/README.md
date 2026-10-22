# Drug Readout standalone web app

**Status: under development.** This is a browser-only version of the Streamlit Drug Readouts app. It is intended for exploration and feedback; verify important interpretations against the original workflow and source data.

## Download and use offline

Download [drug_readouts.html](drug_readouts.html) and open it in a current Chrome, Edge, Firefox, or Safari browser. The HTML contains its styles, scripts, and compressed drug ontology, so no Python installation, web server, or internet connection is needed when using uploaded tables.

For offline analysis, obtain the two result tables from a GNPS2 Feature-Based Molecular Networking (FBMN) job while online, then:

1. Open `drug_readouts.html` locally and choose **Uploaded tables**.
2. Provide `featuretable_reformated.csv` as the quantification table and `merged_results_with_gnps.tsv` as the annotation table.
3. Set the peak-area threshold and, if needed, a regular expression for blank sample IDs. Click **Run analysis**.
4. Review the summaries, top-class chart, UpSet plot, classification overview, and tables. Use the download buttons to save results.

The quantification table needs a `row ID` column and sample columns containing `.mzML` or `.mzXML`. The annotation table needs `#Scan#`, `SpectrumID`, `Compound_Name`, `MQScore`, and `SharedPeaks`. The two tables should come from the same FBMN task.

You can instead choose **Public GNPS2 task**, enter its task ID, and run without downloading tables first. That option, the **Load example** shortcut, and external documentation links require internet access. Uploaded files and analysis results stay in your browser; the task-ID option requests public result files from the selected GNPS2 host.

The page starts with its embedded drug metadata. **Fetch latest metadata** checks the latest version of the [GNPS Drug Library Zenodo record](https://zenodo.org/records/17232042), downloads its two metadata CSVs through Zenodo's browser-readable API, and caches them in this browser using IndexedDB. A cached copy is reused on later visits, including offline visits. Use **Use embedded & clear cache** to return to the built-in dataset. After changing metadata, rerun the analysis; existing results are not silently overwritten. Zenodo metadata may differ from the embedded version and change readouts. The browser cache belongs to this browser/site, not to the HTML file: moving the HTML to another device does not carry the refreshed dataset with it.

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
