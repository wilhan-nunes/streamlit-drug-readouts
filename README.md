
# Streamlit Drug Readouts

This project is a Streamlit-based application designed to process and analyze drug-related data. It allows users to upload input a GNPS2 FBMN task ID to fetch task results, and perform drug readouts from the data.

## Standalone web app (under development)

The browser-based version is [webapp/drug_readouts.html](webapp/drug_readouts.html). It is a self-contained HTML file that can be downloaded and opened locally without Python or a server. Upload the two GNPS2 FBMN result tables to analyze them offline; fetching results by GNPS2 task ID still needs an internet connection. See the [web app README](webapp/README.md) for usage, limitations, and the [changelog](webapp/CHANGELOG.md) for development updates.

The browser app can optionally fetch the latest GNPS Drug Library metadata from Zenodo and cache it in the current browser for later offline use. The downloaded HTML itself keeps its original embedded metadata; see the web app README for the distinction.

When the drug metadata CSVs change, rebuild the standalone page with `python3 webapp/build_drug_readouts.py`; edit [its source HTML](webapp/source/drug_readouts.html) for UI changes.

## Features

- Fetch task result files (quant table and annotation table) using the task ID.
- User-defined intensity thresholds for feature peak areas.
- Visualize and analyze drug readout tables.

## Requirements

- Python 3.8 or higher
- Streamlit
- Required Python libraries (see `requirements.txt`)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/wilhan-nunes/streamlit_drug_readouts.git
   cd streamlit_drug_readouts
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Run the application:
   ```bash
   streamlit run app.py
   ```

## Usage

1. Open the application in your browser (usually at `http://localhost:8501`).
2. Provide the required inputs:
   - Task ID
   - Intensity threshold
3. Click **Run Analysis** to process the data.
4. Visualize and/or download the processed file results.

## File Structure

- `app.py`: Main application file.
- `data/`: Contains metadata files for drug libraries.
- `requirements.txt`: Lists the required Python libraries.
- `README.md`: Project documentation.

## Citation

Zhao H.N., Kvitne K.E., Brungs C., Mohan S., Charron-Lamoureux V. et al.
*A resource to empirically establish drug exposure records directly from
untargeted metabolomics data.* Nature Communications 16, 10600 (2025).
<https://doi.org/10.1038/s41467-025-65993-5>

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Credits

- [Nina Zhao](https://github.com/ninahaoqizhao) for the original R code.

