import streamlit as st

def welcome_page():
    # st.title("💊 Welcome to the Drug Readout App")

    st.markdown("""
    ## Citation
    [Empirically establishing drug exposure records directly from untargeted metabolomics data - bioRxiv](https://www.biorxiv.org/content/10.1101/2024.10.07.617109v1)

    
    ## Purpose
    This app is designed to **contextualize GNPS molecular networking results** into biological metadata using standardized drug-related vocabularies. It enables interpretation of potential **drug exposures** by linking molecular features to reference **pharmacological ontologies**.
    
    ## How It Works
    1. **Enter a GNPS Feature-Based Molecular Networking task ID and set parameters**
    2. The app retrieves and matches compound annotations to entries in a curated drug ontology
    3. Visualize and customize the annotations table and rerun analysis
    4. Results are categorized and summarized at various levels (e.g., drug class, subclass, mechanism)
    5. Outputs are interactive and downloadable for downstream analysis
    
    ## Key Features
    - 🔍 Annotates input feature tables with drug classifications
    - 🧠 Maps detected features to **pharmacological ontologies** across multiple levels
    - 🧬 Leverages curated vocabularies to support biological and clinical interpretation
    
    ## Related Tools
    - [Food Readout App](https://foodreadouts.streamlit.app/): interprets dietary exposure
    - [CMMC Dashboard](https://cmmc-dashboard.gnps2.org/): explores microbiome-derived metabolite annotations
    
    ---
    
    This application is part of the GNPS downstream analysis ecosystem known as **MetaboApps**.
    If you have feedback or suggestions, feel free to reach out to the development team.
    [Checkout other tools](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/toolindex/#gnps2-web-tools)
    """)
