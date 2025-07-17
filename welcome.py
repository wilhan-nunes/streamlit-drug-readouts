import streamlit as st

def welcome_page():
    # st.title("💊 Welcome to the Drug Readout App")

    st.markdown("""
    ### 📖 Citation
    Haoqi Nina Zhao, Kine Eide Kvitne, Corinna Brungs, Siddharth Mohan, Vincent Charron-Lamoureux, *et al.*  
    **Empirically establishing drug exposure records directly from untargeted metabolomics data**  
    *bioRxiv* 2024.10.07.617109; doi: [https://doi.org/10.1101/2024.10.07.617109](https://doi.org/10.1101/2024.10.07.617109)


    
    ### 🧭 Purpose
    This app is designed to **contextualize GNPS molecular networking results** into biological metadata using standardized drug-related vocabularies. It enables interpretation of potential **drug exposures** by linking molecular features to reference **pharmacological ontologies**.
    
    ### 📘 How It Works
    1. **Enter a GNPS Feature-Based Molecular Networking task ID and set parameters**
    2. The app retrieves and matches compound annotations to entries in a curated drug ontology
    3. Visualize and customize the annotations table and rerun analysis
    4. Results are categorized and summarized at various levels (e.g., drug class, subclass, mechanism)
    5. Outputs are interactive and downloadable for downstream analysis
    
    ### 🧩 Key Features
    - 🔍 Annotates input feature tables with drug classifications
    - 🧠 Maps detected features to **pharmacological ontologies** across multiple levels
    - 🧬 Leverages curated vocabularies to support biological and clinical interpretation
    
     ### 🧪 Example Dataset
    To explore the functionality without inputting your own data, use the **"Load example"** checkbox in the sidebar.
    
    ### 🔗 Related Tools
    - [Food Readout App](https://foodreadouts.streamlit.app/): interprets dietary exposure
    - [CMMC Dashboard](https://cmmc-dashboard.gnps2.org/): explores microbiome-derived metabolite annotations
    
    ---    
    """)
    st.info("""
    - This application is part of the GNPS downstream analysis ecosystem known as **MetaboApps**.
    - If you encounter any issues or have suggestions, please reach out to the app maintainers.
    - [Checkout other tools](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/toolindex/#gnps2-web-tools)
    """)
