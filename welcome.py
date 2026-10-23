import streamlit as st

def welcome_page():
    # st.title("💊 Welcome to the Drug Readout App")

    st.markdown("""
    ### 📖 Citation
    Haoqi Nina Zhao, Kine Eide Kvitne, Corinna Brungs, Siddharth Mohan, Vincent Charron-Lamoureux, *et al.*  
    **A resource to empirically establish drug exposure records directly from untargeted metabolomics data**<br>
    *Nature Communications* 16, 10600 (2025); doi: [https://doi.org/10.1038/s41467-025-65993-5](https://doi.org/10.1038/s41467-025-65993-5)


    
    ### 🧭 Purpose
    This app translates GNPS compound annotations into drug exposure information. The app utilizes standardized pharmacologic vocabularies developed for all drug reference spectra on GNPS. By linking molecular features to pharmacological ontologies, it enables rapid identification of drugs present in your data and provides insights into their properties.
    
    ### 📘 How It Works
    1. Download GNPS Drug Library mgf spectra: https://zenodo.org/records/15259192
       (Note: The App still works if you use the default library in FBMN job. However, using the GNPS Drug Library will increase the coverage of your drug detection, especially for drug metabolites).
    2. Run a FBMN job using the GNPS Drug Library MGF file
    3. Enter the GNPS Feature-Based Molecular Networking task ID and set parameters in the DrugReadoutApp
    4. The app retrieves and matches compound annotations to entries in a curated drug ontology
    5. The app also visualizes and allows customization of the annotations table
    
    ### 🧩 Key Features
    - 🔍 Annotates input feature tables with drug classifications
    - 🧠 Maps detected features to **pharmacological ontologies** across multiple levels
    - 🧬 Leverages curated vocabularies to support biological and clinical interpretation
    
     ### 🧪 Example Dataset
    To explore the functionality without inputting your own data, use the **"Load example"** checkbox in the sidebar.
    
    ### 🔗 Related Tools
    - [Food Readout App](https://foodreadouts.gnps2.org/): interprets dietary exposure
    - [CMMC Dashboard](https://cmmc-dashboard.gnps2.org/): explores microbiome-derived metabolite annotations
    
    ---    
    """)
    st.info("""
    - This application is part of the GNPS downstream analysis ecosystem known as **MetaboApps**.
    - If you encounter any issues or have suggestions, please reach out to the app maintainers.
    - [Checkout other tools](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/toolindex/#gnps2-web-tools)
    """)
