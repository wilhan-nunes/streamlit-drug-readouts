import streamlit as st
import uuid  # Added for generating unique identifiers

# Set wide mode as default
st.set_page_config(layout="wide")

from script import *
from utils import fetch_file, highlight_yes

#Example url: http://localhost:8501/?taskid=d6f37a11d90c4f249974280c3fc90108&threshold=1000

# cute badges
BADGE_TASK_ID_ = ":green-badge[Task ID]"

# Streamlit app title
st.title("Drug Readout Analysis")

# defining query params to populate input fields
query_params = st.query_params
gnps_task_id = query_params.get('taskid', '')
threshold = query_params.get('threshold', 100)

# Sidebar inputs
with st.sidebar:
    st.header("Inputs")

    task_id = st.text_input(f"{BADGE_TASK_ID_} FBMN Workflow Task ID (GNPS2)",
                            help="Enter the Task ID from a FBMN Workflow to retrieve the result files.",
                            placeholder='enter task ID...',
                            value=gnps_task_id)
    # slider to set threshold for the number of features
    intensity_thresh = st.number_input("Peak Area Threshold", min_value=100, max_value=50000, value=int(threshold), step=100,
                                 help="Only detections with peak area above this number will be considered.")
    if not task_id:
        st.warning(f"Please enter a {BADGE_TASK_ID_} from a FBMN Workflow to proceed.", )

    run_analysis = st.button("Run Analysis", help="Click to start the analysis with the provided inputs.",
                             use_container_width=True)

# Process files when task ID and sample feature table are provided
if run_analysis:
    try:
        # Generate a unique identifier using uuid
        unique_id = uuid.uuid4().hex

        # Retrieve lib_search using the task ID
        with st.spinner("Downloading Task result files..."):
            quant_file_path = fetch_file(task_id, f"quant_table_{unique_id}.csv", type="quant_table")
            st.markdown(f"Quant table downloaded successfully from task {task_id}! [Download](https://gnps2.org/result?task={task_id}&viewname=quantificationdownload&resultdisplay_type=task)", unsafe_allow_html=True)
            annotation_file_path = fetch_file(
                task_id, f"annotations_{unique_id}.tsv", type="annotation_table"
            )
            st.markdown(
                f"Annotation table downloaded successfully from task {task_id}! [Download](https://gnps2.org/resultfile?task={task_id}&file=nf_output/library/merged_results_with_gnps.tsv)")
            drug_metadata_file = "data/GNPS_Drug_Library_Metadata_Drugs.csv"
            analog_metadata_file = "data/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv"

        # Process data
        with st.spinner("Processing data..."):

            feature_filtered = load_and_filter_features(quant_file_path, intensity_threshold=intensity_thresh)
            annotation_metadata = load_and_merge_annotations(
                annotation_file_path, drug_metadata_file, analog_metadata_file
            )
            feature_annotation = generate_feature_annotation(
                annotation_metadata, feature_filtered
            )
            stratified_df = stratify_by_drug_class(feature_annotation, exclude_analogs=True)
            stratified_df_analogs = stratify_by_drug_class(feature_annotation, exclude_analogs=False)

            #Counting drug class occurrence per sample
            class_count_df = count_drug_class_occurrences(feature_annotation, class_column="pharmacologic_class")
            class_count_df["total_matches"] = class_count_df.sum(axis=1)
            class_count_df_sorted = class_count_df.sort_values("total_matches", ascending=False)


        tab4, tab1, tab2, tab3 = st.tabs(
            ["Summary Statistics", "Feature Annotation", "Drug Detection Tables", "Drug Class Summary"])

        with tab1:
            st.write("Feature Annotation Table")
            st.dataframe(feature_annotation)

        with tab2:
            st.write("Drug Detection Tables")
            with st.expander("Excluding Analogs", expanded=True):
                st.dataframe(stratified_df.style.map(highlight_yes))
            with st.expander("Including Analogs", expanded=False):
                st.dataframe(stratified_df_analogs.style.map(highlight_yes))

        with tab3:
            st.write("Top detected drug classes")
            top_pharm_num = st.number_input("Top classes to display", min_value=1, max_value=20, value=10, key="top_classes",
                                           help="Number of top pharmacologic classes to display based on total detections.")
            top_pharm_classes = class_count_df_sorted.sum(axis=0).nlargest(top_pharm_num)
            st.dataframe(top_pharm_classes.reset_index().rename(columns={'index': 'Pharmacologic Class', 0: 'Total Detections'}))
            class_count_df_sorted.index = class_count_df_sorted.index.str.replace(r'\.(mzML|mzXML) Peak area', '', regex=True)

            st.write("Drug class summary by sample")
            st.dataframe(class_count_df_sorted[["total_matches"] + class_count_df_sorted.columns.tolist()[:-1]].reset_index().rename(columns={'index': 'Sample'}))

        with tab4:
            st.subheader("Summary Statistics")
            antibiotic_count = (stratified_df["antibiotics"] == "Yes").sum()
            antidepressant_count = (stratified_df["antidepressants"] == "Yes").sum()
            sample_count = len(stratified_df)
            antibiotic_pct = (antibiotic_count / sample_count) * 100 if sample_count > 0 else 0
            antidepressant_pct = (antidepressant_count / sample_count) * 100 if sample_count > 0 else 0
            st.markdown(f"**Samples with Antibiotics:** {antibiotic_count} ({antibiotic_pct:.2f}\\%)")
            st.markdown(f"**Samples with Antidepressants:** {antidepressant_count} ({antidepressant_pct:.2f}\\%)")

    except Exception as e:
        st.error(f"An error occurred: {e}")
        # raise
else:
    st.info(
        ":information_source: Please, provide the inputs, then click Run Analysis.")
