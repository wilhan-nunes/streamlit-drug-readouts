import streamlit as st
import uuid  # Added for generating unique identifiers

from script import *
from utils import fetch_file, highlight_yes

#Example url: http://localhost:8501/?taskid=d6f37a11d90c4f249974280c3fc90108&threshold=1000

# cute badges
BADGE_TASK_ID_ = ":green-badge[Task ID]"
UPLOAD_QUANT_TABLE_ = ":blue-badge[Upload Quant Table]"
BADGE_QUANT_TABLE_ = ":orange-badge[Quant Table]"

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
    intensity_thresh = st.slider("Peak Area Threshold", min_value=100, max_value=50000, value=int(threshold), step=100,
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

        tab1, tab2, tab3, tab4 = st.tabs(
            ["Feature Annotation", "Drug Detection (Excluding Analogs)", "Drug Detection (Including Analogs)", "Drug Class Summary"])

        with tab1:
            st.write("Feature Annotation Table")
            st.dataframe(feature_annotation)

        with tab2:
            st.write("Drug detection Table (Excluding Analogs)")
            #Style dataframe to color "yes" and "No"
            st.dataframe(stratified_df.style.map(highlight_yes))


        with tab3:
            st.write("Drug detection Table (Including Analogs)")
            st.dataframe(stratified_df_analogs.style.map(highlight_yes))

        with tab4:
            st.write("Drug Class Summary")
            st.dataframe(class_count_df_sorted)


    except Exception as e:
        st.error(f"An error occurred: {e}")
        # raise
else:
    st.info(
        ":information_source: Please, provide the inputs, then click Run Analysis.")
