import uuid  # Added for generating unique identifiers

import streamlit as st
import re

# Set wide mode as default
st.set_page_config(layout="wide")

if "run_analysis" not in st.session_state:
    st.session_state.run_analysis = False

from script import *
from utils import fetch_file, highlight_yes

# Example url: http://localhost:8501/?taskid=d6f37a11d90c4f249974280c3fc90108&threshold=1000&blank_ids=QC

# cute badges
BADGE_TASK_ID_ = ":green-badge[Task ID]"

# Streamlit app title
st.title("Drug Readout Analysis")

# defining query params to populate input fields
query_params = st.query_params
gnps_task_id = query_params.get('taskid', '')
threshold = query_params.get('threshold', 100)
blank_str = query_params.get('blank_ids', None)

# Sidebar inputs
with st.sidebar:
    st.header("Inputs")

    task_id = st.text_input(f"{BADGE_TASK_ID_} FBMN Workflow Task ID (GNPS2)",
                            help="Enter the Task ID from a FBMN Workflow to retrieve the result files.",
                            placeholder='enter task ID...',
                            value=gnps_task_id)
    # slider to set threshold for the number of features
    intensity_thresh = st.number_input("Peak Area Threshold", min_value=100, max_value=50000, value=int(threshold),
                                       step=100,
                                       help="Only detections with peak area above this number will be considered.")

    blank_ids = st.text_input("Blank IDs (optional)", value=blank_str, placeholder="Example: BLANK|IS|PoolQC|QCmix|SRM",
                              help="Enter substrings to identify blank or control columns, separated by '|'. If given, the table will be filtered to remove these columns from the analysis. If not provided, all columns will be considered.")

    if not task_id:
        st.warning(f"Please enter a {BADGE_TASK_ID_} from a FBMN Workflow to proceed.", )

    run_analysis = st.button("Run Analysis", icon="🏁", help="Click to start the analysis with the provided inputs.",
                             use_container_width=True, key="run_analysis_button")

    if st.button('Restart Session', icon="♻️" ,key='restart_session', use_container_width=True, type='primary'):
        # Reset the session state
        st.session_state.clear()
        st.rerun()

# Process files when task ID and sample feature table are provided
if run_analysis:
    try:
        # Generate a unique identifier using uuid
        unique_id = uuid.uuid4().hex

        # Retrieve lib_search using the task ID
        with st.spinner("Downloading Task result files..."):
            quant_file_path = fetch_file(task_id, f"quant_table_{unique_id}.csv", type="quant_table")
            st.markdown(
                f":white_check_mark: Quant table ({task_id}) [Download](https://gnps2.org/result?task={task_id}&viewname=quantificationdownload&resultdisplay_type=task)",
                unsafe_allow_html=True)
            annotation_file_path = fetch_file(
                task_id, f"annotations_{unique_id}.tsv", type="annotation_table"
            )
            st.markdown(
                f":white_check_mark: Annotation table ({task_id}) [Download](https://gnps2.org/resultfile?task={task_id}&file=nf_output/library/merged_results_with_gnps.tsv)")

            drug_metadata_file = "data/GNPS_Drug_Library_Metadata_Drugs.csv"
            analog_metadata_file = "data/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv"

        # Process data
        with st.spinner("Processing data..."):
            subtract_blanks = True if blank_ids else False
            feature_filtered = load_and_filter_features(quant_file_path, intensity_threshold=intensity_thresh,
                                                        blank_ids=blank_ids, subtract_blanks=subtract_blanks)

            annotation_metadata = load_and_merge_annotations(
                annotation_file_path, drug_metadata_file, analog_metadata_file
            )

            if "feature_annotation_edited" in st.session_state:
                # If the user has edited the feature annotation, use that instead of the original
                feature_annotation = st.session_state.feature_annotation_edited
            else:
                feature_annotation = generate_feature_annotation(
                    annotation_metadata, feature_filtered
                )

            # Perform analysis
            stratified_df = stratify_by_drug_class(feature_annotation, exclude_analogs=True)
            stratified_df_analogs = stratify_by_drug_class(feature_annotation, exclude_analogs=False)

            # Counting drug class occurrence per sample
            class_count_df = count_drug_class_occurrences(feature_annotation, class_column="pharmacologic_class")
            class_count_df["total_matches"] = class_count_df.sum(axis=1)
            class_count_df_sorted = class_count_df.sort_values("total_matches", ascending=False)

            # store all dataframes in session state
            st.session_state.feature_annotation = feature_annotation
            st.session_state.stratified_df = stratified_df
            st.session_state.stratified_df_analogs = stratified_df_analogs
            st.session_state.class_count_df_sorted = class_count_df_sorted

            st.session_state.run_analysis = True

    except Exception as e:
        st.error(f"An error occurred: {e}")
        # raise
else:
    st.info(
        ":information_source: Please, provide the inputs, then click Run Analysis.")

# Display results if analysis has been run
if st.session_state.run_analysis:
    # Summary Statistics Section
    st.header("📊 Summary Statistics")
    stratified_df = st.session_state.get("stratified_df")
    antibiotic_count = (stratified_df["antibiotics"] == "Yes").sum()
    antidepressant_count = (stratified_df["antidepressants"] == "Yes").sum()
    sample_count = len(stratified_df)
    antibiotic_pct = (antibiotic_count / sample_count) * 100 if sample_count > 0 else 0
    antidepressant_pct = (antidepressant_count / sample_count) * 100 if sample_count > 0 else 0

    col1, col2 = st.columns(2)
    with col1:
        st.metric("Samples with Antibiotics", f"{antibiotic_count} ({antibiotic_pct:.2f}%)")
    with col2:
        st.metric("Samples with Antidepressants", f"{antidepressant_count} ({antidepressant_pct:.2f}%)")

    st.divider()

    # Feature Annotation Table Section
    st.header("🔬 Feature Annotation Table")
    st.write("You can edit the table below and then rerun the analysis with your modifications.")

    # Store the original dataframe in session state if not already there
    if "feature_annotation_edited" not in st.session_state:
        st.session_state.feature_annotation_edited = st.session_state.feature_annotation.copy()

    # Use data_editor to allow editing and preserve changes
    edited_df = st.data_editor(
        st.session_state.feature_annotation_edited,
        key="feature_annotation_editor",
        use_container_width=True,
        num_rows="dynamic",
        height=400
    )

    # Update session state with edits
    st.session_state.feature_annotation_edited = edited_df

    # Rerun button
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        rerun_button = st.button("🔄 Rerun Analysis with Edited Data",
                                 use_container_width=True,
                                 type="primary")

    if rerun_button:
        with st.spinner("Reprocessing data with edited annotations..."):
            # Use the edited feature annotation for reanalysis
            feature_annotation = st.session_state.feature_annotation_edited

            # Perform analysis with edited data
            stratified_df = stratify_by_drug_class(feature_annotation, exclude_analogs=True)
            stratified_df_analogs = stratify_by_drug_class(feature_annotation, exclude_analogs=False)

            # Counting drug class occurrence per sample
            class_count_df = count_drug_class_occurrences(feature_annotation, class_column="pharmacologic_class")
            class_count_df["total_matches"] = class_count_df.sum(axis=1)
            class_count_df_sorted = class_count_df.sort_values("total_matches", ascending=False)

            # Update session state with new results
            st.session_state.stratified_df = stratified_df
            st.session_state.stratified_df_analogs = stratified_df_analogs
            st.session_state.class_count_df_sorted = class_count_df_sorted

            st.success("Analysis updated with edited data!")
            st.rerun()

    st.divider()

    # Drug Detection Tables Section
    st.header("🧪 Drug Detection Tables")

    stratified_df = st.session_state.get("stratified_df")
    stratified_df_analogs = st.session_state.get("stratified_df_analogs")

    # optional: clean up sample names
    stratified_df['Sample'] = stratified_df['Sample'].str.replace(r'\.mz[XM]L Peak area', '', regex=True)
    stratified_df_analogs['Sample'] = stratified_df_analogs['Sample'].str.replace(r'\.mz[XM]L Peak area', '', regex=True)

    st.subheader("Excluding Analogs")
    st.dataframe(stratified_df.style.map(highlight_yes), use_container_width=True)

    with st.expander("Show results including analogs"):
        st.dataframe(stratified_df_analogs.style.map(highlight_yes), use_container_width=True)

    st.divider()

    # Drug Class Summary Section

    import matplotlib.pyplot as plt
    from upsetplot import plot, from_indicators

    # Drug Class Summary Section
    st.header("📈 Drug Class Summary")

    class_count_df_sorted = st.session_state.get("class_count_df_sorted")

    # Create tabs for different visualizations
    tab_plot, tab_table = st.tabs(["🔀 UpSet Plot", "📊 Top Classes Table", ])

    with tab_plot:
        st.subheader("Drug Class Co-occurrence Analysis")
        st.write("This UpSet plot shows how different drug classes co-occur across samples. Each bar represents a unique combination of drug classes.")

        # Controls for UpSet plot
        col1, col2 = st.columns(2)
        with col1:
            n_top_classes = st.number_input("Number of Top Classes for UpSet Plot",
                                           min_value=3, max_value=20, value=4,
                                           key="upset_classes_input",
                                           help="Select how many top drug classes to include in the UpSet plot")
        with col2:
            max_samples = st.number_input("Maximum Samples to Include",
                                         min_value=10, max_value=200, value=50,
                                         key="upset_samples_input",
                                         help="Limit the number of samples to avoid overcrowding")

        try:
            # Prepare binary matrix for top classes
            top_classes = class_count_df_sorted.sum(axis=0).nlargest(n_top_classes).index.tolist()
            if "total_matches" in top_classes:
                top_classes.remove("total_matches")

            # Create binary matrix (presence/absence)
            binary_matrix = (class_count_df_sorted[top_classes] > 0).astype(int)
            binary_matrix.index.name = "Sample"

            # Limit number of samples to avoid overcrowding
            limited_matrix = binary_matrix.head(max_samples)

            # Remove samples with no detections in the selected classes
            limited_matrix = limited_matrix[limited_matrix.sum(axis=1) > 0]

            if len(limited_matrix) == 0:
                st.warning("No samples found with detections in the selected drug classes. Try adjusting the parameters.")
            else:
                # Create UpSet plot data
                upset_data = from_indicators(top_classes, limited_matrix.astype(bool))

                # Create the plot
                fig = plt.figure(figsize=(12, 8))
                plot(upset_data,
                     subset_size='count',
                     show_counts=True,
                     fig=fig,
                     sort_by='cardinality',
                     element_size=40)

                plt.suptitle(f'Drug Class Co-occurrence Analysis\n({len(limited_matrix)} samples, {len(top_classes)} drug classes)',
                            fontsize=14, y=0.98)

                st.pyplot(fig)

                # Add interpretation help
                with st.expander("How to interpret this UpSet plot"):
                    st.write("""
                    - **Horizontal bars (left)**: Show the total number of samples containing each individual drug class
                    - **Vertical bars (bottom)**: Show the number of samples with each specific combination of drug classes
                    - **Connected dots**: Indicate which drug classes are part of each combination
                    - **Larger vertical bars**: Represent more common co-occurrence patterns
                    - **Single dots**: Show samples with only one drug class detected
                    - **Multiple connected dots**: Show samples with multiple drug classes detected together
                    """)

                plt.close()  # Clean up the figure

        except Exception as e:
            st.error(f"Error creating UpSet plot: {str(e)}")
            st.info("This might be due to insufficient data or missing dependencies. Make sure you have the 'upsetplot' package installed.")

    with tab_table:
        st.subheader("Top Detected Drug Classes")
        nlarge = st.number_input("Number of Top Classes to Display", min_value=1, value=10, key="top_classes_input")
        top_pharm_classes = class_count_df_sorted.sum(axis=0).nlargest(nlarge)
        st.dataframe(
            top_pharm_classes.reset_index().rename(columns={'index': 'Pharmacologic Class', 0: 'Total Detections'}),
            use_container_width=True
        )


    st.subheader("Drug Class Summary by Sample")
    # Clean up sample names
    class_count_df_display = class_count_df_sorted.copy()
    class_count_df_display.index = class_count_df_display.index.str.replace(r'\.(mzML|mzXML) Peak area', '', regex=True)
    class_count_df_display = class_count_df_display[["total_matches"] + class_count_df_display.columns.tolist()[:-1]].reset_index().rename(
            columns={'index': 'Sample'})

    st.dataframe(
        class_count_df_display,
        use_container_width=True
    )