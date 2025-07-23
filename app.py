from gnpsdata import workflow_fbmn
from streamlit.components.v1 import html
import streamlit as st
from utils import display_comparison_statistics, load_example, get_git_short_rev

# Set page configuration
page_title = "Drug Readout Analysis"

# TODO: Bump version
app_version = "2025-07-17"
git_hash = get_git_short_rev()
repo_link = "https://github.com/wilhan-nunes/streamlit_drug_readouts"

st.set_page_config(
    page_title=page_title,
    layout="wide",
    page_icon=":pill:",
    menu_items={"About": (f"**App version**: {app_version} | "
                          f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})")},
)

# Add a tracking token
html(
    '<script async defer data-website-id="74bc9983-13c4-4da0-89ae-b78209c13aaf" src="https://analytics.gnps2.org/umami.js"></script>',
    width=0,
    height=0,
)

if "run_analysis" not in st.session_state:
    st.session_state.run_analysis = False

from script import *
from utils import fetch_file, highlight_yes, add_sankey_graph

# Example url: http://localhost:8501/?taskid=d6f37a11d90c4f249974280c3fc90108&threshold=1000&blank_ids=QC

# cute badges
BADGE_TASK_ID_ = ":green-badge[Task ID]"

# Streamlit app title
st.title("Drug Readout Analysis")

# defining query params to populate input fields
query_params = st.query_params
gnps_task_id = query_params.get("taskid", "")
threshold = query_params.get("threshold", 1000)
blank_str = query_params.get("blank_ids", None)

# Sidebar inputs
with st.sidebar:
    st.header("Inputs")
    load_example_data = st.checkbox("Load example",
                                    help="Load example from FBMN task ID d6f37a11d90c4f249974280c3fc90108", value=False,
                                    key='load_example_check')
    task_id = st.text_input(
        f"{BADGE_TASK_ID_} FBMN Workflow Task ID (GNPS2)",
        help="Enter the Task ID from a FBMN Workflow to retrieve the result files.",
        placeholder="enter task ID...",
        value=gnps_task_id if not st.session_state.get('load_example_check') else "d6f37a11d90c4f249974280c3fc90108",
        disabled=(load_example_data == True)
    )
    # slider to set threshold for the number of features
    intensity_thresh = st.number_input(
        "Peak Area Threshold",
        min_value=1E2,
        max_value=5E9,
        value=float(threshold),
        step=1E2,
        help="Only detections with peak area above this number will be considered.",
    )

    blank_ids = st.text_input(
        "Blank IDs (optional)",
        value=blank_str if not st.session_state.get('load_example_check') else "QC",
        placeholder="Example: BLANK|IS|PoolQC|QCmix|SRM",
        help="Enter substrings to identify blank or control columns, separated by '|'. If given, the table will be filtered to remove these columns from the analysis. If not provided, all columns will be considered.",
    )

    if not task_id:
        st.warning(
            f"Please enter a {BADGE_TASK_ID_} from a FBMN Workflow to proceed.",
        )

    if st.button(
            "Run Analysis",
            icon="🏁",
            help="Click to start the analysis with the provided inputs.",
            use_container_width=True,
            key="run_analysis_button",
            disabled=not (task_id or load_example_data)
    ):
        st.session_state.run_analysis = True

    if st.button(
            "Restart Session",
            icon="♻️",
            key="restart_session",
            use_container_width=True,
            type="primary",
    ):
        # Reset the session state
        st.session_state.clear()
        st.session_state['load_example_data'] = False
        st.rerun()

    st.subheader("Contributors")
    st.markdown(
        """
    - [Haoqi (Nina) Zhao PhD](https://scholar.google.com/citations?user=xW9jBO0AAAAJ) - UC San Diego
    - [Wilhan Nunes PhD](https://scholar.google.com/citations?user=4cPVoeIAAAAJ) - UC San Diego
    - [Mingxun Wang PhD](https://www.cs.ucr.edu/~mingxunw/) - UC Riverside
    """
    )

    st.subheader("Documentations and Resources")
    st.markdown("""
    [Feature Based Molecular Networking](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/)
    """
                )
# Process files when task ID and sample feature table are provided
if st.session_state.get('run_analysis_button', False) or st.session_state.get("rerun_analysis", False):

    try:
        if not load_example_data:
            # Retrieve lib_search using the task ID
            with st.spinner("Downloading Task result files..."):
                quant_file_df = workflow_fbmn.get_quantification_dataframe(task_id, gnps2=True)
                annotation_file_df = workflow_fbmn.get_library_match_dataframe(task_id, gnps2=True)

                st.markdown(
                    f"[:material/download: Input Quant table]"
                    f"(https://gnps2.org/result?task={task_id}&viewname=quantificationdownload&resultdisplay_type=task) | "
                    f"[:material/download: Input Annotation table]"
                    f"(https://gnps2.org/resultfile?task={task_id}&file=nf_output/library/merged_results_with_gnps.tsv)",
                )
        else:
            quant_file_df, annotation_file_df = load_example()

        DRUG_METADATA_FILE = "data/GNPS_Drug_Library_Metadata_Drugs.csv"
        ANALOG_METADATA_FILE = (
            "data/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv"
        )

        # Process data
        with st.spinner("Processing data..."):
            subtract_blanks = True if blank_ids else False
            feature_filtered = load_and_filter_features(
                quant_file_df,
                intensity_threshold=intensity_thresh,
                blank_ids=blank_ids,
                subtract_blanks=subtract_blanks,
            )

            annotation_metadata = load_and_merge_annotations(
                annotation_file_df, DRUG_METADATA_FILE, ANALOG_METADATA_FILE
            )
            if not st.session_state.get("rerun_analysis", False):
                feature_annotation = generate_feature_annotation(
                    annotation_metadata, feature_filtered
                )
            else:
                feature_annotation = st.session_state.get('feature_annotation_edited')
                st.success(f'Feature annotation get from session state. Shape: {feature_annotation.shape}')
                st.session_state["rerun_analysis"] = False

            # Perform analysis
            stratified_df = stratify_by_drug_class(
                feature_annotation, exclude_analogs=True, peak_threshold=threshold,
            )
            stratified_df_analogs = stratify_by_drug_class(
                feature_annotation, exclude_analogs=False, peak_threshold=threshold,
            )

            # Counting drug class occurrence per sample
            class_count_df = count_drug_class_occurrences(
                feature_annotation, class_column="pharmacologic_class"
            )
            class_count_df["total_matches"] = class_count_df.sum(axis=1)
            class_count_df_sorted = class_count_df.sort_values(
                "total_matches", ascending=False
            )

            # store all dataframes in session state
            st.session_state.feature_annotation = feature_annotation
            st.session_state.stratified_df = stratified_df
            st.session_state.stratified_df_analogs = stratified_df_analogs
            st.session_state.class_count_df_sorted = class_count_df_sorted


    except Exception as e:
        st.error(f"An error occurred: {e}")
        # raise
else:
    st.info(":information_source: Please, provide the inputs, then click Run Analysis.")

# Display results if analysis has been run
if st.session_state.run_analysis:
    # Summary Statistics Section
    st.header("📊 Drug Detection Summary Statistics")
    with st.spinner("Calculating summary..."):
        stratified_df = st.session_state.get("stratified_df", None)
        sample_count = len(stratified_df)

        # Define specific drug categories to highlight
        specific_categories = {
            "antibiotics": "🦠 Antibiotics",
            "antidepressants": "🧠 Antidepressants",
            "statin": "💊 Statins",
            "PPI": "➕ PPIs (Proton Pump Inhibitors)",
            "antihistamine": "🤧 Antihistamines",
            "antihypertensive": "❤️ Antihypertensives",
            "Alzheimer": "🧠 Alzheimer's Meds",
            "antifungal": "🍄 Antifungals",
            "HIVmed": "🏥 HIV Medications",
        }

        # Calculate counts and percentages for specific categories
        category_stats = {}
        for col_name, display_name in specific_categories.items():
            if col_name in stratified_df.columns:
                count = (stratified_df[col_name] == "Yes").sum()
                pct = (count / sample_count) * 100 if sample_count > 0 else 0
                if pct != 0:
                    category_stats[display_name] = {"count": count, "percentage": pct}
                # sort

        # Display specific drug categories in a grid
        if category_stats:
            st.subheader("🎯 Specific Drug Categories")

            # Create columns dynamically based on number of categories
            num_categories = len(category_stats)
            cols_per_row = 3

            for i in range(0, num_categories, cols_per_row):
                cols = st.columns(min(cols_per_row, num_categories - i))

                for j, (display_name, stats) in enumerate(
                        list(category_stats.items())[i: i + cols_per_row]
                ):
                    with cols[j]:
                        st.metric(
                            display_name,
                            f"{stats['count']} ({stats['percentage']:.1f}%)",
                            help=f"Number of samples containing {display_name.split(' ', 1)[1].lower()}",
                        )

        # Overall summary metrics
        st.subheader("📈 Overall Summary")
        col1, col2, col3 = st.columns(3)

        with col1:
            total_drug_columns = len(
                [col for col in stratified_df.columns if col != "Sample"]
            )
            st.metric("Total Drug Categories", total_drug_columns)

        with col2:
            # Calculate samples with any drug detection
            drug_columns = [col for col in stratified_df.columns if col != "Sample"]
            samples_with_drugs = (
                (stratified_df[drug_columns] == "Yes").any(axis=1).sum()
            )
            drug_detection_pct = (
                (samples_with_drugs / sample_count) * 100 if sample_count > 0 else 0
            )
            st.metric(
                "Samples with Any Drug",
                f"{samples_with_drugs} ({drug_detection_pct:.1f}%)",
            )

        with col3:
            st.metric("Total Samples", sample_count)

    st.divider()

    # Feature Annotation Table Section
    st.header("🔬 Feature Annotation Table")
    st.write(
        "You can edit the table below and then rerun the analysis with your modifications. "
        "[Learn how](https://www.youtube.com/watch?v=6tah69LkfxE&list=TLGGKK4Dnf1gepcwNTA2MjAyNQ)"
    )
    st.warning(
        "***Before editing** the data, please clear all filters. Filters display results based on the original, unedited table.")

    feature_annotation = st.session_state.feature_annotation

    from utils import add_df_and_filtering

    filtered_df = add_df_and_filtering(feature_annotation, "feature_annotation")
    edited_df = st.data_editor(
        filtered_df,
        key="feature_annotation_editor",
        use_container_width=True,
        num_rows="dynamic",
        height=400,
    )

    # Rerun button
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button(
                "🔄 Rerun Analysis with Edited Data",
                use_container_width=True,
                type="primary",
                key="rerun_analysis_button",
        ):
            st.session_state['rerun_analysis'] = True

    with st.expander('Features excluded from analysis by default', icon="⛔️"):
        st.info("Excluded features with chemical source 'Endogenous' and/or 'Food'")
        st.dataframe(feature_annotation[
                         feature_annotation["chemical_source"].str.contains(
                             "Endogenous|Food", case=False, na=False
                         )
                     ]
                     )

    if st.session_state.get('rerun_analysis', False):
        rows_to_remove = st.session_state.feature_annotation_editor.get("deleted_rows", [])
        rows_to_edit = st.session_state.feature_annotation_editor.get("edited_rows", {})
        rows_to_add = st.session_state.feature_annotation_editor.get("added_rows", [])

        with st.spinner("Reprocessing data with edited annotations..."):
            # Use the edited feature annotation for reanalysis
            if "feature_annotation_editor" in st.session_state:
                st.write("using edited feature annotation")
                if not st.session_state.get('feature_annotation_editor'):
                    print('Editor empty')
                # If the user has edited the feature annotation, use that instead of the original
            else:
                st.write("using original feature annotation")

            feature_annotation = st.session_state.feature_annotation

            feature_annotation_edited = feature_annotation.copy()

            # Remove rows safely
            if rows_to_remove:
                feature_annotation_edited = feature_annotation_edited.drop(rows_to_remove, errors="ignore")

            # Edit rows
            for index, values_dict in rows_to_edit.items():
                for col, value in values_dict.items():
                    feature_annotation_edited.at[index, col] = value

            # Add new rows if any
            if rows_to_add:
                import pandas as pd

                feature_annotation_edited = pd.concat(
                    [feature_annotation_edited, pd.DataFrame(rows_to_add)], ignore_index=True
                )

            st.session_state.feature_annotation_edited = feature_annotation_edited

            st.rerun()

    st.divider()

    # Drug Detection Tables Section
    with st.spinner("Processing..."):
        st.header("🧪 Drug Detection")
        stratified_df = st.session_state.get("stratified_df")
        stratified_df_analogs = st.session_state.get("stratified_df_analogs")

        # optional: clean up sample names
        stratified_df["Sample"] = stratified_df["Sample"].str.replace(
            r"\.mz[XM]L Peak area", "", regex=True
        )
        stratified_df_analogs["Sample"] = stratified_df_analogs["Sample"].str.replace(
            r"\.mz[XM]L Peak area", "", regex=True
        )
        comparison_tab, table_tab = st.tabs(['Comparison', 'Tables'])
        with comparison_tab:
            display_comparison_statistics()
        with table_tab:
            st.subheader("Excluding Drug Analogs")
            st.dataframe(stratified_df.style.map(highlight_yes), use_container_width=True)

            with st.expander("Show results including drug analogs"):
                st.dataframe(
                    stratified_df_analogs.style.map(highlight_yes), use_container_width=True
                )

    st.divider()

    # Drug Class Summary Section

    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_indicators

    # Drug Class Summary Section
    st.header("📈 Drug Class Summary")

    class_count_df_sorted = st.session_state.get("class_count_df_sorted")

    # Create tabs for different visualizations
    tab_plot, tab_table = st.tabs(
        [
            "🔀 UpSet Plot",
            "📊 Top Classes Table",
        ]
    )

    with tab_plot:
        st.subheader("Drug Class Co-occurrence Analysis")
        st.write(
            "This UpSet plot shows how different drug classes co-occur across samples. Each bar represents a unique combination of drug classes."
        )

        # Controls for UpSet plot
        col1, col2 = st.columns(2)
        with col1:
            n_top_classes = st.number_input(
                "Number of Top Classes for UpSet Plot",
                min_value=1,
                value=4,
                key="upset_classes_input",
                help="Select how many top drug classes to include in the UpSet plot",
            )
        with col2:
            max_samples = st.number_input(
                "Maximum Samples to Include",
                min_value=1,
                value=50,
                key="upset_samples_input",
                help="Limit the number of samples to avoid overcrowding",
            )

        try:
            # Prepare binary matrix for top classes
            top_classes = (
                class_count_df_sorted.sum(axis=0)
                .nlargest(n_top_classes + 1)
                .index.tolist()
            )
            top_classes.remove("total_matches")
            # Create binary matrix (presence/absence)
            binary_matrix = (class_count_df_sorted[top_classes] > 0).astype(int)
            binary_matrix.index.name = "Sample"

            # Limit number of samples to avoid overcrowding
            limited_matrix = binary_matrix.head(max_samples)

            # Remove samples with no detections in the selected classes
            limited_matrix = limited_matrix[limited_matrix.sum(axis=1) > 0]

            if len(limited_matrix) == 0:
                st.warning(
                    "No samples found with detections in the selected drug classes. Try adjusting the parameters."
                )
            else:
                # Create UpSet plot data
                upset_data = from_indicators(top_classes, limited_matrix.astype(bool))

                upset_fig, ax = plt.subplots(figsize=(10, 6))
                ax.set_axis_off()
                UpSet(
                    upset_data,
                    subset_size="count",
                    sort_by="cardinality",
                    show_counts=True,
                ).plot(upset_fig)
                upset_fig.suptitle(
                    f"UpSet Plot for top {n_top_classes} classes and top {max_samples} samples",
                    y=1.05,
                )
                for ax_ in upset_fig.axes:
                    ax_.grid(axis="x", visible=False)

                # Convert the figure to SVG and return as a string
                import io

                buf = io.StringIO()
                upset_fig.savefig(buf, format="svg", bbox_inches="tight")
                svg = buf.getvalue()
                buf.close()
                plt.close(upset_fig)

                _, upset_col, _ = st.columns([1, 6, 1])
                with upset_col:
                    st.image(svg, use_container_width=False)
                    st.download_button(
                        label=":material/download: Download as SVG",
                        data=svg,
                        file_name="upset_plot.svg",
                        mime="image/svg+xml",
                        key="upset_plot_download",
                    )
                    # Add interpretation help
                    with st.expander("How to interpret this UpSet plot"):
                        st.write(
                            """
                        - **Horizontal bars (left)**: Show the total number of samples containing each individual drug class
                        - **Vertical bars (bottom)**: Show the number of samples with each specific combination of drug classes
                        - **Connected dots**: Indicate which drug classes are part of each combination
                        - **Larger vertical bars**: Represent more common co-occurrence patterns
                        - **Single dots**: Show samples with only one drug class detected
                        - **Multiple connected dots**: Show samples with multiple drug classes detected together
                        """
                        )

        except Exception as e:
            st.error(f"Error creating UpSet plot: {str(e)}")
            st.info(
                "This might be due to insufficient data or missing dependencies. Make sure you have the 'upsetplot' package installed."
            )

    with tab_table:
        st.subheader("Top Detected Drug Classes")
        nlarge = st.number_input(
            "Number of Top Classes to Display",
            min_value=1,
            value=10,
            key="top_classes_input",
        )
        top_pharm_classes = class_count_df_sorted.sum(axis=0).nlargest(nlarge)
        st.dataframe(
            top_pharm_classes.reset_index().rename(
                columns={"index": "Pharmacologic Class", 0: "Total Detections"}
            ),
            use_container_width=True,
        )

        st.subheader("Drug Class Summary by Sample")
        # Clean up sample names
        class_count_df_display = class_count_df_sorted.copy()
        class_count_df_display.index = class_count_df_display.index.str.replace(
            r"\.(mzML|mzXML) Peak area", "", regex=True
        )
        class_count_df_display = (
            class_count_df_display[
                ["total_matches"] + class_count_df_display.columns.tolist()[:-1]
                ]
            .reset_index()
            .rename(columns={"index": "Sample"})
        )

        st.dataframe(class_count_df_display, use_container_width=True)

    with st.spinner("Generating Sankey plot..."):
        st.markdown("---")
        add_sankey_graph()

else:
    from welcome import welcome_page

    welcome_page()
