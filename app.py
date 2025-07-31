from gnpsdata import workflow_fbmn
from streamlit.components.v1 import html
import streamlit as st
from utils import display_comparison_statistics, load_example, get_git_short_rev
from dataclasses import dataclass
from typing import Optional
import pandas as pd
from script import *
import warnings

warnings.filterwarnings('ignore', category=FutureWarning, module='upsetplot')

# Set page configuration
page_title = "Drug Readout Analysis"
app_version = "2025-07-31"
git_hash = get_git_short_rev()
repo_link = "https://github.com/wilhan-nunes/streamlit_drug_readouts"

st.set_page_config(
    page_title=page_title,
    layout="wide",
    page_icon=":pill:",
    menu_items={"About": (f"**App version**: {app_version} | "
                          f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})")},
)

# Add tracking token
html(
    '<script async defer data-website-id="74bc9983-13c4-4da0-89ae-b78209c13aaf" src="https://analytics.gnps2.org/umami.js"></script>',
    width=0,
    height=0,
)


@dataclass
class AnalysisData:
    """Data class to manage analysis results and reduce session state calls"""
    feature_annotation: Optional[pd.DataFrame] = None
    stratified_df: Optional[pd.DataFrame] = None
    stratified_df_analogs: Optional[pd.DataFrame] = None
    class_count_df: Optional[pd.DataFrame] = None
    class_count_df_analog: Optional[pd.DataFrame] = None
    default_excluded_features: Optional[pd.DataFrame] = None

    def save_to_session(self):
        """Save all data to session state"""
        for field_name, field_value in self.__dict__.items():
            if field_value is not None:
                st.session_state[field_name] = field_value

    @classmethod
    def load_from_session(cls):
        """Load data from session state"""
        data = cls()
        for field_name in data.__dict__.keys():
            if field_name in st.session_state:
                setattr(data, field_name, st.session_state[field_name])
        return data

    @classmethod
    def data_summary(cls):
        data = cls()
        for field_name in data.__dict__.keys():
            if field_name in st.session_state:
                print(field_name, "(shape):", st.session_state[field_name].shape)
            else:
                print(field_name, "not initialized")


def setup_sidebar():
    """Setup sidebar inputs and return configuration values"""
    with st.sidebar:
        st.header("Inputs")

        # Get query params
        query_params = st.query_params
        gnps_task_id = query_params.get("taskid", "")
        threshold = query_params.get("threshold", 1000)
        blank_str = query_params.get("blank_ids", None)

        load_example_data = st.checkbox(
            "Load example",
            help="Load example from FBMN task ID d6f37a11d90c4f249974280c3fc90108",
            value=False,
            key='load_example_check'
        )

        task_id = st.text_input(
            f":green-badge[Task ID] FBMN Workflow Task ID (GNPS2)",
            help="Enter the Task ID from a FBMN Workflow to retrieve the result files.",
            placeholder="enter task ID...",
            value=gnps_task_id if not st.session_state.get(
                'load_example_check') else "d6f37a11d90c4f249974280c3fc90108",
            disabled=(load_example_data == True)
        )

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
            st.warning("Please enter a Task ID from a FBMN Workflow to proceed.")

        run_analysis = st.button(
            "Run Analysis",
            icon="🏁",
            help="Click to start the analysis with the provided inputs.",
            use_container_width=True,
            key="run_analysis_button",
            disabled=not (task_id or load_example_data)
        )

        if st.button(
                "Restart Session",
                icon="♻️",
                key="restart_session",
                use_container_width=True,
                type="primary",
        ):
            st.session_state.clear()
            st.session_state['load_example_data'] = False
            st.rerun()

        # Contributors and documentation sections
        st.subheader("Contributors")
        st.markdown("""
        - [Haoqi (Nina) Zhao PhD](https://scholar.google.com/citations?user=xW9jBO0AAAAJ) - UC San Diego
        - [Wilhan Nunes PhD](https://scholar.google.com/citations?user=4cPVoeIAAAAJ) - UC San Diego
        - [Mingxun Wang PhD](https://www.cs.ucr.edu/~mingxunw/) - UC Riverside
        """)

        st.subheader("Documentations and Resources")
        st.markdown(
            "[Feature Based Molecular Networking](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/)")

    return {
        'load_example_data': load_example_data,
        'task_id': task_id,
        'intensity_thresh': intensity_thresh,
        'blank_ids': blank_ids,
        'run_analysis': run_analysis
    }


def load_data(config):
    """Load and process initial data"""
    from utils import fbmn_quant_download_wrapper, fbmn_lib_download_wrapper
    if not config['load_example_data']:
        with st.spinner("Downloading Task result files..."):
            st.session_state.quant_file_df = fbmn_quant_download_wrapper(config['task_id'])
            st.session_state.annotation_file_df = fbmn_lib_download_wrapper(config['task_id'])

            st.markdown(
                f"[:material/download: Input Quant table]"
                f"(https://gnps2.org/result?task={config['task_id']}&viewname=quantificationdownload&resultdisplay_type=task) | "
                f"[:material/download: Input Annotation table]"
                f"(https://gnps2.org/resultfile?task={config['task_id']}&file=nf_output/library/merged_results_with_gnps.tsv)",
            )
    else:
        st.session_state.quant_file_df, st.session_state.annotation_file_df = load_example()


def process_analysis_data(quant_file_df, annotation_file_df, config, data: AnalysisData):
    """Process the main analysis data"""
    print('[process_analysis_data] triggered...')
    _drug_metadata_file = "data/GNPS_Drug_Library_Metadata_Drugs.csv"
    _analog_metadata_file = "data/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv"

    with st.spinner("Processing data..."):
        subtract_blanks = True if config['blank_ids'] else False
        feature_filtered = load_and_filter_features(
            quant_file_df,
            intensity_threshold=config['intensity_thresh'],
            blank_ids=config['blank_ids'],
            subtract_blanks=subtract_blanks,
        )

        _annotation_metadata = load_and_merge_annotations(
            annotation_file_df, _drug_metadata_file, _analog_metadata_file
        )

        if not st.session_state.get("rerun_analysis", False):
            print('[process_analysis_data] generating feature_annotation...')
            _feature_annotation, excluded_features = generate_feature_annotation(_annotation_metadata, feature_filtered)
        else:
            print('[process_analysis_data] Rerun - reading feature_annotation from session...')
            _feature_annotation = data.feature_annotation
            excluded_features = data.default_excluded_features
            st.success(f'Feature annotation retrieved from session state. Shape: {_feature_annotation.shape}')
            st.session_state["rerun_analysis"] = False

        # Perform analysis
        _stratified_df = stratify_by_drug_class(
            _feature_annotation, exclude_analogs=True, peak_threshold=config['intensity_thresh'],
        )
        _stratified_df_analogs = stratify_by_drug_class(
            _feature_annotation, exclude_analogs=False, peak_threshold=config['intensity_thresh'],
        )

        # Count drug class occurrences
        _class_count_df, _class_count_df_analog = count_drug_class_occurrences(
            _feature_annotation, class_column="pharmacologic_class"
        )
        _class_count_df["total_matches"] = _class_count_df.sum(axis=1)
        _class_count_df_analog["total_matches"] = _class_count_df.sum(axis=1)

        # Update data object
        data.feature_annotation = _feature_annotation
        data.stratified_df = _stratified_df
        data.stratified_df_analogs = _stratified_df_analogs
        data.class_count_df = _class_count_df
        data.class_count_df_analog = _class_count_df_analog
        data.default_excluded_features = excluded_features
        data.save_to_session()
        print('[process_analysis_data]  saved to session...')


def display_summary_statistics(data: AnalysisData):
    """Display drug detection summary statistics"""
    st.header("📊 Drug Detection Summary Statistics")

    with st.spinner("Calculating summary..."):
        sample_count = len(data.stratified_df)

        # Define specific drug categories
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

        # Calculate category statistics
        category_stats = {}
        for col_name, display_name in specific_categories.items():
            if col_name in data.stratified_df.columns:
                count = (data.stratified_df[col_name] == "Yes").sum()
                pct = (count / sample_count) * 100 if sample_count > 0 else 0
                if pct != 0:
                    category_stats[display_name] = {"count": count, "percentage": pct}

        # Display specific drug categories
        if category_stats:
            # st.subheader("🎯 Specific Drug Categories")
            cols_per_row = 3
            num_categories = len(category_stats)

            for i in range(0, num_categories, cols_per_row):
                cols = st.columns(cols_per_row)
                row_items = list(category_stats.items())[i: i + cols_per_row]
                for j in range(cols_per_row):
                    if j < len(row_items):
                        display_name, stats = row_items[j]
                        with cols[j]:
                            st.metric(
                                display_name,
                                f"{stats['count']} ({stats['percentage']:.1f}%)",
                                help=f"Number of samples containing {display_name.split(' ', 1)[1].lower()}",
                            )
                    else:
                        with cols[j]:
                            st.empty()

        # Overall summary metrics
        st.subheader("📈 Overall Summary")
        col1, col2, col3 = st.columns(3)

        with col1:
            total_drug_columns = len([col for col in data.stratified_df.columns if col != "Sample"])
            st.metric("Total Drug Categories", total_drug_columns)

        with col2:
            drug_columns = [col for col in data.stratified_df.columns if col != "Sample"]
            samples_with_drugs = (data.stratified_df[drug_columns] == "Yes").any(axis=1).sum()
            drug_detection_pct = (samples_with_drugs / sample_count) * 100 if sample_count > 0 else 0
            st.metric("Samples with Any Drug", f"{samples_with_drugs} ({drug_detection_pct:.1f}%)")

        with col3:
            st.metric("Total Samples", sample_count)


def display_feature_annotation_table(data: AnalysisData):
    """Display and handle feature annotation table editing"""
    from utils import add_df_and_filtering

    st.header("🔬 Feature Annotation Table")
    st.write(
        "You can edit the table below and then rerun the analysis with your modifications. "
        "[Learn how](https://www.youtube.com/watch?v=6tah69LkfxE&list=TLGGKK4Dnf1gepcwNTA2MjAyNQ)"
    )
    st.warning(
        "***Before editing** the data, please clear all filters. Filters display results based on the original, unedited table."
    )

    filtered_df = add_df_and_filtering(data.feature_annotation, "feature_annotation_filtered")
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

    # Show excluded features
    with st.expander('Features excluded from analysis by default', icon="⛔️"):
        st.info("Excluded features with chemical source 'Background', 'confidence', 'Endogenous' and/or 'Food'")
        excluded_features = data.default_excluded_features
        st.dataframe(excluded_features)

    # Handle reanalysis with edited data
    if st.session_state.get('rerun_analysis', False):
        data.feature_annotation = edited_df
        data.save_to_session()
        st.rerun()


def display_drug_detection_tables(data: AnalysisData):
    """Display drug detection tables"""
    from utils import highlight_yes

    st.header("🧪 Drug Detection")

    # Clean up sample names
    stratified_df_clean = data.stratified_df.copy()
    stratified_df_analogs_clean = data.stratified_df_analogs.copy()

    stratified_df_clean["Sample"] = stratified_df_clean["Sample"].str.replace(r"\.mz[XM]L Peak area", "", regex=True)
    stratified_df_analogs_clean["Sample"] = stratified_df_analogs_clean["Sample"].str.replace(r"\.mz[XM]L Peak area",
                                                                                              "", regex=True)

    comparison_tab, table_tab = st.tabs(['Comparison', 'Tables'])

    with comparison_tab:
        display_comparison_statistics(data)

    with table_tab:
        st.subheader("Excluding Drug Analogs")
        st.dataframe(stratified_df_clean.style.map(highlight_yes), use_container_width=True)

        with st.expander("Show results including drug analogs"):
            st.dataframe(stratified_df_analogs_clean.style.map(highlight_yes), use_container_width=True)


def display_drug_class_summary(data: AnalysisData):
    """Display drug class summary with UpSet plot and tables"""
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_indicators

    st.header("📈 Drug Class Summary")

    tab_plot, tab_table = st.tabs(["🔀 UpSet Plot", "📊 Top Classes Table"])

    with tab_plot:
        st.subheader("Drug Class Co-occurrence Analysis")
        st.write(
            "This UpSet plot shows how different drug classes co-occur across samples. Each bar represents a unique combination of drug classes.")

        # Controls for UpSet plot
        col1, col2, col3 = st.columns([1, 2, 2])
        with col1:
            upset_analog_inclusion = st.radio("Drug Analogs", options=['Include', 'Exclude'], index=1)
            upset_class_count = data.class_count_df_analog if upset_analog_inclusion == "Include" else data.class_count_df
            upset_class_count = upset_class_count.sort_values("total_matches", ascending=False)
        with col2:
            n_top_classes = st.number_input(
                f"Number of Top Classes for UpSet Plot :green-badge[{len(upset_class_count.columns) - 1} total]",
                min_value=1,
                value=4,
                key="upset_classes_input",
                help="Select how many top drug classes to include in the UpSet plot",
            )
        with col3:
            sample_count = len(data.stratified_df)
            max_samples = st.number_input(
                f"Maximum Samples to Include :green-badge[{sample_count} total]",
                min_value=1,
                value=sample_count,
                key="upset_samples_input",
                help="Limit the number of samples to avoid overcrowding",
            )

        upset_fig, msg = create_upset_plot(upset_class_count, n_top_classes, max_samples, upset_analog_inclusion)
        display_upset_plot(upset_fig, msg)

    with tab_table:
        display_drug_class_tables(upset_class_count, upset_analog_inclusion)


@st.cache_data
def create_upset_plot(upset_class_count, n_top_classes, max_samples, upset_analog_inclusion):
    """Create UpSet plot and return figure and metadata"""
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_indicators

    try:
        # Prepare binary matrix for top classes
        top_classes = upset_class_count.sum(axis=0).nlargest(n_top_classes + 1).index.tolist()
        top_classes.remove("total_matches")

        # Create binary matrix (presence/absence)
        binary_matrix = (upset_class_count[top_classes] > 0).astype(int)
        binary_matrix.index.name = "Sample"

        # Limit number of samples
        limited_matrix = binary_matrix.head(max_samples)
        limited_matrix = limited_matrix[limited_matrix.sum(axis=1) > 0]

        if len(limited_matrix) == 0:
            return None, "No samples found with detections in the selected drug classes. Try adjusting the parameters."

        # Create UpSet plot
        upset_data = from_indicators(top_classes, limited_matrix.astype(bool))

        upset_fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_axis_off()
        UpSet(upset_data, subset_size="count", sort_by="cardinality", show_counts=True).plot(upset_fig)
        upset_fig.suptitle(
            f"UpSet Plot for top {n_top_classes} classes and {max_samples} of {len(binary_matrix)} samples\n(Analogs {upset_analog_inclusion}d)",
            y=1.05,
        )

        for ax_ in upset_fig.axes:
            ax_.grid(axis="x", visible=False)

        return upset_fig, None

    except Exception as e:
        return None, f"Error creating UpSet plot: {str(e)}"


def display_upset_plot(upset_fig, error_message=None):
    """Display UpSet plot in Streamlit interface"""
    import matplotlib.pyplot as plt
    import io

    if error_message:
        if "No samples found" in error_message:
            st.warning(error_message)
        else:
            st.error(error_message)
            st.info(
                "This might be due to insufficient data or missing dependencies. Make sure you have the 'upsetplot' package installed.")
        return

    if upset_fig is None:
        st.error("No plot to display")
        return

    try:
        # Convert to SVG
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

            with st.expander("How to interpret this UpSet plot"):
                st.write("""
                - **Horizontal bars (left)**: Show the total number of samples containing each individual drug class
                - **Vertical bars (bottom)**: Show the number of samples with each specific combination of drug classes
                - **Connected dots**: Indicate which drug classes are part of each combination
                - **Larger vertical bars**: Represent more common co-occurrence patterns
                - **Single dots**: Show samples with only one drug class detected
                - **Multiple connected dots**: Show samples with multiple drug classes detected together
                """)

    except Exception as e:
        st.error(f"Error displaying plot: {str(e)}")


def display_drug_class_tables(upset_class_count, upset_analog_inclusion):
    """Display drug class summary tables"""
    st.subheader(f"Top Detected Drug Classes (Analogs {upset_analog_inclusion}d)")

    nlarge = st.number_input("Number of Top Classes to Display", min_value=1, value=10, key="top_classes_input")
    top_pharm_classes = ((upset_class_count > 0).astype(int)).sum(axis=0).nlargest(nlarge)

    st.dataframe(
        top_pharm_classes.reset_index().rename(
            columns={"class_group": "Pharmacologic Class", 0: "Number of samples containing this class"}
        ),
        use_container_width=True,
    )

    st.subheader("Drug Class Summary by Sample")
    class_count_df_display = upset_class_count.copy()
    class_count_df_display.index = class_count_df_display.index.str.replace(r"\.(mzML|mzXML) Peak area", "", regex=True)
    class_count_df_display = (
        class_count_df_display[["total_matches"] + class_count_df_display.columns.tolist()[:-1]]
        .reset_index()
        .rename(columns={"index": "Sample"})
    )
    st.dataframe(class_count_df_display, use_container_width=True)


###############################
########## MAIN FLOW ##########
###############################

# Initialize session state
if "run_analysis" not in st.session_state:
    st.session_state.run_analysis = False

from utils import fetch_file, highlight_yes, add_sankey_graph

st.title("Drug Readout Analysis")

# Setup sidebar and get configuration
config = setup_sidebar()

# Run analysis if requested
if config['run_analysis'] or st.session_state.get("rerun_analysis", False):
    print('[main] Run triggered, processing data')
    try:
        # Load data
        if not st.session_state.get("rerun_analysis"):
            print('[main] First run - Reading tables from files and downloading')
            load_data(config)
            data = AnalysisData()
        else:
            data = AnalysisData.load_from_session()

        # Process analysis
        process_analysis_data(st.session_state.quant_file_df, st.session_state.annotation_file_df, config, data)

        st.session_state.run_analysis = True

    except Exception as e:
        st.error(f"An error occurred: {e}")
        raise
else:
    st.info(":information_source: Please, provide the inputs, then click Run Analysis.")

# Display results if analysis has been run
if st.session_state.run_analysis:
    print('[main] Loading data from session')
    data = AnalysisData.load_from_session()
    ## Debug
    # data.data_summary()
    # data.feature_annotation.to_csv('debug_feature_annotation.csv')
    # st.session_state.feature_annotation.to_csv('ss_debug_feature_annotation.csv')

    print('[main] Rendering visualizations')
    display_summary_statistics(data)
    st.divider()

    display_feature_annotation_table(data)
    st.divider()

    display_drug_detection_tables(data)
    st.divider()

    display_drug_class_summary(data)

    # Add Sankey graph
    with st.spinner("Generating Sankey plot..."):
        st.markdown("---")
        add_sankey_graph(data.feature_annotation)

else:
    print('Welcome Page Loaded')
    from welcome import welcome_page

    welcome_page()
