import pickle
import time

from streamlit.components.v1 import html
import streamlit as st
from utils import display_comparison_statistics, get_git_short_rev
from dataclasses import dataclass
from typing import Optional
import plotly.express as px
import pandas as pd
from script import *
import warnings

warnings.filterwarnings('ignore', category=FutureWarning, module='upsetplot')

page_title = "Drug Readout Analysis"
app_version = "2025-08-17"
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

# HNRC dataset example
EXAMPLE_1 = {
    "description": "HNRC (HIV) cohort samples",
    "task_id": "a322acf7936c4f91a41fd2f267d9b613",
    "blank_ids": None,
    "threshold": 50000,
}

EXAMPLE_2 = {
    "description": "Alzheimer's Disease Brain samples",
    "task_id": "4d99fc25d84143bdbbf2dd07bf044e5e",
    "blank_ids": "blank",
    "threshold": 1000,
}


@dataclass
class AnalysisData:
    """Data class to manage analysis results and reduce session state calls"""
    feature_annotation: Optional[pd.DataFrame] = None
    stratified_df: Optional[pd.DataFrame] = None
    stratified_df_analogs: Optional[pd.DataFrame] = None
    class_count_df: Optional[pd.DataFrame] = None
    class_count_df_analog: Optional[pd.DataFrame] = None
    default_excluded_features: Optional[pd.DataFrame] = None
    class_compound_dict: Optional[dict] = None
    class_compound_dict_analog: Optional[dict] = None

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
    def load_from_file(cls, file_path):
        """Load data from a pickle file and recreate with current class definition"""
        with open(file_path, 'rb') as f:
            old_data = pickle.load(f)
        
        # Recreate with current class to ensure methods are available
        data = cls(
            feature_annotation=getattr(old_data, 'feature_annotation', None),
            stratified_df=getattr(old_data, 'stratified_df', None),
            stratified_df_analogs=getattr(old_data, 'stratified_df_analogs', None),
            class_count_df=getattr(old_data, 'class_count_df', None),
            class_count_df_analog=getattr(old_data, 'class_count_df_analog', None),
            default_excluded_features=getattr(old_data, 'default_excluded_features', None),
            class_compound_dict=getattr(old_data, 'class_compound_dict', None),
            class_compound_dict_analog=getattr(old_data, 'class_compound_dict_analog', None)
        )
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
        gnps_task_id = query_params.get("task_id", "")
        threshold = query_params.get("threshold", 1000)
        blank_str = query_params.get("blank_ids", None)
        if blank_str:
            blank_str = blank_str.strip()

        load_example = st.checkbox(
            "Load Example",
            help=F"Load analysis data for FBMN task ID",
            value=False,
            key='load_example_check',
            on_change=lambda: st.session_state.pop('run_analysis', None)
        )

        if load_example:
            selected_example = st.selectbox(
                "Select Example Task ID",
                options=["Example 1", "Example 2"],
                format_func=lambda x: EXAMPLE_1['description'] if x == "Example 1" else EXAMPLE_2['description'],
                on_change=lambda: st.session_state.pop('run_analysis', None)
            )
            example = EXAMPLE_1 if selected_example == "Example 1" else EXAMPLE_2
            gnps_task_id = example["task_id"]
            threshold = example["threshold"]
            blank_str = example["blank_ids"]

        task_id = st.text_input(
            f":green-badge[Task ID] FBMN Workflow Task ID (GNPS2)",
            help="Enter the Task ID from a FBMN Workflow to retrieve the result files.",
            placeholder="enter task ID...",
            value=gnps_task_id,
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
            value=blank_str,
            placeholder="Example: BLANK|IS|PoolQC|QCmix|SRM",
            help="Enter substrings to identify blank or control columns, separated by '|'. If given, the table will be filtered to remove these columns from the analysis. If not provided, all columns will be considered.",
        )

        blank_ids = blank_ids.strip() if blank_ids else None

        if not task_id:
            st.warning("Please enter a Task ID from a FBMN Workflow to proceed.")

        run_analysis = st.button(
            "Run Analysis",
            icon="🏁",
            help="Click to start the analysis with the provided inputs.",
            width='stretch',
            key="run_analysis_button",
            disabled=not task_id
        )

        if st.button(
                "Restart Session",
                icon="♻️",
                key="restart_session",
                width='stretch',
                type="primary"
        ):
            st.session_state.clear()
            st.session_state.update({'load_example_check': False})
            st.rerun()

        # Contributors and documentation sections
        st.subheader("Contributors")
        st.markdown("""
        - [Haoqi (Nina) Zhao PhD](https://scholar.google.com/citations?user=xW9jBO0AAAAJ) - UC San Diego
        - [Kine Eide Kvitne PhD](https://scholar.google.com/citations?user=0nP4ie4AAAAJ) - UC San Diego
        - [Wilhan Nunes PhD](https://scholar.google.com/citations?user=4cPVoeIAAAAJ) - UC San Diego
        - [Mingxun Wang PhD](https://www.cs.ucr.edu/~mingxunw/) - UC Riverside
        """)

        st.subheader("Documentations and Resources")
        st.markdown(
            "[Feature Based Molecular Networking](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/)")

    return {
        'load_precomputed_example': load_example,
        'task_id': task_id,
        'intensity_thresh': intensity_thresh,
        'blank_ids': blank_ids,
        'run_analysis': run_analysis
    }


def load_data(task_id):
    """Load and process initial data"""
    from utils import fbmn_quant_download_wrapper, fbmn_lib_download_wrapper
    with st.spinner("Downloading Task result files..."):
        st.session_state.quant_file_df = fbmn_quant_download_wrapper(task_id)
        st.session_state.annotation_file_df = fbmn_lib_download_wrapper(task_id)

        st.markdown(
            f"[:material/download: Input Quant table]"
            f"(https://gnps2.org/result?task={task_id}&viewname=quantificationdownload&resultdisplay_type=task) | "
            f"[:material/download: Input Annotation table]"
            f"(https://gnps2.org/resultfile?task={task_id}&file=nf_output/library/merged_results_with_gnps.tsv)",
        )


def process_analysis_data(quant_file_df, annotation_file_df, config, data: AnalysisData):
    """Process the main analysis data and return updated data object"""
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
            _feature_annotation, _excluded_features = generate_feature_annotation(_annotation_metadata,
                                                                                  feature_filtered)
        else:
            print('[process_analysis_data] Rerun - reading feature_annotation from session...')
            _feature_annotation = data.feature_annotation
            _excluded_features = data.default_excluded_features
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
        _class_count_df, _class_count_df_analog, _class_compounds_dict, _class_compounds_dict_analog = count_drug_class_occurrences(
            _feature_annotation, class_column="pharmacologic_class"
        )
        _class_count_df["total_matches"] = _class_count_df.sum(axis=1)
        _class_count_df_analog["total_matches"] = _class_count_df.sum(axis=1)

        # Update data object with computed values
        updates = {
            'feature_annotation': _feature_annotation,
            'stratified_df': _stratified_df,
            'stratified_df_analogs': _stratified_df_analogs,
            'class_count_df': _class_count_df,
            'class_count_df_analog': _class_count_df_analog,
            'default_excluded_features': _excluded_features,
            'class_compound_dict': _class_compounds_dict,
            'class_compound_dict_analog': _class_compounds_dict_analog
        }
        
        for field, value in updates.items():
            setattr(data, field, value)
        data.save_to_session()

        print('[process_analysis_data] saved to session...')
        return data


@st.cache_data
def cached_process_analysis(quant_file_df, annotation_file_df, config, data):
    """Cached wrapper for process_analysis_data to improve performance on repeated runs"""
    return process_analysis_data(quant_file_df, annotation_file_df, config, data)

def run_analysis_wrapper(task_id, data):
    """Run analysis and return updated data object"""
    if task_id in [EXAMPLE_1["task_id"], EXAMPLE_2["task_id"]]:
        return cached_process_analysis(
            st.session_state.quant_file_df,
            st.session_state.annotation_file_df,
            config,
            data
        )
    else:
        return process_analysis_data(
            st.session_state.quant_file_df,
            st.session_state.annotation_file_df,
            config,
            data,
        )


@st.fragment
def display_summary_statistics(data: AnalysisData):
    """Display drug detection summary statistics"""
    st.header("📊 Drug Detection Summary Statistics")
    
    # Safety check: ensure data is properly loaded
    if data.stratified_df is None:
        st.error("Data not properly loaded. Please rerun the analysis.")
        return
    
    tab1, tab2 = st.tabs(["📈 Summary Metrics", "📊 Top Pharmacological Classes"])
    with tab1:
        with st.spinner("Calculating summary..."):
            sample_count = len(data.stratified_df)

            # Define specific drug categories
            specific_categories = {
                "antibiotics": ["🦠 Antibiotics",
                                "Drugs with pharmacologic_class containing 'microbial', 'bacterial', or 'tetracycline'"],
                "antidepressants": ["🧠 Antidepressants", "Drugs with therapeutic_indication containing 'depression'"],
                "statin": ["💊 Statins",
                           "Drugs with pharmacologic_class containing 'statin' or 'HMG-CoA reductase inhibitor'"],
                "PPI": ["➕ PPIs (Proton Pump Inhibitors)",
                        "Drugs with pharmacologic_class containing 'proton pump inhibitor'"],
                "antihistamine": ["🤧 Antihistamines", "Drugs with pharmacologic_class containing 'histamine-1'"],
                "antihypertensive": ["❤️ Antihypertensives",
                                     "Drugs with therapeutic_indication containing 'hypertension' AND therapeutic_area containing 'cardiology'"],
                "Alzheimer": ["🧠 Alzheimer's Meds", "Drugs with therapeutic_indication containing 'Alzheimer'"],
                "antifungal": ["🍄 Antifungals",
                               "Drugs with pharmacologic_class containing 'antifungal' OR therapeutic_indication containing 'fungal infection'"],
                "HIVmed": ["🏥 HIV Medications",
                           "Drugs with therapeutic_indication containing 'HIV' OR therapeutic_indication containing 'atazanavir'"],
            }

            # Calculate counts and percentages for specific categories
            category_stats = {}
            for col_name, display_name in specific_categories.items():
                if col_name in data.stratified_df.columns:
                    count = (data.stratified_df[col_name] == "Yes").sum()
                    pct = (count / sample_count) * 100 if sample_count > 0 else 0
                    category_stats[display_name[0]] = {"count": count, "percentage": pct,
                                                       "description": display_name[1]}
                    # sort

            # Display specific drug categories in a grid
            if category_stats:
                st.subheader("🎯 Specific Drug Categories", help="Classes with a long dash (—) indicate no detection.")

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
                                f"{stats['count']} ({stats['percentage']:.1f}%)" if stats['percentage'] > 0 else None,
                                help=stats['description'],
                            )

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
    with tab2:
        col1, col2 = st.columns([1, 2])
        with col1:
            upset_analog_inclusion = st.radio("Drug Analogs", options=['Include', 'Exclude'], index=1, horizontal=True)
        with col2:
            nlarge = st.number_input("Number of Top Classes to Display", min_value=1, value=20,
                                     key="top_classes_summary")
        upset_class_count = data.class_count_df_analog if upset_analog_inclusion == "Include" else data.class_count_df
        upset_class_count = upset_class_count.sort_values("total_matches", ascending=False)
        top_pharm_classes = ((upset_class_count > 0).astype(int)).sum(axis=0).nlargest(nlarge)
        total_matches = top_pharm_classes.iloc[0]
        top_pharm_classes.reset_index().rename(
            columns={"class_group": "Pharmacologic Class", 0: "Number of samples containing this class"}
        )
        percentages = (top_pharm_classes / total_matches * 100).round(1)

        compounds_dict = data.class_compound_dict_analog if upset_analog_inclusion == "Include" else data.class_compound_dict
        compounds_list = [compounds_dict.get(pharm_class, []) for pharm_class in top_pharm_classes.index]
        fig = px.bar(
            x=top_pharm_classes.values,
            y=top_pharm_classes.index,
            orientation='h',
            title=f"Top {nlarge} Pharmacologic Classes by Sample Count",
            labels={"x": "Number of Samples", "y": "Pharmacologic Class"},
            hover_data={"compounds": compounds_list},
        )

        fig.update_traces(
            hovertemplate='<b>%{y}</b><br>' +
                          'Sample Count: %{x}<br>' +
                          '<b>Compounds</b>: %{customdata}<extra></extra>'
        )
        # Add percentage text inside bars (centered)
        for i, (count, pct) in enumerate(zip(top_pharm_classes.values, percentages.values)):
            fig.add_annotation(
                x=count / 2,
                y=i,
                text=f"{pct}%",
                showarrow=False,
                font=dict(color="white", size=11, family="Arial")
            )

        fig.update_layout(height=600, margin=dict(l=200), yaxis=dict(autorange="reversed"))
        st.plotly_chart(fig, use_container_width=True)

        download_df = ((upset_class_count > 0).astype(int)).sum(axis=0).rename('counts').to_frame()
        download_df['compounds'] = [";".join(compounds_dict.get(pharm_class, [])) for pharm_class in download_df.index]

        st.download_button(
            label=":material/download: Download data",
            data=download_df.to_csv(sep="\t"),
            file_name="top_pharm_classes.tsv",
            mime="text/tab-separated-values",
            key="top_parhm_class_table_download",
        )


@st.fragment
def display_feature_annotation_table(data: AnalysisData):
    """Display and handle feature annotation table editing"""
    from utils import add_df_and_filtering, conditional_highlighter_low_confidence

    st.header("🔬 Feature Annotation Table")
    col1, col2, col3 = st.columns(3)
    with col1:
        with st.expander(
            "**Edit & Rerun**"):
            st.markdown(
            "You can edit the table below and then rerun the analysis with your modifications.\n\n"
            "[:material/help: Learn how](https://www.youtube.com/watch?v=6tah69LkfxE&list=TLGGKK4Dnf1gepcwNTA2MjAyNQ)"
        )
    with col2:
        with st.expander(
            "**Before Editing**"):
            st.markdown(
            "Please clear all filters before editing the data."
        )
    with col3:
        with st.expander(
            "**Report Issues**", 
            icon=":material/report:"):
            st.markdown(
            f"[Report an annotation issue]({repo_link}/issues/new?assignees=&labels=bug&template=bug_report.md&title=Feature+Annotation+Issue)"
        )
    
    data_df = data.load_from_session().feature_annotation.copy()

    filtered_df = add_df_and_filtering(data_df, "feature_annotation_filtered")

    # Apply styling only if the filtered dataframe is small enough
    edited_df = st.data_editor(
        conditional_highlighter_low_confidence(filtered_df),
        key="feature_annotation_editor",
        width='stretch',
        num_rows="dynamic",
        height=400,
        disabled=["CosineScore", "MatchedPeaks"]
    )

     # Check size before styling the display dataframe
    
    low_confidence_count = ((data.feature_annotation['CosineScore'].astype(float) < 0.9) & (data.feature_annotation['MatchedPeaks'].astype(int) <= 2)).sum()
    
    if data.feature_annotation.size <= 262144:
        st.warning(
            f":red-badge[**{low_confidence_count} Low confidence annotations**]: Features with cosine score < 0.9, matched peaks <= 2 are highlighted in red.")
    else:
        st.warning(f":red-badge[**{low_confidence_count} Low confidence annotations**]: Features with cosine score < 0.9, matched peaks <= 2 are **not highlighted** (dataframe too large). Please, inspect manually.")


    # Rerun button
    col1, col2 = st.columns([1, 1])
    with col1:
        if st.button(
            "Rerun Removing Low Confidence",
            width='stretch',
            icon=':material/replay:',
            disabled=(low_confidence_count == 0),
            help="No low confidence annotations to remove." if low_confidence_count == 0 else "Click to remove low confidence annotations (cosine score < 0.9 and matched peaks <= 2) and rerun the analysis."
            ):
            data.feature_annotation = data.feature_annotation[~((data.feature_annotation['CosineScore'].astype(float) < 0.9) & (data.feature_annotation['MatchedPeaks'].astype(int) <= 2))]
            data.save_to_session()
            st.session_state['rerun_analysis'] = True
            st.rerun()
    with col2:
        contains_filter = st.session_state.get('feature_annotation_filtered_filter_count', 0) > 0
        if st.button(
                "🔄 Rerun Analysis with Edited Data",
                width='stretch',
                type="primary",
                key="rerun_analysis_button",
                disabled=(st.session_state.get('feature_annotation_filtered_filter_count', 0) > 0),
                help="Click to rerun the analysis using the edited feature annotation table." if not contains_filter else "Please :red-badge[clear all filters before editing] and rerunning the analysis."
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
    from utils import conditional_highlighter_yes

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
        #Conditionally highlight "Yes" values in the tables (if df is small enough)
        st.subheader("Excluding Drug Analogs")
        st.dataframe(conditional_highlighter_yes(stratified_df_clean), width='stretch')

        with st.expander("Show results including drug analogs"):
            st.dataframe(conditional_highlighter_yes(stratified_df_analogs_clean), width='stretch')


@st.fragment
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


@st.fragment
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
            st.image(svg, width='content')
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
    st.subheader("Drug Class Summary by Sample")
    class_count_df_display = upset_class_count.copy()
    class_count_df_display.index = class_count_df_display.index.str.replace(r"\.(mzML|mzXML) Peak area", "", regex=True)
    class_count_df_display = (
        class_count_df_display[["total_matches"] + class_count_df_display.columns.tolist()[:-1]]
        .reset_index()
        .rename(columns={"index": "Sample"})
    )
    st.dataframe(class_count_df_display, width='stretch')


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
            load_data(config['task_id'])

            data = AnalysisData()
            # Process analysis
            data = run_analysis_wrapper(config["task_id"], data)
            data.save_to_session()
        else:
            print('[main] Rerun - Reading data from session state')
            data = AnalysisData.load_from_session()
            # Process analysis with updated data
            data = run_analysis_wrapper(config["task_id"], data)

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

    with st.spinner("Rendering Feature Annotation Table..."):
        time.sleep(0.5)
    display_feature_annotation_table(data)
    st.divider()


    with st.spinner("Rendering Drug Detection Tables..."):
        time.sleep(0.5)
    display_drug_detection_tables(data)
    st.divider()


    with st.spinner(("Rendering Drug Class Summary Table...")):
        time.sleep(0.5)
    display_drug_class_summary(data)

    # Add Sankey graph
    st.markdown("---")

    with st.spinner("Rendering Sankey Plot..."):
        time.sleep(0.5)
    add_sankey_graph(data.feature_annotation)
    print('[main] All visualizations rendered')

else:
    with st.spinner("Rendering Welcome page..."):
        time.sleep(0.5)
    print('Welcome Page Loaded')
    from welcome import welcome_page

    welcome_page()
