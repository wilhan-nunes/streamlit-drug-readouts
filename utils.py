import os
from typing import Literal

import pandas as pd
import plotly.graph_objects as go
import requests
import streamlit as st

import colorsys

def generate_colors(n, base_color=(0.55, 0.85, 0.9)):
    # base_color is in HSV, e.g., light blue
    colors = []
    for i in range(n):
        hue = (base_color[0] + i / n) % 1.0
        rgb = colorsys.hsv_to_rgb(hue, base_color[1], base_color[2])
        rgba = f'rgba({int(rgb[0] * 255)}, {int(rgb[1] * 255)}, {int(rgb[2] * 255)}, 0.8)'
        colors.append(rgba)
    return colors

@st.cache_data
def fetch_file(
        task_id: str, file_name: str, type: Literal["quant_table", "annotation_table"]
) -> str:
    """
    Fetches a file from a given task ID and loads it into a pandas DataFrame.

    :param task_id: The task ID to construct the file URL.
    :param file_name: The name of the file to fetch. Must be one of the predefined options.
    :param type: The type of file to fetch. Must be one of "quant_table" or "library_search_table".
    :returns: The path to the downloaded file.
    """
    if type == "annotation_table":
        input_url = f"https://gnps2.org/resultfile?task={task_id}&file=nf_output/library/merged_results_with_gnps.tsv"
    elif type == "quant_table":
        input_url = f"https://gnps2.org/result?task={task_id}&viewname=quantificationdownload&resultdisplay_type=task"
    response = requests.get(input_url)
    response.raise_for_status()  # Raise an error for failed requests
    os.makedirs("input_tables", exist_ok=True)
    output_file_path = f"input_tables/{file_name}"

    with open(output_file_path, "w") as f:
        f.write(response.text)

    return output_file_path


def highlight_yes(val):
    if val == "Yes":
        return "background-color: lightgreen; color: black;"
    return ""


def create_source_to_thera_sankey(feature_annotation: pd.DataFrame, top_n_pharm: int = 10):
    """
    Creates a Sankey diagram showing the flow from chemical source to therapeutical area use.
    The flow values represent the number of samples (columns with "Peak" in name) that have
    non-zero values for each source-therapeutical area combination.

    :param feature_annotation: DataFrame with feature annotations
    :param top_n_pharm: Number of top therapeutic area to include
    :return: Plotly figure object
    """
    # Prepare the data
    df = feature_annotation.copy()

    # Get all columns with "Peak" in the name (these are the sample columns)
    peak_columns = [col for col in df.columns if 'Peak' in col]

    if not peak_columns:
        st.warning("No columns with 'Peak' found in the data.")
        return None

    # Clean up chemical_source - take first value when "|" is present
    df['source_clean'] = df['chemical_source'].fillna('Unknown').apply(
        lambda x: x.split('|')[0] if pd.notna(x) and '|' in str(x) else str(x)
    )
    df['thera_clean'] = df['therapeutic_area'].fillna('Unknown')

    # Filter out rows where both are Unknown or NA
    df = df[(df['source_clean'] != 'Unknown') | (df['thera_clean'] != 'Unknown')]
    df = df[df['thera_clean'] != 'no_match']

    # Expand therapeutic_area (handle multiple values separated by |)
    df_expanded = df.copy()
    #TODO: Maybe this is not necessary?
    df_expanded['thera_clean'] = df_expanded['thera_clean'].astype(str).str.split('|')
    df_expanded = df_expanded.explode('thera_clean')
    df_expanded = df_expanded[df_expanded['thera_clean'].notna()]
    df_expanded = df_expanded[df_expanded['thera_clean'] != 'Unknown']
    df_expanded = df_expanded[df_expanded['thera_clean'] != 'NA']

    if df_expanded.empty:
        st.warning("No valid data found for Sankey diagram.")
        return None

    # Calculate sample counts for each therapeutic area to determine top N
    thera_sample_counts = {}
    for therapeutical_use in df_expanded['thera_clean'].unique():
        thera_data = df_expanded[df_expanded['thera_clean'] == therapeutical_use]
        # Count samples (columns) where at least one feature of this class has > 0 value
        sample_count = 0
        for col in peak_columns:
            if (thera_data[col] > 0).any():
                sample_count += 1
        thera_sample_counts[therapeutical_use] = sample_count

    # Get top N therapeutic area based on sample counts
    top_thera_use = sorted(thera_sample_counts.items(), key=lambda x: x[1], reverse=True)[:top_n_pharm]
    top_thera_use = [thera for thera, _ in top_thera_use]

    # Filter to only include top therapeutic area
    df_filtered = df_expanded[df_expanded['thera_clean'].isin(top_thera_use)]

    if df_filtered.empty:
        st.warning("No data available after filtering for top therapeutic areas.")
        return None

    # Create source-target pairs and count samples with non-zero values
    flow_data = []
    for source in df_filtered['source_clean'].unique():
        for thera in df_filtered['thera_clean'].unique():
            # Get data for this source-pharm combination
            subset = df_filtered[
                (df_filtered['source_clean'] == source) &
                (df_filtered['thera_clean'] == thera)
                ]

            if not subset.empty:
                # Count samples where this combination has non-zero values
                sample_count = 0
                for col in peak_columns:
                    if (subset[col] > 0).any():
                        sample_count += 1

                if sample_count > 0:
                    flow_data.append({
                        'source_clean': source,
                        'thera_clean': thera,
                        'count': sample_count
                    })

    flow_data = pd.DataFrame(flow_data)

    if flow_data.empty:
        st.warning("No valid flows found for Sankey diagram.")
        return None

    # Create unique lists of sources and targets
    sources = flow_data['source_clean'].unique().tolist()
    targets = flow_data['thera_clean'].unique().tolist()

    # Create node labels and indices
    all_nodes = sources + targets
    node_indices = {node: i for i, node in enumerate(all_nodes)}

    # Create source and target indices for the flows
    source_indices = [node_indices[source] for source in flow_data['source_clean']]
    target_indices = [node_indices[target] for target in flow_data['thera_clean']]
    values = flow_data['count'].tolist()

    # Define colors for different chemical sources
    source_colors = {
        'Medical': 'rgba(31, 119, 180, 0.8)',
        'Drug_analog': 'rgba(255, 127, 14, 0.8)',
        'Food': 'rgba(44, 160, 44, 0.8)',
        'Endogenous': 'rgba(214, 39, 40, 0.8)',
        'Drug metabolite': 'rgba(148, 103, 189, 0.8)',
        'Background/QAQC': 'rgba(140, 86, 75, 0.8)',
        'Personal Care': 'rgba(227, 119, 194, 0.8)',
        'Industrial': 'rgba(127, 127, 127, 0.8)',
        'Low_confidence': 'rgba(188, 189, 34, 0.8)',
        'Unknown': 'rgba(23, 190, 207, 0.8)'
    }

    pharm_colors = generate_colors(len(targets))
    pharm_color_map = {pharm: pharm_colors[i] for i, pharm in enumerate(targets)}

    node_colors = []
    for node in all_nodes:
        if node in sources:
            node_colors.append(source_colors.get(node, 'rgba(128, 128, 128, 0.8)'))
        else:
            node_colors.append(pharm_color_map.get(node, 'rgba(173, 216, 230, 0.8)'))
    # Assign colors to nodes
    # node_colors = []
    # for node in all_nodes:
    #     if node in sources:
    #         node_colors.append(source_colors.get(node, 'rgba(128, 128, 128, 0.8)'))
    #     else:
    #         # For therapeutic area, use a lighter version of blue
    #         node_colors.append('rgba(173, 216, 230, 0.8)')

    # Create the Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color=node_colors
        ),
        link=dict(
            source=source_indices,
            target=target_indices,
            value=values,
            color='rgba(128, 128, 128, 0.4)'
        )
    )])

    fig.update_layout(
        title_text=f"Chemical Source to Top {top_n_pharm} Therapeutic Areas<br><sub>Flow width = Number of samples with detections</sub>",
        font_size=14,
        height=600,
        margin=dict(l=40, r=40, t=40, b=40)
    )

    return fig


def add_sankey_graph():
    """
    Adds the Sankey diagram section to the Streamlit app.
    This function should be called from app.py where session state data is available.
    """
    if not st.session_state.get('run_analysis', False):
        return

    st.header("🌊 Chemical Source and Therapeutic area overview")

    # Get the feature annotation data from session state
    feature_annotation = st.session_state.get('feature_annotation')

    if feature_annotation is None:
        st.warning("No feature annotation data available. Please run the analysis first.")
        return

    # User input for number of top therapeutic area
    col1, col2 = st.columns([1, 3])

    with col1:
        top_n_pharm = st.number_input(
            "Top N Therapeutic Areas",
            min_value=5,
            max_value=25,
            value=10,
            step=1,
            help="Select how many top therapeutic areas to include in the Sankey diagram"
        )

    with col2:
        st.write("")  # Empty space for alignment

    # Create and display the Sankey diagram
    with st.spinner("Generating Sankey diagram..."):
        try:
            fig = create_source_to_thera_sankey(feature_annotation, top_n_pharm)

            if fig is not None:
                st.plotly_chart(fig, use_container_width=True)

                # Add interpretation help
                with st.expander("How to interpret this Sankey diagram"):
                    st.write("""
                    - **Left side (Chemical Sources)**: Shows the different sources of detected compounds
                    - **Right side (Therapeutic Area)**: Shows the top therapeutic areas found
                    - **Flow width**: The thickness of each flow represents the **number of samples** where this source-class combination has non-zero detections
                    - **Colors**: Different chemical sources are color-coded for easy identification

                    **Chemical Source Categories:**
                    - **Medical**: Pharmaceutical compounds
                    - **Drug_analog**: Structural analogs of known drugs
                    - **Food**: Food-derived compounds
                    - **Endogenous**: Naturally occurring in the body
                    - **Drug metabolite**: Metabolic products of drugs
                    - And others as detected in your data
                    """)

            else:
                st.info("Unable to generate Sankey diagram. This might be due to insufficient or incompatible data.")

        except Exception as e:
            st.error(f"Error generating Sankey diagram: {str(e)}")
            st.info("Please check your data format and try again.")


if __name__ == "__main__":
    # Example usage
    task_id = "34e2b1b692444bf6ae37e71dd137c300"
    file_name = "merged_results.tsv"

    lib_search_df = fetch_file(task_id, file_name)
    print(lib_search_df)  # Display the first few rows of the DataFrame
