import colorsys
import os
from typing import Literal

import pandas as pd
import plotly.graph_objects as go
import requests
import streamlit as st


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


def create_sankey_plot(feature_annotation: pd.DataFrame, top_areas: int = 5, top_class: int = 10):
    """
    Create a Sankey plot

    Parameters:
    feature_annotation (pd.DataFrame): Input dataframe with drug detection data
    """

    # Clean the name_parent_compound column (equivalent to str_trim)
    feature_annotation['name_parent_compound'] = feature_annotation['name_parent_compound'].str.strip()

    # Select columns containing name_parent_compound and file extensions
    peak_area_cols = [col for col in feature_annotation.columns if
                      'name_parent_compound' in col or
                      '.mzML' in col or
                      'mzXML' in col]
    exo_drug_peak_area = feature_annotation[peak_area_cols]

    # Group by name_parent_compound and sum
    exo_drug_group = exo_drug_peak_area.groupby('name_parent_compound').sum()

    # Calculate occurrences (rowSums > 0)
    exo_drug_detection = pd.DataFrame({
        'occurrence': (exo_drug_group > 0).sum(axis=1)
    })
    exo_drug_detection.reset_index(inplace=True)

    # Join with therapeutic_area and pharmacologic_class
    merge_cols = ['name_parent_compound', 'therapeutic_area', 'pharmacologic_class']
    exo_drug_info = feature_annotation[merge_cols].drop_duplicates().groupby('name_parent_compound').first().reset_index()

    exo_drug_detection = exo_drug_detection.merge(exo_drug_info, on='name_parent_compound', how='left')

    # Handle NA values
    exo_drug_detection = exo_drug_detection.fillna("NA")
    exo_drug_detection.loc[exo_drug_detection['therapeutic_area'] == "NA", 'therapeutic_area'] = "unspecified_area"
    exo_drug_detection.loc[
        exo_drug_detection['pharmacologic_class'] == "NA", 'pharmacologic_class'] = "unspecified_class"

    # Filter out unspecified/NA entries
    pattern = r'NA|no_match|unspecified'
    exo_drug_detection_filtered = exo_drug_detection[
        (~exo_drug_detection['therapeutic_area'].str.contains(pattern, case=False, na=False)) &
        (~exo_drug_detection['pharmacologic_class'].str.contains(pattern, case=False, na=False))
        ]

    # First level: therapeutic_area -> pharmacologic_class
    counts_first = exo_drug_detection_filtered.groupby(['therapeutic_area', 'pharmacologic_class'])[
        'occurrence'].sum().reset_index()
    counts_first.columns = ['source', 'target', 'value']

    # Get top 5 therapeutic areas by total value
    top_sources_first = counts_first.groupby('source')['value'].sum().sort_values(ascending=False).head(top_areas)
    counts_first_top5 = counts_first[counts_first['source'].isin(top_sources_first.index)]

    # Second level: pharmacologic_class -> name_parent_compound
    counts_second = exo_drug_detection_filtered.groupby(['pharmacologic_class', 'name_parent_compound'])[
        'occurrence'].sum().reset_index()
    counts_second.columns = ['source', 'target', 'value']

    # Filter to only include pharmacologic_classes from first level
    counts_second_filtered = counts_second[counts_second['source'].isin(counts_first_top5['target'])]

    # Get top selected pharmacologic classes by total value
    top_sources_second = counts_second_filtered.groupby('source')['value'].sum().sort_values(ascending=False).head(top_class)
    counts_second_top10 = counts_second_filtered[counts_second_filtered['source'].isin(top_sources_second.index)]

    # Combine both levels
    counts = pd.concat([counts_first_top5, counts_second_top10], ignore_index=True)

    # Create links dataframe
    links = pd.DataFrame({
        'source': counts['source'],
        'target': counts['target'],
        'value': counts['value']
    })

    # Create nodes dataframe
    all_nodes = pd.concat([links['source'], links['target']]).unique()
    nodes = pd.DataFrame({'name': all_nodes})

    # Convert names to indices for Plotly
    node_dict = {name: i for i, name in enumerate(nodes['name'])}
    links['source_id'] = links['source'].map(node_dict)
    links['target_id'] = links['target'].map(node_dict)

    node_colors = generate_colors(top_areas) + generate_colors(top_class) + generate_colors(len(nodes['name']) - top_class + top_areas)

    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=10,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes['name'],
            color=node_colors,
        ),
        link=dict(
            source=links['source_id'],
            target=links['target_id'],
            value=links['value']
        )
    )])

    fig.update_layout(
        title_text="Drug Detection Sankey Diagram",
        font_size=14,
        font_family="Arial",
        height=1000,
        font=dict(color="white"),
        margin=dict(l=40, r=200, t=40, b=40),
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
            value=5,
            step=1,
            help="Select how many top therapeutic areas to include in the Sankey diagram"
        )

    with col2:
        top_n_class = st.number_input(
            "Top N Pharmacologic Classes",
            min_value=5,
            max_value=25,
            value=10,
            step=1,
            help="Select how many top pharmacologic classes to include in the Sankey diagram"
        )

    # Create and display the Sankey diagram
    with st.spinner("Generating Sankey diagram..."):
        try:
            fig = create_sankey_plot(feature_annotation, top_areas=top_n_pharm, top_class=top_n_class)

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
