import colorsys
import os
from typing import Literal, List

import pandas as pd
import plotly.graph_objects as go
import requests
import streamlit as st
import plotly.express as px


def display_comparison_statistics():
    """Display comparison between with/without analogs if both are available."""

    stratified_df = st.session_state.get("stratified_df")
    stratified_df_analogs = st.session_state.get("stratified_df_analogs")

    if stratified_df is None or stratified_df_analogs is None:
        st.info("Run analysis with both analog options to see comparison.")
        return

    # Compare detection rates
    specific_categories = ["antibiotics", "antidepressants", "antihistamine",
                           "antihypertensive", "Alzheimer", "antifungal", "HIVmed"]

    all_categories = [name for name in stratified_df.columns.to_list() if name != 'Sample']
    selected_categories = st.multiselect(
        "Select categories to compare:",
        options=all_categories,
        default=["antibiotics", "antidepressants", "antihistamine",
                 "antihypertensive", "Alzheimer", "antifungal", "HIVmed"]
    )

    comparison_data = []
    sample_count = len(stratified_df)

    for category in selected_categories:
        if category in stratified_df.columns and category in stratified_df_analogs.columns:
            count_without = (stratified_df[category] == "Yes").sum()
            count_with = (stratified_df_analogs[category] == "Yes").sum()

            comparison_data.append({
                "Category": category.replace("_", " ").title(),
                "Parent Drugs Only": count_without,
                "With Analogs": count_with,
                "Difference": count_with - count_without,
                "% Increase": ((count_with - count_without) / max(count_without, 1)) * 100
            })

    if comparison_data:
        comparison_df = pd.DataFrame(comparison_data)

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("📊 Detection Comparison")
            st.dataframe(comparison_df, use_container_width=True, hide_index=True)

        with col2:
            st.subheader("📈 Impact of Including Analogs")
            fig_comparison = px.bar(
                comparison_df,
                x="Category",
                y=["Parent Drugs Only", "With Analogs"],
                title="Drug Detection: Parent Drugs vs With Analogs",
                barmode="group"
            )
            fig_comparison.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig_comparison, use_container_width=True)


def generate_colors(n, base_color=(0.55, 0.85, 0.9)):
    # base_color is in HSV, e.g., light blue
    colors = []
    for i in range(n):
        hue = (base_color[0] + i / n) % 1.0
        rgb = colorsys.hsv_to_rgb(hue, base_color[1], base_color[2])
        rgba = (
            f"rgba({int(rgb[0] * 255)}, {int(rgb[1] * 255)}, {int(rgb[2] * 255)}, 0.8)"
        )
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

def load_example():
    quant_file_df = pd.read_csv('data/example_quant_file_d6f37a11d90c4f249974280c3fc90108.csv')
    annotation_file_df = pd.read_csv('data/example_annotation_filed6f37a11d90c4f249974280c3fc90108.tsv', sep='\t')

    return quant_file_df, annotation_file_df


def highlight_yes(val):
    if val == "Yes":
        return "background-color: lightgreen; color: black;"
    return ""


def create_sankey_plot(
    feature_annotation: pd.DataFrame, top_areas: int = 5, top_class: int = 10
):
    """
    Create a Sankey plot

    Parameters:
    feature_annotation (pd.DataFrame): Input dataframe with drug detection data
    """

    # Clean the name_parent_compound column (equivalent to str_trim)
    feature_annotation["name_parent_compound"] = feature_annotation[
        "name_parent_compound"
    ].str.strip()

    # Select columns containing name_parent_compound and file extensions
    peak_area_cols = [
        col
        for col in feature_annotation.columns
        if "name_parent_compound" in col or ".mzML" in col or "mzXML" in col
    ]
    exo_drug_peak_area = feature_annotation[peak_area_cols]

    # Group by name_parent_compound and sum
    exo_drug_group = exo_drug_peak_area.groupby("name_parent_compound").sum()

    # Calculate occurrences (rowSums > 0)
    exo_drug_detection = pd.DataFrame({"occurrence": (exo_drug_group > 0).sum(axis=1)})
    exo_drug_detection.reset_index(inplace=True)

    # Join with therapeutic_area and pharmacologic_class
    merge_cols = ["name_parent_compound", "therapeutic_area", "pharmacologic_class"]
    exo_drug_info = (
        feature_annotation[merge_cols]
        .drop_duplicates()
        .groupby("name_parent_compound")
        .first()
        .reset_index()
    )

    exo_drug_detection = exo_drug_detection.merge(
        exo_drug_info, on="name_parent_compound", how="left"
    )

    # Handle NA values
    exo_drug_detection = exo_drug_detection.fillna("NA")
    exo_drug_detection.loc[
        exo_drug_detection["therapeutic_area"] == "NA", "therapeutic_area"
    ] = "unspecified_area"
    exo_drug_detection.loc[
        exo_drug_detection["pharmacologic_class"] == "NA", "pharmacologic_class"
    ] = "unspecified_class"

    # Filter out unspecified/NA entries
    pattern = r"NA|no_match|unspecified"
    exo_drug_detection_filtered = exo_drug_detection[
        (
            ~exo_drug_detection["therapeutic_area"].str.contains(
                pattern, case=False, na=False
            )
        )
        & (
            ~exo_drug_detection["pharmacologic_class"].str.contains(
                pattern, case=False, na=False
            )
        )
    ]

    # First level: therapeutic_area -> pharmacologic_class
    counts_first = (
        exo_drug_detection_filtered.groupby(
            ["therapeutic_area", "pharmacologic_class"]
        )["occurrence"]
        .sum()
        .reset_index()
    )
    counts_first.columns = ["source", "target", "value"]

    # Get top 5 therapeutic areas by total value
    top_sources_first = (
        counts_first.groupby("source")["value"]
        .sum()
        .sort_values(ascending=False)
        .head(top_areas)
    )
    counts_first_top5 = counts_first[
        counts_first["source"].isin(top_sources_first.index)
    ]

    # Second level: pharmacologic_class -> name_parent_compound
    counts_second = (
        exo_drug_detection_filtered.groupby(
            ["pharmacologic_class", "name_parent_compound"]
        )["occurrence"]
        .sum()
        .reset_index()
    )
    counts_second.columns = ["source", "target", "value"]

    # Filter to only include pharmacologic_classes from first level
    counts_second_filtered = counts_second[
        counts_second["source"].isin(counts_first_top5["target"])
    ]

    # Get top selected pharmacologic classes by total value
    top_sources_second = (
        counts_second_filtered.groupby("source")["value"]
        .sum()
        .sort_values(ascending=False)
        .head(top_class)
    )
    counts_second_top10 = counts_second_filtered[
        counts_second_filtered["source"].isin(top_sources_second.index)
    ]

    # Combine both levels
    counts = pd.concat([counts_first_top5, counts_second_top10], ignore_index=True)

    # Create links dataframe
    links = pd.DataFrame(
        {
            "source": counts["source"],
            "target": counts["target"],
            "value": counts["value"],
        }
    )

    # Create nodes dataframe
    all_nodes = pd.concat([links["source"], links["target"]]).unique()
    nodes = pd.DataFrame({"name": all_nodes})

    # Convert names to indices for Plotly
    node_dict = {name: i for i, name in enumerate(nodes["name"])}
    links["source_id"] = links["source"].map(node_dict)
    links["target_id"] = links["target"].map(node_dict)

    node_colors = (
        generate_colors(top_areas)
        + generate_colors(top_class)
        + generate_colors(len(nodes["name"]) - top_class + top_areas)
    )

    # Create Sankey diagram
    sankey_trace = go.Sankey(
        node=dict(
            pad=10,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=[name.replace("|", "<br>") for name in nodes["name"]],
            color=node_colors,
            align="left",
        ),
        link=dict(
            source=links["source_id"], target=links["target_id"], value=links["value"]
        ),
        textfont=dict(size=18, color="black", shadow="0px -0px 2px white"),
    )

    fig = go.Figure(data=[sankey_trace])

    fig.update_layout(
        title_text="Drug Detection Sankey Diagram",
        font_size=14,
        font_family="Arial",
        height=1000,
        font=dict(color="white"),
        margin=dict(l=40, r=40, t=40, b=40),
    )

    return fig


def add_sankey_graph():
    """
    Adds the Sankey diagram section to the Streamlit app.
    This function should be called from app.py where session state data is available.
    """
    if not st.session_state.get("run_analysis", False):
        return

    st.header("🌊 Chemical Source and Therapeutic area overview")

    # Get the feature annotation data from session state
    feature_annotation = st.session_state.get("feature_annotation")

    if feature_annotation is None:
        st.warning(
            "No feature annotation data available. Please run the analysis first."
        )
        return

    # User input for number of top therapeutic area
    col1, col2 = st.columns([1, 1])

    with col1:
        top_n_areas = st.number_input(
            "Top Therapeutic Areas to show:",
            min_value=1,
            max_value=25,
            value=5,
            step=1,
            help="Select how many top therapeutic areas to include in the Sankey diagram",
        )

    with col2:
        top_n_class = st.number_input(
            "Top Pharmacologic Classes to show:",
            min_value=1,
            max_value=25,
            value=10,
            step=1,
            help="Select how many top pharmacologic classes to include in the Sankey diagram",
        )

    # Create and display the Sankey diagram
    with st.spinner("Generating Sankey diagram..."):
        try:
            fig = create_sankey_plot(
                feature_annotation, top_areas=top_n_areas, top_class=top_n_class
            )

            if fig is not None:
                st.plotly_chart(fig, use_container_width=True)
                svg_bytes = fig.to_image(format="svg")
                st.download_button(
                    label=":material/download: Download Plot as SVG",
                    data=svg_bytes,
                    file_name=f"sankey_plot.svg",
                    mime="image/svg+xml",  # Set the MIME type to SVG
                    key='sankey_plot_download'
                )

                # Add interpretation help
                with st.expander("How to interpret this Sankey diagram"):
                    st.write(
                        """
                    - **Left nodes (Pharmacological Areas)**: Shows the different sources of detected compounds
                    - **Center nodes (Therapeutic Area)**: Shows the therapeutic areas
                    - **Right nodes (Compound name)**: Shows the compound names
                    - **Flow width**: The thickness of each flow represents the **number counts** for the connecting nodes
                    """
                    )

            else:
                st.info(
                    "Unable to generate Sankey diagram. This might be due to insufficient or incompatible data."
                )

        except Exception as e:
            st.error(f"Error generating Sankey diagram: {str(e)}")
            st.info("Please check your data format and try again.")


def add_df_and_filtering(df, key_prefix:str, default_cols: List = None):
    # Session state for tracking number of filters
    if f"{key_prefix}_filter_count" not in st.session_state:
        st.session_state[f"{key_prefix}_filter_count"] = 0

    add_col, remove_col, _, _ = st.columns(4)
    with add_col:
        # Button to add more filter fields
        if st.button("➕ Add Filter Field", use_container_width=True, key=f"{key_prefix}_add_btn"):
            st.session_state[f"{key_prefix}_filter_count"] += 1
    with remove_col:
        if st.button("➖ Remove Filter Field", use_container_width=True, key=f"{key_prefix}_rmv_btn"):
            st.session_state[f"{key_prefix}_filter_count"] -= 1

    filtered_df = df.copy()

    # Generate filter fields
    for i in range(st.session_state[f"{key_prefix}_filter_count"]):
        col1, col2 = st.columns([1, 2])
        with col1:
            if i == 0:
                st.markdown("**Filter Column**")
            selected_col = st.selectbox(
                f"Column {i+1}", filtered_df.columns, key=f"{key_prefix}_col_select_{i}"
            )
        with col2:
            if i == 0:
                st.markdown("**Search String**")
            search_term = st.text_input(
                f"Contains (Column {i+1})", key=f"{key_prefix}_search_input_{i}"
            )

        if selected_col and search_term:
            filtered_df = filtered_df[filtered_df[selected_col].str.contains(search_term, case=False, na=False)]

    # Show result
    st.markdown("### 🔎 Filtered Results")
    st.write(f"Total results: {len(filtered_df)}")
    all_cols = filtered_df.columns
    if default_cols:
        with st.expander('Cols to show'):
            cols_to_show = st.multiselect("Columns to show", options=all_cols, default=default_cols,
                                          label_visibility='collapsed')
    else:
        cols_to_show = all_cols

    return filtered_df[cols_to_show]


if __name__ == "__main__":
    # Example usage
    task_id = "34e2b1b692444bf6ae37e71dd137c300"
    file_name = "merged_results.tsv"

    lib_search_df = fetch_file(task_id, file_name)
    print(lib_search_df)  # Display the first few rows of the DataFrame
