from typing import List

import pandas as pd


def load_and_filter_features(
    features_data: str | pd.DataFrame,
    blank_ids: List[str] | None = None,
    intensity_threshold: int = 10000,
    file_types: List[str] = ["mzML", "mzXML"],
    subtract_blanks: bool = False,
) -> pd.DataFrame:
    """
    Load and filter features from a CSV file based on intensity threshold and file types.
    :param features_data: Path to the CSV file OR dataframe containing features
    :param intensity_threshold: Basically only trust detections with peak area above this number
    :param file_types: List of possible file types to filter the features.
    :param blank_ids: List of substrings to identify blank or control columns.
    :return: pd.Dataframe: Filtered features intensity table.
    """
    if isinstance(features_data, str):
        feature = pd.read_csv(features_data)
    elif isinstance(features_data, pd.DataFrame):
        feature = features_data.copy()

    feature_filtered = feature.set_index("row ID").filter(regex="|".join(file_types))
    feature_filtered[feature_filtered < intensity_threshold] = 0
    if subtract_blanks:
        # Optionally perform blank subtraction if needed
        # Identify blank and sample columns based on their names
        feature_blank = feature_filtered.loc[
            :, feature_filtered.columns.str.contains("|".join(blank_ids))
        ]
        feature_sample = feature_filtered.loc[
            :, ~feature_filtered.columns.str.contains("|".join(blank_ids))
        ]

        # Ensure row indices match
        assert all(feature_blank.index == feature_sample.index)

        # Calculate average blank and sample values per feature
        blk_check = pd.DataFrame(
            {
                "avg_blank": feature_blank.mean(axis=1),
                "avg_sample": feature_sample.mean(axis=1),
            }
        )
        blk_check["ratio_sample_blank"] = (blk_check["avg_sample"] + 1) / (
            blk_check["avg_blank"] + 1
        )
        blk_check["bin"] = (blk_check["ratio_sample_blank"] > 5).astype(int)

        # Keep only features with a high enough sample/blank ratio
        feature_blk_removed = feature_sample[blk_check["bin"] == 1]
        feature_filtered = feature_blk_removed
    return feature_filtered


def load_and_merge_annotations(
    fbmn_annotation_data: str | pd.DataFrame, druglib_file: str, analoglib_file: str
) -> pd.DataFrame:
    """
    Loads and merges annotations with drug library metadata.
    :param fbmn_annotation_data: Path to the annotation TSV file.
    :param druglib_file: Path to the parent drug metadata CSV.
    :param analoglib_file: Path to the drug analog metadata CSV.
    :return: pd.DataFrame: Annotation table merged with metadata.
    """
    if isinstance(fbmn_annotation_data, str):
        annotation = pd.read_csv(fbmn_annotation_data, sep="\t")
    elif isinstance(fbmn_annotation_data, pd.DataFrame):
        annotation = fbmn_annotation_data

    annotation_filtered = annotation[["#Scan#", "SpectrumID", "Compound_Name"]]
    annotation_filtered.columns = ["FeatureID", "SpectrumID", "Compound_Name"]

    druglib = pd.read_csv(
        druglib_file,
        usecols=[
            "gnps_libid",
            "name_parent_compound",
            "chemical_source",
            "pharmacologic_class",
            "therapeutic_area",
            "therapeutic_indication",
        ],
    )

    analoglib = pd.read_csv(
        analoglib_file,
        usecols=[
            "analog_libid",
            "name_connected_compound",
            "chemical_source",
            "pharmacologic_class",
            "therapeutic_area",
            "therapeutic_indication",
        ],
    ).rename(
        columns={
            "analog_libid": "gnps_libid",
            "name_connected_compound": "name_parent_compound",
        }
    )

    all_metadata = pd.concat([druglib, analoglib], ignore_index=True)
    annotation_metadata = annotation_filtered.merge(
        all_metadata, left_on="SpectrumID", right_on="gnps_libid", how="left"
    ).drop("gnps_libid", axis=1)
    return annotation_metadata


def generate_feature_annotation(
    annotation_metadata: pd.DataFrame, feature_filtered: pd.DataFrame
) -> pd.DataFrame:
    """
    Combines filtered feature table with annotated metadata.
    :param annotation_metadata: Merged annotation metadata.
    :param feature_filtered: Intensity-filtered feature table.
    :return: pd.DataFrame: Feature annotation with intensities.
    """
    feature_filtered = feature_filtered.reset_index()
    feature_filtered = feature_filtered.rename(columns={"row ID": "FeatureID"})
    feature_filtered["FeatureID"] = feature_filtered["FeatureID"].astype(str)
    annotation_metadata["FeatureID"] = annotation_metadata["FeatureID"].astype(str)
    merged_feature_metadata = annotation_metadata.merge(feature_filtered, on="FeatureID", how="inner")
    return merged_feature_metadata


def stratify_by_drug_class(
    feature_annotation: pd.DataFrame,
    column_pattern: str = r".*\.mzML|.*\.mzXML",
    exclude_analogs: bool = False,
) -> pd.DataFrame:
    """
    Stratifies samples based on drug class or therapeutic area.

    :param feature_annotation: DataFrame with feature annotations and intensities.
    :param column_pattern: Regex pattern to select sample columns.
    :param exclude_analogs: Whether to exclude analogs.
    :return: pd.DataFrame: Transposed binary DataFrame indicating drug class presence by sample.
    """
    sample_columns = feature_annotation.filter(regex=column_pattern).columns

    if exclude_analogs:
        exo_drug = feature_annotation[
            ~feature_annotation["chemical_source"].str.contains(
                "Background|confidence|Endogenous|Food|analog", case=False, na=False
            )
        ]
    else:
        exo_drug = feature_annotation[
            ~feature_annotation["chemical_source"].str.contains(
                "Background|confidence|Endogenous|Food", case=False, na=False
            )
        ]

    # Fill NaN values
    exo_drug = exo_drug.copy()
    exo_drug.loc[:, "therapeutic_area"] = exo_drug["therapeutic_area"].fillna("NA")
    exo_drug.loc[:, "pharmacologic_class"] = exo_drug["pharmacologic_class"].fillna("NA")
    exo_drug.loc[:, "therapeutic_indication"] = exo_drug["therapeutic_indication"].fillna("NA")

    def aggregate_by_column(exo_df, group_col, prefix=None):
        expanded = exo_df.copy()
        expanded[group_col] = expanded[group_col].astype(str).str.split("|")
        expanded = expanded.explode(group_col)
        expanded = expanded[~expanded[group_col].isin(["NA", "no_match"])]
        if expanded.empty:
            # Return empty DataFrame with correct structure if no data
            empty_df = pd.DataFrame(columns=sample_columns)
            if prefix:
                empty_df.index = [prefix] if isinstance(prefix, str) else prefix
            return empty_df
        aggregate = expanded.groupby(group_col)[sample_columns].sum()
        if prefix:
            aggregate.index = [prefix] * len(aggregate) if isinstance(prefix, str) else prefix
        return aggregate

    # === Specific Drug Categories ===

    # 1. Antidepressants
    depression = exo_drug[
        exo_drug["therapeutic_indication"].str.contains(
            "depression", case=False, na=False
        )
    ]
    agg_depression = depression[sample_columns].sum().to_frame().T
    agg_depression.index = ["antidepressants"]

    # 2. Antibiotics
    abx = exo_drug[
        exo_drug["pharmacologic_class"].str.contains(
            "microbial|bacterial|tetracycline", case=False, na=False
        )
    ]
    agg_abx = abx[sample_columns].sum().to_frame().T
    agg_abx.index = ["antibiotics"]

    # 3. PPIs (Proton Pump Inhibitors)
    ppi = exo_drug[
        exo_drug["pharmacologic_class"].str.contains(
            "proton pump inhibitor", case=False, na=False
        )
    ]
    agg_ppi = ppi[sample_columns].sum().to_frame().T
    agg_ppi.index = ["PPI"]


    # 4. Statins
    statin = exo_drug[
        exo_drug["pharmacologic_class"].str.contains(
            "statin|HMG-CoA reductase inhibitor", case=False, na=False
        )
    ]
    agg_statin = statin[sample_columns].sum().to_frame().T
    agg_statin.index = ["statin"]

    # 5. Antihistamines
    antihistamine = exo_drug[
        exo_drug["pharmacologic_class"].str.contains(
            "histamine-1", case=False, na=False
        )
    ]
    agg_antihistamine = antihistamine[sample_columns].sum().to_frame().T
    agg_antihistamine.index = ["antihistamine"]

    # 6. Antihypertensives
    antihypertensives = exo_drug[
        exo_drug["therapeutic_indication"].str.contains(
            "hypertension", case=False, na=False
        )
    ]
    agg_antihypertensives = antihypertensives[sample_columns].sum().to_frame().T
    agg_antihypertensives.index = ["antihypertensive"]

    # 7. Alzheimer's medications
    alzheimer = exo_drug[
        exo_drug["therapeutic_indication"].str.contains(
            "Alzheimer", case=False, na=False
        )
    ]
    agg_alzheimer = alzheimer[sample_columns].sum().to_frame().T
    agg_alzheimer.index = ["Alzheimer"]

    # 8. Antifungals
    antifungal = exo_drug[
        (exo_drug["pharmacologic_class"].str.contains(
            "antifungal", case=False, na=False
        )) |
        (exo_drug["therapeutic_indication"].str.contains(
            "fungal infection", case=False, na=False
        ))
    ]
    agg_antifungal = antifungal[sample_columns].sum().to_frame().T
    agg_antifungal.index = ["antifungal"]

    # 9. HIV medications
    hiv_med = exo_drug[
        (exo_drug["therapeutic_indication"].str.contains(
            "HIV", case=False, na=False
        )) |
        (exo_drug["therapeutic_indication"].str.contains(
            "atazanavir", case=False, na=False
        ))
    ]
    agg_hiv_med = hiv_med[sample_columns].sum().to_frame().T
    agg_hiv_med.index = ["HIVmed"]

    # === General Categories ===

    # 10. Therapeutic Area aggregation
    exo_drug_agg_thera = aggregate_by_column(exo_drug, "therapeutic_area")

    # 11. Pharmacologic Class (FDA) aggregation
    exo_drug_agg_fda = aggregate_by_column(exo_drug, "pharmacologic_class")

    # === Combine All Categories ===

    # Collect all aggregated DataFrames
    all_dfs = [
        agg_abx,
        agg_depression,
        agg_ppi,
        agg_statin,
        agg_antihistamine,
        agg_antihypertensives,
        agg_alzheimer,
        agg_antifungal,
        agg_hiv_med,
        exo_drug_agg_thera,
        exo_drug_agg_fda
    ]

    # Filter out empty DataFrames
    non_empty_dfs = [df for df in all_dfs if not df.empty]

    if not non_empty_dfs:
        # Return empty DataFrame with correct structure if no data
        empty_result = pd.DataFrame(columns=["Sample"])
        return empty_result

    # Final aggregation and binarization
    final_df = pd.concat(non_empty_dfs, sort=False)

    # Convert to binary (True/False) then to Yes/No
    final_df_TF = final_df > 0
    final_df_TF = final_df_TF.T.map(lambda x: "Yes" if x else "No")
    final_df_TF = final_df_TF.reset_index().rename(columns={"index": "Sample"})

    return final_df_TF


def count_drug_class_occurrences(
    feature_annotation: pd.DataFrame,
    class_column: str = "pharmacologic_class",
    file_extensions: List[str] = ["mzML", "mzXML"],
) -> pd.DataFrame:
    """
    Counts the number of times each drug class appears in each sample.
    :param feature_annotation: Annotated feature intensity data.
    :param class_column: Column name for drug class (e.g., 'pharmacologic_class').
    :param file_extensions: File extensions used to select intensity columns.
    :return: pd.DataFrame: DataFrame with sample-wise counts for each class.
    """
    pattern = "|".join(file_extensions)
    df = feature_annotation.copy()
    df[class_column] = df[class_column].fillna("NA")
    df = df[df[class_column] != "NA"]
    df = df[df[class_column] != "no_match"]
    df = df.assign(class_group=df[class_column].str.split("\\|")).explode("class_group")
    df = df[df["class_group"].notna()]
    df_class = df[["class_group"] + df.filter(regex=pattern).columns.tolist()]
    df_class_binary = df_class.copy()
    df_class_binary[df_class_binary.columns[1:]] = (
        df_class_binary[df_class_binary.columns[1:]].gt(0).astype(int)
    )
    return df_class_binary.groupby("class_group").sum().T


if __name__ == "__main__":
    # --- User-Defined Parameters Section ---
    # This section allows the user to modify parameters for running the script as a standalone file.
    from gnpsdata import workflow_fbmn

    ## Setup file paths and task IDs
    task_id = "d6f37a11d90c4f249974280c3fc90108"
    quant_file_path = workflow_fbmn.get_quantification_dataframe(task_id, gnps2=True)
    annotation_file_path = workflow_fbmn.get_library_match_dataframe(task_id, gnps2=True)


    drug_metadata_file = "data/GNPS_Drug_Library_Metadata_Drugs.csv"
    analog_metadata_file = "data/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv"

    # Load and process data
    blank_ids = ["QC"]  # This should be set to the actual blank IDs if available
    feature_filtered = load_and_filter_features(
        quant_file_path, blank_ids, subtract_blanks=False
    )

    annotation_metadata = load_and_merge_annotations(
        annotation_file_path, drug_metadata_file, analog_metadata_file
    )
    feature_annotation = generate_feature_annotation(
        annotation_metadata, feature_filtered
    )
    stratified_df = stratify_by_drug_class(feature_annotation, exclude_analogs=True)
    stratified_df_analogs = stratify_by_drug_class(
        feature_annotation, exclude_analogs=False
    )

    # adding a summary of drug class occurrence per sample
    class_count_df = count_drug_class_occurrences(
        feature_annotation, class_column="pharmacologic_class"
    )
    class_count_df["total_matches"] = class_count_df.sum(axis=1)
    class_count_df_sorted = class_count_df.sort_values("total_matches", ascending=False)

    # Save the stratified DataFrame to a CSV file
    stratified_df.to_csv("stratified.tsv", sep="\t", index=False)
    stratified_df_analogs.to_csv("stratified_analogs.tsv", sep="\t", index=False)
    feature_annotation.to_csv("feature_annotation.tsv", sep="\t", index=False)
