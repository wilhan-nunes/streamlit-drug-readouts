from typing import List

import pandas as pd


def load_and_filter_features(
    file_path: str,
    blank_ids: List[str] | None = None,
    intensity_threshold: int = 10000,
    file_types: List[str] = ["mzML", "mzXML"],
    subtract_blanks: bool = False,
) -> pd.DataFrame:
    """
    Load and filter features from a CSV file based on intensity threshold and file types.
    :param file_path: Path to the CSV file containing features.
    :param intensity_threshold: Basically only trust detections with peak area above this number
    :param file_types: List of possible file types to filter the features.
    :param blank_ids: List of substrings to identify blank or control columns.
    :return: pd.Dataframe: Filtered features intensity table.
    """
    feature = pd.read_csv(file_path)
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
    annotation_file: str, druglib_file: str, analoglib_file: str
) -> pd.DataFrame:
    """
    Loads and merges annotations with drug library metadata.
    :param annotation_file: Path to the annotation TSV file.
    :param druglib_file: Path to the parent drug metadata CSV.
    :param analoglib_file: Path to the drug analog metadata CSV.
    :return: pd.DataFrame: Annotation table merged with metadata.
    """
    annotation = pd.read_csv(annotation_file, sep="\t")
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
    return annotation_metadata.merge(feature_filtered, on="FeatureID", how="inner")


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
                "Background|confidence|Endogenous|Food|analog", na=False
            )
        ]
    else:
        exo_drug = feature_annotation[
            ~feature_annotation["chemical_source"].str.contains(
                "Background|confidence|Endogenous|Food", na=False
            )
        ]

    exo_drug.loc[:, "therapeutic_area"] = exo_drug["therapeutic_area"].fillna("NA")
    exo_drug.loc[:, "pharmacologic_class"] = exo_drug["pharmacologic_class"].fillna(
        "NA"
    )

    def aggregate_by_column(exo_df, group_col, prefix=None):
        expanded = exo_df.copy()
        expanded[group_col] = expanded[group_col].astype(str).str.split("|")
        expanded = expanded.explode(group_col)
        expanded = expanded[~expanded[group_col].isin(["NA", "no_match"])]
        aggregate = expanded.groupby(group_col)[sample_columns].sum()
        if prefix:
            aggregate.index = [prefix] * len(aggregate)
        return aggregate

    # Therapeutic Area
    exo_drug_agg_thera = aggregate_by_column(exo_drug, "therapeutic_area")

    # Pharmacologic Class (FDA)
    exo_drug_agg_fda = aggregate_by_column(exo_drug, "pharmacologic_class")

    # Antidepressants
    depression = exo_drug[
        exo_drug["therapeutic_indication"].str.contains(
            "depression", case=False, na=False
        )
    ]
    agg_depression = depression[sample_columns].sum().to_frame().T
    agg_depression.index = ["antidepressants"]

    # Antibiotics
    abx = exo_drug[
        exo_drug["pharmacologic_class"].str.contains(
            "microbial|bacterial|tetracycline", case=False, na=False
        )
    ]
    agg_abx = abx[sample_columns].sum().to_frame().T
    agg_abx.index = ["antibiotics"]

    # Final aggregation and binarization
    final_df = pd.concat(
        [agg_abx, agg_depression, exo_drug_agg_thera, exo_drug_agg_fda]
    )
    final_df_TF = final_df > 0
    final_df_TF = final_df_TF.T.map(lambda x: "Yes" if x == True else "No")
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
    from utils import fetch_file

    ## Setup file paths and task IDs
    task_id = "d6f37a11d90c4f249974280c3fc90108"
    quant_file_path = fetch_file(task_id, "quant_table.csv", type="quant_table")
    annotation_file_path = fetch_file(
        task_id, "annotations.tsv", type="annotation_table"
    )
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
