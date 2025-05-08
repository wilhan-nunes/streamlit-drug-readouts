import pandas as pd
from typing import List


def load_and_filter_features(
    file_path: str,
    intensity_threshold: int = 10000,
    file_types: List[str] = ["mzML", "mzXML"],
) -> pd.DataFrame:
    """
    Load and filter features from a CSV file based on intensity threshold and file types.
    :param file_path: Path to the CSV file containing features.
    :param intensity_threshold: Basically only trust detections with peak area above this number
    :param file_types: List of possible file types to filter the features.
    :return: pd.Dataframe: Filtered features intensity table.
    """
    feature = pd.read_csv(file_path)
    feature_filtered = feature.set_index("row ID").filter(regex="|".join(file_types))
    feature_filtered[feature_filtered < intensity_threshold] = 0
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

    druglib = pd.read_csv(druglib_file, usecols=[
        "gnps_libid",
        "name_parent_compound",
        "chemical_source",
        "pharmacologic_class",
        "therapeutic_area",
        "therapeutic_indication",
    ])

    analoglib = pd.read_csv(analoglib_file, usecols=[
        "analog_libid",
        "name_connected_compound",
        "chemical_source",
        "pharmacologic_class",
        "therapeutic_area",
        "therapeutic_indication",
    ]).rename(
        columns={
            "analog_libid": "gnps_libid",
            "name_connected_compound": "name_parent_compound",
        }
    )

    all_metadata = pd.concat([druglib, analoglib], ignore_index=True)
    annotation_metadata = annotation_filtered.merge(
        all_metadata, left_on="SpectrumID", right_on="gnps_libid", how="left"
    ).drop('gnps_libid', axis=1)
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
    feature_filtered = feature_filtered.rename(
        columns={"row ID": "FeatureID"}
    )
    feature_filtered['FeatureID'] = feature_filtered['FeatureID'].astype(str)
    annotation_metadata["FeatureID"] = annotation_metadata["FeatureID"].astype(str)
    return annotation_metadata.merge(feature_filtered, on="FeatureID", how="inner")


def stratify_by_drug_class(
    feature_annotation: pd.DataFrame,
    column_pattern: str = r".*\.mzML|.*\.mzXML",
    exclude_analogs: bool = False
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
        exo_drug = feature_annotation[~feature_annotation['chemical_source'].str.contains("Background|confidence|Endogenous|Food|analog", na=False)]
    else:
        exo_drug = feature_annotation[~feature_annotation['chemical_source'].str.contains("Background|confidence|Endogenous|Food", na=False)]

    exo_drug['therapeutic_area'] = exo_drug['therapeutic_area'].fillna("NA")
    exo_drug['pharmacologic_class'] = exo_drug['pharmacologic_class'].fillna("NA")

    def aggregate_by_column(exo_df, group_col, prefix=None):
        expanded = exo_df.copy()
        expanded[group_col] = expanded[group_col].astype(str).str.split('|')
        expanded = expanded.explode(group_col)
        expanded = expanded[~expanded[group_col].isin(['NA', 'no_match'])]
        aggregate = expanded.groupby(group_col)[sample_columns].sum()
        if prefix:
            aggregate.index = [prefix] * len(aggregate)
        return aggregate

    # Therapeutic Area
    exo_drug_agg_thera = aggregate_by_column(exo_drug, 'therapeutic_area')

    # Pharmacologic Class (FDA)
    exo_drug_agg_fda = aggregate_by_column(exo_drug, 'pharmacologic_class')

    # Antidepressants
    depression = exo_drug[exo_drug['therapeutic_indication'].str.contains("depression", case=False, na=False)]
    agg_depression = depression[sample_columns].sum().to_frame().T
    agg_depression.index = ['antidepressants']

    # Antibiotics
    abx = exo_drug[exo_drug['pharmacologic_class'].str.contains("microbial|bacterial|tetracycline", case=False, na=False)]
    agg_abx = abx[sample_columns].sum().to_frame().T
    agg_abx.index = ['antibiotics']

    # Final aggregation and binarization
    final_df = pd.concat([agg_abx, agg_depression, exo_drug_agg_thera, exo_drug_agg_fda])
    final_df_TF = final_df > 0
    final_df_TF = final_df_TF.T
    final_df_TF = final_df_TF.reset_index().rename(columns={"index": "Sample"})
    return final_df_TF

if __name__ == "__main__":
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

    feature_filtered = load_and_filter_features(quant_file_path)
    annotation_metadata = load_and_merge_annotations(
        annotation_file_path, drug_metadata_file, analog_metadata_file
    )
    feature_annotation = generate_feature_annotation(
        annotation_metadata, feature_filtered
    )
    stratified_df = stratify_by_drug_class(feature_annotation, exclude_analogs=True)
    stratified_df_analogs = stratify_by_drug_class(feature_annotation, exclude_analogs=False)

    # Save the stratified DataFrame to a CSV file
    stratified_df.to_csv('new_stratified.tsv', sep='\t', index=False)
    stratified_df_analogs.to_csv('new_stratified_analogs.tsv', sep='\t', index=False)
    feature_annotation.to_csv('feature_annotation.tsv', sep='\t', index=False)

