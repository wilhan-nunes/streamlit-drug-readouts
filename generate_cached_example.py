"""
Script to download and cache data for task 4d99fc25d84143bdbbf2dd07bf044e5e
This should be run once to generate the cached pickle file and example data files.
"""
import pickle
from gnpsdata import workflow_fbmn
from script import *
import pandas as pd
from dataclasses import dataclass
from typing import Optional

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

def main():
    task_id = "4d99fc25d84143bdbbf2dd07bf044e5e"
    threshold = 10000  # Based on script.py
    blank_ids = "blank"
    
    print(f"Downloading data for task {task_id}...")
    
    # Download files
    quant_file_df = workflow_fbmn.get_quantification_dataframe(task_id, gnps2=True)
    annotation_file_df = workflow_fbmn.get_library_match_dataframe(task_id, gnps2=True)
    
    print("Saving raw data files...")
    quant_file_df.to_csv(f'data/examples/example_quant_file_{task_id}.tsv', sep='\t', index=False)
    annotation_file_df.to_csv(f'data/examples/example_annotation_file_{task_id}.tsv', sep='\t', index=False)
    
    print("Processing analysis data...")
    
    # Load metadata files
    _drug_metadata_file = "data/GNPS_Drug_Library_Metadata_Drugs.csv"
    _analog_metadata_file = "data/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv"
    
    # Filter features
    subtract_blanks = True if blank_ids else False
    feature_filtered = load_and_filter_features(
        quant_file_df,
        intensity_threshold=threshold,
        blank_ids=blank_ids,
        subtract_blanks=subtract_blanks,
    )
    
    # Load and merge annotations
    _annotation_metadata = load_and_merge_annotations(
        annotation_file_df, _drug_metadata_file, _analog_metadata_file
    )
    
    # Generate feature annotation
    _feature_annotation, _excluded_features = generate_feature_annotation(
        _annotation_metadata, feature_filtered
    )
    
    # Perform analysis
    _stratified_df = stratify_by_drug_class(
        _feature_annotation, exclude_analogs=True, peak_threshold=threshold,
    )
    _stratified_df_analogs = stratify_by_drug_class(
        _feature_annotation, exclude_analogs=False, peak_threshold=threshold,
    )
    
    # Count drug class occurrences
    _class_count_df, _class_count_df_analog, _class_compounds_dict, _class_compounds_dict_analog = count_drug_class_occurrences(
        _feature_annotation, class_column="pharmacologic_class"
    )
    _class_count_df["total_matches"] = _class_count_df.sum(axis=1)
    _class_count_df_analog["total_matches"] = _class_count_df.sum(axis=1)
    
    # Create data object
    data = AnalysisData(
        feature_annotation=_feature_annotation,
        stratified_df=_stratified_df,
        stratified_df_analogs=_stratified_df_analogs,
        class_count_df=_class_count_df,
        class_count_df_analog=_class_count_df_analog,
        default_excluded_features=_excluded_features,
        class_compound_dict=_class_compounds_dict,
        class_compound_dict_analog=_class_compounds_dict_analog
    )
    
    print("Saving cached analysis data...")
    pickle.dump(data, open(f'./data/examples/processed_analysis_data_{task_id}.pkl', 'wb'))
    
    print(f"✅ Successfully generated cache files for task {task_id}")
    print(f"   - Raw data files: data/examples/example_quant_file_{task_id}.tsv, data/examples/example_annotation_file_{task_id}.tsv")
    print(f"   - Cached analysis: data/examples/processed_analysis_data_{task_id}.pkl")
    print(f"\nAnalysis Summary:")
    print(f"   - Feature annotations: {len(_feature_annotation)}")
    print(f"   - Samples: {len(_stratified_df)}")
    print(f"   - Excluded features: {len(_excluded_features)}")
    print(f"   - Threshold used: {threshold}")

if __name__ == "__main__":
    main()
