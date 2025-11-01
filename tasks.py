"""
Celery tasks for heavy computation in drug readout analysis.
These tasks are delegated to Celery workers to avoid blocking the Streamlit UI.
"""
import ast
import pandas as pd
from celery import Celery
import os

# Connect to Redis (broker) and store results back in Redis
celery_app = Celery(
    "streamlit-drug-readouts_tasks",
    broker='redis://streamlit-drug-readouts-redis',
    backend='redis://streamlit-drug-readouts-redis'
)

celery_app.conf.update(
    result_expires=900,              # 0.25 hour expiration for results (prevents Redis bloat)
)

@celery_app.task()
def heartbeat_task():
    return "Post-mn worker is alive."

@celery_app.task(bind=True, name='tasks.download_data')
def download_data_task(self, task_id):
    """
    Download quantification and annotation data from GNPS.
    
    Args:
        task_id: GNPS task ID
        
    Returns:
        dict: Serialized dataframes
    """
    from gnpsdata import workflow_fbmn
    
    self.update_state(state='PROGRESS', meta={'current': 0, 'total': 2, 'status': 'Downloading quantification data...'})
    quant_file_df = workflow_fbmn.get_quantification_dataframe(task_id, gnps2=True)
    
    self.update_state(state='PROGRESS', meta={'current': 1, 'total': 2, 'status': 'Downloading annotation data...'})
    annotation_file_df = workflow_fbmn.get_library_match_dataframe(task_id, gnps2=True)
    
    return {
        'quant_file_df': quant_file_df.to_dict('tight'),
        'annotation_file_df': annotation_file_df.to_dict('tight')
    }


@celery_app.task(bind=True, name='tasks.process_analysis')
def process_analysis_task(self, quant_file_data, annotation_file_data, config):
    """
    Process the main analysis data - the most computationally intensive operation.
    
    Args:
        quant_file_data: Quantification dataframe as dict
        annotation_file_data: Annotation dataframe as dict
        config: Configuration dictionary with intensity_thresh, blank_ids, etc.
        
    Returns:
        dict: Serialized AnalysisData object
    """
    from script import (
        load_and_filter_features,
        load_and_merge_annotations,
        generate_feature_annotation,
        stratify_by_drug_class,
        count_drug_class_occurrences
    )
    
    # Deserialize dataframes
    quant_file_df = pd.DataFrame(**quant_file_data)
    annotation_file_df = pd.DataFrame(**annotation_file_data)
    
    _drug_metadata_file = "data/GNPS_Drug_Library_Metadata_Drugs.csv"
    _analog_metadata_file = "data/GNPS_Drug_Library_Metadata_Drug_Analogs_Updated.csv"
    
    # Step 1: Filter features
    self.update_state(state='PROGRESS', meta={'current': 1, 'total': 6, 'status': 'Filtering features...'})
    subtract_blanks = True if config.get('blank_ids') else False
    feature_filtered = load_and_filter_features(
        quant_file_df,
        intensity_threshold=config['intensity_thresh'],
        blank_ids=config.get('blank_ids'),
        subtract_blanks=subtract_blanks,
    )
    
    # Step 2: Load and merge annotations
    self.update_state(state='PROGRESS', meta={'current': 2, 'total': 6, 'status': 'Merging annotations...'})
    _annotation_metadata = load_and_merge_annotations(
        annotation_file_df, _drug_metadata_file, _analog_metadata_file
    )
    
    # Step 3: Generate feature annotation
    self.update_state(state='PROGRESS', meta={'current': 3, 'total': 6, 'status': 'Generating feature annotations...'})
    _feature_annotation, _excluded_features = generate_feature_annotation(
        _annotation_metadata, feature_filtered
    )
    
    # Step 4: Stratify by drug class
    self.update_state(state='PROGRESS', meta={'current': 4, 'total': 6, 'status': 'Stratifying by drug class...'})
    _stratified_df = stratify_by_drug_class(
        _feature_annotation, exclude_analogs=True, peak_threshold=config['intensity_thresh'],
    )
    _stratified_df_analogs = stratify_by_drug_class(
        _feature_annotation, exclude_analogs=False, peak_threshold=config['intensity_thresh'],
    )
    
    # Step 5: Count drug class occurrences
    self.update_state(state='PROGRESS', meta={'current': 5, 'total': 6, 'status': 'Counting drug class occurrences...'})
    _class_count_df, _class_count_df_analog, _class_compounds_dict, _class_compounds_dict_analog = count_drug_class_occurrences(
        _feature_annotation, class_column="pharmacologic_class"
    )
    _class_count_df["total_matches"] = _class_count_df.sum(axis=1)
    _class_count_df_analog["total_matches"] = _class_count_df.sum(axis=1)
    
    # Step 6: Prepare result
    self.update_state(state='PROGRESS', meta={'current': 6, 'total': 6, 'status': 'Finalizing results...'})
    
    # Serialize all data
    result = {
        'feature_annotation': _feature_annotation.to_dict('tight'),
        'stratified_df': _stratified_df.to_dict('tight'),
        'stratified_df_analogs': _stratified_df_analogs.to_dict('tight'),
        'class_count_df': _class_count_df.to_dict('tight'),
        'class_count_df_analog': _class_count_df_analog.to_dict('tight'),
        'default_excluded_features': _excluded_features.to_dict('tight'),
        'class_compound_dict': _class_compounds_dict,
        'class_compound_dict_analog': _class_compounds_dict_analog
    }
    
    return result


@celery_app.task(bind=True, name='tasks.reprocess_analysis')
def reprocess_analysis_task(self, feature_annotation_data, config):
    """
    Reprocess analysis with edited feature annotation table.
    
    Args:
        feature_annotation_data: Edited feature annotation dataframe as dict
        config: Configuration dictionary
        
    Returns:
        dict: Serialized AnalysisData object
    """
    from script import (
        stratify_by_drug_class,
        count_drug_class_occurrences
    )
    
    # Deserialize dataframes
    _feature_annotation = pd.DataFrame(**feature_annotation_data)
    
    # Step 1: Stratify by drug class
    self.update_state(state='PROGRESS', meta={'current': 1, 'total': 3, 'status': 'Re-stratifying by drug class...'})
    _stratified_df = stratify_by_drug_class(
        _feature_annotation, exclude_analogs=True, peak_threshold=config['intensity_thresh'],
    )
    _stratified_df_analogs = stratify_by_drug_class(
        _feature_annotation, exclude_analogs=False, peak_threshold=config['intensity_thresh'],
    )
    
    # Step 2: Count drug class occurrences
    self.update_state(state='PROGRESS', meta={'current': 2, 'total': 3, 'status': 'Recounting drug class occurrences...'})
    _class_count_df, _class_count_df_analog, _class_compounds_dict, _class_compounds_dict_analog = count_drug_class_occurrences(
        _feature_annotation, class_column="pharmacologic_class"
    )
    _class_count_df["total_matches"] = _class_count_df.sum(axis=1)
    _class_count_df_analog["total_matches"] = _class_count_df.sum(axis=1)
    
    # Step 3: Prepare result
    self.update_state(state='PROGRESS', meta={'current': 3, 'total': 3, 'status': 'Finalizing reanalysis...'})
    
    result = {
        'stratified_df': _stratified_df.to_dict('tight'),
        'stratified_df_analogs': _stratified_df_analogs.to_dict('tight'),
        'class_count_df': _class_count_df.to_dict('tight'),
        'class_count_df_analog': _class_count_df_analog.to_dict('tight'),
        'class_compound_dict': _class_compounds_dict,
        'class_compound_dict_analog': _class_compounds_dict_analog
    }
    
    return result