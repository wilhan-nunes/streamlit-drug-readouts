"""
Helper functions to interact with Celery tasks from Streamlit.
Provides async task management and status monitoring.
"""
import streamlit as st
import pandas as pd
from celery.result import AsyncResult
from tasks import celery_app, download_data_task, process_analysis_task, reprocess_analysis_task
import time


def submit_download_task(task_id):
    """Submit download task to Celery and return task ID."""
    task = download_data_task.delay(task_id)
    return task.id


def submit_process_analysis_task(quant_file_df, annotation_file_df, config):
    """Submit analysis processing task to Celery and return task ID."""
    # Serialize dataframes
    quant_data = quant_file_df.to_dict('tight')
    annotation_data = annotation_file_df.to_dict('tight')
    
    task = process_analysis_task.delay(quant_data, annotation_data, config)
    return task.id


def submit_reprocess_analysis_task(feature_annotation_df, config):
    """Submit reanalysis task to Celery and return task ID."""
    feature_data = feature_annotation_df.to_dict('tight')
    
    task = reprocess_analysis_task.delay(feature_data, config)
    return task.id


def get_task_status(task_id):
    """
    Get the status of a Celery task.
    
    Returns:
        dict: Task status information with keys:
            - state: PENDING, STARTED, PROGRESS, SUCCESS, FAILURE
            - meta: Additional metadata (progress info, error messages, etc.)
            - result: Task result if completed
    """
    task = AsyncResult(task_id, app=celery_app)
    
    if task.state == 'PENDING':
        return {
            'state': task.state,
            'meta': {'status': 'Task is waiting to be processed...'},
            'result': None
        }
    elif task.state == 'PROGRESS':
        return {
            'state': task.state,
            'meta': task.info,
            'result': None
        }
    elif task.state == 'SUCCESS':
        return {
            'state': task.state,
            'meta': {'status': 'Task completed successfully!'},
            'result': task.result
        }
    elif task.state == 'FAILURE':
        return {
            'state': task.state,
            'meta': {
                'status': 'Task failed',
                'error': str(task.info)
            },
            'result': None
        }
    else:
        return {
            'state': task.state,
            'meta': task.info if task.info else {'status': f'Task state: {task.state}'},
            'result': None
        }


def wait_for_task(task_id):
    """
    Wait for a Celery task to complete with simple polling.
    
    Args:
        task_id: Celery task ID
        
    Returns:
        dict: Task result or None if failed
    """
    while True:
        status = get_task_status(task_id)
        state = status['state']
        
        if state == 'SUCCESS':
            return status['result']
        elif state == 'FAILURE':
            error_msg = status['meta'].get('error', 'Unknown error')
            raise Exception(f"Task failed: {error_msg}")
        
        time.sleep(1)  # Poll every second


def deserialize_analysis_data(result_dict):
    """
    Convert serialized task result back to pandas DataFrames.
    
    Args:
        result_dict: Dictionary with serialized dataframes
        
    Returns:
        dict: Dictionary with deserialized pandas DataFrames
    """
    from dataclasses import dataclass
    from typing import Optional
    
    deserialized = {}
    
    # Deserialize dataframes
    for key in ['feature_annotation', 'stratified_df', 'stratified_df_analogs', 
                'class_count_df', 'class_count_df_analog', 'default_excluded_features']:
        if key in result_dict and result_dict[key]:
            deserialized[key] = pd.DataFrame(**result_dict[key])
    
    # Keep dictionaries as-is
    for key in ['class_compound_dict', 'class_compound_dict_analog']:
        if key in result_dict:
            deserialized[key] = result_dict[key]
    
    return deserialized


def check_celery_connection():
    """
    Check if Celery is available and workers are running.
    
    Returns:
        bool: True if Celery is available, False otherwise
    """
    try:
        # Try to inspect active workers
        inspect = celery_app.control.inspect()
        active_workers = inspect.active()
        
        if active_workers:
            return True
        else:
            return False
    except Exception as e:
        print(f"Celery connection check failed: {e}")
        return False
