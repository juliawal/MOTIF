"""
Top-level functions.
"""

import pandas as pd
import os
import tempfile
from distributed import Client, LocalCluster
from src.motif_grnboost.core import create_graph, SGBM_KWARGS, RF_KWARGS, EARLY_STOP_WINDOW_LENGTH


def grnboost2(expression_data,
              gene_names=None,
              tf_names='all',
              client_or_address='local',
              early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
              limit=None,
              seed=None,
              verbose=False,
              target_genes='all',
              output_file=None):
    """
    Launch arboreto with [GRNBoost2] profile.
    """
    result = diy(expression_data=expression_data, regressor_type='GBM', regressor_kwargs=SGBM_KWARGS,
               gene_names=gene_names, tf_names=tf_names, client_or_address=client_or_address,
               early_stop_window_length=early_stop_window_length, limit=limit, seed=seed, 
               verbose=verbose, target_genes=target_genes)
    
    if output_file:
        _safe_save_result(result, output_file, verbose)
    return result


def genie3(expression_data,
           gene_names=None,
           tf_names='all',
           client_or_address='local',
           limit=None,
           seed=None,
           verbose=False,
           target_genes='all',
           output_file=None):
    """
    Launch arboreto with [GENIE3] profile.
    """
    result = diy(expression_data=expression_data, regressor_type='RF', regressor_kwargs=RF_KWARGS,
               gene_names=gene_names, tf_names=tf_names, client_or_address=client_or_address,
               limit=limit, seed=seed, verbose=verbose, target_genes=target_genes)
    
    if output_file:
        _safe_save_result(result, output_file, verbose)
    return result


def _safe_save_result(result, output_file, verbose=False):
    """Safely save results to file using temporary directory."""
    try:
        tmpdir = tempfile.gettempdir()
        temp_file = os.path.join(tmpdir, f"temp_{os.path.basename(output_file)}")
        
        if verbose:
            print(f"Temporärer Speicherort: {temp_file}")
            print(f"Endgültiger Speicherort: {output_file}")
        
        # First save to temp location
        result.to_csv(temp_file, sep="\t", index=False)
        
        # Atomic file move
        os.replace(temp_file, output_file)
        
        if verbose:
            print(f"Ergebnis erfolgreich gespeichert nach: {output_file}")
    except Exception as e:
        print(f"Fehler beim Speichern der Ergebnisse: {str(e)}")
        # Cleanup temp file if it exists
        if os.path.exists(temp_file):
            os.remove(temp_file)
        raise


def diy(expression_data,
        regressor_type,
        regressor_kwargs,
        gene_names=None,
        tf_names='all',
        client_or_address='local',
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        limit=None,
        seed=None,
        verbose=False,
        target_genes='all'):
    """
    Main computation function.
    """
    if verbose:
        print('preparing dask client')

    client, shutdown_callback = _prepare_client(client_or_address)

    try:
        if verbose:
            print('parsing input')

        expression_matrix, gene_names, tf_names = _prepare_input(expression_data, gene_names, tf_names)

        if verbose:
            print('creating dask graph')

        graph = create_graph(expression_matrix,
                           gene_names,
                           tf_names,
                           client=client,
                           regressor_type=regressor_type,
                           regressor_kwargs=regressor_kwargs,
                           early_stop_window_length=early_stop_window_length,
                           limit=limit,
                           seed=seed,
                           target_genes=target_genes)

        if verbose:
            print(f'{graph.npartitions} partitions')
            print('computing dask graph')

        result = client.compute(graph, sync=True).sort_values(by='importance', ascending=False)
        
        return result

    finally:
        shutdown_callback(verbose)
        if verbose:
            print('finished')


def _prepare_client(client_or_address):
    """
    Prepare dask client with memory-safe settings.
    """
    if client_or_address is None or str(client_or_address).lower() == 'local':
        local_cluster = LocalCluster(
            n_workers=1,                      # Reduzierte Worker für weniger Overhead
            threads_per_worker=1,              # Vermeidet Thread-Konflikte
            memory_limit="4GB",                # Strengeres Memory-Limit
            memory_target_fraction=0.5,       # Früheres Spilling
            memory_spill_fraction=0.7,        # Aggressiveres Spilling
            diagnostics_port=None,
            local_directory="/tmp/dask-tmp",  # Fallback auf /tmp
            silence_logs=40                   # Weniger Logging
        )
        client = Client(local_cluster)

        def close_client_and_local_cluster(verbose=False):
            if verbose:
                print('shutting down client and local cluster')
            client.close()
            local_cluster.close()

        return client, close_client_and_local_cluster

    elif isinstance(client_or_address, str) and client_or_address.lower() != 'local':
        client = Client(
            client_or_address,
            memory_limit="4GB",
            memory_target_fraction=0.5,
            memory_spill_fraction=0.7,
            local_directory="/tmp/dask-tmp",
            silence_logs=40
        )

        def close_client(verbose=False):
            if verbose:
                print('shutting down client')
            client.close()

        return client, close_client

    elif isinstance(client_or_address, Client):
        def close_dummy(verbose=False):
            if verbose:
                print('not shutting down externally created client')
        return client_or_address, close_dummy

    else:
        raise ValueError(f"Invalid client specified {str(client_or_address)}")


def _prepare_input(expression_data,
                   gene_names,
                   tf_names):
    """Prepare input data."""
    if isinstance(expression_data, pd.DataFrame):
        expression_matrix = expression_data.to_numpy()
        gene_names = list(expression_data.columns)
    else:
        expression_matrix = expression_data
        assert expression_matrix.shape[1] == len(gene_names)

    if tf_names is None or tf_names == 'all':
        tf_names = gene_names
    else:
        if len(tf_names) == 0:
            raise ValueError('Specified tf_names is empty')
        if not set(gene_names).intersection(set(tf_names)):
            raise ValueError('Intersection of gene_names and tf_names is empty.')

    return expression_matrix, gene_names, tf_names