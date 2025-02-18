import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from loguru import logger

# dict for convertion from nucleotide to one-hot encoding
extended_mapping = {
    'A': [1.0, 0.0, 0.0, 0.0],  # Adenine
    'T': [0.0, 1.0, 0.0, 0.0],  # Thymine
    'G': [0.0, 0.0, 1.0, 0.0],  # Guanine
    'C': [0.0, 0.0, 0.0, 1.0],  # Cytosine
    'R': [0.5, 0.0, 0.5, 0.0],  # A or G (purine)
    'Y': [0.0, 0.5, 0.0, 0.5],  # C or T (pyrimidine)
    'S': [0.0, 0.0, 0.5, 0.5],  # G or C
    'W': [0.5, 0.5, 0.0, 0.0],  # A or T
    'K': [0.0, 0.5, 0.5, 0.0],  # G or T
    'M': [0.5, 0.0, 0.0, 0.5],  # A or C
    'B': [0.0, 0.33, 0.33, 0.33],  # C, G, or T
    'D': [0.33, 0.33, 0.33, 0.0],  # A, G, or T
    'H': [0.33, 0.33, 0.0, 0.33],  # A, C, or T
    'V': [0.33, 0.0, 0.33, 0.33],  # A, C, or G
    'N': [0.25, 0.25, 0.25, 0.25],  # Any base (A, T, G, or C)
}


def seq2vec(seq: str) -> np.ndarray:
    """
    Convert a nucleotide sequence to a 4-dimensional one-hot encoded matrix,
    supporting an extended alphabet with ambiguous nucleotides.

    Args:
        seq (str): Input nucleotide sequence (e.g., "ATGRY").

    Returns:
        np.ndarray: One-hot encoded matrix of shape (L, 4), where L is the sequence length.
    """
    # Create an empty matrix of shape (L, 4)
    one_hot = np.zeros((len(seq), 4))
    
    # Fill the matrix
    for i, nucleotide in enumerate(seq):
        if nucleotide.upper() in extended_mapping:
            one_hot[i] = extended_mapping[nucleotide.upper()]
    
    return one_hot


def qual2vec(quality_str: str, offset: int = 33) -> np.ndarray:
    """
    Convert a quality string to a vector of quality scores.

    Args:
        quality (str): Input quality string (e.g., "!#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]_`abcdefghijklmnopqrstuvwxyz{|}~").
        offset (int): The ASCII offset (33 for Phred+33, 64 for Phred+64).

    Returns:
        np.ndarray: Vector of quality scores of shape (L,).
    """
    return np.frombuffer(quality_str.encode(), dtype=np.uint8) - offset


def calculate_log_average_quality(quals: np.ndarray, s_axis: int = 0) -> np.ndarray:
    """
    Calculate the logarithmic average quality score using NumPy.

    Args:
        quals (np.ndarray): Input array of integer quality scores.

    Returns:
        float: The computed logarithmic average quality score.
    """
    # Compute 10^(q / -10) for each quality score
    exp_quals = 10 ** (quals / -10)
    # Compute the average of the transformed scores
    avg_exp_quals = np.mean(exp_quals, axis=s_axis)
    # Compute the final result: -10 * log10(average)
    result = -10 * np.log10(avg_exp_quals)
    
    return result


def calculate_mertics(
    record, 
    kmer_length: int = 15, 
    window_size: int = 100, 
    window_overlap: int = 50, 
    max_outer_downgrade: float = 0.5, 
    min_inner_complexity: float = 0.5, 
    max_inner_downgrade: float = 0.66, 
    desired_q_value: float = 20, 
    min_q_change: float = 0.01):
    """
    Calculate metrics for a given record.

    Args:
        record (dict): A dictionary containing the following keys:
            - "id": The ID of the read.
            - "sequence": The sequence of the read.
    """
    read_id = record["id"]
    
    # последовательность прочтения, записанная как вектор
    seq = seq2vec(record["sequence"])
    logger.debug(f"{read_id} 'seq' shape: {seq.shape}")
    # последовательность качества прочтения, записанная как вектор
    qual = qual2vec(record["quality"])
    logger.debug(f"{read_id} 'qual' shape: {qual.shape}")
    # последовательность k-мер вектора последовательности
    vector_kmers_repr = sliding_window_view(seq, window_shape=(kmer_length), axis=0)
    logger.debug(f"{read_id} vector of kmers shape: {vector_kmers_repr.shape}")
    # вектор качества каждого нуклеотида в каждой к-мере
    vector_kmers_qual = sliding_window_view(qual, window_shape=kmer_length)
    logger.debug(f"{read_id} vector of kmers quality shape: {vector_kmers_qual.shape}")
    # вектор среднего качества к-меры
    vector_kmers_avg_qual = calculate_log_average_quality(vector_kmers_qual, s_axis=1)
    logger.debug(f"{read_id} vector of kmers avg quality shape: {vector_kmers_avg_qual.shape}")
    # вектор отличия текущей следующей к-меры от текущей к-меры
    vector_kmers_difference = np.diff(vector_kmers_repr, axis=0)
    logger.debug(f"{read_id} vector of kmer outer difference: {vector_kmers_difference.shape}")
    # вектор отличия следующего нуклеотида от предыдущего в пределах к-меры
    vector_kmers_inner_diff = np.diff(vector_kmers_repr, axis=2, )
    logger.debug(f"{read_id} vector of kmer inner difference: {vector_kmers_inner_diff.shape}")
    
    # норма вектора внешней разницы
    norm_kmers_difference = np.linalg.norm(vector_kmers_difference, axis=(1, 2), ord=1)
    logger.debug(f"{read_id} norm of outer difference shape: {norm_kmers_difference.shape}")
    # норма сложности каждой следующей к-меры, то сколько раз в пределах к-меры следующий нуклеотид отличается от предыдущего
    norm_kmer_complexity = np.count_nonzero(abs(vector_kmers_inner_diff).sum(axis=1), axis=1) / kmer_length
    logger.debug(f"{read_id} norm of inner difference shape: {norm_kmer_complexity.shape}")
    # производная изменения качества
    vector_kmers_avg_qual_change = np.diff(vector_kmers_avg_qual)
    logger.debug(f"{read_id} vector of kmer quality change: {vector_kmers_avg_qual_change.shape}")
    
    failed_regions = list()
    
    idx = 0
    while idx < vector_kmers_avg_qual_change.size:
        window_left = idx 
        window_right = idx + window_size
    
        outer = norm_kmers_difference[window_left:window_right]
        inner = norm_kmer_complexity[:-1][window_left:window_right]
        qual = vector_kmers_avg_qual[:-1][window_left:window_right]
        qual_ch = vector_kmers_avg_qual_change[window_left:window_right]
    
        # если изменения к-мер пропали
        decision_outer = any(outer <= max_outer_downgrade)
        # если сложность падает до минимально разрешенного размера или перепад сложности более чем на 0.66
        decision_inner = any(inner <= min_inner_complexity) or max(inner) - min(inner) >= max_inner_downgrade
        # ищем, где качество к-меры падает ниже требуемого значения
        low_q_arr_condition = qual <= desired_q_value
        # после чего выкидываем такой участок, если более половины окна меньше порогового значения
        decision_qual = len(qual[low_q_arr_condition]) >= 0.5 * (window_right - window_left)
        # условие на слабую скорость изменения такого качества в этом участке
        decision_qual_ch = any(abs(qual_ch[low_q_arr_condition]) <= min_q_change)
    
        if (decision_outer and decision_inner) or (decision_qual and decision_qual_ch):
            logger.trace(f"{read_id} window {[window_left, window_right]} passed with {[decision_outer, decision_inner, decision_qual, decision_qual_ch]}")
            
            if len(failed_regions) > 0 and failed_regions[-1][1] >= window_left:
                failed_regions[-1][1] = window_right
            else:
                failed_regions.append([window_left, window_right])
        
        idx = window_right - window_overlap
    
    logger.debug(f"{read_id} failed regions: {failed_regions}")
    logger.debug(f"{read_id} failed regions size: {len(failed_regions)}")

    result = {
        "id": read_id,
        "metrics": {
            "norm_kmers_difference": norm_kmers_difference,
            "norm_kmer_complexity": norm_kmer_complexity,
            "vector_kmers_avg_qual": vector_kmers_avg_qual,
            "vector_kmers_avg_qual_change": vector_kmers_avg_qual_change,
            "failed_regions": failed_regions
        },
        'length': len(record['sequence'])
    }
    return result


def get_wanted_regions(sequence_length, dont_want_list):
    """
    Get regions of the sequence that are not in the `dont_want_list`.

    Args:
        sequence_length (int): Total length of the sequence.
        dont_want_list (list): List of [start, end] ranges to exclude.

    Returns:
        list: List of [start, end] ranges for wanted regions.
    """
    # Sort the excluded ranges by start position
    dont_want_list.sort(key=lambda x: x[0])
    
    # Initialize the list of wanted regions
    wanted_regions = []
    prev_end = 0  # Start from the beginning of the sequence
    
    # Iterate through the excluded ranges and find the wanted regions
    for start, end in dont_want_list:
        if start > prev_end:
            # Add the region before the excluded range
            wanted_regions.append([prev_end, start - 1])
        prev_end = max(prev_end, end + 1)  # Update the end of the last excluded range
    
    # Add the region after the last excluded range
    if prev_end < sequence_length:
        wanted_regions.append([prev_end, sequence_length - 1])
    
    return wanted_regions
