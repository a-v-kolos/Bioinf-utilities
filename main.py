from modules.filter_fastq_modules import * 
from modules.dna_rna_modules import *


def run_dna_rna_tools(*args):
    
    if len(args) < 1:
        return None
    
    *seqs, operation = args  

    if len(seqs) == 0: 
        return None

    if operation == 'is_nucleic_acid': 
        results = [is_nucleic_acid(seq) for seq in seqs] 
    elif operation == 'transcribe': 
        results = [transcribe(seq) for seq in seqs]   
    elif operation == 'reverse':      
        results = [reverse(seq) for seq in seqs]  
    elif operation == 'complement':      
        results = [complement(seq) for seq in seqs]      
    elif operation == 'reverse_complement':  
        results = [reverse_complement(seq) for seq in seqs]  
    else:
        return None
        
    return results[0] if len(results) == 1 else results


def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    """
    Filter fastq sequences by GC content, length and quality
    
    Args:
    seqs: dictionary sequence_name: (seq, quality)
    gc_bounds: GC content range in default: (0, 100)
    length_bounds: length range in default: (0, 2**32)
    quality_threshold: minimum average quality by default: 0
    
    Returns:
    dictionary with filtered sequences
    """
    filtered_seqs = {}
    
    for seq_name, (seq, quality) in seqs.items():
        seq_length = len(seq)
        if not check(seq_length, length_bounds):
            continue
        
        gc_content = calculate_gc_content(seq)
        if not check(gc_content, gc_bounds):
            continue
        
        avg_quality = calculate_average_quality(quality)
        if avg_quality < quality_threshold:
            continue
        
        filtered_seqs[seq_name] = (seq, quality)
    
    return filtered_seqs
