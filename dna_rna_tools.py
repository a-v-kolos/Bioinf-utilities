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


def is_nucleic_acid(seq):
        
    dna_nucleotides = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
    rna_nucleotides = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
    
    is_dna = all(char in dna_nucleotides for char in seq)
    
    is_rna = all(char in rna_nucleotides for char in seq)
    
    has_t = any(char in {'T', 't'} for char in seq)
    has_u = any(char in {'U', 'u'} for char in seq)
    
    if has_t and has_u:
        return False
    
    return is_dna or is_rna


def transcribe(seq):
    if not is_nucleic_acid(seq):
        return None
        
    transcription_map = str.maketrans('Tt', 'Uu')
    
    return seq.translate(transcription_map)


def reverse(seq):
    if not is_nucleic_acid(seq):
        print("Error: Check what's wrong with sequence or sequences")
        return None
    
    return seq[::-1]


def complement(seq):
    
    if not is_nucleic_acid(seq):
        return None
        
    dna_complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    rna_complement = str.maketrans('AUCGaucg', 'UAGCuagc')
    
    has_t = any(char in {'T', 't'} for char in seq)
    has_u = any(char in {'U', 'u'} for char in seq)
    
    if has_t:
        return seq.translate(dna_complement)
    elif has_u:
        return seq.translate(rna_complement)
    else:
        return seq.translate(dna_complement)


def reverse_complement(seq):
    if not is_nucleic_acid(seq):
        return None
    
    return complement(reverse(seq))
