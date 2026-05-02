import argparse
import json
import sys

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


def parse_bounds(value):
    """Parse a bound argument: either a single number or 'min,max'."""
    parts = value.split(',')
    if len(parts) == 1:
        return float(parts[0])
    elif len(parts) == 2:
        return (float(parts[0]), float(parts[1]))
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid bounds '{value}'. Use a single number or 'min,max'."
        )
 
 
def cmd_dna_rna(args):
    result = run_dna_rna_tools(*args.seqs, args.operation)
    if result is None:
        print("Error: invalid sequences or operation.", file=sys.stderr)
        sys.exit(1)
    if isinstance(result, list):
        for item in result:
            print(item)
    else:
        print(result)

def cmd_filter_fastq(args):
    """
    Expects sequences passed as repeating triplets via --seq:
        --seq NAME SEQUENCE QUALITY
    Example:
        --seq read1 ATGC !!!! --seq read2 GCTA ????
    """
    seqs = {}
    if args.seq:
        for triplet in args.seq:
            if len(triplet) != 3:
                print(
                    f"Error: --seq requires exactly 3 values (name seq quality), got: {triplet}",
                    file=sys.stderr,
                )
                sys.exit(1)
            name, sequence, quality = triplet
            seqs[name] = (sequence, quality)
 
    gc_bounds = parse_bounds(args.gc_bounds)
    length_bounds = parse_bounds(args.length_bounds)
 
    result = filter_fastq(
        seqs,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=args.quality_threshold,
    )
 
    if not result:
        print("No sequences passed the filters.")
    else:
        for name, (seq, quality) in result.items():
            print(f"{name}\t{seq}\t{quality}")
 
 
def build_parser():
    parser = argparse.ArgumentParser(
        prog="main.py",
        description="Bioinformatics toolkit: DNA/RNA tools and FASTQ filtering.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
 

    dna_rna_parser = subparsers.add_parser(
        "dna-rna",
        help="Perform DNA/RNA operations on one or more sequences.",
    )
    dna_rna_parser.add_argument(
        "operation",
        choices=["is_nucleic_acid", "transcribe", "reverse", "complement", "reverse_complement"],
        help="Operation to apply to each sequence.",
    )
    dna_rna_parser.add_argument(
        "seqs",
        nargs="+",
        metavar="SEQ",
        help="One or more nucleotide sequences.",
    )
    dna_rna_parser.set_defaults(func=cmd_dna_rna)
 

    fq_parser = subparsers.add_parser(
        "filter-fastq",
        help="Filter FASTQ sequences by GC content, length, and quality.",
    )
    fq_parser.add_argument(
        "--seq",
        nargs=3,
        metavar=("NAME", "SEQUENCE", "QUALITY"),
        action="append",
        help="A sequence entry: name, nucleotide sequence, quality string. "
             "Repeat for multiple sequences.",
    )
    fq_parser.add_argument(
        "--gc-bounds",
        default="0,100",
        metavar="MIN,MAX",
        help="GC content range as 'min,max' (default: 0,100). "
             "Pass a single number for an upper bound only.",
    )
    fq_parser.add_argument(
        "--length-bounds",
        default=f"0,{2**32}",
        metavar="MIN,MAX",
        help="Sequence length range as 'min,max' (default: 0,4294967296). "
             "Pass a single number for an upper bound only.",
    )
    fq_parser.add_argument(
        "--quality-threshold",
        type=float,
        default=0,
        metavar="THRESHOLD",
        help="Minimum average Phred33 quality score (default: 0).",
    )
    fq_parser.set_defaults(func=cmd_filter_fastq)
 
    return parser
 
 
if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
