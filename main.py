import argparse
import json
import sys

from modules.filter_fastq_modules import * 
from modules.dna_rna_modules import *


def setup_logging(log_file: str = "bio_toolkit.log") -> logging.Logger:
    """Configure logging to both file and stderr."""
    logger = logging.getLogger("bio_toolkit")
    logger.setLevel(logging.DEBUG)
 
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
 
    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
 
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging.WARNING)
    console_handler.setFormatter(formatter)
 
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
 
    return logger
 
 
logger = setup_logging()




def run_dna_rna_tools(*args):
    
    if len(args) < 1:
        return None
    
    *seqs, operation = args  

    if len(seqs) == 0: 
        return None

    logger.info("Running '%s' on %d sequence(s): %s", operation, len(seqs), seqs)
 
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
        logger.error("Unknown operation: '%s'", operation)
        return None
 
    logger.info("'%s' completed. Results: %s", operation, results)
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
logger.info(
        "Starting filter_fastq: %d sequence(s), gc_bounds=%s, length_bounds=%s, quality_threshold=%s",
        len(seqs), gc_bounds, length_bounds, quality_threshold,
    )
 
    filtered_seqs = {}
 
    for seq_name, (seq, quality) in seqs.items():
        seq_length = len(seq)
        if not check(seq_length, length_bounds):
            logger.debug("'%s' rejected: length %d outside bounds %s", seq_name, seq_length, length_bounds)
            continue
 
        gc_content = calculate_gc_content(seq)
        if not check(gc_content, gc_bounds):
            logger.debug("'%s' rejected: GC content %.1f%% outside bounds %s", seq_name, gc_content, gc_bounds)
            continue
 
        avg_quality = calculate_average_quality(quality)
        if avg_quality < quality_threshold:
            logger.debug("'%s' rejected: avg quality %.2f below threshold %.2f", seq_name, avg_quality, quality_threshold)
            continue
 
        filtered_seqs[seq_name] = (seq, quality)
 
    logger.info(
        "filter_fastq done: %d/%d sequence(s) passed.",
        len(filtered_seqs), len(seqs),
    )
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
        logger.error(
            "dna-rna command failed: operation='%s', seqs=%s",
            args.operation, args.seqs,
        )
        print("Error: invalid sequences or operation.", file=sys.stderr)
        sys.exit(1)
    if isinstance(result, list):
        for item in result:
            print(item)
    else:
        print(result)
 
 
def cmd_filter_fastq(args):
    seqs = {}
    if args.seq:
        for triplet in args.seq:
            if len(triplet) != 3:
                logger.error("--seq received %d value(s) instead of 3: %s", len(triplet), triplet)
                print(
                    f"Error: --seq requires exactly 3 values (name seq quality), got: {triplet}",
                    file=sys.stderr,
                )
                sys.exit(1)
            name, sequence, quality = triplet
            seqs[name] = (sequence, quality)
 
    if not seqs:
        logger.error("filter-fastq called with no sequences provided.")
        print("Error: no sequences provided. Use --seq NAME SEQUENCE QUALITY.", file=sys.stderr)
        sys.exit(1)
 
    try:
        gc_bounds = parse_bounds(args.gc_bounds)
        length_bounds = parse_bounds(args.length_bounds)
    except argparse.ArgumentTypeError as e:
        logger.error("Failed to parse bounds: %s", e)
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
 
    result = filter_fastq(
        seqs,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=args.quality_threshold,
    )
 
    if not result:
        logger.info("No sequences passed the filters.")
        print("No sequences passed the filters.")
    else:
        for name, (seq, quality) in result.items():
            print(f"{name}\t{seq}\t{quality}")
 
def build_parser():
    parser = argparse.ArgumentParser(
        prog="main.py",
        description="Bioinformatics toolkit: DNA/RNA tools and FASTQ filtering.",
    )
    parser.add_argument(
        "--log-file",
        default="bio_toolkit.log",
        metavar="PATH",
        help="Path to the log file (default: bio_toolkit.log).",
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
        help="GC content range as 'min,max' (default: 0,100).",
    )
    fq_parser.add_argument(
        "--length-bounds",
        default=f"0,{2**32}",
        metavar="MIN,MAX",
        help="Sequence length range as 'min,max' (default: 0,4294967296).",
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
 
    # Re-configure logging if user specified a custom log file
    if args.log_file != "bio_toolkit.log":
        for handler in logger.handlers[:]:
            if isinstance(handler, logging.FileHandler):
                logger.removeHandler(handler)
        new_file_handler = logging.FileHandler(args.log_file, encoding="utf-8")
        new_file_handler.setLevel(logging.DEBUG)
        new_file_handler.setFormatter(logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        ))
        logger.addHandler(new_file_handler)
 
    logger.info("Command started: %s", " ".join(sys.argv))
    args.func(args)
    logger.info("Command finished: %s", args.command)
