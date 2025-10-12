# Bioinf-utilities

A set of tools for working with biological sequences, including DNA/RNA analysis and FASTQ file filtering.

## Installation

1. Clone the repository:
```bash
git clone git@github.com:a-v-kolos/Bioinf-utilities.git
```

2. Make sure you have Python 3.6 or higher installed.

## Supported operations:


1. The main function run_dna_rna_tools() takes sequences and an operation:

is_nucleic_acid - nucleic acid validation

transcribe - DNA→RNA transcription

reverse - sequence reversal

complement - complementary sequence

reverse_complement - reverse-complementary sequence


2. The filter_fastq() function filters sequences by GC content, length, and quality:

gc_bounds: GC content range in % (default: 0-100)

length_bounds: sequence length range (default: 0-2³²)

quality_threshold: minimum average read quality (phred33)

## Modules
dna_rna_modules.py
Contains functions for nucleic acid manipulation:

is_nucleic_acid(seq) - checks if sequence is DNA or RNA

transcribe(seq) - transcribes DNA to RNA

reverse(seq) - returns reversed sequence

complement(seq) - returns complementary sequence

reverse_complement(seq) - returns reverse-complementary sequence

filter_fastq_modules.py
Contains helper functions for FASTQ filtering:

calculate_gc_content(seq) - calculates GC content in %

calculate_average_quality(quality_string) - calculates average read quality

check(value, bounds) - checks if value is within specified bounds

bio_files_processor.py
Contains utility functions for processing bioinformatics file formats:

parse_blast_output(input_file: str, output_file: str) -> None - extracts the best protein matches from BLAST results and saves them to a sorted list

convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str | None = None) -> None - converts a FASTA file with multiline sequences to a format where each sequence is contained in a single line


