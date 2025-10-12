def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str | None = None) -> None:
    """
    Converts a FASTA file with multiline sequences to a format with single-line sequences.
    
    Args:
        input_fasta: Path to the input FASTA file
        output_fasta: Path to the output FASTA file. If not provided, it will be generated automatically.
    """
    
    if output_fasta is None:
        if input_fasta.endswith('.fasta') or input_fasta.endswith('.fa'):
            base_name = input_fasta.rsplit('.', 1)[0]
            output_fasta = f"{base_name}_oneline.fasta"
        else:
            output_fasta = f"{input_fasta}_oneline.fasta"
    
    try:
        with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
            current_header: str | None = None
            current_sequence: list[str] = []
            
            for line in infile:
                line = line.strip()
                
                if line.startswith('>'): 
                   
                    if current_header is not None:
                        
                        complete_sequence = ''.join(current_sequence)
                        outfile.write(f"{current_header}\n{complete_sequence}\n")
                    
                    current_header = line
                    current_sequence = []
                else:
                    current_sequence.append(line)
            
            if current_header is not None:
                complete_sequence = ''.join(current_sequence)
                outfile.write(f"{current_header}\n{complete_sequence}\n")
                
        print(f"File successfully converted: {output_fasta}")
        
    except FileNotFoundError:
        print(f"Error: File {input_fasta} not found")
    except Exception as e:
        print(f"An error occurred: {e}")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parse BLAST results and extract the best match descriptions for each query.
    
    Args:
        input_file: Path to the input file with BLAST results
        output_file: Path to the output file for saving results
    """
    best_matches = []
    current_query = None
    in_significant_alignments = False
    first_match_found = False
    
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('Query #'):
                current_query = line
                in_significant_alignments = False
                first_match_found = False
            
            elif line.startswith('Sequences producing significant alignments:'):
                in_significant_alignments = True
                first_match_found = False
            
            elif in_significant_alignments and line.startswith('Description'):
                continue
            
            elif in_significant_alignments and line.startswith('Alignments:'):
                in_significant_alignments = False
                first_match_found = False
            
            elif in_significant_alignments and line and not first_match_found:
                description = line.split('    ')[0].strip()
                if description and not description.startswith('Scientific') and not description.startswith('---'):
                    best_matches.append(description)
                    first_match_found = True
    
    best_matches.sort()
    
    with open(output_file, 'w', encoding='utf-8') as f:
        for match in best_matches:
            f.write(match + '\n')
