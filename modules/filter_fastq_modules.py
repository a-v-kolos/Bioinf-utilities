def calculate_gc_content(seq):
    """
    Calculate GC content of a sequence in percentage
    """
    if not seq:
        print("Error: something wrong with the sequence")
        return None
    
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    return (gc_count / len(seq)) * 100


def calculate_average_quality(quality_string):
    """
    Calculate average read quality in phred33 format
    """
    if not quality_string:
        print("Error: something wrong with the string")
        return None
    
    total_quality = sum(ord(char) - 33 for char in quality_string)
    return total_quality / len(quality_string)


def check(value, bounds):
    """
    Check if value is within specified bounds
    """
    if isinstance(bounds, (int, float)):
        return 0 <= value <= bounds
    else:
        lower_bound, upper_bound = bounds
        return lower_bound <= value <= upper_bound
