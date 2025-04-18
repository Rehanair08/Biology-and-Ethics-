def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_length = len(sequence)
    if total_length == 0:
        return 0
    return (gc_count / total_length) * 100

def calculate_tm(sequence):
    """Estimate melting temperature (Tm) using a basic formula."""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    at_count = sequence.count('A') + sequence.count('T')
    length = len(sequence)
    # Wallace rule for short sequences (<50 bp): Tm = 2(A+T) + 4(G+C)
    if length < 50:
        tm = 2 * at_count + 4 * gc_count
    else:
        # Simplified nearest-neighbor approximation for longer sequences
        tm = 64.9 + 41 * (gc_count - 16.4) / length
    return tm

def check_repetitive_patterns(sequence, min_repeat_length=6):
    """Check for repetitive patterns in the sequence."""
    sequence = sequence.upper()
    length = len(sequence)
    repeat_score = 0
    for i in range(length - min_repeat_length + 1):
        substring = sequence[i:i + min_repeat_length]
        occurrences = sequence.count(substring)
        if occurrences > 1:
            repeat_score += occurrences - 1
    return repeat_score > 5  # True if highly repetitive

def check_palindromes(sequence, min_length=8):
    """Check for palindromic sequences that could form cruciforms."""
    sequence = sequence.upper()
    length = len(sequence)
    for i in range(length - min_length + 1):
        subseq = sequence[i:i + min_length]
        # Reverse complement
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        rev_comp = ''.join(complement[base] for base in subseq[::-1])
        if subseq == rev_comp:
            return True  # Found a palindrome
    return False

def check_dinucleotide_bias(sequence):
    """Check for extreme CpG or AT content."""
    sequence = sequence.upper()
    length = len(sequence)
    cpg_count = sequence.count('CG')
    at_count = sequence.count('AT') + sequence.count('TA')
    cpg_ratio = cpg_count / (length - 1) if length > 1 else 0
    at_ratio = at_count / (length - 1) if length > 1 else 0
    return cpg_ratio > 0.2 or at_ratio > 0.4  # Arbitrary thresholds for instability

def check_hairpin_potential(sequence, stem_length=4):
    """Check for potential hairpin-forming regions."""
    sequence = sequence.upper()
    length = len(sequence)
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i in range(length - 2 * stem_length - 1):
        stem1 = sequence[i:i + stem_length]
        for j in range(i + stem_length + 1, length - stem_length + 1):
            stem2 = sequence[j:j + stem_length]
            rev_comp_stem2 = ''.join(complement[base] for base in stem2[::-1])
            if stem1 == rev_comp_stem2:
                return True  # Hairpin possible
    return False

def predict_stability(sequence):
    """Predict if a DNA sequence is stable or unstable with all conditions."""
    sequence = sequence.upper()
    
    # Validate input
    valid_bases = set('ATCG')
    if not all(base in valid_bases for base in sequence) or len(sequence) == 0:
        return "Invalid sequence: Only A, T, C, G are allowed."

    # Calculate basic properties
    gc_content = calculate_gc_content(sequence)
    tm = calculate_tm(sequence)
    length = len(sequence)
    is_repetitive = check_repetitive_patterns(sequence)
    has_palindromes = check_palindromes(sequence)
    has_dinucleotide_bias = check_dinucleotide_bias(sequence)
    has_hairpin_potential = check_hairpin_potential(sequence)
    
    # Stability rules and explanations
    reasons = []
    instability_factors = 0
    
    # GC content
    if 40 <= gc_content <= 60:
        reasons.append(f"GC content ({gc_content:.1f}%) is in stable range (40-60%).")
    else:
        reasons.append(f"GC content ({gc_content:.1f}%) is outside stable range (40-60%).")
        instability_factors += 1
    
    # Melting temperature
    if 50 <= tm <= 80:
        reasons.append(f"Melting temperature ({tm:.1f}°C) is in stable range (50-80°C).")
    else:
        reasons.append(f"Melting temperature ({tm:.1f}°C) is outside stable range (50-80°C).")
        instability_factors += 1
    
    # Repetition
    if is_repetitive:
        reasons.append("High repetition detected, indicating potential instability.")
        instability_factors += 1
    else:
        reasons.append("No significant repetition detected.")
    
    # Palindromes
    if has_palindromes:
        reasons.append("Palindromic sequence detected, may form cruciforms (instability).")
        instability_factors += 1
    else:
        reasons.append("No significant palindromes detected.")
    
    # Dinucleotide bias
    if has_dinucleotide_bias:
        reasons.append("Extreme CpG or AT bias detected, suggesting instability.")
        instability_factors += 1
    else:
        reasons.append("No extreme dinucleotide bias detected.")
    
    # Hairpin potential
    if has_hairpin_potential:
        reasons.append("Hairpin-forming potential detected, indicating instability.")
        instability_factors += 1
    else:
        reasons.append("No hairpin-forming potential detected.")
    
    # Length
    if length < 20:
        reasons.append("Sequence is very short (<20 bases), may be less stable.")
        instability_factors += 1
    else:
        reasons.append(f"Sequence length ({length} bases) is sufficient.")
    
    # Final prediction
    if instability_factors == 0:
        stability = "Stable"
    elif instability_factors <= 2:
        stability = "Moderately Stable"
    else:
        stability = "Unstable"
    
    explanation = "\n".join(reasons)
    return f"Prediction: {stability}\nExplanation:\n{explanation}"

# Get user input and predict stability directly
if __name__ == "__main__":
    sequence = input("Enter a DNA sequence (A, T, C, G only): ").strip()
    result = predict_stability(sequence)
    print("\n" + result)
    