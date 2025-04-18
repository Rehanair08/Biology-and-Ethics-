def find_all_mutations(original, mutated):
    return [(i, o, m) for i, (o, m) in enumerate(zip(original, mutated)) if o != m]

def find_pam_site(sequence):
    """Finds the first valid NGG PAM site."""
    for i in range(len(sequence) - 2):
        triplet = sequence[i:i+3]
        if triplet[1:] == 'GG':
            return i
    return None

def cas9_cleave_site(pam_position):
    return pam_position - 3

def correct_sequence(original, mutated):
    corrected = list(mutated)
    for i in range(len(original)):
        if original[i] != mutated[i]:
            corrected[i] = original[i]
    return ''.join(corrected)

def similarity_percentage(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return round(matches / len(seq1) * 100, 2)

# Sequences
original_seq = "ATGCTAGCTAGGCTACGTAGCTAGGATCGTAGGCTAACGTAGCTAGCTAG"
mutated_seq  = "ATGCTAGTTAGGCTTCGTAGTTAGGATGGTACGCTAGCCTAGATAGTTAG"

# Step 1 - Find mutations
mutations = find_all_mutations(original_seq, mutated_seq)

# Step 2 - Find PAM
pam_pos = find_pam_site(mutated_seq)

# Output
print("="*60)
print("           CRISPR-Cas9 Gene Editing Simulation        ")
print("="*60)

print("\n>> Step 1 - Input Sequences:")
print("Original Sequence :", original_seq)
print("Mutated Sequence  :", mutated_seq)

if mutations:
    print("\n>> Detected Mutations:")
    for pos, o, m in mutations:
        print(f"  - Position {pos+1}: {m} → {o}")
else:
    print(" No mutations detected.")

if pam_pos is not None:
    cleavage_pos = cas9_cleave_site(pam_pos)
    pam_seq = mutated_seq[pam_pos:pam_pos+3]

    print("\n>> Step 2 - PAM Site & Cas9 Cleavage:")
    print(f" PAM site found at position {pam_pos+1} ({pam_seq})")
    print(f" Cas9 Cleaves at position {cleavage_pos+1}")

    # Visual with cleavage
    visual_seq = list(mutated_seq)
    if 0 <= cleavage_pos < len(visual_seq):
        visual_seq[cleavage_pos] = f"[{visual_seq[cleavage_pos]}]"
    print("Visual:\n" + ''.join(visual_seq))

    print("\n>> Step 3 - Homology-Directed Repair (HDR):")
    for pos, o, m in mutations:
        print(f"Corrected base at position {pos + 1} ({m} → {o})")

    # Correct mutations
    edited_seq = correct_sequence(original_seq, mutated_seq)
    similarity = similarity_percentage(original_seq, edited_seq)

    print("\n" + "="*60)
    print("                   Final Edited Sequence              ")
    print("="*60)
    print("Edited Sequence:", edited_seq)

    print("\n" + "="*60)
    print("                    Similarity Check                  ")
    print("="*60)
    print(f" Similarity to Original: {similarity}%")
    if similarity == 100.0:
        print(" All mismatches corrected. Sequence successfully restored!")
    else:
        print(" Some mismatches remain.")
else:
    print("\n No valid NGG PAM site found. Cas9 editing not possible.")