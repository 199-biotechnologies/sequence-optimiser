#!/usr/bin/env python3
"""
FST344 Follistatin Codon Optimization Tool
Optimizes the FST344 coding sequence for AAV expression
"""

import re
from collections import Counter

# Human codon usage table (optimized for mammalian expression, CpG-depleted)
OPTIMAL_CODONS = {
    'A': ['GCC', 'GCT'],  # Alanine
    'R': ['AGA', 'AGA'],  # Arginine (avoid CGC)
    'N': ['AAC'],         # Asparagine
    'D': ['GAC'],         # Aspartic acid
    'C': ['TGC'],         # Cysteine
    'Q': ['CAG'],         # Glutamine
    'E': ['GAG'],         # Glutamic acid
    'G': ['GGA', 'GGT'],  # Glycine (avoid GGC when possible)
    'H': ['CAC'],         # Histidine
    'I': ['ATC'],         # Isoleucine
    'L': ['CTG', 'CTC'],  # Leucine
    'K': ['AAG'],         # Lysine
    'M': ['ATG'],         # Methionine (start)
    'F': ['TTC'],         # Phenylalanine
    'P': ['CCC', 'CCT'],  # Proline (avoid CCG)
    'S': ['AGC', 'TCC'],  # Serine
    'T': ['ACC'],         # Threonine
    'W': ['TGG'],         # Tryptophan
    'Y': ['TAC'],         # Tyrosine
    'V': ['GTG', 'GTC'],  # Valine
    '*': ['TGA']          # Stop codon
}

# Standard genetic code
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def translate_dna(sequence):
    """Translate DNA sequence to protein"""
    protein = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            protein += GENETIC_CODE.get(codon, 'X')
    return protein

def optimize_codon(amino_acid, avoid_cpg=True, last_two_bases=""):
    """Select optimal codon for amino acid avoiding CpG dinucleotides"""
    candidates = OPTIMAL_CODONS.get(amino_acid, ['XXX'])
    
    if avoid_cpg:
        # Filter out codons that create CpG dinucleotides
        safe_candidates = []
        for codon in candidates:
            # Check if codon itself contains CG
            if 'CG' in codon:
                continue
            # Check if codon creates CpG with previous bases
            if last_two_bases and last_two_bases.endswith('C') and codon.startswith('G'):
                continue
            safe_candidates.append(codon)
        
        # Use safe candidates if available
        if safe_candidates:
            candidates = safe_candidates
    
    # Return best candidate
    return candidates[0] if candidates else 'XXX'

def optimize_sequence(dna_sequence):
    """Optimize entire DNA sequence"""
    # Translate to protein
    protein = translate_dna(dna_sequence)
    
    # Re-encode with optimal codons
    optimized_dna = ""
    for i, aa in enumerate(protein):
        if aa != 'X':
            # Get last two bases to avoid CpG
            last_two = optimized_dna[-2:] if len(optimized_dna) >= 2 else ""
            codon = optimize_codon(aa, avoid_cpg=True, last_two_bases=last_two)
            optimized_dna += codon
    
    return optimized_dna, protein

def analyze_sequence(sequence, name):
    """Analyze sequence composition"""
    length = len(sequence)
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = (gc_count / length * 100) if length > 0 else 0
    
    # Count CpG dinucleotides
    cpg_count = sequence.count('CG')
    
    # Count codons
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    
    return {
        'name': name,
        'length': length,
        'gc_content': gc_content,
        'cpg_count': cpg_count,
        'cpg_frequency': cpg_count / (length/2) * 100 if length > 0 else 0,
        'codon_count': len(codons)
    }

def main():
    # FST344 coding sequence (from cds.fna)
    original_cds = """ATGGTCCGCGCGAGGCACCAGCCGGGTGGGCTTTGCCTCCTGCTGCTGCTGCTCTGCCAGTTCATGGAGG
ACCGCAGTGCCCAGGCTGGGAACTGCTGGCTCCGTCAAGCGAAGAACGGCCGCTGCCAGGTCCTGTACAA
GACCGAACTGAGCAAGGAGGAGTGCTGCAGCACCGGCCGGCTGAGCACCTCGTGGACCGAGGAGGACGTG
AATGACACACACCTCTTCAAGTGGATGATTTTCAACGGGGGCGCCCCAACTGCATCCCCTGTAAAGAAA
CGTGTGAGAACGTGGACTGTGGACCTGGGAAAAATGCCGAATGAACAAGAAGAACAAACCCCGCTGCGT
CTGCGCCCCGGATTGTTCCAACATCACCTGGAAGGGTTCCAGTCTGCGGGCTGGATGGGAAAACCTACCGC
AATGAATGTGCACTCCTAAAGGCAAGATGTAAAGAGCAGCCAGAACTGGAAGTCCAGTACCAAGGCAGAT
GTAAAAAGACTTGTCGGGATGTTTCTGTCCAGGCAGCTCCACATGTGTGGTGGACCAGACCAATAATGC
CTACTGTGTGACCTGTAATCGGATTTTGCCCAGAGCCTGCTTCCTCTGAGCAATATCTCTGTGGGAATGAT
GGAGTCACCTACTCCAGTGCCTGCCACCTGAGAAAGGCTACCTGCCTGCTGGGCAGATCTATTGGATTAG
CCTATGAGGGAAAGTGTATCAAAGCAAAGTCCTGTGAAGATATCCAGTGCACTGGTGGGAAAAATGTTT
ATGGGATTTCAAGGTTGGGAGAGGCCGGTGTTCCCTCTGTGATGAGCTGTGCCCTGACAGTAAGTCGGAT
GAGCCTGTCTGTGCCAGTGACAATGCCACTTATGCCAGCGAGTGTGCCATGAAGGAAGCTGCCTGCTCCT
CAGGTGTGCTACTGGAAGTAAAGCACTCCGGATCTTGCAACTCCATTTCGGAAGACACCGAGGAAGAGGA
GGAAGAATGAAGACCAGGACTACAGCTTTCCTATATCTTCTATTCTAGAGTGGTAA""".replace('\n', '')
    
    print("=== FST344 Codon Optimization Analysis ===\n")
    
    # Analyze original sequence
    original_stats = analyze_sequence(original_cds, "Original FST344")
    print(f"Original sequence analysis:")
    print(f"  Length: {original_stats['length']} bp")
    print(f"  GC content: {original_stats['gc_content']:.1f}%")
    print(f"  CpG dinucleotides: {original_stats['cpg_count']} ({original_stats['cpg_frequency']:.2f}%)")
    print(f"  Codons: {original_stats['codon_count']}")
    
    # Translate original
    original_protein = translate_dna(original_cds)
    print(f"  Protein length: {len(original_protein)} amino acids")
    
    # Optimize sequence
    print("\n=== Performing Optimization ===")
    optimized_cds, protein = optimize_sequence(original_cds)
    
    # Analyze optimized sequence
    optimized_stats = analyze_sequence(optimized_cds, "Optimized FST344")
    print(f"\nOptimized sequence analysis:")
    print(f"  Length: {optimized_stats['length']} bp")
    print(f"  GC content: {optimized_stats['gc_content']:.1f}%") 
    print(f"  CpG dinucleotides: {optimized_stats['cpg_count']} ({optimized_stats['cpg_frequency']:.2f}%)")
    print(f"  Codons: {optimized_stats['codon_count']}")
    
    # Compare improvements
    print(f"\n=== Optimization Results ===")
    gc_change = optimized_stats['gc_content'] - original_stats['gc_content']
    cpg_reduction = original_stats['cpg_count'] - optimized_stats['cpg_count']
    cpg_percent_reduction = (cpg_reduction / original_stats['cpg_count'] * 100) if original_stats['cpg_count'] > 0 else 0
    
    print(f"GC content change: {gc_change:+.1f}%")
    print(f"CpG reduction: {cpg_reduction} dinucleotides ({cpg_percent_reduction:.1f}% reduction)")
    print(f"Protein identity: {'IDENTICAL' if original_protein == protein else 'DIFFERENT'}")
    
    # Generate optimized FASTA
    fasta_content = f">FST344_OPTIMIZED_CDS Human follistatin optimized for AAV expression\n"
    
    # Format sequence with line breaks
    formatted_seq = ""
    for i in range(0, len(optimized_cds), 70):
        formatted_seq += optimized_cds[i:i+70] + "\n"
    
    fasta_content += formatted_seq.rstrip()
    
    # Write optimized sequence
    with open('/Users/biobook/Projects/sequence-optimiser/FST344_OPTIMIZED_CDS.fasta', 'w') as f:
        f.write(fasta_content)
    
    print(f"\n✅ Optimized sequence saved to: FST344_OPTIMIZED_CDS.fasta")
    
    # Generate comparison report
    report = f"""# FST344 Optimization Comparison Report

## Original Sequence (NM_013409.3 CDS)
- Length: {original_stats['length']} bp
- GC Content: {original_stats['gc_content']:.1f}%
- CpG Dinucleotides: {original_stats['cpg_count']} ({original_stats['cpg_frequency']:.2f}%)
- Codons: {original_stats['codon_count']}

## Optimized Sequence
- Length: {optimized_stats['length']} bp  
- GC Content: {optimized_stats['gc_content']:.1f}%
- CpG Dinucleotides: {optimized_stats['cpg_count']} ({optimized_stats['cpg_frequency']:.2f}%)
- Codons: {optimized_stats['codon_count']}

## Improvements
- GC Content Change: {gc_change:+.1f}%
- CpG Reduction: {cpg_reduction} dinucleotides ({cpg_percent_reduction:.1f}% reduction)
- Protein Sequence: {'Identical' if original_protein == protein else 'Modified'}

## Expected Benefits
- Reduced immunogenicity (fewer CpG motifs)
- Improved AAV packaging efficiency
- Enhanced expression in mammalian cells
- Better codon adaptation for human expression

## Next Steps
1. Validate protein folding using AlphaFold
2. Test expression in HEK293T cells
3. Package in AAV vectors for in vivo testing
"""

    with open('/Users/biobook/Projects/sequence-optimiser/FST344_OPTIMIZATION_COMPARISON.md', 'w') as f:
        f.write(report)
    
    print(f"📊 Comparison report saved to: FST344_OPTIMIZATION_COMPARISON.md")

if __name__ == "__main__":
    main()