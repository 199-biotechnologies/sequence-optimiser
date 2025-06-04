#!/usr/bin/env python3
"""
FST344 Structural Validation Tool
Analyzes protein structure, binding sites, and validates optimization
"""

import re
from collections import defaultdict

# Standard genetic code for translation
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

def parse_pdb_sequence(pdb_file):
    """Extract protein sequence from PDB file"""
    sequence = ""
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('SEQRES'):
                    # Parse amino acids from SEQRES lines
                    parts = line.split()
                    for i in range(4, len(parts)):
                        if parts[i] in aa_map:
                            sequence += aa_map[parts[i]]
    except FileNotFoundError:
        print(f"PDB file {pdb_file} not found")
        return ""
    
    return sequence

def extract_confidence_scores(pdb_file):
    """Extract AlphaFold confidence scores from PDB file"""
    confidence_scores = []
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[13:15] == 'CA':  # C-alpha atoms
                    # B-factor column contains confidence score
                    confidence = float(line[60:66].strip())
                    confidence_scores.append(confidence)
    except FileNotFoundError:
        print(f"PDB file {pdb_file} not found")
        return []
    
    return confidence_scores

def analyze_functional_domains():
    """Define known functional domains in follistatin"""
    domains = {
        'signal_peptide': (1, 29),
        'nd1_domain': (30, 120),     # N-terminal domain 1
        'nd2_domain': (121, 205),    # N-terminal domain 2  
        'nd3_domain': (206, 290),    # N-terminal domain 3
        'c_terminal': (291, 344),    # C-terminal domain
        'activin_binding_1': (58, 85),    # Key activin binding region 1
        'activin_binding_2': (145, 175),  # Key activin binding region 2
        'activin_binding_3': (220, 250),  # Key activin binding region 3
    }
    
    critical_residues = {
        # Critical cysteines for disulfide bonds
        'disulfide_cysteines': [39, 43, 56, 69, 90, 104, 120, 135, 148, 162, 177, 191, 205, 218, 231, 245, 260, 274],
        # Known activin contact residues (based on literature)
        'activin_contacts': [65, 67, 83, 85, 152, 155, 172, 174, 225, 227, 245, 247],
    }
    
    return domains, critical_residues

def validate_structure_integrity(original_seq, optimized_seq, confidence_scores):
    """Validate that optimization doesn't affect critical structural elements"""
    
    if len(original_seq) != len(optimized_seq):
        return {"error": "Sequence lengths don't match"}
    
    domains, critical_residues = analyze_functional_domains()
    
    validation_results = {
        'identical_sequences': original_seq == optimized_seq,
        'length_match': len(original_seq) == len(optimized_seq),
        'domain_analysis': {},
        'critical_residues_intact': True,
        'confidence_analysis': {},
        'recommendations': []
    }
    
    # Check each functional domain
    for domain_name, (start, end) in domains.items():
        start_idx = start - 1  # Convert to 0-based indexing
        end_idx = end
        
        if end_idx <= len(original_seq):
            original_domain = original_seq[start_idx:end_idx]
            optimized_domain = optimized_seq[start_idx:end_idx]
            domain_identical = original_domain == optimized_domain
            
            validation_results['domain_analysis'][domain_name] = {
                'positions': f"{start}-{end}",
                'identical': domain_identical,
                'original': original_domain,
                'optimized': optimized_domain
            }
            
            if not domain_identical:
                validation_results['critical_residues_intact'] = False
    
    # Check critical residues
    for residue_type, positions in critical_residues.items():
        mismatches = []
        for pos in positions:
            if pos <= len(original_seq):
                if original_seq[pos-1] != optimized_seq[pos-1]:
                    mismatches.append(pos)
        
        validation_results[f'{residue_type}_intact'] = len(mismatches) == 0
        if mismatches:
            validation_results[f'{residue_type}_mismatches'] = mismatches
    
    # Analyze AlphaFold confidence scores
    if confidence_scores:
        avg_confidence = sum(confidence_scores) / len(confidence_scores)
        high_confidence_regions = [i for i, score in enumerate(confidence_scores) if score > 70]
        low_confidence_regions = [i for i, score in enumerate(confidence_scores) if score < 50]
        
        validation_results['confidence_analysis'] = {
            'average_confidence': avg_confidence,
            'high_confidence_count': len(high_confidence_regions),
            'low_confidence_count': len(low_confidence_regions),
            'high_confidence_regions': high_confidence_regions[:10],  # First 10
            'low_confidence_regions': low_confidence_regions[:10]     # First 10
        }
    
    # Generate recommendations
    if validation_results['identical_sequences']:
        validation_results['recommendations'].append("✅ SAFE: Protein sequences are identical - optimization preserved all amino acids")
    else:
        validation_results['recommendations'].append("❌ ERROR: Protein sequences differ - optimization changed amino acids!")
    
    if validation_results.get('disulfide_cysteines_intact', True):
        validation_results['recommendations'].append("✅ SAFE: All disulfide bond cysteines preserved")
    else:
        validation_results['recommendations'].append("❌ CRITICAL: Disulfide bond cysteines altered - protein folding compromised")
    
    if validation_results.get('activin_contacts_intact', True):
        validation_results['recommendations'].append("✅ SAFE: Activin binding residues preserved")
    else:
        validation_results['recommendations'].append("⚠️ WARNING: Activin binding residues altered - test binding affinity")
    
    return validation_results

def generate_pymol_visualization_script(pdb_file, critical_residues):
    """Generate PyMOL script to visualize critical regions"""
    
    script = f"""# PyMOL script for FST344 structural analysis
# Load structure
load {pdb_file}, fst344

# Basic representation
hide all
show cartoon, fst344
color lightblue, fst344

# Highlight functional domains
select nd1_domain, resi 30-120
select nd2_domain, resi 121-205  
select nd3_domain, resi 206-290
select c_terminal, resi 291-344

color yellow, nd1_domain
color orange, nd2_domain
color red, nd3_domain
color purple, c_terminal

# Highlight activin binding sites
select binding_site_1, resi 58-85
select binding_site_2, resi 145-175
select binding_site_3, resi 220-250

color green, binding_site_1
color green, binding_site_2  
color green, binding_site_3
show sticks, binding_site_1 or binding_site_2 or binding_site_3

# Highlight critical cysteines
select disulfide_cysteines, resi 39+43+56+69+90+104+120+135+148+162+177+191+205+218+231+245+260+274
color yellow, disulfide_cysteines
show spheres, disulfide_cysteines

# Set view
orient
zoom

# Labels
label disulfide_cysteines and name CA, "%s%s" % (resn,resi)

print "FST344 structure loaded with functional domains highlighted"
print "Green: Activin binding sites"
print "Yellow spheres: Disulfide cysteines"
print "Color domains: ND1(yellow), ND2(orange), ND3(red), C-term(purple)"
"""
    
    return script

def main():
    print("=== FST344 Structural Validation Analysis ===\n")
    
    # Load sequences
    print("Loading sequences...")
    
    # Original CDS
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
    
    # Load optimized sequence
    try:
        with open('/Users/biobook/Projects/sequence-optimiser/FST344_OPTIMIZED_ADVANCED.fasta', 'r') as f:
            content = f.read()
            optimized_cds = ''.join(content.split('\n')[1:])  # Skip header line
    except FileNotFoundError:
        print("Error: Optimized sequence file not found")
        return
    
    # Translate sequences
    original_protein = translate_dna(original_cds)
    optimized_protein = translate_dna(optimized_cds)
    
    print(f"Original protein length: {len(original_protein)} AA")
    print(f"Optimized protein length: {len(optimized_protein)} AA")
    
    # Load AlphaFold structure
    print("\nAnalyzing AlphaFold structure...")
    alphafold_sequence = parse_pdb_sequence('/Users/biobook/Projects/sequence-optimiser/alphafold_fst_structure.pdb')
    confidence_scores = extract_confidence_scores('/Users/biobook/Projects/sequence-optimiser/alphafold_fst_structure.pdb')
    
    print(f"AlphaFold sequence length: {len(alphafold_sequence)} AA")
    print(f"Confidence scores extracted: {len(confidence_scores)} residues")
    
    # Validate sequences match AlphaFold
    if alphafold_sequence:
        alphafold_match = original_protein == alphafold_sequence
        print(f"Original sequence matches AlphaFold: {alphafold_match}")
        if not alphafold_match:
            print("WARNING: Original protein sequence doesn't match AlphaFold structure!")
    
    # Perform structural validation
    print("\n=== Structural Validation Results ===")
    validation = validate_structure_integrity(original_protein, optimized_protein, confidence_scores)
    
    # Print results
    print(f"Sequences identical: {validation['identical_sequences']}")
    print(f"Length match: {validation['length_match']}")
    print(f"Critical residues intact: {validation['critical_residues_intact']}")
    
    if 'confidence_analysis' in validation:
        conf = validation['confidence_analysis']
        print(f"Average AlphaFold confidence: {conf['average_confidence']:.1f}")
        print(f"High confidence regions (>70): {conf['high_confidence_count']}")
        print(f"Low confidence regions (<50): {conf['low_confidence_count']}")
    
    # Print domain analysis
    print("\n=== Functional Domain Analysis ===")
    for domain, info in validation['domain_analysis'].items():
        status = "✅" if info['identical'] else "❌"
        print(f"{status} {domain}: {info['positions']} - {'Identical' if info['identical'] else 'CHANGED'}")
    
    # Print recommendations
    print("\n=== Safety Assessment ===")
    for rec in validation['recommendations']:
        print(rec)
    
    # Generate PyMOL script
    print("\n=== Generating Visualization Files ===")
    domains, critical_residues = analyze_functional_domains()
    pymol_script = generate_pymol_visualization_script('alphafold_fst_structure.pdb', critical_residues)
    
    with open('/Users/biobook/Projects/sequence-optimiser/fst344_pymol_visualization.py', 'w') as f:
        f.write(pymol_script)
    
    print("✅ PyMOL visualization script saved: fst344_pymol_visualization.py")
    
    # Generate comprehensive structural report
    report = f"""# FST344 Structural Validation Report

## Summary
Comprehensive structural analysis of optimized FST344 sequence using AlphaFold structure validation.

## Sequence Validation
- **Original protein**: {len(original_protein)} amino acids
- **Optimized protein**: {len(optimized_protein)} amino acids  
- **Sequences identical**: {validation['identical_sequences']}
- **AlphaFold match**: {alphafold_sequence == original_protein if alphafold_sequence else 'N/A'}

## AlphaFold Confidence Analysis
"""
    
    if 'confidence_analysis' in validation:
        conf = validation['confidence_analysis']
        report += f"""- **Average confidence**: {conf['average_confidence']:.1f}/100
- **High confidence regions (>70)**: {conf['high_confidence_count']}/{len(confidence_scores)}
- **Low confidence regions (<50)**: {conf['low_confidence_count']}/{len(confidence_scores)}
- **Overall assessment**: {'EXCELLENT' if conf['average_confidence'] > 70 else 'GOOD' if conf['average_confidence'] > 50 else 'POOR'}
"""
    
    report += f"""
## Functional Domain Integrity
"""
    
    for domain, info in validation['domain_analysis'].items():
        status = "✅ PRESERVED" if info['identical'] else "❌ ALTERED"
        report += f"- **{domain}** ({info['positions']}): {status}\n"
    
    report += f"""
## Critical Residue Analysis
- **Disulfide cysteines**: {'✅ All preserved' if validation.get('disulfide_cysteines_intact', True) else '❌ ALTERED - CRITICAL ISSUE'}
- **Activin binding sites**: {'✅ All preserved' if validation.get('activin_contacts_intact', True) else '⚠️ Some altered - test binding'}

## Safety Assessment
"""
    
    for rec in validation['recommendations']:
        report += f"- {rec}\n"
    
    report += f"""
## Experimental Validation Required
1. **Protein folding assays**: Circular dichroism, dynamic light scattering
2. **Activin binding assays**: Surface plasmon resonance, competitive binding
3. **Functional assays**: Activin neutralization, FSH inhibition
4. **Structural validation**: NMR, X-ray crystallography (if critical)

## Files Generated
- `fst344_pymol_visualization.py` - PyMOL script for structure visualization
- `alphafold_fst_structure.pdb` - AlphaFold structure file
- This validation report

## Conclusion
"""
    
    if validation['identical_sequences']:
        report += "**SAFE TO PROCEED**: Optimization preserved all amino acids and functional domains. No structural concerns identified."
    else:
        report += "**REQUIRES REVIEW**: Optimization altered amino acid sequence. Structural validation needed before proceeding."
    
    with open('/Users/biobook/Projects/sequence-optimiser/FST344_STRUCTURAL_VALIDATION_REPORT.md', 'w') as f:
        f.write(report)
    
    print("📊 Structural validation report saved: FST344_STRUCTURAL_VALIDATION_REPORT.md")
    
    # Save detailed validation data
    import json
    with open('/Users/biobook/Projects/sequence-optimiser/fst344_validation_data.json', 'w') as f:
        json.dump(validation, f, indent=2)
    
    print("📊 Detailed validation data saved: fst344_validation_data.json")

if __name__ == "__main__":
    main()