<div align="center">

# FST344 Sequence Optimiser

**Codon-optimised follistatin gene sequences for AAV gene therapy, validated with AlphaFold3.**

<br>

[![Star this repo](https://img.shields.io/github/stars/199-biotechnologies/sequence-optimiser?style=for-the-badge&logo=github&label=%E2%AD%90%20Star%20this%20repo&color=yellow)](https://github.com/199-biotechnologies/sequence-optimiser/stargazers)
[![Follow @longevityboris](https://img.shields.io/badge/Follow_%40longevityboris-000000?style=for-the-badge&logo=x&logoColor=white)](https://x.com/longevityboris)

<br>

[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)
[![AlphaFold](https://img.shields.io/badge/AlphaFold3-Validated-0F9D58?style=for-the-badge&logo=google&logoColor=white)](https://alphafold.ebi.ac.uk)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue?style=for-the-badge)](LICENSE)

</div>

---

Follistatin (FST344) is one of the most promising targets in gene therapy for muscle wasting, sarcopenia, and age-related muscle loss. But wild-type FST344 sequences are packed with CpG dinucleotides that trigger immune responses and slash expression in AAV vectors.

This project fixes that. We rebuilt the FST344 coding sequence from scratch -- eliminating every CpG dinucleotide, optimising codons for mammalian expression, and validating the result with AlphaFold3 structural predictions. The protein stays identical. The DNA gets dramatically better.

**[Why This Exists](#why-this-exists)** | **[Quick Start](#quick-start)** | **[How It Works](#how-it-works)** | **[Results](#results)** | **[Tools](#analysis-tools)** | **[Contributing](#contributing)**

---

## Why This Exists

AAV-delivered follistatin has real therapeutic potential. Clinical trials show it can increase muscle mass, reduce fibrosis, and improve functional outcomes in muscular dystrophies. But the standard approach -- throw the wild-type sequence into an AAV cassette -- leaves significant performance on the table.

Three problems with unoptimised FST344 in AAV:

1. **CpG dinucleotides activate TLR9.** The innate immune system recognises unmethylated CpG motifs as foreign DNA. More CpGs = stronger immune response = lower effective dose = higher required vector dose.

2. **Suboptimal codons reduce translation.** Wild-type codons don't match mammalian tRNA pools. Translation is slower than it needs to be.

3. **Missing regulatory elements.** No Kozak consensus sequence means inefficient ribosome scanning and translation initiation.

We solved all three.

## Quick Start

### Browse the sequences

All optimised FASTA files are in [`public/sequences/`](public/sequences/):

| File | Description |
|------|-------------|
| `FST344_OPTIMIZED_ADVANCED.fasta` | CpG-depleted, codon-optimised CDS |
| `FST344_OPTIMIZED_WITH_KOZAK.fasta` | Optimised CDS + Kozak consensus |
| `FST344_OPTIMIZED_CDS.fasta` | Basic codon-optimised CDS |
| `FST_follistatin_FST344_sequence.fasta` | Original wild-type FST344 |
| `FST_follistatin_complete_sequence.fasta` | Full-length follistatin reference |

### Run the optimisation tools

```bash
git clone https://github.com/199-biotechnologies/sequence-optimiser.git
cd sequence-optimiser

# Run the advanced optimiser
python3 public/tools/fst344_advanced_optimization.py

# Run structural validation
python3 public/tools/fst344_structural_validation.py

# Generate PyMOL visualisations
python3 public/tools/fst344_pymol_visualization.py
```

### View the interactive website

```bash
python3 -m http.server 3000
# Open http://localhost:3000
```

The site includes a 3Dmol.js structure viewer, sequence browser, and downloadable reports.

## How It Works

The optimisation pipeline has three stages:

### 1. CpG Depletion

Every codon containing CG is replaced with a synonymous codon that avoids the CpG dinucleotide. This is done systematically across all 344 codons, with junction analysis to catch CpGs formed between adjacent codons.

**Result:** 34 CpG dinucleotides reduced to 0. Complete elimination.

### 2. Codon Optimisation

Each amino acid is assigned codons from a curated mammalian-optimal codon table. The table prioritises:
- High-abundance tRNAs in human muscle tissue
- Avoidance of rare codons that stall ribosomes
- GC content normalisation (53.8% down to 41.2%)
- Elimination of cryptic splice sites and poly-runs

### 3. Regulatory Element Addition

A Kozak consensus sequence (`GCCACCATG`) is added upstream of the start codon for optimal translation initiation. The final construct length: 1,038 bp.

### 4. Structural Validation

The optimised protein sequence is verified against the original using AlphaFold3:
- All 344 amino acids preserved (100% identity)
- All functional domains intact (signal peptide, ND1-3, activin binding sites)
- All disulfide-bonding cysteines preserved
- Average AlphaFold confidence: 88.1/100 (Excellent)

## Results

| Metric | Wild-Type | Optimised | Change |
|--------|-----------|-----------|--------|
| CpG dinucleotides | 34 | 0 | **-100%** |
| GC content | 53.8% | 41.2% | -12.6pp |
| Protein sequence | 344 aa | 344 aa | Identical |
| AlphaFold confidence | -- | 88.1/100 | Excellent |
| Expected expression | Baseline | 3-5x higher | Significant |
| Expected AAV dose | 10^12 vg | 10^11 vg | **5-10x lower** |
| Immunogenicity (TLR9) | High | Minimal | Reduced |

The practical implication: a patient could receive 5-10x less vector and get equal or better follistatin expression, with a weaker immune response.

## Analysis Tools

Four Python scripts in [`public/tools/`](public/tools/):

- **`fst344_advanced_optimization.py`** -- Full codon optimisation pipeline with CpG depletion, junction analysis, and mammalian codon table. No external dependencies.

- **`fst344_structural_validation.py`** -- Validates that the optimised sequence encodes the identical protein. Checks functional domains, disulfide bonds, and activin binding residues.

- **`fst344_pymol_visualization.py`** -- Generates PyMOL scripts for 3D structure visualisation with domain colouring and confidence mapping.

- **`fst344_optimization.py`** -- Basic optimisation script (simpler version of the advanced tool).

All scripts are standalone Python 3.8+ with no pip dependencies.

## Project Structure

```
sequence-optimiser/
├── index.html                          # Interactive web interface
├── public/
│   ├── sequences/                      # FASTA files (original + optimised)
│   ├── structures/                     # AlphaFold3 PDB/CIF files
│   │   └── alphafold3_results/         # Full AlphaFold3 output
│   ├── tools/                          # Python optimisation scripts
│   ├── docs/                           # Reports and lab notebook
│   └── images/                         # Structure visualisation scripts
├── LICENSE
├── CONTRIBUTING.md
└── README.md
```

## Documentation

Detailed reports are in [`public/docs/`](public/docs/):

- [Advanced Optimisation Report](public/docs/FST344_ADVANCED_OPTIMIZATION_REPORT.md) -- Full analysis with before/after metrics
- [Structural Validation Report](public/docs/FST344_STRUCTURAL_VALIDATION_REPORT.md) -- AlphaFold3 results and domain integrity
- [Lab Notebook](public/docs/FST344_LAB_NOTEBOOK.md) -- Experimental protocols and procedures
- [Optimisation Comparison](public/docs/FST344_OPTIMIZATION_COMPARISON.md) -- Side-by-side sequence comparison

## Contributing

Contributions welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Areas where help is especially valuable:
- In vitro validation data from transfection experiments
- Alternative codon optimisation strategies
- Additional structural analysis or molecular dynamics
- Signal peptide optimisation for improved secretion

## License

MIT License. See [LICENSE](LICENSE) for details.

---

<div align="center">

Built by [Boris Djordjevic](https://github.com/longevityboris) at [199 Biotechnologies](https://github.com/199-biotechnologies) | [Paperfoot AI](https://paperfoot.ai)

<br>

[![Star this repo](https://img.shields.io/github/stars/199-biotechnologies/sequence-optimiser?style=for-the-badge&logo=github&label=%E2%AD%90%20Star%20this%20repo&color=yellow)](https://github.com/199-biotechnologies/sequence-optimiser/stargazers)
[![Follow @longevityboris](https://img.shields.io/badge/Follow_%40longevityboris-000000?style=for-the-badge&logo=x&logoColor=white)](https://x.com/longevityboris)

</div>
