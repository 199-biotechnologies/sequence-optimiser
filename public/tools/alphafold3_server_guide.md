# AlphaFold3 Server Usage Guide for FST344 Validation

## Overview
AlphaFold3 Server (alphafoldserver.com) allows modeling of protein complexes, protein-DNA/RNA interactions, and small molecule binding. This is perfect for validating our FST344 optimization.

## Current Status
- **Remaining jobs**: 30/30 daily quota
- **Token limit**: 5,000 per job
- **Capabilities**: Multi-molecular predictions

## How to Use AlphaFold3 Server for FST344

### Step 1: Access the Server
1. Go to https://alphafoldserver.com
2. Sign in with Google account (free, non-commercial use)
3. You get 30 jobs per day, refreshed daily

### Step 2: Set Up FST344 Prediction Job

#### Basic FST344 Structure Prediction:
1. **Entity Type**: Select "Protein"
2. **Copies**: 1
3. **Sequence Input**: Paste our optimized FST344 protein sequence

```
MVRAHQPGGLCLLLLLLCQFMEDRSAQAGNCWLRQAKNGRCRQLYKTELSKECCSTGRLSTSWTEEDVNDNTLFKWMIFNGAGPNCIPCKEETCENVDCGPGKKCRMMNKKNKPRCVCAPDCSNITAWKGPVCGLDGKTYRNECALLKAROCKEQELEVRQGRCKKCTCRDVFCPGSSTCVVDQTNNAYCCVTCNRICPEPASSEQYLCGNDGVTYSSACHLRKATCLLGRSIGLAEGKCIKAASCEDIRRCTGGKKCLDVKGRGRCSLDCELCPDSKSDEPVCASDNATYASECAMKEAACSSGGVLLVEVKHSGSCSISEDTEEEEDQDYSFPISSILEW
```

#### Advanced: FST344-Activin Complex Prediction:
1. **Add Entity**: Click "+ Add entity"
2. **Entity 1**: FST344 (our optimized sequence)
3. **Entity 2**: Activin A (if studying binding)
   - Activin A sequence: Add activin beta A chain
4. This will model how they interact together

### Step 3: Job Configuration
1. **Job Name**: "FST344_Optimized_Structure_Validation"
2. **Description**: "Validation of codon-optimized FST344 structure"
3. Click "Continue and preview job"

### Step 4: Review and Submit
1. Review the sequences
2. Check token usage (should be <1000 for single protein)
3. Click "Save job" to submit

### Step 5: Results Analysis
Results typically available in 1-10 minutes:
1. **Confidence scores**: Compare with original AlphaFold2 (AF-P19883-F1)
2. **3D structure**: Download PDB file
3. **Interface quality**: If modeling complexes

## Specific Use Cases for FST344 Project

### 1. Structure Validation
**Purpose**: Confirm optimized sequence folds identically
**Input**: Optimized FST344 protein sequence
**Expected**: High confidence, identical fold to AF-P19883-F1

### 2. Activin Binding Complex
**Purpose**: Validate binding site preservation
**Input**: FST344 + Activin A/B chains
**Expected**: Proper binding interface formation

### 3. Signal Peptide Analysis
**Purpose**: Model signal peptide cleavage
**Input**: Full FST344 with signal peptide + signal peptidase
**Expected**: Proper cleavage site prediction

### 4. Heparin Binding (if relevant)
**Purpose**: Model heparin interaction
**Input**: FST344 + heparin oligosaccharide
**Expected**: Binding site validation

## Token Calculation
- **Single protein (344 AA)**: ~344 tokens
- **Protein complex**: ~500-1500 tokens
- **With small molecules**: +50-200 tokens

## Comparison Strategy
1. **Run original FST344** (from AF-P19883-F1 sequence)
2. **Run optimized FST344** (our codon-optimized version)
3. **Compare structures** using PyMOL or ChimeraX
4. **RMSD calculation** to quantify differences

## Expected Results
- **Confidence**: >90% for well-folded regions
- **RMSD**: <1Å between original and optimized
- **Binding sites**: Preserved in complexes
- **Overall fold**: Identical to AF-P19883-F1

## Troubleshooting
- **Low confidence**: May indicate folding issues
- **High RMSD**: Could suggest optimization problems
- **Failed binding**: Might indicate lost function

## Integration with Our Project
1. **Download PDB files** → Add to `structures/` folder
2. **Compare with AF-P19883-F1** → Document differences
3. **Update validation report** → Include AF3 results
4. **Visualization** → Add to 3D viewer

## Next Steps After Results
1. **Structure comparison** using PyMOL alignment
2. **Confidence score analysis** 
3. **Binding site validation** if complexes modeled
4. **Documentation update** with findings

This will provide the most advanced validation available for our optimization!