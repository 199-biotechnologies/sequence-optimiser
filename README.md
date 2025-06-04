# FST344 Follistatin Optimization Project Website

## Overview
Interactive website showcasing the complete FST344 follistatin sequence optimization project for AAV gene therapy applications.

## Features
- **Interactive 3D Structure Viewer** using 3Dmol.js
- **Complete Project Documentation** with downloadable files
- **Optimization Results** with before/after comparisons
- **Validation Reports** and structural analysis
- **Responsive Design** optimized for all devices

## Project Highlights
- ✅ **100% CpG Elimination** (34 → 0 dinucleotides)
- ✅ **5-10x Dosage Reduction** expected (10¹² → 10¹¹ vg)
- ✅ **3-5x Expression Increase** predicted
- ✅ **Identical Protein Sequence** (344 amino acids preserved)
- ✅ **AlphaFold Validated** (88.1% confidence)

## Deployment

### Local Development
```bash
cd website
python3 -m http.server 3000
# Visit http://localhost:3000
```

### Vercel Deployment
```bash
cd website
vercel --prod
```

## File Structure
```
website/
├── index.html              # Main website
├── public/
│   ├── sequences/          # FASTA files and optimized sequences
│   ├── structures/         # PDB structure files
│   ├── tools/             # Python optimization scripts
│   └── docs/              # Markdown documentation
├── package.json
├── vercel.json
└── README.md
```

## Technologies Used
- **3Dmol.js** - Interactive molecular visualization
- **Plotly.js** - Data visualization and charts
- **Vanilla JavaScript** - Interactive features
- **CSS Grid/Flexbox** - Responsive layout
- **Vercel** - Static site hosting

## Documentation Access
All documentation files are accessible via direct links:
- Lab Notebook: `/public/docs/FST344_LAB_NOTEBOOK.md`
- Optimization Report: `/public/docs/FST344_ADVANCED_OPTIMIZATION_REPORT.md`
- Structural Validation: `/public/docs/FST344_STRUCTURAL_VALIDATION_REPORT.md`
- Literature Review: `/public/docs/literature_review_follistatin_aav_optimization.md`

## Downloads Available
- Optimized FASTA sequences
- AlphaFold PDB structure file
- Python optimization scripts
- Complete analysis data

## Browser Support
- Chrome/Edge (recommended)
- Firefox
- Safari
- Mobile browsers

## License
MIT License - See project documentation for details.