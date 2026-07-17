# Protein Cluster Analysis and Disorder Prediction Pipeline

## Overview

This project provides a comprehensive bioinformatics pipeline for analyzing conserved protein clusters across species and predicting intrinsically disordered regions. The workflow integrates MongoDB for cluster management, IUPred2A for disorder prediction, and phylogenetic analysis tools.

## Project Structure

```
project/
├── src/                          # Source code directory
│   ├── script1.ipynb             # Get conserved clusters (intra-species)
│   ├── script2.ipynb             # Get conserved clusters (inter-species)
│   ├── script3.ipynb             # Protein disorder prediction analysis
│   ├── script4.py                # Amino acid composition analysis
│	└── script5.ipynb             # Create annotated FASTA files with CCC proteins and GO annotations
├── data/                          # Original data directory
│   ├── original_data.json        # Raw dataset
└── results/                       # Final output and results
    ├── IUPRED2A/                 # Disorder prediction results
    │   ├── NCCC_inter_prot_id.result
    │   ├── NCCC_intra_protein_id.result
    │   ├── PCCC_inter_protein_id.result
    │   └── PCCC_intra_protein_id.result
    ├── IQTree/                   # Phylogenetic analysis results
        ├── NCCC_inter_proteins.afa.iqtree
		├── NCCC_intra_protein.aln.iqtree
        ├── PCCC_inter_proteins.afa.iqtree
        └── protein_PCCC_intra.afa.iqtree
```

## Prerequisites

### Required Software
- Python 3.8 or higher
- MongoDB (running on localhost:27017)
- IUPred2A (for disorder predictions)
- MAFFT (for sequence alignment)
- IQTree (for phylogenetic analysis)

### Python Dependencies
```bash
pip install pymongo pandas matplotlib
```

## Script Descriptions

### 1. Script 1: Intra-Species Conserved Clusters (`script1.ipynb`)

**Purpose:** Identifies conserved protein clusters within the same species.

**Database Requirements:**
- MongoDB database: `LPM_w15_ncc` (nuclear proteome)
- Output database: `Conserved_cluster_intra`

**Key Operations:**
- Iterates through species collections
- Identifies proteins with shared cluster IDs
- Groups proteins by cluster within each species
- Stores results with species and cluster information

**Input:** MongoDB collections containing protein cluster data  
**Output:** Conserved cluster documents organized by species

---

### 2. Script 2: Inter-Species Conserved Clusters (`script2.ipynb`)

**Purpose:** Identifies conserved protein clusters across different species.

**Database Requirements:**
- Input databases:
  - `GAM_w15_ncc` (Mitochondrial - GAM)
  - `GANM_w25_ncc` (Nuclear - GANM)
  - `LPM_w15_ncc` (Nuclear - LPM)
  - `LPNM_w20_ncc` (Nuclear - LPNM)
  - `OPM_w30_ncc` (Mitochondrial - OPM)
  - `OPNM_w20_ncc` (Nuclear - OPNM)
- Output database: `Conserved_cluster_inter`

**Key Operations:**
- Cross-references clusters between different databases
- Identifies clusters conserved across multiple species
- Stores cluster-to-species relationships
- Handles case-insensitive species matching

**Input:** Multiple MongoDB databases with species-specific collections  
**Output:** Inter-species conserved cluster documents

---

### 3. Script 3: Protein Disorder Prediction Analysis (`script3.ipynb`)

**Purpose:** Analyzes intrinsically disordered regions in proteins using IUPred2A predictions.

**Prerequisites:** Script 5 should be run first to generate required FASTA files.

**Input Data:**
- IUPred2A output file: `NCCC_intra_protein_id.result` (from IUPred2A predictions)
- FASTA file: `NCCC_intra_pos.fasta` (with position information, created by Script 5)

**Processing Steps:**
1. **Parse IUPred2A Output:**
   - Extracts disorder scores (0-1, higher = more disordered)
   - Extracts anchor scores (0-1, anchoring residues)
   - Associates scores with residue positions

2. **Parse FASTA File:**
   - Extracts protein identifiers
   - Extracts sequence coordinates
   - Reconstructs protein names from header information

3. **Generate Visualizations:**
   - Creates line plots for each protein
   - Plots IUPRED_SCORE (blue line)
   - Plots ANCHOR_SCORE (green line)
   - Highlights regions of interest (orange shaded areas)
   - Includes disorder threshold line (y=0.5, red dashed)

**Output:**
- `fasta.csv` - Parsed FASTA data
- `bar_plots/` directory containing PNG plots (1000 DPI, publication quality)
  - Files named: `{protein_name}.png`

**Key Parameters:**
- Disorder threshold: 0.5
- Plot resolution: 1000 DPI
- Font: Arial, 16pt

---

### 4. Script 4: Amino Acid Composition Analysis (`script4.py`)

**Purpose:** Calculates amino acid composition for conserved clusters across all species.

**Prerequisites:** Scripts 1 & 2 must be run first to populate conserved cluster databases.

**Database Requirements:**
- Input databases: All 6 proteome databases
- Conserved cluster database: `Conserved_cluster_inter`

**Processing:**
- Retrieves sequences for each conserved cluster
- Merges sequences from all species
- Calculates percentage composition for each amino acid
- Generates statistics on cluster distribution

**Output:** 
- `amino_acid_composition.csv` containing:
  - Cluster ID
  - Number of sequences
  - Total length
  - Percentage of each amino acid (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)

**Standard Amino Acids Analyzed:**
```
ACDEFGHIKLMNPQRSTVWY
```
---

### 5. Script 5: CCC Protein Extraction and FASTA File Creation (`script5.ipynb`)

**Purpose:** Extracts conserved charged cluster (CCC) proteins from MongoDB databases and creates annotated FASTA files enriched with Gene Ontology (GO) molecular function annotations.

**Prerequisites:** Scripts 1 and 2 must be run first to populate the conserved cluster databases.

**Database Requirements:**
- Proteome databases:
  - `GANM` (Gene - Nuclear/Mitochondrial)
  - `LPM` (Legume - Mitochondrial)
  - `LPNM` (Legume - Nuclear/Mitochondrial)
- CCC cluster databases:
  - `GANM_w25_ncc`, `LPM_w15_ncc`, `LPNM_w20_ncc` (Nuclear charged clusters)
  - `GANM_w25_pccc`, `LPM_w15_pccc`, `LPNM_w20_pccc` (Phosphorylation charged clusters)
- Conserved cluster databases (populated by Scripts 1 & 2):
  - `Conserved_cluster_inter` (inter-species clusters)
  - `Conserved_cluster_intra` (intra-species clusters)

**GO Annotation Files Required:**
- `GO_summary_NCCC_inter_proteins.xlsx` — Inter-species nuclear CCC annotations
- `GO_summary_NCCC_intra_proteins.xlsx` — Intra-species nuclear CCC annotations
- `GO_summary_PCCC_inter_proteins.xlsx` — Inter-species phosphorylation CCC annotations
- `GO_summary_PCCC_intra_proteins.xlsx` — Intra-species phosphorylation CCC annotations

**Key Operations:**
1. **Database Connection:** Establishes connections to proteome and CCC cluster databases
2. **GO Annotation Loading:** Reads molecular function annotations from Excel files
3. **Protein Extraction:** Retrieves CCC proteins across four categories:
   - NCCC inter-species (charged clusters conserved across species)
   - NCCC intra-species (charged clusters within species)
   - PCCC inter-species (phosphorylation-related charged clusters across species)
   - PCCC intra-species (phosphorylation-related charged clusters within species)
4. **Header Enrichment:** Builds comprehensive FASTA headers with:
   - Database and species family information (GANM/LPM/LPNM_Genus_species)
   - Protein ID
   - Molecular function annotations from GO
   - Protein description and name
5. **FASTA File Generation:** Writes annotated sequences to output files

**FASTA Header Format:**
```
>GANM_Genus_species_ProteinID | MF:Function1;Function2 | Desc:Description | Name:ProteinName
MKTSVISLALLTLAGCQAKVEQAER...
```

**Output Files:**
- `NCCC_inter_proteins_annotated.fasta` — Inter-species NCCC proteins with sequences
- `NCCC_intra_proteins_annotated.fasta` — Intra-species NCCC proteins with sequences
- `PCCC_inter_proteins_annotated.fasta` — Inter-species PCCC proteins with sequences
- `PCCC_intra_proteins_annotated.fasta` — Intra-species PCCC proteins with sequences

**Key Functions:**
- `load_go_annotations(excel_path)` — Loads GO molecular function annotations from Excel
- `get_species_list(proteome_db)` — Retrieves all species from a database
- `find_protein_by_id(proteome_dbs, protein_id)` — Searches for protein across all databases
- `extract_ccc_proteins(collection, proteome_dbs, category_type)` — Extracts CCC proteins with statistics
- `build_fasta_header()` — Creates enriched FASTA header with annotations
- `write_fasta_file()` — Writes annotated sequences to FASTA file
- `print_stats()` — Displays statistical summary

**Important Notes:**
- This script MUST be run AFTER Scripts 1 & 2 (uses their database output)
- Update GO annotation file paths in cells 5-8 to match your file locations
- Requires valid MongoDB connections with populated cluster databases
- Generates publication-ready, annotated FASTA sequences for downstream analysis
- Produces statistics useful for methods sections in publications

---

## Running the Pipeline

### Step 1: MongoDB Setup
Ensure MongoDB is running:
```bash
mongod --dbpath /path/to/data/directory
```

### Step 2: Run Cluster Analysis
```bash
# Intra-species clusters
jupyter notebook src/script1.ipynb

# Inter-species clusters
jupyter notebook src/script2.ipynb
```

### Step 3: Create Annotated FASTA Files
```bash
jupyter notebook src/script5.ipynb
```
This generates the annotated FASTA files required by downstream analyses.

### Step 4: Run IUPred2A Predictions
Generate disorder predictions on the FASTA files created by Script 5:
```bash
iupred2a NCCC_intra_proteins_annotated.fasta long > results/NCCC_intra_protein_id.result
iupred2a NCCC_inter_proteins_annotated.fasta long > results/NCCC_inter_protein_id.result
```

### Step 5: Generate Disorder Analysis Visualizations
```bash
jupyter notebook src/script3.ipynb
```

### Step 6: Phylogenetic Analysis with IQTree
Run multiple sequence alignment and phylogenetic tree construction on the annotated FASTA files:
```bash
# For each protein category, first align sequences
mafft NCCC_inter_proteins_annotated.fasta > NCCC_inter_proteins.afa

# Then construct phylogenetic trees
iqtree -s NCCC_inter_proteins.afa -m MFP -bb 1000 -alrt 1000
iqtree -s NCCC_intra_proteins.afa -m MFP -bb 1000 -alrt 1000
iqtree -s PCCC_inter_proteins.afa -m MFP -bb 1000 -alrt 1000
iqtree -s PCCC_intra_proteins.afa -m MFP -bb 1000 -alrt 1000
```

### Step 7: Amino Acid Composition Analysis
```bash
python src/script4.py
```

## Expected Outputs

### From Script 5 (FASTA File Creation):
The following annotated FASTA files are generated:
```
results/
├── NCCC_inter_proteins_annotated.fasta
├── NCCC_intra_proteins_annotated.fasta
├── PCCC_inter_proteins_annotated.fasta
└── PCCC_intra_proteins_annotated.fasta
```

### From IUPred2A Analysis:
Disorder prediction output files:
```
results/IUPRED2A/
├── NCCC_inter_prot_id.result (6.9 MB)
├── NCCC_intra_protein_id.result (11.4 MB)
├── PCCC_inter_protein_id.result (882 KB)
└── PCCC_intra_protein_id.result (1.2 MB)
```
Each `.result` file contains per-residue disorder and anchor scores for the corresponding protein set.

### From Script 3 (Disorder Prediction Visualization):
```
Each plot shows disorder prediction profiles with highlighted regions of interest and disorder threshold lines.
```
### From IQTree Phylogenetic Analysis:
Phylogenetic tree files for each protein category:
```
results/IQTree/
├── NCCC_inter_proteins.afa.iqtree (320 KB)
├── NCCC_intra_protein.aln.iqtree (536 KB)
├── PCCC_inter_proteins.afa.iqtree (57 KB)
└── protein_PCCC_intra.afa.iqtree (59 KB)
```
These files contain the phylogenetic relationships and evolutionary analyses for each CCC protein category.

### From Script 4 (Amino Acid Composition):
```csv
cluster,n_sequences,total_length,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
...
```
Summary statistics for amino acid composition across all conserved clusters.

## Data Requirements

### For Disorder Prediction (Script 3):
- IUPred2A output format:
  ```
  >protein_name
  # Position  AA  IUPred_Score  Anchor_Score
  1  M  0.123  0.045
  2  T  0.234  0.056
  ...
  ```

- FASTA format with position information:
  ```
  >protein_id_cluster_start_position
  MVTLSPADKTNVIRAAQNCYDITPEeisevkdkskqvvvkgksklf...
  ```

## Troubleshooting

### MongoDB Connection Issues
- Verify MongoDB is running on localhost:27017
- Check connection string in scripts
- Ensure database and collection names match your setup

### IUPred2A File Not Found
- Verify file paths in script3.ipynb
- Check file exists at specified location
- Update file paths for your system

## Author Notes

This pipeline was developed for analyzing conserved protein clusters with emphasis on:
- Cross-species conservation patterns
- Disorder prediction across multiple proteomes
- Amino acid composition analysis
- Publication-quality visualizations

## Related Publications

This pipeline was developed to support the analysis presented in:

**Nebli, et al. (2026).** "Conservation and functional significance of charged clusters in mitochondria-located proteins of green plants." *Journal of Molecular Evolution.*

The pipeline facilitates the identification of conserved protein clusters across species and their characterization through disorder prediction and compositional analysis.

## License

See LICENSE file in project root.

## Citation

If you use this pipeline in your research, please cite the related publication:

Nebli, S., Rebai, A. & Ayadi, I. Sequence Conservation and Functional Significance of Charged Clusters in Mitochondria-Located Proteins of Green Plants. J Mol Evol (2026). https://doi.org/10.1007/s00239-026-10335-2
```
---
Additionally, please acknowledge the following tools and libraries:
- **IUPred2A:** Mészáros et al. (2018) for intrinsic disorder prediction
- **MongoDB:** https://www.mongodb.com/
- **Matplotlib:** https://matplotlib.org/
- **Pandas:** https://pandas.pydata.org/

---

**Last Updated:** July 17, 2026  
**Python Version:** 3.8+  
**Status:** Active Development
