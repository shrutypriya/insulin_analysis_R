# Bioinformatics-R: Sequence-to-Structure Analysis of Insulin

A beginner-friendly computational biology project built in R to bridge the gap between 1D sequence data and 3D protein structures. This workflow automates data retrieval from major biological repositories and performs structural flexibility analysis.

## Project Overview
Using **Insulin** as a model protein, this project demonstrates how to programmatically interact with global biological databases to download, process, and visualize structural data. 



### Key Features:
* **Automated NCBI Fetching:** Downloads cross-species protein sequences using accession numbers.
* **PDB Parsing:** Imports and extracts metadata (chains, residues, heteroatoms) from `.pdb` files.
* **Flexibility Analysis:** Computes and plots residue-level B-factors to evaluate structural rigidity.
* **Interactive 3D Visuals:** Generates rotatable 3D ribbon models directly in the R environment.

---

## 📁 Repository Structure
```text
├── data/                  # Automated download directory for .fasta and .pdb files
├── 01_data_fetching.R     # Complete executable R script 
└── README.md              # Project documentation
