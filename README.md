# TFG_2026 — Bioinformatics Scripts

This repository contains the bioinformatics scripts developed for the analysis of CRISPR spacer acquisition data in Collector Plasmid experiments, as part of the Bachelor's Thesis **"Optimization of Collector Plasmids, a tool to monitor horizontal gene transfer in bacterial communities"** (Universidad de Cantabria, 2025–2026). Scripts process raw sequencing output, perform spacer extraction and alignment, reconstruct conjugative trajectories, and analyze spacer mutations. Some scripts were developed by Andrea Fernández as part of earlier work in the Llosa laboratory; adapted versions and new scripts are indicated accordingly.

---

## Scripts

### `Spacer_extraction.py`
This script processes raw FASTQ files from CRISPR array sequencing to build a non-redundant library of spacers. Developed by Andrea Fernández.

---

### `SAM_to_TSV.py`
Python script that converts a SAM alignment file into a TSV matrix, automatically grouping CRISPR spacers by their belonging array, calculating array lengths, and merging multi-mapping hits (/) of equivalent alignment quality based on their mismatch scores.

---

### `Chord_diagram_adapted.R`
This R script reconstructs the movement trajectories of the Collector Plasmid between different bacterial hosts. Compared to the original (Andrea Fernández) genus-level script, the adapted version narrows the taxonomic resolution to five specific target species and introduces strict filtering logic ('walls and tunnels') to eliminate ambiguous mappings and bypass self-plasmid sequences before calculating trajectories.

---

### `Mutation_count.py`
This script screens a FASTQ file to locate a target CRISPR spacer sequence with a customizable Levenshtein distance, outputting a detailed report of identified sequence variants, their exact mutation patterns, and overall frequencies.

---

### `Mutation_table.py`
This script screens a FASTQ file for a specific spacer reference using sliding-window Levenshtein matching and maps the exact genomic positions of any mutations to generate an Excel file (.xlsx) of position-by-position variability.

---

### `Spacer_extraction_adapted.py`
This adapted script introduces a Levenshtein-distance screening step to identify arrays containing a specific terminal spacer sequence, while preserving the original repeat-matching and data-filtering structure.

---

### `SAM_to_TSV_adapted.py`
Based on the original `SAM_to_TSV.py`, this adapted script integrates an additional sequence-screening step using Levenshtein distance to check reads against two specific target sequences of interest, automatically prefixing matching outputs with an `@` symbol flag.

---

### `Sankey_diagrams.R`
This R script generates Sankey diagrams to visualize the downstream spacer acquisition trajectories inferred from sequencing data after a specific target spacer. It represents ordered host-to-host transitions and stops propagation upon reaching unknown spacers or ambiguous multi-species mappings. Putative intra-species transfers are considered, although they cannot be unambiguously distinguished from double acquisitions within the same cell.

---

### `LCA_algorithm.py`
When analyzing acquisitions from complex communities (e.g., wastewater), this script resolves multi-species hits from BLASTn searches.LCA Logic: It evaluates hits sharing the highest bitscore and assigns the taxonomy to the Lowest Common Ancestor (LCA). Developed by Andrea Fernández.
