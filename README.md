# Cross-Kingdom miRNA Target Predictor

## 📖 Introduction
This repository provides an exploratory computational pre-screening framework designed specifically for predicting cross-kingdom regulation, particularly **plant-derived exosomal miRNAs (ELNVs) targeting gut bacterial mRNAs**. 

Unlike traditional miRNA prediction tools designed for eukaryotic systems (which rely on 3' UTR binding and Argonaut complex mechanisms), this pipeline is optimized for the unique constraints of cross-kingdom interactions.

## 🧬 Methodology
The pipeline identifies putative target candidates through a highly robust statistical framework:
* **Background Modeling**: Constructs a local 1st-order Markov Model to account for bacterial genomic GC bias and dinucleotide transition frequencies.
* **Statistical Testing**: Employs a **Poisson Process** approximation to evaluate the probability of rare seed motif occurrences.
* **Steric Constraints**: Strictly enforces `overlap = FALSE` during motif scanning to better satisfy the independent event assumption and approximate the physical steric hindrance of ribonucleoprotein complexes.
* **Multiple Testing**: Incorporates Benjamini-Hochberg (BH) FDR correction across the complete valid searchable universe.

## 🚀 Quick Start (Toy Dataset)
We provide a lightweight toy dataset (`test_target_regions.fasta`) for quick pipeline validation.
1. Download the `predict_targets.R` and `test_target_regions.fasta` files.
2. Ensure you have the `Biostrings` package installed in R:
   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install("Biostrings")
   ```
3. Run the script. The entire pipeline executes in seconds and outputs a .csv matrix containing the computed $\lambda$ expectations, Poisson P-values, and FDR correlations.

## Citation & Acknowledgements
The conceptual foundation of identifying plant-derived exosomal miRNA targets in gut bacteria was inspired by Teng et al., 2018 (Cell Host & Microbe).
This repository implements a highly refactored and statistically robust algorithm, featuring mathematical boundary continuity (underflow protection) and global multiple testing universe completeness.

## License
This project is licensed under the MIT License. You are free to use, modify, and distribute this software, provided that proper academic citation is given.

## 📬 Contact & Feedback
If you encounter any bugs, have questions about the methodology, or want to discuss potential collaborations, please feel free to open an **Issue** here on GitHub or contact me directly at: 
**pengxu0907@gmail.com** 
