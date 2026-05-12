# Long-distance dispersal drives global tropical distributions in a widespread moth lineage (Lepidoptera: Limacodidae)

This repository accompanies the project:

> **"Long-distance dispersal drives global tropical distributions in a widespread moth lineage (Lepidoptera: Limacodidae)"**

by **Taberer et al.**  
Contact: `tabitha.taberer@biology.ox.ac.uk` / `tabrtab@aol.com`

The project explores the evolutionary history of the *Parasa* species complex through the construction of a global phylogeny primarily based on museum specimens, followed by biogeographical analyses. The phylogeny is based on complete and fragmented BUSCO genes extracted from *de novo* genome assemblies generated from Illumina whole-genome sequencing (WGS) reads.

---

# Overview

The goal of this pipeline is to process Illumina WGS reads into a time-calibrated phylogeny suitable for downstream biogeographical analyses.

The workflow includes:

- Raw read quality control
- Read trimming
- De novo genome assembly
- BUSCO gene extraction
- Alignment cleaning and trimming
- Concatenated and gene-tree phylogenetic inference
- Species tree estimation
- Time calibration
- Historical biogeographical analysis

The pipeline was run on a Unix-based system using a SLURM job scheduler for computationally intensive steps.

---

# Repository Structure

```text
├── pipeline/
├── data/
│   └── Raw_reads/
└── README.md
```

---

# Pipeline

The `pipeline/` directory contains the scripts used throughout the phylogenomic and biogeographical workflow.

- `1_FastQC.sh` -> A simple quality check of raw reads using FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc)
- `2_QCwithFastP.sh` -> Another quality checker tool, FastP (Chen et al., 2018), that allows trimming of read ends
- `3_SPAdes_quast_loop.sh` -> de novo genome assembly with SPAdes (Bankevich et al., 2012), followed by quality check with quast (Gurevich et al., 2013)
- `4_Rename_SPAdes_scaffolds_for_buscophy.sh` -> Renames assemblies for next step
- `5_Buscophy.sh` -> Extracts complete and fragmented Lepidopteran* BUSCO genes from the genome assembly using buscophy (https://gitlab.leibniz-lib.de/smartin/buscophy)
- `6_OliInSeq_aa.sh` -> Removes outliers from amino acid alignments with OliInSeq (https://github.com/cmayer/OliInSeq)
- `7_OliInSeq_nt.sh` -> Same as above but for nucleotide alignments
- `8_trimAl_aa.sh` -> trimAl (Capella-Gutiérrez et al., 2009) to remove columns with less than 50% data in amino acid alignments
- `9_trimAl_nt.sh` -> As above but for nucleotide alignments
- `10_AMAS.sh` -> Creates summary statistics for each alignment with AMAS (Borowiec, 2016), and concatenates the alignments together
- `11_IQ-Tree_MF.sh` -> IQ-TREE (Minh et al., 2020) model finder on concatenated alignment
- `12_IQ-Tree_TS.sh` -> IQ-TREE tree search on concatenated alignment
- `13_IQ-Tree_GT.sh` -> IQ-TREE tree searches on individual nucleotide alignments
- `14_ASTRAL.sh` -> Creates consensus tree with ASTRAL-IV (Zhang et al., 2025) from gene trees
- `15_Concordance.sh` -> Tests for concordance between concatenated tree and consensus tree
- `16_Concordance_vector.R` -> As above
- `17_Change_labels.R` -> As above
- `18_Concordance_table.R` -> As above
- `19_LSD2.sh` -> Time calibrates the phylogeny with LSD2 (https://github.com/tothuhien/lsd2) based on user-specified calibration points
- `20_Constrained_bs_for_treePL.sh` -> Creates bootstrap replicates for running treePL (Smith et al., 2012)
- `21_Root_bootstrap_replicates.R` -> Roots bootstrap trees from above
- `22_Stage1_treePL.sh` -> treePL for time calibration
- `23_Stage2_treePL.sh` -> treePL for time calibration
- `24_Stage3_treePL.sh` -> treePL for time calibration
- `25_BioGeoBears.R` -> Ancestral range reconstruction analysis with BioGeoBears (Matzke, 2013)
- `26_BioGeoBears_BSM.R` -> Biogeographical stochastic mapping with BioGeoBears

*Can be specified to your taxonomic group of choice (see https://busco.ezlab.org/list_of_lineages.html)

---

# Data

The `data/` directory contains template configuration files and sample lists required to run the workflow.

- calibrations.txt -> Required for `19_LSD2.sh`
- config_STAGE1.txt -> config file required for `22_Stage1_treePL.sh`
- config_STAGE2.txt -> config file required for `23_Stage2_treePL.sh`
- config_STAGE3.txt -> config file required for `24_Stage3_treePL.sh`
- dispersal_multipliers_fn -> Required for `25_BioGeoBears.R` and `26_BioGeoBears_BSM.R`
- outgroup.txt -> Required for `19_LSD2.sh`
- Parasa_geog.data -> Required for `25_BioGeoBears.R` and `26_BioGeoBears_BSM.R`
- sample_listFastP.txt -> Sample list required for `2_QCwithFastP.sh`
- sample_listFastQC.txt -> Sample list required for `1_FastQC.sh`
- sample_listGenetreealignments.txt -> Sample list required for `13_IQ-Tree_GT.sh`
- sample_listOliInSeq_aa.txt -> Sample list required for `6_OliInSeq_aa.sh`
- sample_listOliInSeq_nt.txt -> Sample list required for `7_OliInSeq_nt.sh`
- sample_listSPAdes_scaffolds.txt -> Sample list required for `4_Rename_SPAdes_scaffolds_for_buscophy.sh`
- sample_listSPAdes.txt -> Sample list required for `3_SPAdes_quast_loop.sh`
- sample_listtrimAl_aa.txt -> Sample list required for `8_trimAl_aa.sh`
- sample_listtrimAl_nt.txt -> Sample list required for `9_trimAl_nt.sh`
- times.txt -> Required for `25_BioGeoBears.R` and `26_BioGeoBears_BSM.R`

---

# Raw Reads

The `Raw_reads/` directory contains paired-end Illumina sequencing reads for each sample.

```text
data/
└── Raw_reads/
    ├── T1_013/
    ├── T1_018/
    └── ...
```

These files were generated by Novogene.

Due to file size limitations, raw genomic data are hosted externally with NCBI Sequence Read Archive under the accession number PRJNA1462647.

---

# Software Requirements

The following software and tools were used throughout the workflow, with specific versions found in the script files and manuscript:

- FastQC
- fastp
- SPAdes
- QUAST
- buscophy
- OliInSeq
- trimAl
- Geneious
- AMAS
- IQ-TREE 2
- ASTRAL-IV
- LSD2
- treePL
- BioGeoBEARS
- R

---

# System Requirements

- Unix/Linux operating system
- SLURM job scheduler recommended for parallel analyses
- Sufficient RAM and storage for genome assembly and phylogenomic analyses

---

# References

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., ... Pyshkin, A. V. (2012). SPAdes: A new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology*, 19, 455–477.

Borowiec, M. L. (2016). AMAS: A fast tool for alignment manipulation and computing of summary statistics. *PeerJ*, 4, e1660.

Capella-Gutiérrez, S., Silla-Martínez, J. M., & Gabaldón, T. (2009). trimAl: A tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*, 25, 1972–1973.

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34, i884–i890.

Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. *Bioinformatics*, 29, 1072–1075.

Matzke, N. J. (2013). BioGeoBEARS: BioGeography with Bayesian (and Likelihood) Evolutionary Analysis in R Scripts [Computer software]. *University of California, Berkeley*.  
https://CRAN.R-project.org/package=BioGeoBEARS

Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution*, 37, 1530–1534.

Smith, S. A., & O’Meara, B. C. (2012). treePL: Divergence time estimation using penalized likelihood for large phylogenies. *Bioinformatics*, 28, 2689–2690.

Zhang, C., Nielsen, R., & Mirarab, S. (2025). ASTER: A package for large-scale phylogenomic reconstructions. *Molecular Biology and Evolution*, 42, msaf172.
