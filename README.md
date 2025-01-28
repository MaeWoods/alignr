# alignr 0.1.0: An R package that produces a full length human T cell receptor by calling fragments and genes from a 10X genomics single cell annotation JSON file.

The purpose of this micro tool is to compartmentalize single cell RNA sequencing steps used to manufacture transgenic receptors.

# Rationale for initial package development

The sofware was initially developed to identify receptor sequences from single cell RNA sequencing coupled with amplification of the receptor sequence. Specifically Version 0.1.0 has been used to quantify cells that are activated through the binding of their T cell receptors (TCRs) to cognate antigens. Therefore, the tool provides a method to identify full length TCRs with specificity directed toward a target antigen, that can then be tested for use as transgenes. However the aim is to develop functionality for the analysis of a broad range of sequencing of polyclonal data.

**Introduction to TCR sequencing**

T cells contribute to our immune system by recognizing targets on infected cells via the expression of the TCR receptor. The TCR sequence is important to study because T cell specificity depends on the **CDR3** region of the TCR and this varys in sequence between different people and single cells, equipping T cells with capacity to mount a response against a broad range of infected targets. T cells can be grouped into clonotypes that share a common CDR3 beta chain and this way, used to estimate the frequency of target specific T cells in the blood. Harnessing this heterogeneity in sequence between T cells for a quantitative analysis of adaptive immunology has broad applicability in immunology and immunotherapy because the clearing of infection and cancer depends on availability of immune cells (including T cells) with capacity to mount a response. 
Immunosequencing is a PCR-based based method that exploits the capacity of high-throughput sequencing technology to characterize tens of thousands of TCR CDR3 chains simultaneously and **aocseq** has been developed to analyse and annotate this data.

Specifically, this is a software package of statistical tools that can be used to trace, analyse, annotate and query clonotypes subject to amplification or reduction in frequency following antigen stimulation or between experimental conditions. **aocseq** has initially been applied to Virus specific T cells (VSTs) because these clinical blood products contain non specific bystander T cells in addition to potent virus specific clonotypes. However, the tool can be adapted to model other barcoded time series frequency data and the accompanying vignette demonstrates how the tool could be used to track clonotypes *in vivo*, using Adaptive's ImmuneAccess database. 

# Installation and running the software: 

install_github("MaeWoods/aocseq");

library("alignr")

See the **vignettes** for further help loading and running the software.

**Vignettes**

* To read in gene expression data and analyze clonotypes:
[Getting started with aocseq](./html/vignette.svg)
