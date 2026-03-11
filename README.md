During the course of my project several groups explored the use of scRNA-seq to analyse gene
expression within the MGE of embryonic mice, providing new insights into the heterogeneity
within this region and to the processes controlling neurogenesis (Mi et al. 2018) (Lee et al.
2022). During development, the MGE is patterned into ventral and dorsal domains through the
actions of secreted growth factors, such as SHH and WNT (Mi et al. 2018). The formation and
proportion of interneuron subtypes rely on the establishment of specific domains within the
MGE, with the ventral region primarily producing PV interneurons and the dorsal region
primarily producing SST interneurons. Advancements in next-generation sequencing provide
new avenues for insight and understanding of MGE patterning and to the molecular basis of
interneuron diversity. In particular, scRNA-seq allows for unsupervised genome-wide analysis of
gene expression at a cellular level (Chen, G, Ning & Shi 2019; Tang, F et al. 2009).
The workflow for scRNA-seq requires biological samples to be processed to produce single cells
which are each provided with unique barcodes (Picelli 2017). cDNA is created from the cells
using reverse transcription on single-stranded RNA, which is then amplified (Haque et al. 2017).
Sequencing is performed on the transcripts, generating raw data that is analysed through
bioinformatics to assign the transcriptome for individual cells (Haque et al. 2017). The
bioinformatic workflow consists of well-established processes, including quality control,
normalisation, and dimensional reduction using techniques such as uniform manifold
approximation and projection (UMAP) or t-distributed stochastic neighbor embedding (tSNE)
(Satija et al. 2015; Stuart et al. 2019). Furthermore, the data must be deconvoluted into clusters 
of similarity to assess gene expression within and across cell types, and across the entire dataset.
As many clustering methods presume equal quality of cells analysed, they do not consider
biological and technical factors that can alter cell classification (Andrews, TS & Hemberg 2018).
Exceptional consideration must be taken when identifying clusters for accurate scRNA-seq
analysis. While scRNAseq provides excellent data for hypothesis generation, both in vitro and in
vivo analyses are crucial to validate the findings (Su-Feher et al. 2022).
The earliest scRNA-seq datasets used the entire MGE in attempt to deconvolute the complexity
of neural precursors, finding both proliferating neural precursors and immature neurons which
both contained sub-populations. Although lineage of these subpopulations is connected to
GABAergic interneurons, they were unable to connect progenitor subpopulations to mature
interneurons (Chen, YJ et al. 2017). Using physical microdissection of ventral and dorsal
domains of the MGE, which preferentially give rise to PV and SST interneurons, Mi et al (218)
were able to identify distinct transcriptional profiles within ventral and dorsal domains that they
propose represent subtypes of PV and SST interneuron precursors (Mi et al. 2018). Excitingly,
they identified 14-3-3ζ as one of the differentially expressed genes within one of the PV
subtypes, which, if real, could help to explain the mild deficiency of PV interneurons seen in 14-
3-3ζ-/- mice. The most recent sc-RNAseq data set from GFP lineage traced Nestin-dVenus mice
found significant transcriptional heterogeneity between lateral, medial and caudal ganglionic
eminences and within individual spatial subdomains of the ganglionic eminences (Lee et al.
2022). Furthermore, they also found differential gene expression between VZ cells at E12.5 and
E14.5. However, in contrast to Mi et al (2018) they were unable to find any interneuron
precursor subtypes. Within this chapter, I utilised existing public scRNA-seq datasets in 
combination with experimental methods to explore a potential role for 14-3-3ζ in specific subdomains within the MGE during neurogenesis.
