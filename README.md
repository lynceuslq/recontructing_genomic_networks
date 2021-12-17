# recontructing_genomicnetworks 
To study the evolutionary history of bacteriophages (particularly those assembled from metagenomes) by reconstructing networks


An example of the gene-sharing network of 21 viral clusters from GPD, visualized with Cytoscape (http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/)
![vcset3_vcontact2 vc gpd c1 ntw](https://user-images.githubusercontent.com/55744039/125021832-4ea18080-e073-11eb-98ad-326c6f987e77.png)

This is a prototype method to detect and estimate relative abundance for metagenome-assembled phage genomes classified with de nov clustering methods (eg. vConTACT2), here, GPD (http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/) was used as the database for the pipelines.

# Manual




# Literature outlines

Outlines: the identification of core proteins for genus-level clusters of phages and a workflow of abundance estimation of viruses with core proteins

Summary
A new approach is proposed in this research to automatically identify core proteins for subfamily- or genus-level clusters of gut phages, by which relative abundance of phage clusters can be profiled within highly reduced computational time, compared to genome alignments.
With profiling marker proteins, impacts due to homologous bacterial genomes with phages on the accuracy of abundance estimation were solved.
 The new method can help us investigate the interaction of phages and bacteria in guts, and consequently, the association between gut phages and health status in an efficient and systematic manner. 
The research may also bring ideas to unravel the evolution of phage, an unsolved mystery complicated by their interaction with hosts and massive gene interchanges. 

Background
1. Current studies concerning novel phage genomes obtained from meta genomes (human guts and other systems): 
Limitations in phage diversity in Refseq and other resources: metagenome-assembled phages expand our knowledge in gut phageomes greatly
Limited information to assign metagenome-assembled phages to known lineages (eg. ICTV liaisons)

2. The importance of quantifying ecosystem in guts
Disease-association studies: gut microbiome is known to contribute to multiple diseases and disorders (IBD, tumour development etc.)
Providing information for the interaction and evolution of gut microbiota 
Quantification tools have been well-developed to estimate and compare bacterial compositions in meta genomes, which boost the development of metagenomic studies for the last decade, yet, there is no standard tool to compare phage abundance, therefore, their contribution in human health remains unclear.

3. Challenges in the estimation of phage abundance in meta genomes
Current commonly used methods: aligning reads against phage genomes (genomes from phage databases for assembled from the same meta genome)
Problems: time-consuming & bacterial contamination
Current methods to estimate microbial abundance in meta genomes 


4. Aims of the study 
Identifying core proteins for de novo phage clusters and using them to estimate relative abundance of gut phages
Identifying protein clusters and viral clusters 
Selection on core proteins -> proteins for abundance estimation
Profiling PE reads against selected proteins
Assessing accuracy 

Methodology

1. Obtaining phage genome and proteome databases and gut bacteria databases
Description of GPD (num. genomes, clusters, resource diversity)
Genome clusters included in the study and why
Using MGnify as gut bacterial resource

2. Reconstructing viral clusters and identifying core proteins
Vcontact2 workflow and adjustments in its parameters
Identifying phage clusters at subfamily or genus level
Defining core proteins for subfamily- or genus-level phages
Functional conservation of phage proteomes (?) 

3. Selecting protein clusters for abundance estimation
Estimation of intra-cluster protein homology
Criteria of primary selection
Why and how to perform cleaving on protein alignment for secondary selection
Removing inter-homologs PCs (id > 75 && alignment length > 30 aa)
and limiting the number of PCs for each VC (?)

4. Estimating abundance of phage clusters
Making interleaved reads (if PE reads)
Generating alignments to selected PCs with Diamond
Calculation coverage rate and coverage depth with bedtools
Filtering proteins for abundance estimation
Calculating abundance estimates and SD for each cluster with negative log likelihood (or mean values)

5. Designing artificial metagenomes
Artificial metagnomes simulated with insilicoseq to resemble human gut samples
Scenarios: abundance distribution models, phage species diversity, host compositions and diversity, etc
For most cases, 3 replications were generated 


6. Assessing abundance estimates by marker protein clusters and genomes
Parameters to compare (Pearson’s regression with true phage abundance) 
Abundance estimation with phage genomes
Reduction in computational time


Results

1. Overview of VCs and core-PCs identified by the pipelines, and their relation to the original literature
Total amount of genomes in GPD140000 -> 85471 genomes in clusters with >= 3 genomes
85471 input -> 68251 phages clustered with vcontact2 (3533 VCs)
62149 phages from 2918 viral clusters (42472 PCs from primary selection with Jaccard index >= 0.7)
43571 phages (2319 VC) and 19241 PCs were selected based on primary selection standards
7321 PCs were added by secondary selection standards
Removing inter-homo PCs (886: protein-pairwise homology, 1904: protein-genome homology)
2825 VCs (56509 genomes) 

2. Choosing suitable parameters for abundant estimation with marker proteins
Comparing the outcomes of TP and FP alignments in diamond blast table
Selecting different sets of alignment identity and e-value thresholds and comparing rates of TP and PF after filtering
Decision: two sets of parameters: “more sensitive” and “more specific”
 
3. Assessing the accuracy of phage abundance estimated by marker proteins 
For phage-only meta genomes (across all abundance models and diversity levels_, the accuracy (resembled by Pearson’s regression) of estimates by marker proteins is comparable to the accuracy of genome estimates (coefficient ~ 60 to 80%)
For meta genomes integrated with bacterial hosts, the accuracy of marker protein estimates is independent of bacterial compositions, whereas the accuracy of genome estimates drops in the presence of hosts
The accuracy of genome estimates significantly goes down with the fraction of host genomes in metagenomes, the accuracy of protein estimates exceeds genome estimates when total host abundance is more than 65%. 

4. Assessing the enhancement in computational performance by estimating abundance with marker proteins
Comparing computational time between aligning reads to marker proteins and genomes (differing reads counts)


Discussion
1. A conclusion 

2. Limitations in the method 
Limited knowledge in the classification metagenome-assembled phages and generate more sensible viral clusters
Unable to assign more phages to genome clusters due to limitations in computing capacity -> need better tools to perform de novo clustering on viral genomes
The lack of representative proteins when using other de novo genome cluster methods (GPD genome clusters), which requires less computational resources but worse clustering results. 

3. Uncertainties in quantifying the composition of microbiota 
Metagenome-assembled phages: incompleteness in genomes, mistakenly assembled, host assignments, etc
Completeness in database: a definite diversity of gut phages, or more will be observed if we sample more? 

4. Future perspectives
Integrated phage databases
Improvements in clustering methods to allow assigning phages with higher divergence to genome clusters ( clustering rates) 
Improvements in the standards of quantifying phage evolution (whole-genome diversity or gene-sharing networks)


Reference

