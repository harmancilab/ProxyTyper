# Basics and Setup 

We lay out general description and setup of the entities involved in imputation outsourcing. We also introduce the basic idea behind using proxy-panels in genotype imputation. Following directories provide extensive examples of using mechanisms to build panel processing pipelines for privacy-aware genotype imputation.

# Basic Setup of Genotype Imputation Outsourcing

Genotype Imputation is an important step in many genomics data analysis pipelines, e.g., GWAS. In a nutshell, a cohort of individuals are genotyped using, for example, an array platform and a sparse genotype matrix is generated (as a VCF file). To impute the other variants, researchers who hold the sparsely "typed"* variant genotypes can use imputation software to "fill-in" the genotypes for remaining set of variants, which we refer to as untyped variants. We refer to the researchers who hold the typed variant genotype panel as the *querying site* to facilitate the outsourcing scenario below.

Imputation may involve usage of large panels that may not always be accessible to research institutes. In addition, imputation is a computationally intensive procedure. To accommodate these considerations, imputation process may be outsourced to a *reference site* (e.g., Michigan Imputation Server), who are geared for executing the imputation software. The reference site In this setup, the query site uploads their sparse panel to the reference site. Reference site uses a large local imputation panel (e.g., TOPMed, 1000 Genomes, Haplotype Reference Consortium) This, however, requires two factors:

1) Query site must ensure that the typed variant genotype panel can be uploaded to the imputation server. This may not be always feasible.
2) The reference panel that is used for imputation at the reference site is often not publicly available and is owned by entities that are separate from the imputation server. Thus, the owners of the reference panel may lose control over the reference panel dataset when it is used in plain form at the reference site. Thus, we may have a third *hidden* entity, who is the actual owner of the reference panel. This generally complicates the concerns about privacy when reference panel is used and managed by the imputation service (what if the imputation servers are hit by ransomware attacks?).  

*ProxyTyper*'s main goal is to make sure that these concerns are minimized while ensuring that imputation can still be practically performed. The basic idea is to hash, encode, and anonymize reference and query panels to generate *proxy-panels* such that they do not leak any information about the original variants and the genotypes. Each time ProxyTyper generates a proxy-panel, it uses randomness in the operations (e.g., hashing, resampling, random permutes, and noise addition). This is essential to provide the privacy since proxy-panels should not be easily inverted to original panels. Thus, the proxy panels generated from same panels will not be similar in most properties, e.g., allele frequencies, number of variants, variant positions, and genetic maps, when they are proxized independently.

# Panel Protection Mechanisms Basics

ProxyTyper makes use of numerous panel processing operations (i.e., mechanisms) to build protocols for genotype imputation:

<ol>
<li> **Panel Resampling**: Resampling is used for shuffling the haplotypes on phased panels to create "mosaic" of the original panel. *ProxyTyper* can resample panels to any desired number of subjects. The resampling can slightly change allele frequency distribution of the variants due to the random nature of sampling of the haplotypes. Panel resampling does not rely on the whether a variant is typed or untyped. 
<li> **Typed variant genotype augmentation**: Augmentation aims to increase the number of typed variants and randomly place them in vicinity of original typed variants.
<li> **Typed variant genotype hashing** Given a phased genotype panel, the alleles on each haplotype is hashed (or sometimes we call this "proxized") by binary combination of the surrounding variants. The variants that are used for hashing and their effect on hashing are randomly selected by ProxyTyper and stored in a *proxization model*.
<li> **Typed variant permutation** 
<li> **Untyped variant decomposition** Variants that are untyped constitute variants whose genotypes are imputed by the reference site. These variants are dominated by low frequency variants, especially in large panels. The untyped variant decomposition mechanism decomposes the high frequency untyped variants into new *proxy-variants*, which aims at hiding  
<li> **Coordinate and Genetic Map Anonymization** Most attacks (Homer t-test, Bustamante beacon, Srirankaraman LRT-test) on genotype panels require matching between a reference panel and the pool panel. Thus, anonymizing the coordinates and genetic maps ensures that the coordinates of variants cannot be immediately mapped to the . Of note, BEAGLE does not make use of varaint coordinates if the genetic map is provided.
</ol>
* *Throughout this documentation we may inadvertantly refer to typed variants as "tags" and untyped variants as "targets", we hope this does not create confusion.*

# Model-based Mechanisms Require Initializing Random Proxization Models

Among these mechanisms listed above; augmentation, hashing, permutation, anonymization, and decomposition mechanisms make use of random *proxization models* that are built by ProxyTyper. These models are shared between the query site and the data owners, who encode their plaintext panels and generate the proxized panels. If the models do not match, the proxized panels are not *compatible* and cannot be used for imputation. Also, the knowledge of proxization models can enable any other entity to decode the proxy-panels. Thus, the proxization models should be kept confidential, similar to symmetric encryption keys in encryption algorithms.

# ProxyTyper is a Genotype Panel Processing Pipeline

Each mechanism processes and modifies a certain property of the input panel and *maps* them to a new domain. The new panel can then be input to other mechanisms. As such, we can feed the output from these mechanisms into each other to build new genotype proxization protocols. 

As an example let's assume that the query site has a panel that contains 20,000 untyped variants and 10,000 subjects. We can use the following pipeline to generate a randomly proxized panel:
1. Resample the panel up-to 40,000 subjects:
 [20,000 variantx10,000 subjects]=> Resample => [20,000 variantx40,000 subjects]
2. Augment the variants with 0.5 probability at each typed variant:
a) Build augmentation model with 0.5 probability.
b) [20,000x40,000] => Augment(0.5) => [28,000x40,000]
Here, we assume the augmentation added 8,000 new typed variants to the panel.
3. Hash-proxize the variant genotypes using 11-variant vicinity of each variant.
a) Build a hash-proxization model using 11-vicinity of each variant.
b) [28,000x40,000] => Hash-Proxize(11 vicinity) => [28,000x40,000]
This panel should have very different allele-frequency characteristics than the original panel. The local haplotype-frequency of the new panel is also substantially different from the previous panel.
4. Permute-proxize the typed variants in [28,000x40,000] panel.
5. Anonymize coordinates and genetic maps.

The last two steps permute the typed variants locally (keeping genotypes same) and anonymizes the coordinates. The final panel can be used for imputation.
Obviously, data owner site needs to process the same typed variants so that they are "compatible" and can be used for imputation. This means that the proxization models need to be shared between query site and the data owners. In addition, the data owners can proxize the untyped variants using the decomposition mechanism. ProxyTyper stores the models as bed files and binary files and can be compressed to a small manageable size to transfer between data owner and query sites.

## Protocols are Panel Mapping/Mutating Pipelines 

The main point of this introductory text is to motivate that ProxyTyper's mechanisms can be viewed as a genotype panel processing (mapping/mutating/augmenting/decomposing) pipeline to map original panels to proxy panels. The implementation of a protocol using ProxyTyper is thus a flow of calls to the corresponding mechanisms that map the panels in a sequential order. For the model-based mechanisms, ProxyTyper ensures that the model is compatible with the input panel and gives error messages if they cannot be used together. 

# Multithreading

Typed variant augmentation and variant decomposition ProxyTyper makes use of multithreading extensively. This is very important to decrease the run time while processing large panels.

# Executables 

ProxyTyper is composed executables. The main script named *ProxyTyper.sh* is the front-end and implements each mechanism and the executable file *ProxyTyper_Release* implements the underlying functionalities that are used for generating the proxized panels. *ProxyTyper.sh* contains extensive error checks and input validation steps and should be the starting point for any user who wants to use existing protocols or design new imputation protocols using the protection mechanisms offered by ProxyTyper. It should be seamless to integrate *ProxyTyper.sh* into new pipelines.

*ProxyTyper.ini* contains numerous options that are used for setting up the parameters needed by the different mechanisms. We will extensively review these scripts and implemented options in later examples. After building ProxyTyper, it is necessary that the executable under *bin/ProxyTyper* is included in the *$PATH*.

---
