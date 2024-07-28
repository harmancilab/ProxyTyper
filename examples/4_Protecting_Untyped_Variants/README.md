# Protecting Untyped Variants by Decompositions

This folder contains examples of the mechanisms for protecting the untyped variants using decomposition of the untyped variants.

In the previous examples, we focused on protecting the typed variants. Here we focus on the untyped variant protection. In protocols, these mechanisms should be combined to protect all variants.

The untyped variants reside on the reference panel (not the query or the imputation server). As such, the reference site(s) must develop policies about sharing them with other sites. Given that most untyped variants are ultra-rare (most likely singletons or doubletons, etc) even the knowledge of their coordinates may pose privacy risks.

## How does ProxyTyper Protect Untyped Variants?
To protect the untyped variants, ProxyTyper adopts a "hiding in the crowd" approach: For each variant above a certain alternate allele frequnecy, ProxyTyper flips a coin to select the variant (or skip it) and, if selected, it partitions the haplotypes that harbor the alternate allele into 2 distinct sets. The partitioned haplotypes are used to set the genotypes of two new proxy variants. This way, ProxyTyper expands the number of untyped variants. The proxy variants have, by definition, lower alternate allele frequency. This expands the allele frequency spectrum of the untyped variants and hides original untyped variants. 

The reference site can run decomposition as many times as they do so desire. However, each decomposition will create double the number of untyped variants above the selected allele frequency threshold and may impose extensive computational burden on the imputation server. ProxyTyper performs one round of decomposition by default.

## Only Reference Site Performs the Untyped Variant Decomposition (Not the imputation server)
Given a reference panel, following command performs the untyped variant decomposition and writes decomposition file at reference site:
```
REFERENCE_panel_id=REFERENCE_PANEL
./ProxyTyper.sh -decompose_variants TYPED_VARIANTS.bed ${REFERENCE_panel_id} ${REFERENCE_panel_id}_decomposed
```
When this command completes, we get a new panel named REFERENCE_PANEL_decomposed and an interval file that describes the decomposition procedure. This interval file is intentionally not copied to the model directory (i.e., *PROXIZATION_MODEL*) because it is a private model that can only be accessed by the reference site. This file should not be shared with the query site or the imputation server.

Counting the number of variants in the decomposed panel indicates the number of proxy untyped variants added to the panel:
```
wc -l ${REFERENCE_panel_id}_variants.bed
wc -l ${REFERENCE_panel_id}_decomposed_variants.bed
```

## Parameters that Govern Untyped Variant Decomposition
Following parameters in the configuration file (*PROXYTYPER.ini*) govern the untyped variant decomposition:
```
variant_decomp_min_AAF=0.1  : Minimum alternate allele frequency for the untyped variants to be decomposed. All untyped variants with alternate allele frequency above this value will be decomposed. Decreasing to a low value will create a very large number of proxy variants and increase computational burden (a lot). 
shuffle_decomp_vars=1       : This flag indicates that we will shuffle the decomposed variants to a random location within the typed-typed variant block that they reside in. Don't change this.
USE_MT_RECOMPOSITION=1      : This flag indicates ProxyTyper to use multithreading in decomposition. Don't modify it.
```

## What happens to the expanded set of untyped variants after they are decomposed?
The reference site sends the decomposed panel to the imputation server (after all mechanisms applied, with the typed variants). Imputation server uses this as the reference panel in BEAGLE and imputes the query panel (this panel should also be "proxized" by the same mechanisms). The imputed results are sent back to the query site. This file is a VCF file with the allele probabilities for each of the untyped variant.

The imputed untyped variants are sent to the query site. *At this point, the reference site must choose a policy to reveal which untyped variants will be revealed to the query site.* Based on this policy, the reference site can choose to filter the interval file to only share the decomposition patterns of certain variants. 

Given the BEAGLE imputed VCF file, we use following command to convert the VCF file into a "recomposed" panel with the original untyped variants from the reference panel:
```
BEAGLE_VCF_FILE=imputed_query_proxy_panel.vcf.gz
QUERY_PANEL_ID=query_panel
UNTYPED_DECOMPOSING_INTERVAL_FILE=untyped_decomposing_patterns.interval
./ProxyTyper.sh -recompose_BEAGLE_VCF_variants ${BEAGLE_VCF_FILE} \
${QUERY_PANEL_ID}_subjects.list \
${UNTYPED_DECOMPOSING_INTERVAL_FILE} \
QUERY_RESAMP_IMPUTED
```
After completion, this command processes the allele frequencies and recomposes them using the decomposition interval file to generate the final allele frequencies of the original untyped variants in the reference panel. (Note that we did not explore this step yet since it pertains to running BEAGLE.)

---

The demonstration of decomposition is implemented in the bash script named *Decompose_Untyped_Variants.sh* in this folder.