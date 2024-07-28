# Resampling Phased Panels

Resampling refers to generating mosaic panels by sampling phased genomes in a panel. 

Resampling panels is useful for increasing the haplotype diversity by generating mosaic panels, which enables query and reference panel to protect chromosome long haplotype sequences at the individual level. If we shared the chromosome long haplotypes, even very small leakage of each variant will enable easy re-identification of individuals. Note that even with re-sampling, even some simple statistics (and more complex attacks) can reveal participation of an individual in the resampled panel. That is why it is important to follow resampling with other mechanisms, e.g, permute and allele hashing, augmentation, to protect variant coordinates and identities.

Using ProxyTyper, the query and reference panels can be resampled to higher sample sizes. While this procedure may create artificial haplotypes, it increases the haplotype diversity.

## How does ProxyTyper resample haplotypes?
ProxyTyper uses a Markov chain-based approach to randomly select a haplotype at "anchor" sites. The anchor sites refer to variants where we can perform a haplotype switch. For large panels with millions of variants, choosing anchor sites has very low impact on the quality of resampled panels and decreases the computational burden on ProxyTyper. By default, ProxyTyper selects variants that are at least 0.005 cM away from each other as anchor sites. Even with this small distance, the computation becomes much more computationally efficient. 

## How is Resampling related between Query and Reference Panels?
They are not related at all. Resampling is done at each site locally and independent of other sites. Sites do not have to expect any input from other site to perform resampling. This is why resampling is almost always one of the first steps in the protocols. 

When query panel is phased, it is useful to increase its size by resampling it and protect chromosome-long haplotypes. Of note, this means when it is sent to imputation server, the imputed panel that is returned from the server cannot be immediately mapped to query site's original subjects. As we will see later, ProxyTyper protocols include a local re-imputation step using the imputed panel returned from the resampled query panel.

## How can I run Resampling?
After importing a VCF file (see previous example), we are ready to resample it. In ProxyTyper, resampling is performed with *-resample_panel* option:
```
PANEL_ID=vcf_import_panel
./ProxyTyper.sh -resample_panel ${PANEL_ID} 5000 ${PANEL_ID}_resampled
```
In this command, we start from a panel whose identifier is "vcf_import_panel" and we saved the resampled panel as *${PANEL_ID}_resampled*. You can easily check the number of subjects and variants in this resampled panel:
```
wc -l ${PANEL_ID}_resampled_variants.bed
wc -l ${PANEL_ID}_resampled_subject.list
```

After resampling, one simple diagnostic is to calculate the alternate allele frequencies:
```
./ProxyTyper.sh -calculate_panel_AAF ${PANEL_ID} ${PANEL_ID}_AAF.bed
```
This command writes the alternate allele frequencies for each variant to a BED formatted file. The 5th column of the file contains the allele frequency in the resampled panel. This can be compared with the original panel's frequencies (these should already be in the variant identifier columns, since we asked for this while importing the VCF file.)

## Parameters that Govern Resampling
Following parameters in ProxyTyper's configuration file (PROXYTYPER.ini) govern the resampling procedure. It is generally not necessary to modify these.
```
haplocoded_flag=1                       : Save resampled panels as phased. 
pre_mosaic_Ne=0.125                     : This tunes the strength of recombination. Higher N_e means more recombinations in resampling. Between two variants separated by genetic distance of d centiMorgans, prob_recomb=(1-exp(-4*pre_mosaic_Ne*d))/N_haps
pre_mosaic_allele_eps=0.0               : Allele probability, this is not used in any protocol. If you want uniform IID errors to be introduced, set this.
pre_mosaic_n_threads=${n_threads}       : # of threads to be used for resampling.
min_cM_delta_per_anchoring_vars=0.001   : Mininum genetic distance between the anchors used for resampling. We do not change this in any experiment.
save_recomb_patterns=0                  : This option is not used, keep it set to 0.
max_l_seg_n_bps=10000000                : Maximum length of any segment in base pairs. This is to ensure that ProxyTyper does not randomly sample very long segments by chance.
max_l_seg_cM=0                          : Maximum resampled segment length in cMs. 0 means ignore this cutoff (note we have a cutoff in bps)
max_l_seg_nvars=0                       : Maximum resampled segment length in number of variants. 0 means ignore this cutoff (note we have a cutoff in bps)
```

---

We include an example of re-sampling with several options to explore the resampled panel in a script named *Resampling_Panels.sh*. Make sure *ProxyTyper.sh* and *PROXYTYPER.ini* are copied to this folder, and the *ProxyTyper_Release* is in PATH (See installation/ for more information) before running *Resampling_Panels.sh*.