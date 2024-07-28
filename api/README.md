# ProxyTyper.sh API

We list the options that are used by *ProxyTyper.sh*. 

---

#### ./ProxyTyper.sh -import_VCF [VCF file path] [Is VCF phased? (0/1)] [Add Ref/Alt allele-AF Info to name (0/1)] [Output prefix]
Import a VCF file and save it in a panel. VCF file is expected to be compressed by gzip. Phased VCF files should set the corresponding flag to 1. This option can also add the ref/alt variant and allele frequency information to the variant identifiers, if requested. The output panel is saved to the prefix.

---

#### ./ProxyTyper.sh -setup_BEAGLE_genetic_maps [4-column PLINK formatted (chr rsid cM posn) genetic map file URL]
Downloads and sets up BEAGLE genetic maps in the default genetic maps directory defined in the configuration file *PROXTYPER.ini*. By default, this directory is named "*genetic_maps*".

You can use the default URL, https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip for hg19.

---

#### ./ProxyTyper.sh -resample_panel [Panel prefix] [Resample size] [Output prefix]"
Resamples an imported panel upto the resample size and saves it into a panel using the output prefix.

---

#### ./ProxyTyper.sh -calculate_panel_AAF [Panel prefix] [Alternate allele frequency BED file]"
Calculates the alternate allele frequency of variants in a panel and saves it to a BED file.

---

#### ./ProxyTyper.sh -uniquefy_panel_variants [Panel prefix] [Output prefix]"
Removes duplicated variants from a panel and saves the unique panel to the output prefix.

---

#### ./ProxyTyper.sh -subset_panel_variants [Panel prefix] [Variants BED file] [Output prefix]"
Selects subset of variants in the BED file from a panel and saves to output prefix.

---

#### ./ProxyTyper.sh -subset_panel_subjects [Panel prefix] [Subject list file] [Output prefix]"
Selects subset of subjects in the list (one subject identifier in each row) in a panel and saves to output prefix.

---

#### ./ProxyTyper.sh -generate_tag_proxizer_model [Tag variants BED file]
Generates tag hashing model using the tag variants in the BED file. Model is saved under ${MODELS_DIR}, defined in *PROXYTYPER.ini* file.

---

#### ./ProxyTyper.sh -random_phase_panel [Panel prefix] [Output prefix]
Randomly phases an unphased panel and saves it to the output prefix. 

---

#### ./ProxyTyper.sh -proxize_tag_variants [Panel prefix] [Output prefix]
Calculates hash for the phased panel using the hashing model (*-generate_tag_proxizer_model* option) and saves it to the output prefix.

---

#### ./ProxyTyper.sh -generate_tag_permuter_model [Tag variants BED]
Generates a tag permutation model for the typed variants in the BED File.

---

#### ./ProxyTyper.sh -permute_proxize_tag_variants [Panel prefix] [Output prefix]
Permutes the tag variants in the panel and saves the panel to output prefix. The permutation model must be generated before this option can be called

---

#### ./ProxyTyper.sh -generate_coordinate_anonymizer_model [tag variants BED (query site variants)] [tag+target variants BED (reference site variants)]
Generates the coordinates and genetic map anonymizing model for the typed and typed+target variant sets. Anonymized genetic map for the typed variants is saved and the anonymized coordiantes are saved under ${MODELS_DIR}.

---

#### ./ProxyTyper.sh -anonymize_coordinates [Panel prefix] [Output prefix]
Anonymizes the coordinates of all variants for a given panel. The coordinates anonymizing model must be generated before calling this option using "*-generate_coordinate_anonymizer_model*". Output panel is saved to the output prefix.  All variants in the panel must have a corresponding entry in the anonymized coordinate mapping model.

#### ./ProxyTyper.sh -deanonymize_coordinates [Panel prefix] [Output prefix]
De-anonymizes the coordinates of all variants for a given panel.

#### ./ProxyTyper.sh -generate_tag_augmenter_model [Tag variants BED] [All panel variants BED]
Build a typed variant augmentation model for the tags variants and the set of all panel variants. Model is saved under ${MODELS_DIR} defined in *PROXYTYPER.ini* file.

#### ./ProxyTyper.sh -augment_tag_variants [Panel prefix] [Output prefix]
Augments typed variants in the panel and saves the new panel to output prefix. Untyped variants are copied as they were.

#### ./ProxyTyper.sh -copy_panel [Panel prefix] [Output prefix]
Copies a panel and saves it to output prefix.

#### ./ProxyTyper.sh -intersect [Source BED] [Destination BED] [Region selector (Reg1/Reg2/Reg12)] [Overlap BED]
Intersects source BED file (Reg1) with destination BED (Reg2) file. Output region selector describes what should be saved from the intersection.

#### ./ProxyTyper.sh -exclude [Source BED] [Exclusion regions BED] [Output BED]
Excludes the exclusion regions from the source regions and saves to output BED file.

#### ./ProxyTyper.sh -decompose_variants [Tag variants BED] [Tag/Target Reference Panel prefix] [Output prefix]
Decomposes the untyped variants in the typed+untyped variants panel and saves to output prefix. Typed variants are copied as intact. The decomposition model is saved 

#### ./ProxyTyper.sh -recompose_BEAGLE_VCF_variants [BEAGLE VCF file] [BEAGLE imputed tag sample subject id's] [Untyped variant decomposition interval file] [Output prefix]
Takes a BEAGLE's imputed VCF.gz file from a decomposed panel and saves to a recomposed panel.

#### ./ProxyTyper.sh -run_BEAGLE [Query typed variant panel prefix] [Reference typed/untyped variant panel prefix] [BEAGLE map file, e.g., ${ANONYMIZED_COORDINATES_DIR}/20_BEAGLE.map] [Impute targets? (true/false)] [Phased query tag? (0/1)] [BEAGLE directory]
Sets up the necessary files for running BEAGLE using query panel's typed variants, reference panel's typed+untyped panel, the BEAGLE map file. BEAGLE is run to save the results to the BEAGLE directory.

#### ./ProxyTyper.sh -haplo2geno [HC Panel prefix] [GC Output prefix]
Converts a phased panel to an unphased panel and saves to output prefix.

#### ./ProxyTyper.sh -summarize_accuracy [prxytypr accuracy list file] [Min MAF] [Max MAF]
Summarizes average genotype-level R2 for the variants in the MAF range defined by min_MAF and max_MAF arguments from accuracy file generated by "*-genotype_R2*" option and write to the screen.

#### ./ProxyTyper.sh -genotype_R2 [Target variants BED file] [Imputed genotypes matbed] [Known genotypes matbed] [Output file]
Calculates accuracy statistics by comparing the imputed and known genotypes panel for the typed variants in the BED file and saves to output file. 

#### ./ProxyTyper.sh -delete_panel [Panel prefix]
Permanently deletes a panel with prefix.

#### ./ProxyTyper.sh -sort_panel [Panel prefix] [Output prefix]
Sorts the variants in a panel and saves to the output prefix.

#### ./ProxyTyper.sh -setup_BEAGLE_files [Haplocoded query tag genotypes matbed] [Query sample ids] [Haplocoded reference tag+target genotype matbed] [Reference sample ids] [Use phased(1)/unphased(0) genotypes] [Output directory]
This option sets up BEAGLE VCF files for running BEAGLE in output directory.

#### ./ProxyTyper.sh -export_VCF [Panel identifier] [Assembly ID (e.g., hg19)] [Phased output? (0/1)] [Output VCF.gz file path]
Exports a panel to a VCF.gz file using the assembly identifier. Users can choose to save the panel as phased or unphased with the flag.




