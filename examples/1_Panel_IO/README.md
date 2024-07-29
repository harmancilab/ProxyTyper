# Panel I/O (Panel_IO.sh)

Panel input and output is important to read external datasets and saving results. We do recognize that some of the operations can be performed with existing tools, e.g., tabix and vcf/bcftools. However, building new tools relying on these requires many more dependencies that turned out complex in a non-self managed environment and we chose to implement these options for completeness. 

We follow an example for importing, subsetting and saving the 1000 Genomes Project VCF file for chromosome 20. Before running, make sure to test the script:
```
./ProxyTyper.sh
```

This should list the options from ProxyTyper and quit. If not, there should be an error message. The error should describe the reason for the issue (executable not found in PATH, or configuration file is missing.) *installation/* folder contains the instructions to install ProxyTyper before first run.

## Download 1kG Panel VCF File
```
wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

## Import VCF to ProxyTyper panel matrix format:
When we have a raw VCF file, we first "import" the VCF file and save it in the native genotype format used by ProxyTyper:
```
all_REF_panel_ID=KG_chr20
is_panel_phased=1
update_variant_ID=1
./ProxyTyper.sh -import_VCF ${ref_VCF} ${is_panel_phased} ${update_variant_ID} ${all_REF_panel_ID}
```
Here, we first set the panel identifier. Next, we indicate that the panel is phased. Finally, we set a flag that tells the script to update variant identifiers with ref/alt allele and allele frequency information. After we VCF import is complete, you can see that ProxyTyper created 3 files whose names are specially formatted:
```
KG_chr_20_variants.bed          : The variants in the panel.
KG_chr_20_subjects.list         : Subjects list extracted from the VCF file.
KG_chr_20_genotypes.matrix.gz   : Genotypes in a plain matrix.
```
These three files should not be modified manually after they are created. We use ProxyTyper to process and proxize these panels as we will see in the later examples.

## ProxyTyper Uses Panel Identifiers to Easily Manage Flow of Genotype Panels in Proxizing Protocols
ProxyTyper's protocols makes extensive use of identifier-based accession to the panels. We found this to be the easiest to manage the panels when protocols have many steps. The basic idea is to process each panel, save the result to a panel with an updated identifier, then update the current panel identifier to the new panel identifier. This way, we can easily track the names of panels as they flow through the protocols. Although ProxyTyper writes and reads many files underneath this process, the code simply processes and updates the panel identifiers.  

##  Extract unique variants and Sort the Panel:
The panel we imported contains some variants that share locations. We should filter these out:
```
./ProxyTyper.sh -uniquefy_panel_variants ${all_REF_panel_ID} ${all_REF_panel_ID}_unique
all_REF_panel_ID=${all_REF_panel_ID}_unique
```
Note that we update the panel identifier after we process the panel. This makes it easier to follow the pipeline protocol since we do not have to give a new name to every one of the panels. If there is no use for the previously generated panels, these can be deleted before panel identifiers are updated. 

```
./ProxyTyper.sh -sort_panel ${all_REF_panel_ID} ${all_REF_panel_ID}_sorted
all_REF_panel_ID=${all_REF_panel_ID}_sorted
```
Similar to before, we sort the variants and update the panel identifier.

## Subset the Variants and Subjects
Subset the variants in the region 20:20,000,000-30,000,000:
```
awk {'if($2>20000000 && $3<30000000){print $0}'} ${all_REF_panel_ID}_variants.bed > subset_variants.bed

./ProxyTyper.sh -subset_panel_variants ${all_REF_panel_ID} subset_variants.bed ${all_REF_panel_ID}_variant_subset
all_REF_panel_ID=${all_REF_panel_ID}_variant_subset
```

Next, select random 1000 subjects and save:
```
shuf ${all_REF_panel_ID}_subjects.list | head -n 1000 > subset_subjects.list

# extract the subset panel.
./ProxyTyper.sh -subset_panel_subjects ${all_REF_panel_ID} subset_subjects.list ${all_REF_panel_ID}_subject_subset

# Update the panel identifier in the protocol.
all_REF_panel_ID=${all_REF_panel_ID}_subject_subset
```
Similar to before, we subset the variants and update the panel identifier in both subsetting operations.

## Export the Processed Panel to a VCF file:
We export the panel to a VCF file using the panel identifier:
```
./ProxyTyper.sh -export_VCF ${all_REF_panel_ID} hg19 1 ${all_REF_panel_ID}.vcf.gz
```
This write a vcf file for the panel with identifier stored in *${all_REF_panel_ID}*. You can use bgzip, etc on this panel to process it further with external tools.

---

These commands are implemented in the script named *Panel_IO.sh* that users can review and test. This script also includes numerous error checks at each step to make sure errors are caught. We observed that these checks are very useful while testing long protocols and we strongly encourage users to check for successful completion of the commands in pipelines.