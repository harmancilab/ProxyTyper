# Anonymizing Coordinates and Genetic Maps

Anonymizing coordinates and genetic maps is essential to ensure that none of the typed or untyped variants can be immediately located among the variant set uploaded to the imputation server. 

ProxyTyper performs coordinate anonymization uniformly placing the variants within a predefined chromosome length.

For anonymizing the genetic maps, ProxyTyper adds noise to the original genetic distances (in centiMorgans) and saves this file. Of note, this genetic map still maintains the original length since noise is much (much) smaller than the total genetic length of the chromosomes and it is fairly straightforward to estimate which chromosome each variant resides on. *Q: Then why do you add noise to genetic maps??* The noise addition aims to introduce local ambiguity compared to the publicly available genetic maps up to a certain level of difference. This makes it more challenging to map the genetic distances directly to genomic coordinates by cross referencing public genetic maps to anonymized genetic maps.

Anonymization steps are model-based, as we saw earlier typed variant protection mechanisms. To build the anonymization model, we use following command:
```
./ProxyTyper.sh -generate_coordinate_anonymizer_model TYPED_VARIANTS.bed TYPED_AND_UNTYPED_VARIANTS.bed
```
Upon completion, this command saves the original-to-anonymized coordinate mapping BED files and the anonymized genetic maps. Of note, ProxyTyper uses PLINK formatted genetic maps (4-column formatted as "chrom rsid cM posn"). The anonymized genetic maps are saved only for the typed variants.

After anonymization models are generated, we can anonymize the coordinates of a panel using following command:
```
PANEL_ID=vcf_imported_panel
./ProxyTyper.sh -anonymize_coordinates ${PANEL_ID} ${PANEL_ID}_anon_coords
```
you can now view the variant coordinates (vcf_imported_panel_anon_coords_variants.bed) and see that the coordinates are mapped between 1 and length of the anonymized chromosome.

To de-anonymize and map the coordinates back to the original coordinates, we use following:
```
PANEL_ID=vcf_imported_panel
./ProxyTyper.sh -deanonymize_coordinates ${PANEL_ID}_anon_coords ${PANEL_ID}_anon_coords_deanonym
```
## Parameters that govern Coordinate and Genetic Map Aonymization:
In PROXYTPER.ini file, following entries define the parameters used in coordinate and genetic map anonymization:
```
coord_anon_cM_noise_SD=0.05 : Standard Deviation of the noise (in cMs) added to the original genetic map for each variant.
l_anonym_chrom=100000000    : The length of the anonymized chromosome in base pairs. 
```
**REMARK**: After a panel is mapped to anonymized coordinates, all further operations should use the same coordinate system. For instance, if the sites run imputation, the genetic maps corresponding to the appropriate coordinates system should be use while running BEAGLE.

**REMARK**: *l_anonym_chrom* should not be made too large or too small. If it is too small, mechanisms that increase typed and untyped variant number may not be effective. If it is too large compared to the number of variants, BEAGLE may complain about integer value overflow issue. In general, 100 megabase is a reasonable length for all chromosomes whose length are below 100 megabase and 250 megabase is reasonable for longer chromosomes.

---

You can use *Coordinate_Genetic_Map_Anonymization.sh* script as a demonstration of the commands we discussed above.