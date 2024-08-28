# Running BEAGLE Imputation Server

BEAGLE is a highly popular and well maintained genotype imputation method that performs phasing and imputation in one package <span style="font-size:10px;">(yes, we do believe it is pretty awesome.)</span>. ProxyTyper uses BEAGLE as the default imputation method. ProxyTyper implements the options to run BEAGLE automatically with one command. This process should, in principle, take place at the imputation server. When resampling is used at the query site, it is necessary to perform a local imputation at the query site to estimate the imputed genotypes for the query subjects. 

Before running this example, make sure there is a java interpreter installed on the system:
```
java --version
```
This command should run and return the version that is used with java.

Given the proxized reference and query panels (with all mechanisms that we discussed in previous examples applied), we simulate running BEAGLE on the imputation server with the following command:
```
QUERY_PANEL_ID=QUERY_proxy
REFERENCE_PANEL_ID=REFERENCE_proxy
anonymized_genetic_map=anon_genetic_maps/20_BEAGLE.map
do_impute=true
phased_query=0  # Query is always sent as unphased to the imputation server, but it does not have to be. Set this to 1 to perform phased imputation with BEAGLE (it is much faster).
./ProxyTyper.sh -run_BEAGLE ${QUERY_PANEL_ID} ${REFERENCE_PANEL_ID} ${anonymized_genetic_map} ${do_impute} ${phased_query} PROXY_IMPUTED_BEAGLE_DIR
```
This command should setup the necessary files and run BEAGLE.

After BEAGLE is run, we most likely have to recompose the untyped variants and deanonymize their coordinates:
```
# We have seen these commands before:
./ProxyTyper.sh -recompose_BEAGLE_VCF_variants PROXY_IMPUTED_BEAGLE_DIR/imputed.op.vcf.gz ${QUERY_PANEL_ID}_subjects.list ${REFERENCE_PANEL_ID}.interval ${QUERY_PANEL_ID}_imputed
./ProxyTyper.sh -deanonymize_coordinates ${QUERY_PANEL_ID}_imputed ${QUERY_PANEL_ID}_imputed_deanonymized
```

We can now calculate the accuracy of the untyped variants using following command:
```
./ProxyTyper.sh -genotype_R2 UNTYPED_VARIANTS.bed ${QUERY_PANEL_ID}_imputed_deanonymized QUERY_PANEL_KNOWN_GENOTYPES imputation_accuracy.txt
```
In this command, UNTYPED_VARIANTS.bed contains only the untyped variants that we will score and QUERY_PANEL_KNOWN_GENOTYPES is the panel identifier for the known genotypes of the untyped variants for the query panel sample. The accuracy statistics are written to *imputation_accuracy.txt*.

Users can summarize the accuracy statistics:
```
# Get the average R2 for variants with 0.005>MAF>0.001
./ProxyTyper.sh -summarize_accuracy imputation_accuracy.txt 0.001 0.005
```

Following options are used in PROXYTYPER.ini configuration file to setup BEAGLE:
```
BEAGLE_JAR=beagle.13Mar20.38e.jar   : BEAGLE's jar java archive file.
JAVA_MEM_FLAG=Xmx400000m            : The memory option for running java, set to an appropriate size to give enough memory to java (400000m -> 400 gigabytes)
```

**REMARK:** while we chose to use BEAGLE as the imputation tool, other methods (e.g., eagle+Minimac4) should be suitable to be used with proxized panels.

---

You can use *Running_BEAGLE_imputation.sh* script as a demonstration of the commands we discussed above.