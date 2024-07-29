# Imputation Protocols

This folder contains 4 protocols (2 for phased query and 2 for unphased query) that can be run after setting up data.

Following panels are expected to be present for these scripts to run
```
chr_id=20
PROTOCOL_REF_ID=${chr_id}_REF_tags_targets
PROTOCOL_QUERY_ID=${chr_id}_QUERY_tags
PROTOCOL_QUERY_targets=${chr_id}_QUERY_targets
```

After setting up the reference and query panels, the scripts can be run, for example, as following:
```
cd Phased_Panel_Imputation_Protocol
chmod 755 Phased_Protocol*.sh
N_REF_RESAMPLE_SIZE=10000
N_QUERY_RESAMPLE_SIZE=2000
./Phased_Protocol1.sh ${chr_id} ${N_REF_RESAMPLE_SIZE} ${N_QUERY_RESAMPLE_SIZE}
./Phased_Protocol2.sh ${chr_id} ${N_REF_RESAMPLE_SIZE} ${N_QUERY_RESAMPLE_SIZE}
```
Unphased protocols do not perform query resampling and use only reference panel resampling size:
```
cd Unphased_Panel_Imputation_Protocol
chmod 755 Unphased_Protocol*.sh
N_REF_RESAMPLE_SIZE=10000
./Unphased_Protocol1.sh ${chr_id} ${N_REF_RESAMPLE_SIZE}
./Unphased_Protocol1.sh ${chr_id} ${N_REF_RESAMPLE_SIZE}
```

**REMARK:** If it is not clear how to generate the appropriate panels, you can use the data setup script in Example 7, *example/7_Example_Imputation_Pipelines/Setup_Query_Ref_KG_data.sh* to set up these panels from 1000 Genomes Project data, which automatically downloads and sets up the necessary reference and query panels in appropriate names (also runs BEAGLE).

---

We recommend with getting familiarized with running ProxyTyper's main script (*ProxyTyper.sh*) by reviewing the *examples/* directory before running these scripts (you don't have to, obviously).
