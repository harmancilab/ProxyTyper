all: ProxyTyper 

CC = g++
exec_name = bin/ProxyTyper_Release
lib_flags = -lz -lpthread -lgsl -lgslcblas

comp_flags = -c -O3 -Wall

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
src/main.o \
src/prxytypr_proxytyper.o \
src/prxytypr_vcf_io_utils.o \
src/prxytypr_haplotype_resampling_utils.o \
src/prxytypr_haplotype_resampling_utils_optimized.o \
src/prxytypr_unphased_panel_resampling_utils.o \
src/prxytypr_typed_variant_proxization_utils.o \
src/prxytypr_typed_variant_proxization_utils_MT.o \
src/prxytypr_untyped_variant_proxization_utils.o \
src/prxytypr_reidentification_attack_utils.o \
src/prxytypr_variant_matrix_decomposition_utils.o \
src/prxytypr_variant_proxy_attack_per_kmer_tracking_utils.o \
src/prxytypr_ansi_cli.o \
src/prxytypr_utils.o \
src/prxytypr_ansi_thread.o \
src/prxytypr_ansi_mutex.o \
src/prxytypr_config.o \
src/prxytypr_xlog_math.o \
src/prxytypr_matrix_linalg_utils.o \
src/prxytypr_ansi_string.o \
src/prxytypr_annot_region_tools.o \
src/prxytypr_variation_tools.o \
src/prxytypr_imputation_utils.o \
src/prxytypr_genome_sequence_tools.o \
src/prxytypr_nucleotide.o \
src/prxytypr_histogram.o \
src/prxytypr_nomenclature.o \
src/prxytypr_genomics_coords.o \
src/prxytypr_rng.o \
src/prxytypr_seed_manager.o 

ProxyTyper: ${objs}
	@echo Linking executable $@
	@${CC} ${objs} ${lib_flags} -o ${exec_name} -O3

clean:
	@echo Cleaning..
	@rm -f ${objs} ${exec_name} 
