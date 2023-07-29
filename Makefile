all: ProxyTyper 

CC = g++
exec_name = bin/ProxyTyper_Release
lib_flags = -lz -lpthread -lgsl -lgslcblas
LIB_DIR = src

comp_flags = -c -O3 -Wall

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/prxytypr_main.o \
${LIB_DIR}/prxytypr_proxytyper.o \
${LIB_DIR}/prxytypr_haplotype_resampling_utils.o \
${LIB_DIR}/prxytypr_typed_variant_proxization_utils.o \
${LIB_DIR}/prxytypr_typed_variant_proxization_utils_MT.o \
${LIB_DIR}/prxytypr_untyped_variant_proxization_utils.o \
${LIB_DIR}/prxytypr_reidentification_attack_utils.o \
${LIB_DIR}/prxytypr_variant_proxy_attack_per_kmer_tracking_utils.o \
${LIB_DIR}/prxytypr_x_chisqr.o \
${LIB_DIR}/prxytypr_x_sigmoid.o \
${LIB_DIR}/prxytypr_file_utils.o \
${LIB_DIR}/prxytypr_ansi_thread.o \
${LIB_DIR}/prxytypr_xlog_math.o \
${LIB_DIR}/prxytypr_ansi_string.o \
${LIB_DIR}/prxytypr_exception_obj.o \
${LIB_DIR}/prxytypr_annot_region_tools.o \
${LIB_DIR}/prxytypr_variation_tools.o \
${LIB_DIR}/prxytypr_genome_sequence_tools.o \
${LIB_DIR}/prxytypr_nucleotide.o \
${LIB_DIR}/prxytypr_histogram.o \
${LIB_DIR}/prxytypr_nomenclature.o \
${LIB_DIR}/prxytypr_genomics_coords.o \
${LIB_DIR}/prxytypr_rng.o \
${LIB_DIR}/prxytypr_seed_manager.o 

ProxyTyper: ${objs}
	@echo Linking executable $@
	@${CC} -O3 ${lib_flags} -o ${exec_name} ${objs}

clean:
	@echo Cleaning..
	@rm -f ${objs} ${exec_name} 
