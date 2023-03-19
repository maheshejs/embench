#python script creats ram init files for simulation and hardware
###############################################################
RISCV_PREFIX ?= riscv32-unknown-elf-
ELF_TO_HW_INIT_OPTIONS ?= $(RISCV_PREFIX) 0x80000000 65536
###############################################################

###############################################################
#Embench
#Assumes binaries are in the BENCHMARK_DIR
EMBENCH_DIR=.
EMBENCH_BENCHMARKS =  \
edn \
matmult-int \
nbody \
st \
ud
#aha-mont64 \
#crc32 \
#cubic \
#huffbench \
#md5sum \
#minver \
#nbody \
#nettle-aes \
#nettle-sha256 \
#nsichneu \
#picojpeg \
#primecount \
#qrduino \
#sglib-combined \
#slre \
#st \
#statemate \
#tarfind \
#ud \
#wikisort


#add file path to benchmarks
embench_bins = $(addprefix $(EMBENCH_DIR)/build/bin/, $(EMBENCH_BENCHMARKS))

#embench benchmarks copied into a bin folder to simplify makefile rules
.PHONY: build
build :
	cd $(EMBENCH_DIR);\
	./build_all.py -v --clean --use-vector $(USE_VECTOR) --builddir build --arch riscv32 --chip generic --board ri5cyverilator --cc /home/brumaire/build/gnu-fpu/bin/riscv32-unknown-elf-gcc --cflags="-c -O2 -ffunction-sections -march=rv32imafdv -mabi=ilp32d" --ldflags="-Wl,-gc-sections" --user-libs="-lm"
	mkdir -p $(EMBENCH_DIR)/build/bin
	$(foreach x,$(EMBENCH_BENCHMARKS), mv $(EMBENCH_DIR)/build/src/$(x)/$(x) $(EMBENCH_DIR)/build/bin/$(x);)

.PHONY: run
run :
	spike --isa=rv32imav --varch=vlen:4096,elen:64 ~/build/pk/riscv32-unknown-elf/bin/pk -s $(EMBENCH_DIR)/build/bin/$(BENCHMARK) ; echo $?
	
#Benchmarks built by build_embench
.PHONY : $(embench_bins)

.PHONY: clean
clean : 
	rm -rf build/ logs/

