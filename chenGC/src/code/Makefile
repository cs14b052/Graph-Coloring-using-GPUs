include ../common.mk
INCLUDES += -I../../cub-1.1.1
B40_DIR=./back40computing-read-only
B40C_INC=-I$(B40_DIR) -I$(B40_DIR)/test
EXTRA := $(NVFLAGS) $(CFLAGS) $(INCLUDES) -DITERATIONS=1 -DBLKSIZE=1024
SRC=main.cu

all: chengc

chengc:
	$(NVCC) $(EXTRA) $(SRC) -DVARIANT=8 -o $@
	mv $@ $(BIN)

clean:
	rm -f $(BIN)/data_best 

