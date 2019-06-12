BIN=../bin/

all: baseline baselineOpt1 sirgpervertex

baseline: 
	nvcc -O3 -w -arch=sm_60 graphColoringVariants.cu Baseline.cu -lcurand -lcusparse -o $@
	mv $@ $(BIN)

baselineOpt1:	
	nvcc -O3 -w -arch=sm_60 graphColoringVariants.cu Baseline+Opt1.cu -lcurand -lcusparse -o $@
	mv $@ $(BIN)

sirgpervertex:
	nvcc -O3 -w -arch=sm_60 graphColoringVariants.cu sirgPerVertex.cu -lcurand -lcusparse -o $@
	mv $@ $(BIN)

clean:
	rm -f $(BIN)/baseline $(BIN)/baselineOpt1 $(BIN)/sirgpervertex
