BIN=../bin/

csrcolor: 
	nvcc -O3 -w -arch=sm_60 main.cu nvidia_csrcolor.cu -lcurand -lcusparse -o $@
	mv $@ $(BIN)

clean:
	rm $(BIN)/csrcolor
