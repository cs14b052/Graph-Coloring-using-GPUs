#!/bin/bash

rm -f outputs/results_9.out

echo '% improvement in execution time' >> outputs/results_9.out
./improvExecTime.sh outputs/output_sirg.out outputs/output_sirgpervertex.out >> outputs/results_9.out

echo >> outputs/results_9.out 
echo '% reduction in memory usage' >> outputs/results_9.out
./memory.sh outputs/output_sirg.out outputs/output_sirgpervertex.out >> outputs/results_9.out
