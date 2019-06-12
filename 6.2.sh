#!/bin/bash

rm -f outputs/results_7.out
echo 'Speedup of Opt1 over Baseline' >> outputs/results_7.out
./speedup.sh outputs/output_baselineOpt1.out outputs/output_baseline.out >> outputs/results_7.out

echo  >>  outputs/results_7.out
echo 'Speedup of Opt1 + Opt2 over Baseline' >> outputs/results_7.out
./speedup.sh outputs/output_sirg.out outputs/output_baseline.out >> outputs/results_7.out

rm -f outputs/results_8.out
echo '% memory reduction in memory usage' >> outputs/results_8.out
./memory.sh outputs/output_sirg.out outputs/output_baselineOpt1.out >> outputs/results_8.out
