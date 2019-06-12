#!/bin/bash

rm -f outputs/results_6_a.out
echo 'Speedup of SIRG over ChenGC' >> outputs/results_6_a.out
./speedup.sh outputs/output_sirg.out outputs/output_chengc.out >> outputs/results_6_a.out

rm -f outputs/results_6_b.out
echo 'Speedup of SIRG over csrcolor' >> outputs/results_6_b.out
./speedup.sh outputs/output_sirg.out outputs/output_csrcolor.out >> outputs/results_6_b.out
