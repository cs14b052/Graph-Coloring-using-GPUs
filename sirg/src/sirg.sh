#!/bin/bash
rm -f ../../outputs/output_sirg.out
make
for file in ../../inputs/*.mtx;
do
  ../bin/sirg $file 0 >> ../../outputs/output_sirg.out
done
