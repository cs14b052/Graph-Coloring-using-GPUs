#!/bin/bash
rm -f ../../outputs/output_csrcolor.out
make
for file in ../../inputs/*.mtx;
do
  ../bin/csrcolor $file >> ../../outputs/output_csrcolor.out
done
