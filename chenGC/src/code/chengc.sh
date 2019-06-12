#!/bin/bash
rm -f ../../../outputs/output_chengc.out
make
for file in ../../../inputs/*.mtx;
do
  ../../bin/chengc $file >> ../../../outputs/output_chengc.out
done
