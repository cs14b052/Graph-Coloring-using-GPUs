#!/bin/bash
if [ $# != 1 ]
  then
  echo "Usage: $(basename $0) <baseline | baselineOpt1 | sirgpervertex>"
  exit 1
fi
rm -f ../../outputs/output_$1.out
make $1
for file in ../../inputs/*.mtx;
do
  ../bin/$1 $file 0 >> ../../outputs/output_$1.out
done
