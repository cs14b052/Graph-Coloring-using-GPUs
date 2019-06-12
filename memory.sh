#!/bin/bash

declare -a lines1
while IFS= read -r line
  do
    lines1+=("$line")
  done < "$1"

declare -a lines2
while IFS= read -r line
  do
    lines2+=("$line")
  done < "$2"

declare -a filenames1
declare -a filenames2
declare -a values1
declare -a values2
for line in "${lines1[@]}"
  do
    filenames1+=("$(cut -d' ' -f1 <<<"$line")")
    values1+=("$(cut -d' ' -f4 <<<"$line")")
  done

for line in "${lines2[@]}"
  do
    filenames2+=("$(cut -d' ' -f1 <<<"$line")")
    values2+=("$(cut -d' ' -f4 <<<"$line")")
  done

len1=${#filenames1[@]}
len2=${#filenames2[@]}


if [ $len1 != $len2 ]
  then
  echo "Number of graphs in both the files need to be equal"
  exit 1
fi

for (( i=0; i<$len1; i++ ));
  do
    if [ "${filenames1[$i]}" != "${filenames2[$i]}" ]
      then
      echo "Graph Names need to be in same order in both files"
      exit 1
    fi
  done

##### MEMORY REDUCTION CALCULATION #####

declare -a output
geomean=1
count=0
for (( i=0; i<$len1; i++ ));
  do
   output+=($(bc -l <<< "((${values2[$i]} - ${values1[$i]})*100)/${values2[$i]}"))
   if [ "${output[$i]}" != 0 ]
     then
        geomean=($(bc -l <<< "(($geomean*${output[$i]}))")) 
        count=($(bc -l <<< "$count+1"))
   fi
  done

output+=($(bc -l <<< "e( l($geomean)/$count )"))

for (( i=0; i<$len1; i++ ));
  do
    echo -ne "${filenames1[$i]} "
    printf "%0.2f\n" "${output[$i]}"
  done

echo -ne "GeoMean "
printf "%0.2f\n" "${output[$len1]}"
