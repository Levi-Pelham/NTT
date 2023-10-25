#!/bin/bash
count=0
files=0
total=0;

for file in *.cpp
do
	echo "assembling $file ..."
	g++ -O3 -S -fno-unroll-loops $file -o assembly/"${file%.*}".s 
	let files=files+1
done

echo "Ran $files Files."
   