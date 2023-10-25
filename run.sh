#!/bin/bash
count=0
files=0
total=0;

for file in *.cpp
do
   echo "running $file ..."
   g++ -O3 -fno-unroll-loops $file && ./a.out 
   let files=files+1
done

echo "Ran $files Files."
   