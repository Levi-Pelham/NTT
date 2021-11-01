#!/bin/bash
files=0
count=0

for file in *}.txt
do
   let count=count+1
   sort -no $file $file
done

echo "Sorted $count Files."



