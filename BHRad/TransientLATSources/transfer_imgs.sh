#!/bin/bash

mkdir ./correct_unassociated/

for i in 4 5 6 13 15 16 19 23 25 26 27 29 32 34 44 48 56 61 62 74 80 91 95 96 100 104 108 112 117 119 122 123 126 128 130
do
	mv ./unassociated/*_"$i"_hr* ./correct_unassociated/
done
