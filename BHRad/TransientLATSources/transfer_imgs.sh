#!/bin/bash

mkdir ./unassociated/

for i in 4 5 6 13 15 16 19 23 25 26 27 29 31 33 35 36 46 50 58 63 64 76 82 84 95 99 100 104 108 112 116 121 123 126 127 130 132 134 #141 142
do
	mv ./apjsac072afigset17/*_"$i"_hr* ./unassociated/
done
