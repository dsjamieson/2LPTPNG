#!/bin/bash
i=0
for seed in $( awk '{print $0}' seeds.txt )
do
	sed -i "s/Seed.*/Seed    ${seed}/g" test_potential.param
	for mode in "LC" "EQ" "OR"
	do
		echo $mode $i
		../2LPTNG${mode} test_potential.param
		mv Meshes/potential_NG_${mode}.dat Meshes/${mode}/potential_${mode}_${i}.dat
		if [ $mode == "LC" ]
		then
			mv Meshes/potential_G.dat Meshes/G/potential_G_${i}.dat
		fi
	done
	i=$(( $i + 1 ))
done
