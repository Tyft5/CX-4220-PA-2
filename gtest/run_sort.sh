#!/usr/bin/env bash

echo "Compiling..."
make

MPIRUN=/usr/lib64/openmpi/bin/mpirun

echo "Running..."
n=10000
p=6
c=0
while [ $p -le 48 ]; do
	echo "" > out_p-$p.txt
	while [ $n -le 1000000 ]; do
		while [ $c -lt 5 ]; do
			$MPIRUN -np $p --hostfile $PBS_NODEFILE ./sort -r -n $n >> out_p-$p.txt
			let c=c+1
		done
		let n=n+10000
		let c=0
	done
	let p=p+6
	let n=10000
done

echo "Done."
