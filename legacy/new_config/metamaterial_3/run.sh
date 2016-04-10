#!/bin/bash

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

sed -i "s/^eps_slae = .*/eps_slae = 1e-1/" config.ini
./vfem

for EPS in "5e-2" "1e-2" "5e-3" "1e-3" "5e-4" "1e-4" "5e-5" "1e-5" "5e-6" "1e-6" "5e-7" "1e-7" "5e-8" "1e-8" "5e-9" "1e-9"
do
	sed -i "s/^eps_slae = .*/eps_slae = ${EPS}/" config.ini
	./vfem -continue
done

