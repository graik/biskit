#!/usr/bin/env zsh

./lamboot_if_needed.py

for f in 2*; do
	cd $f
	./start.csh
	cd -
done

cd 3*
./start.csh
cd -
