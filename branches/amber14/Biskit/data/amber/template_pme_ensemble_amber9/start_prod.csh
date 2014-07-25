#!/usr/bin/env zsh

for f in 2*; do
	cd $f
	./start.csh > "start_$f.out" &
	cd ..
done
