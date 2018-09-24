#!/bin/bash

echo "USAGE: "
echo "py2to3.zsh sourceModule.py destModule.py"
echo "run from biskit3 root folder, assumes biskit2 folder"
echo "available for diff one level up (../biskit2/)"

exists() { type -t "$1" > /dev/null 2>&1; }

git mv $1 $2
git commit -m 'move Biskit -> biskit without modification'

2to3 -n -W $2

~/data/py/biskit3/scripts/replace_imports.py $2 -replacefile

if exists meld; then
    meld $2 ../biskit2/$1
elif exists bcomp; then
    bcomp $2 ../biskit2/$1
fi
