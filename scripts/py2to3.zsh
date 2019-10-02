#!/bin/bash

echo "USAGE: "
echo "py2to3.zsh Biskit/sourceModule.py biskit/destModule.py"
echo "run from biskit3 root folder, assumes archive_biskit2 folder"
echo "available for diff one level up (../biskit2/Biskit/...)"

exists() { type -t "$1" > /dev/null 2>&1; }

git mv archive_biskit2/$1 $2
git commit -m 'move Biskit -> biskit without modification'

2to3 -n -W $2

~/data/py/biskit3/scripts/replace_imports.py $2 -replacefile

if exists meld; then
    meld $2 ../archive_biskit2/$1
elif exists bcomp; then
    bcomp $2 ../archive_biskit2/$1
fi
