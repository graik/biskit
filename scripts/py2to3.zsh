!#/usr/bin/env zsh

echo "USAGE: "
echo "py2to3.zsh biskitModule.py"

TARGETFOLDER=~/data/py/biskit3/biskit
SOURCEFOLDER=~/data/py/biskit2/Biskit

2to3-3.5 -n -W -o $TARGETFOLDER $SOURCEFOLDER/$1

~/data/py/biskit3/scripts/replace_imports.py $TARGETFOLDER/$1 -replacefile

meld $TARGETFOLDER/$1 $SOURCEFOLDER/$1
