!#/usr/bin/env zsh

echo "USAGE: "
echo "py2to3.zsh biskitModule.py"

biskit3=~/data/py/biskit3/biskit
Biskit2=~/data/py/biskit3/Biskit
Biskit2Copy=~/data/py/biskit2/Biskit

git mv $Biskit2/$1 $biskit3/$1
git commit -m 'move Biskit -> biskit without modification'

2to3-3.5 -n -W $biskit3/$1

~/data/py/biskit3/scripts/replace_imports.py $biskit3/$1 -replacefile

meld $biskit3/$1 $Biskit2Copy/$1
