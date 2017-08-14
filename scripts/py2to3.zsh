!#/usr/bin/env zsh

echo "USAGE: "
echo "py2to3.zsh sourceModule.py destModule.py"
echo "run from biskit3 root folder, assumes biskit2 folder"
echo "available for diff one level up (../biskit2/)"

git mv $1 $2
git commit -m 'move Biskit -> biskit without modification'

2to3-3.5 -n -W $2

~/data/py/biskit3/scripts/replace_imports.py $2 -replacefile

meld $2 ../biskit2/$1
