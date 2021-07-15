grep -v '#' $1 | sort -k 1,1 -k 2,2n > body
grep '#' $1 > head
cat head body > $1
rm head body