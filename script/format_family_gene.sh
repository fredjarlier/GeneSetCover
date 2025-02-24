cut --complement -s -f1,3,4,5,6,7,8,9,10,11,12,14,15,16 $1  1> $3"/tmp.txt"

awk -F'	' 'BEGIN {OFS="	"} { print $2,$1 }' $3"/tmp.txt" 1> $3"/tmp2.txt"

rm $3"/tmp.txt"

sort -k1 -d -r $3"/tmp2.txt" > $3"/tmp3.txt"

rm $3"/tmp2.txt"

gcc -o $2".o" $2".c"

$2".o" $3"/tmp3.txt" 1>$3/"tmp4.txt"

rm $3"/tmp3.txt"


sort -k2 -n $3"/tmp4.txt" > $3"/tmp5.txt"

rm $3"/tmp4.txt"

awk -F'	' 'BEGIN {OFS="	"} { print $2,$3 }' $3"/tmp5.txt" 1> $3"/family_gene.txt"

rm $3"/tmp5.txt"