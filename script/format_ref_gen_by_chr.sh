

cut --complement -s -f1,4,7,8,9,10,11,12,13,14,15,16 $2 1>$3/tmp.ref.out

nbLine=$(cat $4 | wc -l)

for i in $( seq 1 $nbLine )
do
    chtmp=$(sed -n $i"p" $4)
    cat $3/tmp.ref.out |  grep "	$chtmp	" 1>$3/"B"$chtmp".gtf"
   
done

rm $3/tmp.ref.out


for i in $( seq 1 $nbLine )
do
	chtmp=$(sed -n $i"p" $4)
	sort -k3 -n $3/"B"$chtmp".gtf" > $3/$chtmp".gtf"
	cut --complement -s -f2 $3/$chtmp".gtf" 1>$3/"B"$chtmp".gtf"
	awk -F'	' 'BEGIN {OFS="	"} { print $2,$3,$1 }' $3/"B"$chtmp".gtf" 1>$3/$chtmp".gtf"

	rm $3/"B"$chtmp".gtf"
done

rm $4


