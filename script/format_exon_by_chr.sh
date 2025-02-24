
nbLine=$(cat $2 | wc -l)

for i in $( seq 1 $nbLine )
do
	chtmp=$(sed -n $i"p" $2)
	cat $1 |  grep "	$chtmp	" | cut -f2,9,10,11,16 1>$3/"tmp_exon_"$chtmp".gtf"
done


for i in $( seq 1 $nbLine )
do
	chtmp=$(sed -n $i"p" $2)
	sort -k1 -n $3/"tmp_exon_"$chtmp".gtf" > $3/$chtmp"_exon.gtf"
	rm $3/"tmp_exon_"$chtmp".gtf"
done

