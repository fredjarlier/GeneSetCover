filename="$1"
refgene="$2"

cut -f1 --delimiter='=' $filename 1>$3/ref.cut

ref=$3/ref.cut

cut -f2 --delimiter='=' $filename 1>$3/ch.cut

ch=$3/ch.cut

nbLine=$(cat $ref | wc -l)

for i in $(seq 1 $nbLine)
do
	reftmp=$(sed -n $i"p" $ref)
	chtmp=$(sed -n $i"p" $ch)
	cat $refgene | grep "	$reftmp	" | sed "s/	$reftmp	/	$chtmp	/g" 1>>$3/new_refGene.txt
	
done

cat $3/new_refGene.txt | grep "NM_" |sed 's/NM_//g' 1>$3/new_refGene_only_NM.txt


rm $ref

