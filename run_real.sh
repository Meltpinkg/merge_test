
# extract different svtype from $1 and generate INS/DEL/INV/DUP_$2.vcf.gz
parse_svtype() {
	for svtype in {INS,DEL,INV,DUP}
	do
		grep '#' $1 > head
		grep SVTYPE=$svtype $1 > body
		cat head body > $svtype\_$2.vcf
		rm head body
		bgzip -c $svtype\_$2.vcf > $svtype\_$2.vcf.gz
		tabix $svtype\_$2.vcf.gz
	done
}
sortandindex() {
	grep -v '#' $1 | sort -k 1,1 -k 2,2n > body
	grep '#' $1 > head
	cat head body > $1
	rm head body
	bgzip -c $1 > $1.gz
	tabix $1.gz
}
bench_low() {
    # bench [base.vcf] [compare.vcf]
    echo 'start bench'
    bgzip -c $1 > $1.gz
    tabix $1.gz
    sortandindex $2
    parse_svtype $1 base
	parse_svtype $2 comp
    rm -r cmp_ins/ cmp_del/ cmp_inv/ cmp_dup/ cmp_all/
    truvari bench -b INS_base.vcf.gz -c INS_comp.vcf.gz -o cmp_ins -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
    truvari bench -b DEL_base.vcf.gz -c DEL_comp.vcf.gz -o cmp_del -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
    truvari bench -b INV_base.vcf.gz -c INV_comp.vcf.gz -o cmp_inv -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
    truvari bench -b DUP_base.vcf.gz -c DUP_comp.vcf.gz -o cmp_dup -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
    truvari bench -b $1.gz -c $2.gz -o cmp_all/ -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
    python truvari/ss.py cmp_ins/summary.txt cmp_del/summary.txt cmp_inv/summary.txt cmp_dup/summary.txt cmp_all/summary.txt answer.txt
	echo 'write to answer'
}

bench() {
    # bench [base.vcf] [compare.vcf]
    echo 'start bench'
    bgzip -c $1 > $1.gz
    tabix $1.gz
    sortandindex $2
    truvari bench -b $1.gz -c $2.gz -o cmp_all/ -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	echo 'finish bench'
}

bench $1 $2