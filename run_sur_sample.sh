sortandindex() {
	grep -v '#' $1 | sort -k 1,1 -k 2,2n > body
	grep '#' $1 > head
	cat head body > $1
	rm head body
	bgzip -c $1 > $1.gz
	tabix $1.gz
}
# extract different svtype from $1 and generate INS/DEL/INV/DUP.vcf.gz
parse_svtype() {
	for svtype in {INS,DEL,INV,DUP}
	do
		grep '#' $1 > head
		grep SVTYPE=$svtype $1 > body
		cat head body > $svtype.vcf
		rm head body
		bgzip -c $svtype.vcf > $svtype.vcf.gz
		tabix $svtype.vcf.gz
	done
}

run_sur() {
	# run_sur support coverage
	for i in {CHM1,CHM13,HG00268,HG00514,HG00733,HG01352,HG02059,HG02106,HG02818,HG04217,HX1,NA12878,NA19240,NA19434}
	do
		#echo $i
		/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/$2/$i.fofn 100 $1 1 0 1 30 benchmark/survivor_vcf/$i\_$2.vcf
		grep -v '#' benchmark/survivor_vcf/$i\_$2.vcf | sort -k 1,1 -k 2,2n > body
		grep '#' benchmark/survivor_vcf/$i\_$2.vcf > head
		cat head body > benchmark/survivor_vcf/$i\_$2.vcf
		rm head body
	done
	ls benchmark/survivor_vcf/*_$2.vcf > benchmark/survivor_vcf/$2.fofn

	/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge benchmark/survivor_vcf/$2.fofn 100 1 1 0 1 30 benchmark/survivor_vcf/sur_$2.vcf
	sortandindex benchmark/survivor_vcf/sur_$2.vcf
}

main() {
	# compress_coverage folder
	for i in {5x,10x,15x,20x,30x}
	do
		run_sur 3 $i
	done
	echo 'finish merge'

	# coverage
	rm -r cmp_5x/ cmp_10x/ cmp_15x/ cmp_20x/ cmp_30x/
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_5x.vcf.gz -o cmp_5x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_10x.vcf.gz -o cmp_10x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_15x.vcf.gz -o cmp_15x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_20x.vcf.gz -o cmp_20x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_30x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	# type
	parse_svtype benchmark/survivor_vcf/sur_30x.vcf
	rm -r cmp_ins/ cmp_del/ cmp_inv/ cmp_dup/
	truvari bench -b benchmark/nstd_INS.vcf.gz -c INS.vcf.gz -o cmp_ins -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_DEL.vcf.gz -c DEL.vcf.gz -o cmp_del -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_INV.vcf.gz -c INV.vcf.gz -o cmp_inv -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_DUP.vcf.gz -c DUP.vcf.gz -o cmp_dup -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	#truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_all -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	# length
	rm -r cmp_1/ cmp_2/ cmp_3/ cmp_4/ cmp_5/ cmp_6/
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_1 -p 0 -r 1000 -s 30 --sizemax 99 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_2 -p 0 -r 1000 -s 100 --sizemax 499 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_3 -p 0 -r 1000 -s 500 --sizemax 999 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_4 -p 0 -r 1000 -s 1000 --sizemax 4999 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_5 -p 0 -r 1000 -s 5000 --sizemax 9999 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_6 -p 0 -r 1000 -s 10000 --sizemax 10000000 --multimatch
	#truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_7 -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt \
					cmp_ins/summary.txt cmp_del/summary.txt cmp_inv/summary.txt cmp_dup/summary.txt \
					cmp_1/summary.txt cmp_2/summary.txt cmp_3/summary.txt cmp_4/summary.txt cmp_5/summary.txt cmp_6/summary.txt benchmark/survivor_vcf/answer.txt
	
}
rm benchmark/survivor_vcf/*
main
