sortandindex() {
	grep -v '#' $1 | sort -k 1,1 -k 2,2n > body
	grep '#' $1 > head
	cat head body > $1
	rm head body
	bgzip -c $1 > $1.gz
	tabix $1.gz
}

run_jas() {
	# run_jas support coverage
	#--dup_to_ins genome_file=/home/tjiang/ref/hs37d5.fa
	for i in {CHM1,CHM13,HG00268,HG00514,HG00733,HG01352,HG02059,HG02106,HG02818,HG04217,HX1,NA12878,NA19240,NA19434}
	do
		echo $i
		jasmine file_list=input/$2/$i.fofn out_file=output_sample/jasmine_$i\_$2.vcf min_support=$1
	done
	ls output_sample/jasmine_*_$2.vcf > input/jas_samples_$2.fofn
	wc -l input/jas_samples_$2.fofn

	jasmine file_list=input/jas_samples_$2.fofn out_file=output74/jasmine_samples_$2.vcf --output_genotypes
	sortandindex output419/jasmine_samples_$2.vcf
	#wc -l output419/jasmine_samples_$2.vcf
}

compress_coverage() {
	for i in {5x,10x,15x,20x,30x}
	do
		run_jas 3 $i
	done
	echo 'finish merge'
	rm -r cmp_5x/ cmp_10x/ cmp_15x/ cmp_20x/ cmp_30x/
	truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/jasmine_samples_5x.vcf.gz -o cmp_5x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/jasmine_samples_10x.vcf.gz -o cmp_10x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/jasmine_samples_15x.vcf.gz -o cmp_15x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/jasmine_samples_20x.vcf.gz -o cmp_20x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/jasmine_samples_30x.vcf.gz -o cmp_30x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt answer.txt
	echo 'write to answer'
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

# output svtype benchmark of $1(*.vcf)
compress_type() {
	#compress_type coverage
	run_jas 3 $1
	echo 'finish merge'
	parse_svtype output419/jasmine_samples_$1.vcf
	rm -r cmp_ins/ cmp_del/ cmp_inv/ cmp_dup/ cmp_all/
	truvari bench -b benchmark/nstd_INS.vcf.gz -c INS.vcf.gz -o cmp_ins -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_DEL.vcf.gz -c DEL.vcf.gz -o cmp_del -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_INV.vcf.gz -c INV.vcf.gz -o cmp_inv -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_DUP.vcf.gz -c DUP.vcf.gz -o cmp_dup -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/jasmine_samples_$1.vcf.gz -o cmp_all -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_ins/summary.txt cmp_del/summary.txt cmp_inv/summary.txt cmp_dup/summary.txt cmp_all/summary.txt answer.txt
	echo 'write to answer'
}

#compress_coverage
#compress_type 30x

run_jas 3 30x
#rm cmp_jas
#truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/jasmine_samples_30x.vcf.gz -o cmp_jas -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch