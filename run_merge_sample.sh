run() {
	# run support_num coverage folder
	# run 3 30x 831
	for i in {CHM1,CHM13,HG00268,HG00514,HG00733,HG01352,HG02059,HG02106,HG02818,HG04217,HX1,NA12878,NA19240,NA19434}
	do
		python first_release/cuteSV_merge.py input/$2/$i.fofn $3/$i\_$2.vcf ./ -t 1 -S $1
	done
	ls $3/*_$2.vcf > $3/$2.fofn
	wc -l $3/$2.fofn
	python first_release/cuteSV_merge.py $3/$2.fofn $3/$2.vcf ./ -t 1

	bgzip -c $3/$2.vcf > $3/$2.vcf.gz
	tabix $3/$2.vcf.gz
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
parse_nstd() {
	for svtype in {INS,DEL,INV,DUP}
	do
		grep '#' benchmark/nstd_merge.vcf > head
		grep SVTYPE=$svtype benchmark/nstd_merge.vcf > body
		cat head body > benchmark/nstd_$svtype.vcf
		rm head body
		bgzip -c benchmark/nstd_$svtype.vcf > benchmark/nstd_$svtype.vcf.gz
		tabix benchmark/nstd_$svtype.vcf.gz
	done
}

# output benchmark of $1(*.vcf)
compress_type() {
	#compress_type folder
	echo 'COMPRESS TYPE'
	run_src 3 30x $1
	parse_svtype $1/merged_samples_30x.vcf
	rm -r cmp_ins/ cmp_del/ cmp_inv/ cmp_dup/ cmp_all/
	nohup truvari bench -b benchmark/nstd_INS.vcf.gz -c INS.vcf.gz -o cmp_ins -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_DEL.vcf.gz -c DEL.vcf.gz -o cmp_del -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_INV.vcf.gz -c INV.vcf.gz -o cmp_inv -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_DUP.vcf.gz -c DUP.vcf.gz -o cmp_dup -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_all -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_ins/summary.txt cmp_del/summary.txt cmp_inv/summary.txt cmp_dup/summary.txt cmp_all/summary.txt answer.txt
	#echo 'write to answer'
}
compress_coverage() {
	# compress_coverage folder
	echo 'COMPRESS COVERAGE'
	for i in {5x,10x,15x,20x,30x}
	do
		run_src 3 $i $1
	done
	#echo 'finish merge'
	rm -r cmp_5x/ cmp_10x/ cmp_15x/ cmp_20x/ cmp_30x/
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_5x.vcf.gz -o cmp_5x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_10x.vcf.gz -o cmp_10x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_15x.vcf.gz -o cmp_15x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_20x.vcf.gz -o cmp_20x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_30x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt answer.txt
	#echo 'write to answer'
}
compress_svlength() {
	# compress_svlength folder
	echo 'COMPRESS SVLENGTH'
	run_src 3 30x $1
	#echo 'finish merge'
	rm -r cmp_1/ cmp_2/ cmp_3/ cmp_4/ cmp_5/ cmp_6/ cmp_7/
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_1 -p 0 -r 1000 -s 30 --sizemax 99 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_2 -p 0 -r 1000 -s 100 --sizemax 499 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_3 -p 0 -r 1000 -s 500 --sizemax 999 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_4 -p 0 -r 1000 -s 1000 --sizemax 4999 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_5 -p 0 -r 1000 -s 5000 --sizemax 9999 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_6 -p 0 -r 1000 -s 10000 --sizemax 10000000 --multimatch
	nohup truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/merged_samples_30x.vcf.gz -o cmp_7 -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_1/summary.txt cmp_2/summary.txt cmp_3/summary.txt cmp_4/summary.txt cmp_5/summary.txt cmp_6/summary.txt cmp_7/summary.txt answer.txt
	#echo 'write to answer'
}

main() {
	# compress_coverage folder
	rm -r $1
	mkdir $1
	for i in {5x,10x,15x,20x,30x}
	do
		run 3 $i $1
	done
	echo 'finish merge'

	# coverage
	rm -r cmp_5x/ cmp_10x/ cmp_15x/ cmp_20x/ cmp_30x/
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/5x.vcf.gz -o cmp_5x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/10x.vcf.gz -o cmp_10x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/15x.vcf.gz -o cmp_15x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/20x.vcf.gz -o cmp_20x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/30x.vcf.gz -o cmp_30x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	# type
	parse_svtype $1/30x.vcf
	rm -r cmp_ins/ cmp_del/ cmp_inv/ cmp_dup/
	truvari bench -b benchmark/nstd_INS.vcf.gz -c INS.vcf.gz -o cmp_ins -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_DEL.vcf.gz -c DEL.vcf.gz -o cmp_del -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_INV.vcf.gz -c INV.vcf.gz -o cmp_inv -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_DUP.vcf.gz -c DUP.vcf.gz -o cmp_dup -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	#truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_all -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	# length
	rm -r cmp_1/ cmp_2/ cmp_3/ cmp_4/ cmp_5/ cmp_6/
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/30x.vcf.gz -o cmp_1 -p 0 -r 1000 -s 30 --sizemax 99 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/30x.vcf.gz -o cmp_2 -p 0 -r 1000 -s 100 --sizemax 499 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/30x.vcf.gz -o cmp_3 -p 0 -r 1000 -s 500 --sizemax 999 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/30x.vcf.gz -o cmp_4 -p 0 -r 1000 -s 1000 --sizemax 4999 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/30x.vcf.gz -o cmp_5 -p 0 -r 1000 -s 5000 --sizemax 9999 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/30x.vcf.gz -o cmp_6 -p 0 -r 1000 -s 10000 --sizemax 10000000 --multimatch
	#truvari bench -b benchmark/nstd_merge.vcf.gz -c benchmark/survivor_vcf/sur_30x.vcf.gz -o cmp_7 -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt \
					cmp_ins/summary.txt cmp_del/summary.txt cmp_inv/summary.txt cmp_dup/summary.txt \
					cmp_1/summary.txt cmp_2/summary.txt cmp_3/summary.txt cmp_4/summary.txt cmp_5/summary.txt cmp_6/summary.txt answer.txt
	
}

main $1