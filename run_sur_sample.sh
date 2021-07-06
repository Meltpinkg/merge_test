:<<!
for i in {CHM1,CHM13,HG00268,HG00514,HG00733,HG01352,HG02059,HG02106,HG02818,HG04217,HX1,NA12878,NA19240,NA19434}
do
	/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/$i.fofn 1000 2 1 0 1 30 output_sample/survivor_s2_$i.vcf
	/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/$i.fofn 1000 3 1 0 1 30 output_sample/survivor_s3_$i.vcf
done
ls output_sample/survivor_s2_* > input/sur_s2_samples.fofn
ls output_sample/survivor_s3_* > input/sur_s3_samples.fofn
/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/sur_s2_samples.fofn 1000 1 1 0 1 30 output419/survivor_s2_sample.vcf
/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/sur_s3_samples.fofn 1000 1 1 0 1 30 output419/survivor_s3_sample.vcf

grep -v '#' output419/survivor_s2_sample.vcf | sort -k 1,1 -k 2,2n > body
grep '#' output419/survivor_s2_sample.vcf > head
cat head body > output419/survivor_s2_sample.vcf
grep -v '#' output419/survivor_s3_sample.vcf | sort -k 1,1 -k 2,2n > body
grep '#' output419/survivor_s3_sample.vcf > head
cat head body > output419/survivor_s3_sample.vcf
rm head body
rm -r benchmark/tru_sur_s2 benchmark/tru_sur_s3
echo 'support == 2'
#python src/eval.py benchmark/cmp_sur_s2 chr1_merge.vcf output419/survivor_s2_sample.vcf
bgzip -c output419/survivor_s2_sample.vcf > output419/survivor_s2_sample.vcf.gz
tabix output419/survivor_s2_sample.vcf.gz
truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/survivor_s2_sample.vcf.gz -o benchmark/tru_sur_s2 -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
python src/parse_truvari.py benchmark/tru_sur_s2
echo 'support == 3'
#python src/eval.py benchmark/cmp_sur_s2 chr1_merge.vcf output419/survivor_s3_sample.vcf
bgzip -c output419/survivor_s3_sample.vcf > output419/survivor_s3_sample.vcf.gz
tabix output419/survivor_s3_sample.vcf.gz
truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/survivor_s3_sample.vcf.gz -o benchmark/tru_sur_s3 -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
python src/parse_truvari.py benchmark/tru_sur_s3
echo 'run sur sample finished'
!
sortandindex() {
	grep -v '#' $1 | sort -k 1,1 -k 2,2n > body
	grep '#' $1 > head
	cat head body > $1
	rm head body
	bgzip -c $1 > $1.gz
	tabix $1.gz
}

run_sur() {
	# run_sur support coverage
	for i in {CHM1,CHM13,HG00268,HG00514,HG00733,HG01352,HG02059,HG02106,HG02818,HG04217,HX1,NA12878,NA19240,NA19434}
	do
		echo $i
		/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/$2/$i.fofn 1000 $1 1 0 1 30 output_sample/survivor_$i\_$2.vcf
	done
	ls output_sample/survivor_*_$2.vcf > input/sur_samples_$2.fofn
	wc -l input/sur_samples_$2.fofn

	/data/0/sqcao/tools/SURVIVOR/Debug/SURVIVOR merge input/sur_samples_$2.fofn 1000 1 1 0 1 30 output74/survivor_samples_$2.vcf
	sortandindex output74/survivor_samples_$2.vcf
}

compress_coverage() {
	# compress_coverage folder
:<<!
	for i in {5x,10x,15x,20x,30x}
	do
		run_sur 3 $i
	done
	echo 'finish merge'
!
	rm -r cmp_5x/ cmp_10x/ cmp_15x/ cmp_20x/ cmp_30x/
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/survivor_samples_5x.vcf.gz -o cmp_5x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/survivor_samples_10x.vcf.gz -o cmp_10x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/survivor_samples_15x.vcf.gz -o cmp_15x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/survivor_samples_20x.vcf.gz -o cmp_20x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	truvari bench -b benchmark/nstd_merge.vcf.gz -c $1/survivor_samples_30x.vcf.gz -o cmp_30x -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch
	python truvari/ss.py cmp_5x/summary.txt cmp_10x/summary.txt cmp_15x/summary.txt cmp_20x/summary.txt cmp_30x/summary.txt answer.txt
	echo 'write to answer'
}

compress_coverage output419
#run_sur 3 10x
#rm cmp_sur
#truvari bench -b benchmark/nstd_merge.vcf.gz -c output419/survivor_samples_10x.vcf.gz -o cmp_sur -p 0 -r 1000 -s 30 --sizemax 10000000 --multimatch