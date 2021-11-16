simulated() {

}

real() {
    for i in {10,20,50,100,200,500,1000,2504}
    do
        jasmine file_list=real_test/real_$i.fofn out_file=benchmark/jasmine_vcf/real_$i.vcf min_support=1 min_dist=1000 --output_genotypes threads=
    done
}

#simulated
real