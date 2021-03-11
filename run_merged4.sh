for i in {cutesv,pbsv,sniffles,svim,cutesv0}
do
    /usr/bin/time -v python src/cuteSV_merge.py input/$i.fofn output/merged1_$i.vcf ./ -t 16
done
