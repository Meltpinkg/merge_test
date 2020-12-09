import numpy as np

'''
    node -> list(Record)
    return (center)start, end
'''
def cal_center(node):
    #return len(node) / 2
    r_start = 0
    r_end = 0
    start_list = []
    end_list = []
    for record in node:
        r_start += record.start
        r_end += record.end
        start_list.append(record.start)
        end_list.append(record.end)
    r_start = r_start / len(node)
    r_end = r_end / len(node)
    r_idx = 0
    start_dif = 0x3f3f3f3f
    for idx in range(len(node)):
        if abs(node[idx].start - r_start) < start_dif:
            start_dif = node[idx].start - r_start
            r_idx = idx
    return r_idx, cal_ci(start_list), cal_ci(end_list)


def cal_ci(input_list):
    pos = int(1.96 * np.std(input_list) / len(input_list) ** 0.5)
    return "-%d,%d"%(pos, pos)


def pre_vcf(filenames):
    #'''
    vcf_filenames = []
    vcfgz_filenames = []
    os.system('rm -r temp')
    os.mkdir('temp')
    start_time = time.time()
    with open(filenames, 'r') as f:
        for line in f:
            line = line.strip()
            vcf_filenames.append(line)
            if line[-2: 0] != 'gz':
                cmd = 'bgzip -c ' + line + ' > temp/' + line + '.gz'
                os.system(cmd)
            line = 'temp/' + line + '.gz'
            os.system('tabix -p vcf ' + line)
            vcfgz_filenames.append(line)
    #'''
    #vcf_filenames = ['cuteSV1.vcf', 'cuteSV2.vcf', 'cuteSV3.vcf']
    #vcfgz_filenames = ['temp/cuteSV1.vcf.gz', 'temp/cuteSV2.vcf.gz', 'temp/cuteSV3.vcf.gz']
    chrom_set = set()
    contiginfo = set()
    for vcf_filename in vcf_filenames:  # 默认header的config中存放了所有chrom信息
        vcf_reader = VariantFile(vcf_filename, 'r')
        for contig in vcf_reader.header.contigs:
            chrom_set.add(contig)  # contig 染色体号
    chrom_cnt = []  # 1:2989, 2:3009, 10:1956, X:1094
    '''
    base_cmd1 = 'grep -v \'#\' ' + vcf_filenames[0] + ' | awk -F \'\\t\' \'{print $1}\' | grep -x \''
    base_cmd2 = '\' | wc -l'
    for chrom in chrom_set:
        fd = os.popen(base_cmd1 + chrom + base_cmd2)
        mi = fd.read().strip()
        chrom_cnt.append([chrom, int(mi)])
    chrom_cnt.sort(key = lambda x:x[1], reverse = True)
    '''
    for chrom in chrom_set:
        chrom_cnt.append([chrom])
    #print(chrom_cnt)
    return vcf_filenames, vcfgz_filenames, chrom_set, chrom_cnt, contiginfo