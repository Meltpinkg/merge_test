import os
import sys
from cuteSV_linkedList import parse_svtype
# 1	10469	nssv14474157	C	<DUP>	.	.	DBVARID;SVTYPE=DUP;CIPOS=-100,100;CIEND=-100,100;IMPRECISE;END=23975;SVLEN=13507;EXPERIMENT=1;SAMPLE=NA12878;REGIONID=nsv3320972
def parse_nstd(file):
    file = open(file + '.tmp', 'w')
    sample_ids = ['CHM1', 'CHM13', 'HG00268', 'HG00514', 'HG00733', 'HG01352', 'HG02059', 'HG02106', 'HG02818', 'HG04217', 'HX1', 'NA12878', 'NA19240', 'NA19434']
    for xx in sample_ids:
        file.write(xx + '\t')
    file.write('\n')
    sample_dict = dict() # sample_dict[sv_type][start][end/len][sample] = gt
    svlen_dict = dict()
    for sv_type in ['INS', 'DEL', 'INV', 'DUP']:
        sample_dict[sv_type] = dict()
        svlen_dict[sv_type] = dict()
    with open('nstd_withgt.vcf', 'r') as f:
        idx = 0
        for line in f:
            lineseq = line.strip().split('\t')
            chrom1 = lineseq[0]
            start = int(lineseq[1])
            sv_type = lineseq[7].split('SVTYPE=')[1].split(';')[0]
            if 'DUP' in sv_type:
                sv_type = 'DUP'
            idx += 1
            end = abs(int(lineseq[7].split(';END=')[1].split(';')[0]))
            svlen = abs(int(lineseq[7].split(';SVLEN=')[1].split(';')[0]))
            name = 'variant' + str(idx)
            sample = lineseq[7].split('SAMPLE=')[1].split(';')[0]
            gt = lineseq[7].split('GT=')[1].split(';')[0]
            if start not in sample_dict[sv_type]:
                sample_dict[sv_type][start] = dict()
                svlen_dict[sv_type][start] = dict()
            if svlen not in sample_dict[sv_type][start]:
                sample_dict[sv_type][start][svlen] = dict()
                svlen_dict[sv_type][start][svlen] = end
            sample_dict[sv_type][start][svlen][sample] = gt
    
    for sv_type in sample_dict:
        for start in sample_dict[sv_type]:
            for svlen in sample_dict[sv_type][start]:
                af = 0
                for xx in sample_ids:
                    if xx in sample_dict[sv_type][start][svlen]:
                        if sample_dict[sv_type][start][svlen][xx] == '0/1':
                            af += 1
                        elif sample_dict[sv_type][start][svlen][xx] == '1/1':
                            af += 2
                af /= 28
                file.write('{CHR}\t{POS}\tSVTYPE={SVT};END={END};SVLEN={SVLEN};AF={AF}'.format(
                    CHR = '1',
                    POS = start,
                    SVT = sv_type,
                    END = svlen_dict[sv_type][start][svlen],
                    SVLEN = svlen,
                    AF = round(af, 4)
                ))
                for xx in sample_ids:
                    if xx in sample_dict[sv_type][start][svlen]:
                        file.write('\t' + sample_dict[sv_type][start][svlen][xx])
                    else:
                        file.write('\t0/0')
                file.write('\n')
            
def add_gt_from_laser():
    laser = dict()
    for sample in ['CHM1', 'CHM13', 'HX1', 'HG00268', 'HG00514', 'HG00733', 'HG01352']:
        lst = list()
        with open('/data/0/sqcao/data/VISOR2/laser125.bed') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[4] == '50.0':
                    gt = '0/1'
                else:
                    gt = '1/1'
                lst.append([int(seq[1]), int(seq[2]), gt])
        laser[sample] = lst
    for sample in ['HG02059', 'HG02106', 'HG02818', 'HG04217', 'NA12878', 'NA19240', 'NA19434']:
        lst = list()
        with open('/data/0/sqcao/data/VISOR2/laser' + sample + '.bed') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[4] == '50.0':
                    gt = '0/1'
                else:
                    gt = '1/1'
                lst.append([int(seq[1]), int(seq[2]), gt])
        laser[sample] = lst
    return laser

def add_gt_to_nstd(laser_dict):
    file = open('nstd_withgt.vcf', 'w')
    with open('nstd162.GRCh37.variant_call.vcf', 'r') as f:
        idx = 0
        for line in f:
            if line[0] == '#':
                continue
            lineseq = line.strip().split('\t')
            chr = lineseq[0]
            start = int(lineseq[1])
            if chr != '1':
                break
            sample = lineseq[7].split('SAMPLE=')[1].split(';')[0]
            for item in laser_dict[sample]:
                if item[0] <= start < item[1]:
                    gt = item[2]
            file.write(line.strip() + ';GT=' + gt + '\n')

def parse_nstd_tru():
    file = open('nstd_merge.vcf', 'w')
    svs = set()
    with open('nstd162.GRCh37.variant_call.vcf', 'r') as f:
        idx = 0
        for line in f:
            if line[0] == '#':
                if line[1] == '#':
                    file.write(line)
                else:
                    file.write(line.strip() + '\tFORMAT\tNULL\n')
                continue
            lineseq = line.strip().split('\t')
            chrom1 = lineseq[0]
            if chrom1 != '1':
                break
            sv_type = parse_svtype(lineseq[7].split('SVTYPE=')[1].split(';')[0])
            pos = int(lineseq[1])
            svlen = abs(int(lineseq[7].split(';SVLEN=')[1].split(';')[0]))
            svend = abs(int(lineseq[7].split(';END=')[1].split(';')[0]))
            if idx == 0:
                idx = 1
                file.write(line.strip() + '\tGT\t./.\n')
                svs.add((sv_type, pos, svlen, svend))
                continue
            if (sv_type, pos, svlen, svend) in svs:
                pass
            else:
                file.write(line.strip() + '\tGT\t./.\n')
                svs.add((sv_type, pos, svlen, svend))
    file.close()

if __name__ == '__main__':
    #laser_dict = add_gt_from_laser()
    #add_gt_to_nstd(laser_dict)
    #parse_nstd(sys.argv[1])
    #os.system('sort -k 2,2n ' + sys.argv[1] + '.tmp > ' + sys.argv[1])
    parse_nstd_tru()

