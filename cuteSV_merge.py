from cuteSV_AVLTree import TreeNode, AVLTree, Record
from cuteSV_calculation import cal_center, cal_ci
from cuteSV_linkedList import add_node, print_list, ListNode
from pysam import VariantFile
from multiprocessing import Pool, Manager
import os
import vcf
import sys
import time


def test():
    vcf_reader = VariantFile('temp/cuteSV1.vcf.gz', 'r')
    for record in vcf_reader.fetch('1'):
        r = Record(record, '1')
        print(record)
        print(r.start)
        print(r.end)
        

def solve_chrom(vcf_filenames, chrom):
    tree = AVLTree()
    #print('start solving chrom: ' + chrom)
    time0 = time.time()
    for i in range(len(vcf_filenames)):
        idx = 0
        vcf_reader = VariantFile(vcf_filenames[i], 'r')
        for record in vcf_reader.fetch(chrom):
            idx += 1
            tree.insert(i, Record(record, i))
            '''
            print(idx)
            print(Record(record, i).to_string())
            #if record.start == 142535438 and record.info['SVTYPE'] == 'DUP':
            if idx == 1560:
                print(record)
                node_list = []
                vis_list = []
                tree.inorder(tree.root, node_list, vis_list)
                for node in node_list:
                    for r in node:
                        print(r.to_string(), end='')
                        print('\t')
                    print('\n')
                break
            #'''
        #print(idx)
    #print('start output chrom: ' + chrom)
    time1 = time.time()
    output_chrom(tree, len(vcf_filenames), chrom)
    #print('finish solving chrom: ' + chrom)
    time2 = time.time()
    #print(str(time1 - time0) + ' ' + str(time2 - time1))


'''
##FORMAT=<ID=ID,Number=1,Type=String,Description="Variant ID from input.">
##FORMAT=<ID=RAL,Number=1,Type=String,Description="Reference allele sequence reported from input.">
##FORMAT=<ID=AAL,Number=1,Type=String,Description="Alternative allele sequence reported from input.">
##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinates">
'''
def output_chrom(head, input_cnt, chrom):
    file = open('chrom_ans/' + chrom + '_ans.vcf', 'w')
    #print('start writing to file ' + chrom)
    while head != None:
        can_idx, cipos, ciend = cal_center(head.variant_list)
        can_record = head.variant_list[can_idx]
        supp_id = ''
        for i in range(input_cnt):
            if i in head.vis:
                supp_id += '1'
            else:
                supp_id += '0'
        if can_record.qual == "." or can_record.qual is None:
            filter_lable = "PASS"
        else:
            filter_lable = "PASS" if float(can_record.qual) >= 5.0 else "q5"
        if can_record.type == 'INS':
            sv_len = can_record.end
            sv_end = can_record.start - 1
            info_list = "SUPP={SUPP};SUPP_ID={SUPP_ID};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN}".format(
                    SUPP = len(head.variant_list),
                    SUPP_ID = supp_id,
                    SVTYPE = can_record.type, 
                    SVLEN = sv_len, 
                    END = sv_end, 
                    CIPOS = cipos,
                    CILEN = ciend)
            if can_record.strand != '.':
                info_list.append(';STRAND=' + can_record.strand)
        elif can_record.type == 'DEL':
            sv_len = can_record.start - can_record.end - 1
            sv_end = can_record.end
            info_list = "SUPP={SUPP};SUPP_ID={SUPP_ID};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CIEND={CIEND}".format(
                    SUPP = len(head.variant_list),
                    SUPP_ID = supp_id,
                    SVTYPE = can_record.type, 
                    SVLEN = sv_len, 
                    END = sv_end, 
                    CIPOS = cipos, 
                    CIEND = ciend)
            if can_record.strand != '.':
                info_list.append(';STRAND=' + can_record.strand)
        elif can_record.type == 'INV' or can_record.type == 'DUP':
            sv_len = can_record.end - can_record.start + 1
            sv_end = can_record.end
            info_list = "SUPP={SUPP};SUPP_ID={SUPP_ID};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END}".format(
                    SUPP = len(head.variant_list),
                    SUPP_ID = supp_id,
                    SVTYPE = can_record.type, 
                    SVLEN = sv_len, 
                    END = sv_end, 
                    CIPOS = cipos, 
                    CIEND = ciend)
            if can_record.strand != '.':
                info_list.append(';STRAND=' + can_record.strand)
        elif can_record.type == 'BND':
            info_list = "SUPP={SUPP};SUPP_ID={SUPP_ID};SVTYPE={SVTYPE};CHR2={CHR2};END={END}".format(
                    SUPP = len(head.variant_list),
                    SUPP_ID = supp_id,
                    SVTYPE = can_record.type, 
                    CHR2 = can_record.chrom2, 
                    END = can_record.end)
            
        file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\n".format(
            CHR = can_record.chrom1, 
            POS = can_record.start,
            ID = can_record.name,
            REF = can_record.ref,
            ALT = can_record.alt[0], 
            QUAL = '.' if can_record.qual is None else can_record.qual,
            PASS = filter_lable,
            INFO = info_list, 
            FORMAT = "ID:RAL:AAL",
            ))
        head = head.next
    file.close()


def generate_header(filename, contiginfo):
    file = open(filename, 'w')
    # General header
    file.write('##fileformat=VCFv4.2\n')
    file.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
	#for i in contiginfo:
	#	file.write("##contig=<ID=%s,length=%d>\n"%(i[0], i[1]))

	# Specific header
	# ALT
    file.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n")
    file.write("##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n")
    file.write("##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n")
    file.write("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n")
    file.write("##ALT=<ID=BND,Description=\"Breakend of translocation\">\n")

	# INFO
    file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n")
    file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n")
    file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
    file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
    file.write("##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n")
    # file.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
    file.write("##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n")
    file.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
    file.write("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n")
    file.write("##FILTER=<ID=q5,Description=\"Quality below 5\">\n")
    # FORMAT
    # file.write("\n")
    file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    file.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
    file.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
    file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
    file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")
    #file.write("##CommandLine=\"cuteSV %s\"\n"%(" ".join(argv)))


def pre_vcf(filenames):
    '''
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
    '''
    vcf_filenames = ['cuteSV1.vcf', 'cuteSV2.vcf', 'cuteSV3.vcf']
    vcfgz_filenames = ['temp/cuteSV1.vcf.gz', 'temp/cuteSV2.vcf.gz', 'temp/cuteSV3.vcf.gz']
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
    print(chrom_cnt)
    return vcf_filenames, vcfgz_filenames, chrom_set, chrom_cnt, contiginfo


def ll_solve_chrom(vcf_filenames, chrom, max_dist, max_inspro):
    #print('start ' + chrom)
    start_time = time.time()
    list_head = ListNode(-1, None)
    cur_node = list_head
    for i in range(len(vcf_filenames)):
        idx = 0
        vcf_reader = VariantFile(vcf_filenames[i], 'r')
        for record in vcf_reader.fetch(chrom):
            #print(record,end='')
            idx += 1
            cur_node = add_node(cur_node, i, Record(record, i), max_dist, max_inspro)
            #print(cur_node.to_string())
            #print()
    #print(cur_node.to_string())
    #print_list(list_head.next)
    output_chrom(list_head.next, len(vcf_filenames), chrom)
    print('finish ' + chrom + ' in ' + str(time.time() - start_time))


def main(argv):
    start_time = time.time()
    os.system('rm -r chrom_ans')
    os.mkdir('chrom_ans')
    max_dist = 1000
    max_inspro = 0.7

    vcf_filenames, vcfgz_filenames, chrom_set, chrom_cnt, contiginfo = pre_vcf(argv[0])
    '''
    for chrom in chrom_set:
        ll_solve_chrom(vcfgz_filenames, chrom, 1000, 0.7)
        print(chrom + '\t' + str(time.time() - start_time))
    '''
    pool = Pool(processes = int(argv[2]))
    for iter in chrom_cnt:
        # multi processes
        # print(iter[0])
        #if iter[1] == 0:
        #    continue
        pool.apply_async(ll_solve_chrom, args=(vcfgz_filenames, iter[0], max_dist, max_inspro)) 
    pool.close()
    pool.join()
    print('xx')
    print(time.time() - start_time)
    chrom_cnt.sort(key = lambda x:x[0], reverse = False)
    generate_header(argv[1], contiginfo)
    for i in chrom_cnt:
        #if i[1] == 0:
        #    break
        os.system('cat chrom_ans/' + i[0] + '_ans.vcf >> ' + argv[1])
    print('finish merge in ' + str(time.time() - start_time) + 'seconds')

    #os.system('rm -r temp')
    #'''
    

if __name__ == '__main__':
    #main(sys.argv[1:])  # filenames.txt, threads
    if len(sys.argv) == 2:
        test()
    else:
        main(['filenames.txt', 'sample_merged.vcf', 16])
    #test()


'''
chrom_set = set()
base_cmd = 'cat '
for vcf_filename in vcf_filenames:
    base_cmd += vcf_filename + ' '
base_cmd += '| grep -v \'#\' | awk -F \'\\t\' \'{print $1}\' | uniq | sort | uniq'
print(base_cmd)
'''