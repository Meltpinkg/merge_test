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


def create_index(line, gz_filename, work_dir):
    cmd = 'bgzip -c ' + line + ' > ' + gz_filename
    try:
        os.system(cmd)
        os.system('tabix -p vcf ' + gz_filename)
    except Exception as ee:
        print(ee)

def pre_vcf(filenames, threads, work_dir):
    process_pool = Pool(processes = threads)
    vcf_filenames = []
    vcfgz_filenames = []
    start_time = time.time()
    with open(filenames, 'r') as f:
        for line in f:
            line = line.strip()
            vcf_filenames.append(line)
            filename = line.split('/')[-1]
            if filename[-2: 0] != 'gz':
                filename += '.gz'
                gz_filename = work_dir + line.split('/')[-1] + '.gz'
                vcfgz_filenames.append(gz_filename)
                process_pool.apply_async(create_index, (line, gz_filename, work_dir))
            #os.system('tabix -p vcf ' + line)
            else:
                vcfgz_filenames.append(line)
    print(vcfgz_filenames)
    process_pool.close()
    process_pool.join()
    #vcf_filenames = ['cuteSV1.vcf', 'cuteSV2.vcf', 'cuteSV3.vcf']
    #vcfgz_filenames = ['temp/cuteSV1.vcf.gz', 'temp/cuteSV2.vcf.gz', 'temp/cuteSV3.vcf.gz']
    chrom_set = set()
    contiginfo = set()
    for vcf_filename in vcfgz_filenames:  # 默认header的config中存放了所有chrom信息
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
    #output_chrom(list_head.next, len(vcf_filenames), chrom)

    result = list()
    head = list_head.next
    while head != None:
        candidate_idx, cipos, ciend = cal_center(head.variant_dict)
        candidate_record = head.variant_dict[candidate_idx]  # Record
        result.append([candidate_record.chrom1, candidate_record.start, candidate_record, cipos, ciend, head.variant_dict])
        head = head.next
    print('finish %s in %s, total %dnodes'%(chrom, str(time.time() - start_time), len(result)))
    return result  # [[CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT], [], ...]


'''
    var_list: [[a, b, info], ...]
    return [info] for which a <= breakpoint < b
'''
def div_in(var_list, breakpoint):
    left = 0
    right = len(var_list) - 1
    mid = 0
    info = []
    while left < right:
        mid = (left + right + 1) >> 1
        if var_list[mid][0] <= breakpoint:
            left = mid
        else:
            right = mid - 1
    for i in range(left, -1, -1):
        if breakpoint < var_list[i][1]:
            info.append(var_list[i])
        if 
    return info

'''
    var_list: [[a, b, info], ...]
    return [info] for which [start, end] and [a, b] overlap
'''
def div_overlap(var_list, start, end):
    left = 0
    right = len(var_list) - 1
    mid = 0
    info = []
    while left < right:
        mid = (left + right + 1) >> 1
        if var_list[mid][0] <= start:
            left = mid
        else:
            right = mid - 1
    for i in range(left, -1, -1):
        if start < var_list[i][1]:
            info.append(var_list[i])

    for i in range(left + 1, len(var_list), 1):
        if var_list[i][0] < end:
            info.append(var_list[i])
    return info


def add_anotation(var_list, sv_type, start, end):
    if len(var_list) == 0:
        return []
    left = 0
    right = len(var_list) - 1
    mid = 0
    while left < right:
        mid = (left + right + 1) >> 1
        if var_list[mid][0] <= start:
            left = mid
        else:
            right = mid - 1
    for i in range(left, -1, -1):
        if start < var_list[i][0]: 

    if right > 0 and pos - var_list[right - 1][1] <= 1000:
        for i in range(right - 1, -1, -1):
            if check_same_variant(var_type, var_list[i][2], sv_end):
                read_id_list.add(var_list[i][3])
                search_start = var_list[i][1]
            if i > 0 and var_list[i][1] - var_list[i - 1][1] > bias:
                break
    if var_list[right][1] - pos <= 1000:
        for i in range(right, len(var_list)):
            if check_same_variant(var_type, var_list[i][2], sv_end):  # if abs(var_list[i][2] - sv_end) < 1000:
                read_id_list.add(var_list[i][3])
                search_end = var_list[i][1]
            if i < len(var_list) - 1 and var_list[i + 1][1] - var_list[i][1] > bias:
                break
    return []


def main(argv):
    start_time = time.time()
    max_dist = 1000
    max_inspro = 0.7
    file_info = argv[0]
    file_output = argv[1]
    threads = argv[2]
    work_dir = argv[3]
    if work_dir[-1] != '/':
        work_dir += '/'
    if not os.path.exists(work_dir + 'index'):
        os.mkdir(work_dir + 'index')
    vcf_filenames, vcfgz_filenames, chrom_set, chrom_cnt, contiginfo = pre_vcf(file_info, threads, work_dir + 'index/')
    '''
    for chrom in chrom_set:
        ll_solve_chrom(vcfgz_filenames, chrom, 1000, 0.7)
        print(chrom + '\t' + str(time.time() - start_time))
    '''
    pool = Pool(processes = int(threads))
    result = list()
    for iter in chrom_cnt:
        # multi processes
        # print(iter[0])
        #if iter[1] == 0:
        #    continue
        para = [vcfgz_filenames, iter[0], max_dist, max_inspro]
        result.append(pool.map_async(ll_solve_chrom, para))
    pool.close()
    pool.join()
    print('finish merging', end='')
    print(time.time() - start_time)

    semi_result = list()
    for res in result:
        try:
            semi_result += res.get()[0]
        except:
            pass
    print('semi_result length=%d'%(len(semi_result)))
    semi_result = sorted(semi_result, key = lambda x:(x[0], int(x[1])))

    output_result(semi_result, samples_id, file_output, contiginfo)

    print('finish merge in ' + str(time.time() - start_time) + 'seconds')

    #os.system('rm -r temp')
    #'''
    

if __name__ == '__main__':
    #main(sys.argv[1:])  # filenames.txt, threads
    if len(sys.argv) == 2:
        test()
    else:
        main(['test.info', 'sample_merged25.vcf', 16, './'])
    #test()


'''
chrom_set = set()
base_cmd = 'cat '
for vcf_filename in vcf_filenames:
    base_cmd += vcf_filename + ' '
base_cmd += '| grep -v \'#\' | awk -F \'\\t\' \'{print $1}\' | uniq | sort | uniq'
print(base_cmd)
'''