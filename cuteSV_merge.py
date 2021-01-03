from cuteSV_calculation import cal_center, cal_ci
from cuteSV_linkedList import add_node, print_list, ListNode, Record
from cuteSV_output import output_result
from pysam import VariantFile
from multiprocessing import Pool, Manager
import os
import sys
import time
import argparse


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
    #print(vcfgz_filenames)
    process_pool.close()
    process_pool.join()
    chrom_set = set()
    contiginfo = dict()
    sample_set = set()
    sample_ids = list()
    for vcf_filename in vcfgz_filenames:  # 默认header的config中存放了所有chrom信息
        vcf_reader = VariantFile(vcf_filename, 'r')
        for i in range(len(vcf_reader.header.contigs)):
            chrom_set.add(str(vcf_reader.header.contigs[i].name))
            contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
        sample_id = vcf_reader.header.samples[0]
        temp_sample_id = sample_id
        temp_idx = 0
        while temp_sample_id in sample_set:
            temp_sample_id = sample_id + '_' + str(temp_idx)
            temp_idx += 1
        sample_ids.append(temp_sample_id)
        sample_set.add(temp_sample_id)
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
    return vcfgz_filenames, chrom_set, chrom_cnt, contiginfo, sample_ids


#def ll_solve_chrom(vcf_filenames, chrom, max_dist, max_inspro):
def ll_solve_chrom(para):
    vcf_filenames = para[0]
    chrom = para[1]
    max_dist = para[2]
    max_inspro = para[3]
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
    #print('finish %s in %s, total %dnodes'%(chrom, str(time.time() - start_time), len(result)))
    return result  # [[CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT], [], ...]


def main(argv):
    start_time = time.time()
    max_dist = 1000
    max_inspro = 0.7
    if args.work_dir[-1] != '/':
        args.work_dir += '/'
    if not os.path.exists(args.work_dir + 'index'):
        os.mkdir(args.work_dir + 'index')
    vcfgz_filenames, chrom_set, chrom_cnt, contiginfo, sample_ids = pre_vcf(args.input, args.threads, args.work_dir + 'index/')
    pool = Pool(processes = args.threads)
    result = list()
    '''
    for iter in chrom_cnt: 
        result.append(ll_solve_chrom(vcfgz_filenames, iter[0], max_dist, max_inspro))
    '''
    for iter in chrom_cnt:
        # multi processes
        # print(iter[0])
        #if iter[1] == 0:
        #    continue
        para = [(vcfgz_filenames, iter[0], max_dist, max_inspro)]
        result.append(pool.map_async(ll_solve_chrom, para))
        #result.append(pool.map_async(call, para))
    pool.close()
    pool.join()
    print('finish merging in ', end='')
    print(time.time() - start_time)

    semi_result = list()
    for res in result:
        try:
            semi_result += res.get()[0]
            #semi_result += res
        except:
            pass
    print('semi_result length=%d'%(len(semi_result)))
    #print(semi_result)
    #exit(0)
    semi_result = sorted(semi_result, key = lambda x:(x[0], int(x[1])))

    output_result(semi_result, sample_ids, args.output, contiginfo, args.annotation)

    print('finish in ' + str(time.time() - start_time) + 'seconds')

    #os.system('rm -r temp')
    #'''
    

if __name__ == '__main__':
    #main(sys.argv[1:])  # filenames.txt, threads
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
            metavar = 'input_file',
            type = str,
            help = 'filenames of VCF files to be merged')
    parser.add_argument('output',
            type = str,
            help = 'output VCF format file')
    parser.add_argument('work_dir',
            type = str,
            help = 'work directory for temporary files')
    parser.add_argument('-t', '--threads',
            type = int,
            default = 16,
            help = 'number of threads to use[%(default)s]')
    parser.add_argument('-a', '--annotation',
            type = str,
            default = None,
            help = 'annotation file to add')
    args = parser.parse_args(sys.argv[1:])
    main(args)


'''
chrom_set = set()
base_cmd = 'cat '
for vcf_filename in vcf_filenames:
    base_cmd += vcf_filename + ' '
base_cmd += '| grep -v \'#\' | awk -F \'\\t\' \'{print $1}\' | uniq | sort | uniq'
print(base_cmd)
'''