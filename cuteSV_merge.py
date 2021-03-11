from cuteSV_calculation import cal_center, cal_ci
from cuteSV_linkedList import add_node, print_list, ListNode, Record
from cuteSV_output import output_result, solve_annotation, generate_header
from pysam import VariantFile
from multiprocessing import Pool, Manager
import os
import sys
import time
import argparse
        

def avl_solve_chrom(vcf_filenames, chrom):
    tree = AVLTree()
    #print('start solving chrom: ' + chrom)
    time0 = time.time()
    for i in range(len(vcf_filenames)):
        idx = 0
        vcf_reader = VariantFile(vcf_filenames[i], 'r')
        for record in vcf_reader.fetch(chrom):
            idx += 1
            tree.insert(i, Record(record, i))
    time1 = time.time()
    output_chrom(tree, len(vcf_filenames), chrom)
    time2 = time.time()


def create_index(line, gz_filename):
    try:
        cmd = 'bgzip -c ' + line + ' > ' + gz_filename
        os.system(cmd)
        os.system('tabix -p vcf ' + gz_filename)
    except Exception as ee:
        print(ee)

def pre_vcf(filenames, threads, work_dir):
    print('start indexing...')
    start_time = time.time()
    if args.work_dir[-1] != '/':
        args.work_dir += '/'
    if not os.path.exists(args.work_dir + 'index'):
        os.mkdir(args.work_dir + 'index')
    process_pool = Pool(processes = threads)
    vcf_filenames = []
    vcfgz_filenames = []
    #start_time = time.time()
    with open(filenames, 'r') as f:
        for line in f:
            line = line.strip()
            vcf_filenames.append(line)
            filename = line.split('/')[-1]
            if filename[-3:] != 'vcf':
                print('input error')
                continue
            filename += '.gz'
            gz_filename = work_dir + line.split('/')[-1] + '.gz'
            vcfgz_filenames.append(gz_filename)
            process_pool.apply_async(create_index, (line, gz_filename))
    process_pool.close()
    process_pool.join()
    '''
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
    chrom_cnt = []

    base_cmd1 = 'grep -v \'#\' ' + vcf_filenames[0] + ' | awk -F \'\\t\' \'{print $1}\' | grep -x \''
    base_cmd2 = '\' | wc -l'
    for chrom in chrom_set:
        fd = os.popen(base_cmd1 + chrom + base_cmd2)
        mi = fd.read().strip()
        chrom_cnt.append([chrom, int(mi)])
    chrom_cnt.sort(key = lambda x:x[1], reverse = True)
    for chrom in chrom_set:
        chrom_cnt.append([chrom])
    '''
    print('finish indexing in %s'%(round(time.time() - start_time, 4)))
    #return vcfgz_filenames, chrom_cnt, contiginfo, sample_ids
    return vcfgz_filenames


def parse_vcf(para):
    filename = para[0]
    idx = para[1]
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict()
    # contigs
    contiginfo = dict()
    for i in range(len(vcf_reader.header.contigs)):
        contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
    record_dict['contig'] = contiginfo
    # records
    for record in vcf_reader.fetch():
        chr = record.chrom
        if chr not in record_dict:
            record_dict[chr] = []
        record_dict[chr].append(Record(record, idx))
    # sample_id
    record_dict['sample'] = vcf_reader.header.samples[0]
    #print(filename)
    #print(record_dict.keys())
    return record_dict


def parse_vcfs(filenames, threads):
    print('start parse vcfs...')
    start_time = time.time()
    pool = Pool(processes = threads)
    file_dict = dict()
    idx = 0
    with open(filenames, 'r') as f:
        for filename in f:
            file_dict[idx] = pool.map_async(parse_vcf, [(filename.strip(), idx)])
            idx += 1
    pool.close()
    pool.join()
    #print(file_dict.keys())
    chr_dict = dict()
    sample_temp = ['' for i in range(idx)]
    contiginfo = dict()
    for fileidx in file_dict:
        temp = file_dict[fileidx].get()[0]  #{chr -> List(Record)}
        sample_temp[fileidx] = temp['sample']
        contig_temp = temp['contig']
        temp.pop('sample')
        temp.pop('contig')
        for chr in temp:
            if chr not in chr_dict:
                chr_dict[chr] = dict()
                if chr not in contiginfo:
                    contiginfo[chr] = contig_temp[chr]
            chr_dict[chr][fileidx] = temp[chr]

    for fileidx in file_dict:
        temp = file_dict[fileidx].get()[0]  #{chr -> List(Record)}
        for chr in temp:
            if chr not in chr_dict:
                chr_dict[chr] = dict()
            chr_dict[chr][fileidx] = temp[chr]

    #print(contiginfo)
    print(chr_dict.keys())
    sample_set = set()
    sample_ids = list()
    for sample_id in sample_temp:
        temp_sample_id = sample_id
        temp_idx = 0
        while temp_sample_id in sample_set:
            temp_sample_id = sample_id + '_' + str(temp_idx)
            temp_idx += 1
        sample_ids.append(temp_sample_id)
        sample_set.add(temp_sample_id)
    #print(sample_ids)
    print('finish parsing in %s'%(str(time.time() - start_time)))
    return chr_dict, sample_ids, contiginfo


# all variants on the chrom
def resolve_chrom(filenames, output, chrom, sample_ids, threads, max_dist, max_inspro, seperate_svtype, anno):
    #print('start resolving chrom ' + chrom)
    start_time = time.time()
    pool = Pool(processes = threads)
    record_dict = dict()
    idx = 0
    for filename in filenames:
        record_dict[idx] = pool.map_async(parse_vcf_chrom, [(filename, chrom, idx)])
        idx += 1
    pool.close()
    pool.join()
    if seperate_svtype:
        result = solve_chrom2(record_dict, max_dist, max_inspro, anno)
    else:
        result = solve_chrom1(record_dict, max_dist, max_inspro, anno)
    result = sorted(result, key = lambda x:(x[0], int(x[1])))
    output_result(result, sample_ids, output)
    #print('finish resolving %s, %d in %s'%(chrom, len(result), round(time.time() - start_time, 4)))

def parse_vcf_chrom(para):
    filename = para[0]
    chrom = para[1]
    idx = para[2]
    vcf_reader = VariantFile(filename, 'r')
    record_list = list()
    # records
    for record in vcf_reader.fetch(chrom):
        record_list.append(Record(record, idx))
    return record_list

def solve_chrom1(record_dict, max_dist, max_inspro, anno):
    list_head = ListNode(-1, None)
    cur_node = list_head
    for fileidx in record_dict:
        for record in record_dict[fileidx].get()[0]:
            cur_node = add_node(cur_node, fileidx, record, max_dist, max_inspro)
    result = list()
    head = list_head.next
    while head != None:
        candidate_idx, cipos, ciend = cal_center(head.variant_dict)
        candidate_record = head.variant_dict[candidate_idx]  # Record
        annotation = solve_annotation(candidate_record.type, anno, candidate_record.start, candidate_record.end)  # dict{'gene_id' -> str}
        result.append([candidate_record.chrom1, candidate_record.start, candidate_record, cipos, ciend, head.variant_dict, annotation])
        head = head.next
    return result  # [[CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT, ANNOTATION], [], ...]

def solve_chrom2(record_dict, max_dist, max_inspro, anno):
    list_head = dict()
    cur_node = dict()
    #list_head = ListNode(-1, None)
    #cur_node = list_head
    for fileidx in record_dict:
        record_list = record_dict[fileidx].get()[0]
        for record in record_list: 
            sv_type = record.type
            if sv_type not in list_head:
                list_head[sv_type] = ListNode(-1, None)
                cur_node[sv_type] = list_head[sv_type]
            cur_node[sv_type] = add_node(cur_node[sv_type], fileidx, record, max_dist, max_inspro)
    result = list()
    for sv_type in list_head:
        head = list_head[sv_type].next
        while head != None:
            candidate_idx, cipos, ciend = cal_center(head.variant_dict)
            candidate_record = head.variant_dict[candidate_idx]  # Record
            annotation = solve_annotation(candidate_record.type, anno, candidate_record.start, candidate_record.end)  # dict{'gene_id' -> str}
            result.append([candidate_record.chrom1, candidate_record.start, candidate_record, cipos, ciend, head.variant_dict, annotation])
            head = head.next
    return result  # [[CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT, ANNOTATION], [], ...]

def resolve_contigs(filenames, threads):
    print('start resolving contigs...')
    start_time = time.time()
    sample_temp = list()
    sample_set = set()
    sample_ids = list()
    contiginfo = dict()
    result = list()
    pool = Pool(processes = threads)
    for filename in filenames:
        result.append(pool.map_async(parse_contigs, [(filename)]))
    pool.close()
    pool.join()
    for res in result:
        temp = res.get()[0]
        for contig in temp:
            if contig == 'sample':
                sample_id = temp['sample']
                temp_sample_id = sample_id
                temp_idx = 0
                while temp_sample_id in sample_set:
                    temp_sample_id = sample_id + '_' + str(temp_idx)
                    temp_idx += 1
                sample_ids.append(temp_sample_id)
                sample_set.add(temp_sample_id)
            else:
                if contig not in contiginfo:
                    contiginfo[contig] = temp[contig]
    print('finish resolving contigs in %s'%(round(time.time() - start_time, 4)))
    return sample_ids, contiginfo

def parse_contigs(para):
    filename = para
    contiginfo = dict()
    vcf_reader = VariantFile(filename, 'r')
    for i in range(len(vcf_reader.header.contigs)):
        contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
    contiginfo['sample'] = vcf_reader.header.samples[0]
    return contiginfo

def parse_annotation_file(annotation_file):
    if annotation_file == None:
        return None
    annotation_dict = dict()  # [chrom -> [[a, b, dict()], ...]]
    start_time = time.time()
    with open(annotation_file, 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            chr = seq[0]
            if chr[:3] == "chr":
                chr = chr[3:]
            if chr not in annotation_dict:
                annotation_dict[chr] = []
            info_seq = seq[8].strip().split(';')
            info_dict = dict()
            for info_item in info_seq[:-1]:
                info_item = info_item.strip().split(' ')
                if len(info_item) != 2:
                    print('error:length of gtf_info_item is not 2')
                    continue
                info_dict[info_item[0]] = info_item[1][1:-1]
            annotation_dict[chr].append([int(seq[3]), int(seq[4]), info_dict])    
    print(time.time() - start_time)
    return annotation_dict


def ll_solve_chrom(para):
    chrom = para[1]
    #print('start %s at %s'%(chrom, str(time.time())))
    #vcf_filenames = para[0]
    vcf_files = para[0]
    max_dist = para[2]
    max_inspro = para[3]
    anno = para[4]
    start_time = time.time()
    list_head = ListNode(-1, None)
    cur_node = list_head
    '''
    for i in range(len(vcf_filenames)):
        idx = 0
        vcf_reader = VariantFile(vcf_filenames[i], 'r')
        for record in vcf_reader.fetch(chrom):
            idx += 1
            cur_node = add_node(cur_node, i, Record(record, i), max_dist, max_inspro)
    '''
    for fileidx in vcf_files:
        for record in vcf_files[fileidx]:
            cur_node = add_node(cur_node, fileidx, record, max_dist, max_inspro)

    result = list()
    head = list_head.next
    while head != None:
        candidate_idx, cipos, ciend = cal_center(head.variant_dict)
        candidate_record = head.variant_dict[candidate_idx]  # Record
        annotation = solve_annotation(candidate_record.type, anno, candidate_record.start, candidate_record.end)  # dict{'gene_id' -> str}
        '''
        if candidate_record.type == 'TRA':
            annotation = solve_annotation_tra(annotation_dict[candidate_record.chrom2], candidate_record.end, annotation)
            print(annotation)
        '''
        result.append([candidate_record.chrom1, candidate_record.start, candidate_record, cipos, ciend, head.variant_dict, annotation])
        head = head.next
    #print('finish %s in %s, total %dnodes'%(chrom, str(time.time() - start_time), len(result)))
    #print('finish %s at %s'%(chrom, str(time.time())))
    return result  # [[CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT, ANNOTATION], [], ...]


def main_ctrl(args):
    start_time = time.time()
    max_dist = 1000
    max_inspro = 0.7
    filenames = pre_vcf(args.input, args.threads, args.work_dir + 'index/')
    annotation_dict = parse_annotation_file(args.annotation)
    sample_ids, contiginfo = resolve_contigs(filenames, args.IOthreads)
    print('%d samples find'%(len(sample_ids)))
    file = open(args.output, 'w')
    generate_header(file, contiginfo, sample_ids)
    file.close()
    for contig in contiginfo:
        if annotation_dict != None and contig in annotation_dict:
                anno = annotation_dict[contig]
        else:
            anno = []
        resolve_chrom(filenames, args.output, contig, sample_ids, args.IOthreads, max_dist, max_inspro, args.seperate_svtype, anno)
    print('finish merging in %s'%(round(time.time() - start_time, 4)))
    os.system('rm -r ' + args.work_dir + 'index/')

def main(args):
    start_time = time.time()
    max_dist = 1000
    max_inspro = 0.7
    chr_dict, sample_ids, contiginfo = parse_vcfs(args.input, args.IOthreads)
    '''
    if args.work_dir[-1] != '/':
        args.work_dir += '/'
    if not os.path.exists(args.work_dir + 'index'):
        os.mkdir(args.work_dir + 'index')
    vcfgz_filenames, chrom_cnt, contiginfo, sample_ids = pre_vcf(args.input, args.threads, args.work_dir + 'index/')
    '''
    annotation_dict = parse_annotation_file(args.annotation)
    pool = Pool(processes = args.threads)
    result = list()
    if args.threads == 1:
        for chr in chr_dict:
            if annotation_dict != None and chr in annotation_dict:
                anno = annotation_dict[chr]
            else:
                anno = []
            result.append(ll_solve_chrom([chr_dict[chr], chr, max_dist, max_inspro, anno]))
    else:
        para = []
        for chr in chr_dict:
            if annotation_dict != None and chr in annotation_dict:
                anno = annotation_dict[chr]
            else:
                anno = []
            para.append([chr_dict[chr], chr, max_dist, max_inspro, anno])
        result = pool.map(ll_solve_chrom, para)

    pool.close()
    pool.join()
    print('finish merging in ', end='')
    print(time.time() - start_time)
    semi_result = list()
    for res in result:
        try:
            semi_result += res
            # semi_result += res.get()[0]
        except:
            pass
    print('semi_result length=%d'%(len(semi_result)))
    semi_result = sorted(semi_result, key = lambda x:(x[0], int(x[1])))

    file = open(args.output, 'w')
    generate_header(file, contiginfo, sample_ids)
    file.close()
    output_result(semi_result, sample_ids, args.output)

    print('finish in ' + str(time.time() - start_time) + 'seconds')

    #os.system('rm -r ' + args.work_dir + 'index/')
    #'''
    

if __name__ == '__main__':
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
    parser.add_argument('-tio', '--IOthreads',
            type = int,
            default = 16,
            help = 'number of threads of IO[%(default)s]')
    parser.add_argument('-a', '--annotation',
            type = str,
            default = None,
            help = 'annotation file to add')
    parser.add_argument('--seperate_svtype',
            action="store_true",
            default = False)
    parser.add_argument('--seperate_chrom',
            action="store_true",
            default = False)
    args = parser.parse_args(sys.argv[1:])
    if args.seperate_chrom:
        main_ctrl(args)
    else:
        main(args)
    # /usr/bin/time -v python src/cuteSV_merge.py input/file_SRR.fofn sample_merged25.vcf ./ -t 16 -a hg19.refGene.gtf
    # /usr/bin/time -v python src/cuteSV_merge.py test.info test_merged.vcf ./ -t 16 -a hg19.refGene.gtf
    # /usr/bin/time -v python src/cuteSV_merge.py input/cutesv.fofn merged_cutesv.vcf ./ -t 16
    # /usr/bin/time -v python src/cuteSV_merge.py input/sniffles.fofn merged1_sniffles.vcf ./ -t 16
    # /usr/bin/time -v python src/cuteSV_merge.py input/pbsv.fofn output/merged_pbsv.vcf ./ -t 16

'''
chrom_set = set()
base_cmd = 'cat '
for vcf_filename in vcf_filenames:
    base_cmd += vcf_filename + ' '
base_cmd += '| grep -v \'#\' | awk -F \'\\t\' \'{print $1}\' | uniq | sort | uniq'
print(base_cmd)
'''
'''
for iter in chrom_cnt:
    if annotation_dict != None and iter[0] in annotation_dict:
        anno = annotation_dict[iter[0]]
    else:
        anno = []
    result.append(ll_solve_chrom([vcfgz_filenames, iter[0], max_dist, max_inspro, anno]))
'''
'''
for iter in chrom_cnt:
    # multi processes
    #print(iter[0])
    #if iter[1] == 0:
    #    continue
    if annotation_dict != None and iter[0] in annotation_dict:
        anno = annotation_dict[iter[0]]
    else:
        anno = []
    para = [(vcfgz_filenames, iter[0], max_dist, max_inspro, anno)]
    result.append(pool.map_async(ll_solve_chrom, para))
'''
'''
para = []
for iter in chrom_cnt:
    if annotation_dict != None and iter[0] in annotation_dict:
        anno = annotation_dict[iter[0]]
    else:
        anno = []
    para.append([vcfgz_filenames, iter[0], max_dist, max_inspro, anno])
result = pool.map(ll_solve_chrom, para, chunksize = 1)
'''