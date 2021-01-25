from cuteSV_calculation import cal_center, cal_ci
from cuteSV_linkedList import add_node, print_list, ListNode, Record
from cuteSV_output import output_result, solve_annotation
from pysam import VariantFile
from multiprocessing import Pool, Manager
import os
import sys
import time
import argparse
        

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
    #'''
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
    '''
    #print(chrom_cnt)
    return vcfgz_filenames, chrom_cnt, contiginfo, sample_ids


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
    print(contiginfo)
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
    print('finish %s in %s, total %dnodes'%(chrom, str(time.time() - start_time), len(result)))
    #print('finish %s at %s'%(chrom, str(time.time())))
    return result  # [[CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT, ANNOTATION], [], ...]


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
        result = pool.map(ll_solve_chrom, para, chunksize = 1)
        # result = pool.map(ll_solve_chrom, para)
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

    output_result(semi_result, sample_ids, args.output, contiginfo)

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
    args = parser.parse_args(sys.argv[1:])
    main(args)
    # /usr/bin/time -v python src/cuteSV_merge.py input/file_SRR.fofn sample_merged25.vcf ./ -t 16 -a hg19.refGene.gtf
    # /usr/bin/time -v python src/cuteSV_merge.py test.info test_merged.vcf ./ -t 16 -a hg19.refGene.gtf
    # /usr/bin/time -v python src/cuteSV_merge.py input/cutesv.fofn merged_cutesv.vcf ./ -t 16

'''
chrom_set = set()
base_cmd = 'cat '
for vcf_filename in vcf_filenames:
    base_cmd += vcf_filename + ' '
base_cmd += '| grep -v \'#\' | awk -F \'\\t\' \'{print $1}\' | uniq | sort | uniq'
print(base_cmd)
'''