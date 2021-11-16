#####################################################
#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||  :  |||// \
#                / _||||| -:- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  ___/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#         \  \ `_.   \_ __\ /__ _/   .-` /  /
#     =====`-.____`.___ \_____/___.-`___.-'=====
#                       `=---='
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#               佛祖保佑         永无BUG
#####################################################
from asyncio.events import get_event_loop
import threading
import asyncio
import sys
from queue import Queue
import time
from pysam import VariantFile
from multiprocessing import Pool
from cuteSV_linkedList import *
from heapq import *
from cuteSV_merge import first_sort

# parse vcf and return {chr -> svtype -> list(Record)}
def parse_vcf(para):
    start_time = time.time()
    filename = para[0]
    idx = para[1]
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict()
    # contigs
    contiginfo = dict()
    for i in range(len(vcf_reader.header.contigs)):
        try:
            contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[str(vcf_reader.header.contigs[i].name)] = 0
    record_dict['contig'] = contiginfo
    # records
    for record in vcf_reader.fetch():
        chr = record.chrom
        svtype = 'BND' if record.info['SVTYPE'] == 'TRA' else record.info['SVTYPE']
        if chr not in record_dict:
            record_dict[chr] = dict()
        if svtype not in record_dict[chr]:
            record_dict[chr][svtype] = []
        record_dict[chr][svtype].append(Record(record, idx))
    # sample_id
    record_dict['sample'] = vcf_reader.header.samples[0]
    #print(record_dict.keys())
    #print('finish %d in %s'%(idx, str(time.time() - start_time)))
    return record_dict
def parse_vcfs(filenames, threads):
    print('multi processing')
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
        temp = file_dict[fileidx].get()[0]  #{chr -> svtype -> List(Record)}
        sample_temp[fileidx] = temp['sample']
        contig_temp = temp['contig']
        temp.pop('sample')
        temp.pop('contig')
        for chr in temp:
            if chr not in chr_dict:
                chr_dict[chr] = dict()
            if chr not in contiginfo:
                contiginfo[chr] = contig_temp[chr]
            for svtype in temp[chr]:
                if svtype not in chr_dict[chr]:
                    chr_dict[chr][svtype] = dict()
                chr_dict[chr][svtype][fileidx] = temp[chr][svtype]

    #print(contiginfo)
    #print(chr_dict.keys())
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

def mp():
    start_time = time.time()
    chr_dict, sample_ids, contiginfo = parse_vcfs('real_test/real_200.fofn', 16)
    print('finish in %s'%(round(time.time() - start_time, 6)))
    for chr in chr_dict:
        print(chr)
        print(len(chr_dict[chr][0]))
        print(len(chr_dict[chr][1]))

def job(data, q):
    for i in range(len(data)):
        data[i] = data[i] ** 2
    q.put(data)
def mt():
    start_time = time.time()
    q = Queue()
    threads = []
    data = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    for i in range(3):
        t = threading.Thread(target=job, args=(data[i], q))
        t.start()
        t.join()
    #for thread in threads:
    #    thread.join()
    results = []
    while not q.empty():
        results.append(q.get())
    print(results)

def parse_vcf_thread(filename, idx, que):
    start_time = time.time()
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict()
    # contigs
    contiginfo = dict()
    for i in range(len(vcf_reader.header.contigs)):
        try:
            contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[str(vcf_reader.header.contigs[i].name)] = 0
    record_dict['contig'] = contiginfo
    # records
    for record in vcf_reader.fetch():
        chr = record.chrom
        if chr not in record_dict:
            record_dict[chr] = []
        record_dict[chr].append(Record(record, idx))
    # sample_id
    record_dict['sample'] = vcf_reader.header.samples[0]
    record_dict['idx'] = idx
    #print(record_dict.keys())
    #print('finish %d in %s'%(idx, str(time.time() - start_time)))
    que.put(record_dict)
def parse_vcfs_thread(filenames, threads):
    print('multi threads')
    start_time = time.time()
    que = Queue()
    idx = 0
    with open(filenames, 'r') as f:
        for filename in f:
            #file_dict[idx] = pool.map_async(parse_vcf, [(filename.strip(), idx)])
            t = threading.Thread(target=parse_vcf_thread, args=(filename.strip(), idx, que))
            t.start()
            t.join()
            idx += 1
    #print(file_dict.keys())
    chr_dict = dict()
    sample_temp = ['' for i in range(idx)]
    contiginfo = dict()
    #for fileidx in file_dict:
    while not que.empty():
        temp = que.get()
        fileidx = temp['idx']
        #temp = file_dict[fileidx].get()[0]  #{chr -> List(Record)}
        sample_temp[fileidx] = temp['sample']
        contig_temp = temp['contig']
        temp.pop('sample')
        temp.pop('contig')
        temp.pop('idx')
        for chr in temp:
            if chr not in chr_dict:
                chr_dict[chr] = dict()
                if chr not in contiginfo:
                    contiginfo[chr] = contig_temp[chr]
            chr_dict[chr][fileidx] = temp[chr]

    #print(contiginfo)
    #print(chr_dict.keys())
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

async def parse_vcf_asyncio(filename, idx):
    start_time = time.time()
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict()
    # contigs
    contiginfo = dict()
    for i in range(len(vcf_reader.header.contigs)):
        try:
            contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[str(vcf_reader.header.contigs[i].name)] = 0
    record_dict['contig'] = contiginfo
    # records
    for record in vcf_reader.fetch():
        chr = record.chrom
        if chr not in record_dict:
            record_dict[chr] = []
        record_dict[chr].append(Record(record, idx))
    # sample_id
    record_dict['sample'] = vcf_reader.header.samples[0]
    record_dict['idx'] = idx
    #print(record_dict.keys())
    #print('finish %d in %s'%(idx, str(time.time() - start_time)))
    return record_dict
def parse_vcfs_asyncio(filenames, threads):
    print('asyncio')
    start_time = time.time()
    loop = asyncio.get_event_loop()
    tasks = []
    idx = 0
    with open(filenames, 'r') as f:
        for filename in f:
            #file_dict[idx] = pool.map_async(parse_vcf, [(filename.strip(), idx)])
            #t = threading.Thread(target=parse_vcf_asyncio, args=(filename.strip(), idx, que))
            tasks.append(loop.create_task(parse_vcf_asyncio(filename.strip(), idx)))
            idx += 1
    wait_coro = asyncio.wait(tasks)
    loop.run_until_complete(wait_coro)
    print('finish multi in %s'%(str(time.time() - start_time)))
    chr_dict = dict()
    sample_temp = ['' for i in range(idx)]
    contiginfo = dict()
    for task in tasks:
        temp = task.result()
        fileidx = temp['idx']
        #temp = file_dict[fileidx].get()[0]  #{chr -> List(Record)}
        sample_temp[fileidx] = temp['sample']
        contig_temp = temp['contig']
        temp.pop('sample')
        temp.pop('contig')
        temp.pop('idx')
        for chr in temp:
            if chr not in chr_dict:
                chr_dict[chr] = dict()
                if chr not in contiginfo:
                    contiginfo[chr] = contig_temp[chr]
            chr_dict[chr][fileidx] = temp[chr]

    #print(contiginfo)
    #print(chr_dict.keys())
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

def parse_vcf_single(filename, idx):
    start_time = time.time()
    vcf_reader = VariantFile(filename, 'r')
    record_dict = dict()
    # contigs
    contiginfo = dict()
    for i in range(len(vcf_reader.header.contigs)):
        try:
            contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
        except:
            print('contig length not find')
            contiginfo[str(vcf_reader.header.contigs[i].name)] = 0
    record_dict['contig'] = contiginfo
    # records
    for record in vcf_reader.fetch():
        chr = record.chrom
        if chr not in record_dict:
            record_dict[chr] = []
        record_dict[chr].append(Record(record, idx))
    # sample_id
    record_dict['sample'] = vcf_reader.header.samples[0]
    record_dict['idx'] = idx
    #print(record_dict.keys())
    #print('finish %d in %s'%(idx, str(time.time() - start_time)))
    return record_dict
def parse_vcfs_single(filenames, threads):
    print('single')
    start_time = time.time()
    file_dict = dict()
    idx = 0
    with open(filenames, 'r') as f:
        for filename in f:
            file_dict[idx] = parse_vcf_single(filename.strip(), idx)
            idx += 1
    print('finish multi in %s'%(str(time.time() - start_time)))
    chr_dict = dict()
    sample_temp = ['' for i in range(idx)]
    contiginfo = dict()
    for fileidx in file_dict:
        temp = file_dict[fileidx]
        fileidx = temp['idx']
        #temp = file_dict[fileidx].get()[0]  #{chr -> List(Record)}
        sample_temp[fileidx] = temp['sample']
        contig_temp = temp['contig']
        temp.pop('sample')
        temp.pop('contig')
        temp.pop('idx')
        for chr in temp:
            if chr not in chr_dict:
                chr_dict[chr] = dict()
                if chr not in contiginfo:
                    contiginfo[chr] = contig_temp[chr]
            chr_dict[chr][fileidx] = temp[chr]

    #print(contiginfo)
    #print(chr_dict.keys())
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

def parse_vcf_fuse(para):
    start_time = time.time()
    filenames = para[0]
    idx = int((para[1] - 1) / 10) * 10 + 1
    record_dict_list = []
    for filename in filenames:
        vcf_reader = VariantFile(filename, 'r')        
        record_dict = dict()
        # contigs
        contiginfo = dict()
        for i in range(len(vcf_reader.header.contigs)):
            try:
                contiginfo[str(vcf_reader.header.contigs[i].name)] = int(vcf_reader.header.contigs[i].length)
            except:
                print('contig length not find')
                contiginfo[str(vcf_reader.header.contigs[i].name)] = 0
        record_dict['contig'] = contiginfo
        # records
        for record in vcf_reader.fetch():
            chr = record.chrom
            if chr not in record_dict:
                record_dict[chr] = []
            record_dict[chr].append(Record(record, idx))
        # sample_id
        record_dict['sample'] = vcf_reader.header.samples[0]
        record_dict_list.append(record_dict)
        idx += 1
    print('finish %d in %s'%(idx, str(time.time() - start_time)))
    return record_dict_list
def parse_vcfs_fuse(filenames, threads):
    print('processing + threading')
    start_time = time.time()
    pool = Pool(processes = threads)
    file_dict = dict()
    file_list = []
    idx = 0
    with open(filenames, 'r') as f:
        for filename in f:
            idx += 1
            file_list.append(filename.strip())
            if idx % 10 == 0:
                file_dict[idx] = pool.map_async(parse_vcf_fuse, [(file_list, idx)]) # [idx-9, idx]
                file_list = []
    print(idx)
    if idx % 10 != 0:
        file_dict[idx] = pool.map_async(parse_vcf_fuse, [(file_list, idx)]) # remain
    pool.close()
    pool.join()
    print('join')
    #print(file_dict.keys())
    chr_dict = dict()
    sample_temp = ['' for i in range(idx + 1)]
    contiginfo = dict()
    for idx in file_dict:
        temp_list = file_dict[idx].get()[0] # record_dict_list
        fileidx = int((idx - 1) / 10) * 10 + 1
        for temp in temp_list:  #{chr -> List(Record)}
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
            fileidx += 1

    #print(contiginfo)
    #print(chr_dict.keys())
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

async def fac(num):
    await asyncio.sleep(1)
    print(time.time())
    return num*num
def mm():
    loop = asyncio.get_event_loop()
    #tasks = loop.create_task(coro)
    #tasks = [fac(1), fac(2), fac(3)]
    tasks = []
    for i in range(3):
        tasks.append(loop.create_task(fac(i)))
    wait_coro = asyncio.wait(tasks)
    #loop.run_until_complete(asyncio.wait(tasks))
    loop.run_until_complete(wait_coro)
    for task in tasks:
        print(task.result())
    loop.close()
    #loop.run_until_complete(tasks)
    #print(tasks.result())

if __name__ == '__main__':
    #mp()
    #mt()
    #mm()
    '''
    coro = fac(1)
    loop = asyncio.get_event_loop()
    task = loop.create_task(coro)
    print('运行情况：', task)

    loop.run_until_complete(task)
    print('再看下运行情况：', task)
    loop.close()
    '''
    chr_dict, sample_ids, contiginfo = parse_vcfs(sys.argv[1], int(sys.argv[2]))
    #chr_dict, sample_ids, contiginfo = parse_vcfs_thread(sys.argv[1], int(sys.argv[2]))
    #chr_dict, sample_ids, contiginfo = parse_vcfs_asyncio(sys.argv[1], int(sys.argv[2]))
    #chr_dict, sample_ids, contiginfo = parse_vcfs_single(sys.argv[1], int(sys.argv[2]))
    #chr_dict, sample_ids, contiginfo = parse_vcfs_fuse(sys.argv[1], int(sys.argv[2]))
    cnt1 = cnt2 = 0
    '''for chr in chr_dict:
        print('===  ' + chr)
        for svtype in chr_dict[chr]:
            print('+++  ' + svtype)
            for fileidx in chr_dict[chr][svtype]:
                print('---  ' + str(fileidx))
                for record in chr_dict[chr][svtype][fileidx]:
                    print(record.to_string())'''
    for fileidx in chr_dict['1']['INS']:
        cnt1 += len(chr_dict['1']['INS'][fileidx])
    ans = first_sort(chr_dict['1']['INS'], 1000)
    print(len(ans))
    for node in ans:
        print('===  ' + str(len(node)))
        for idx in node:
            print(idx)
            for xx in node[idx]:
                print(xx.to_string())
            cnt2 += len(node[idx])
    print(cnt1)
    print(cnt2)
'''
start_time = time.time()
q = Queue()
threads = []
data = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
for i in range(3):
    t = threading.Thread(target=job, args=(data[i], q))
    t.start()
    threads.append(t)
for thread in threads:
    thread.join()
results = []
for i in range(3):
    results.append(q.get())
print(results)
'''