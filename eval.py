from pysam import VariantFile
import vcf
import sys

def parse_to_int(sth):
    if sth == None:
        return 0
    elif isinstance(sth, str):
        return int(sth)
    elif isinstance(sth, list):
        return parse_to_int(sth[0])
    elif isinstance(sth, tuple):
        return parse_to_int(sth[0])
    elif isinstance(sth, int):
        return sth
    else:
        return sth

def phase_bnd(alt, chr, pos):
	# print(alt)
	# print(str(alt[0]))
	if alt[0] == ']':
		form = ']]N'
		chr2 = alt.split(':')[0][1:]
		pos2 = int(alt.split(':')[1][:-2])
	elif alt[0] == '[':
		form = '[[N'
		chr2 = alt.split(':')[0][1:]
		pos2 = int(alt.split(':')[1][:-2])
	else:
		# print(type(alt))
		if alt[1] == ']':
			form = 'N]]'
			chr2 = alt.split(':')[0][2:]
			pos2 = int(alt.split(':')[1][:-1])
		else:
			form = 'N[['
			chr2 = alt.split(':')[0][2:]
			pos2 = int(alt.split(':')[1][:-1])

	try:
		if int(chr) <= int(chr2):
			if form == 'N[[':
				form = ']]N'
			if form == ']]N':
				form = 'N[['
			return form, chr, pos, chr2, pos2
		else:
			return form, chr2, pos2, chr, pos
	except:
		return form, chr, pos, chr2, pos2

# end：end[TRA], svlen[OTHERS]
def parse_record(record):
    sample_ids = ['CHM13', 'CHM1', 'HG00268', 'HG00514', 'HG00733', 'HG01352', 'HG02059', 'HG02106', 'HG02818', 'HG04217', 'HX1', 'NA12878', 'NA19240', 'NA19434']
    sv_type = record.info['SVTYPE'][:3]
    if sv_type == 'TRA':
        sv_type = 'BND'
    chrom1 = record.chrom
    start = parse_to_int(record.pos)
    try:
        svlen = record.info['SVLEN']
        try:
            if len(svlen):
                svlen = svlen[0]
        except:
            pass
    except:
        svlen = 0
    try:
        end = record.stop
    except:
        end = start
    if sv_type == 'BND':
        bnd = str(record.alts[0])
        form, chrom1, start, chrom2, end = phase_bnd(bnd, chrom1, start)

    supp_id = record.info['SUPP_VEC']
    supp_vec = list()
    for i in range(len(supp_id)):
        if supp_id[i] == '1':
            supp_vec.append(sample_ids[i])
    return sv_type, chrom1, start, end, abs(svlen), supp_vec


def check_same_variant(sv_type, start1, start2, end1, end2, len1, len2):
    #print('%s, %d, %d, %d, %d'%(sv_type, start1, start2, end1, end2))
    if sv_type == 'INS':
        if abs(start1 - start2) < 1000 and 0.7 < min(len1, len2) / max(len1, len2) <= 1:
            return True
    else:
        if abs(start1 - start2) < 1000 and 0.7 < min(len1, len2) / max(len1, len2) <= 1 and max(start1, start2) <= min(end1, end2):
            return True
    return False


def parse_vcf(filename):
    vcf_reader = VariantFile(filename, 'r')
    sv_dict = dict()
    for record in vcf_reader.fetch():
        sv_type, chrom1, start, end, svlen, supp_id = parse_record(record)
        if sv_type not in sv_dict:
            sv_dict[sv_type] = dict()
        if chrom1 in sv_dict[sv_type]:
            sv_dict[sv_type][chrom1].append([start, end, svlen, supp_id])
        else:
            sv_dict[sv_type][chrom1] = [[start, end, svlen, supp_id]]
    return sv_dict

def parse_chr_merge(filename):
    sample_ids = ['CHM1', 'CHM13', 'HG00268', 'HG00514', 'HG00733', 'HG01352', 'HG02059', 'HG02106', 'HG02818', 'HG04217', 'HX1', 'NA12878', 'NA19240', 'NA19434']
    sv_dict = dict()
    with open(filename, 'r') as f:
        for line in f:
            if line[:3] == 'CHM':
                continue
            lineseq = line.strip().split('\t')
            chrom1 = lineseq[0]
            start = int(lineseq[1])
            sv_type = lineseq[2].split('SVTYPE=')[1].split(';')[0]
            if sv_type not in sv_dict:
                sv_dict[sv_type] = dict()
            svlen = int(lineseq[2].split(';SVLEN=')[1].split(';')[0])
            end = int(lineseq[2].split(';END=')[1].split(';')[0])
            samples = list()
            for i in range(14):
                if lineseq[i + 3] != '0/0':
                    samples.append(sample_ids[i])
            if chrom1 in sv_dict[sv_type]:
                sv_dict[sv_type][chrom1].append([start, end, svlen, samples])
            else:
                sv_dict[sv_type][chrom1] = [[start, end, svlen, samples]]
    return sv_dict

# sample_flag == 1: take sample into account
def solve_two_vcf(sv_dict1, sv_dict2, output, sample_flag):  # sv_dict1是被比较的
    file = open(output, 'w')
    right = 0
    wrong = 0
    vis = set()
    idx = 0
    for sv_type in ['INS', 'DEL', 'INV', 'DUP']:
        if sv_type not in sv_dict1:
            continue
        for chrom in sv_dict1[sv_type].keys():
            for record_list in sv_dict1[sv_type][chrom]:
                '''
                if sv_type == 'DEL':
                    print('==')
                    print(record_list)
                    idx += 1
                '''
                flag = 0
                temp_start = 0
                temp_end = 0
                if chrom in sv_dict2[sv_type]:
                    for temp_list in sv_dict2[sv_type][chrom]:
                        '''
                        if record_list[0] == 64704785 and temp_list[0] == 64704786:
                            print(record_list[1])
                            print(temp_list[1])
                            print(check_same_variant(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1]))
                        '''
                        if temp_list[0] < 500000 and sv_type == 'DEL':
                            #print(temp_list)
                            pass
                        if check_same_variant(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1], record_list[2], temp_list[2]):
                            if sample_flag == 0:
                                file.write('same variant on chr%s (%d, %d), (%d, %d)\n'%(chrom, record_list[0], record_list[1], temp_list[0], temp_list[1]))
                                flag = 1
                                right += 1
                            else:
                                for sample_id in record_list[2]:
                                    if sample_id in temp_list[2]:
                                        right += 1
                            break
                if flag == 0:
                    file.write('wrong on chr%s %s, (%d, %d, %d)\n'%(chrom, sv_type, record_list[0], record_list[1], record_list[2]))
                    if sample_flag == 0:
                        wrong += 1
                    else:
                        wrong += len(record_list[2])
                #if idx == 8:
                #    exit(0)
    vis = set()
    if 'BND' in sv_dict1:
        for chrom in sv_dict1['BND']:
            for record_list in sv_dict1['BND'][chrom]:
                #print(record_list)
                if sample_flag == 0:
                    wrong += 1
                else:
                    wrong += len(record_list[2])
    
    '''
    for sv_type in ['TRA']:
        for chrom in sv_dict1[sv_type].keys():
            for chrom2 in sv_dict1[sv_type][chrom].keys():
                if chrom2 not in sv_dict2[sv_type][chrom]:
                    for record_list in sv_dict1[sv_type][chrom][chrom2]:
                        file.write('wrong on chr%s %s, TRA, (%d, %d)'%(chrom, chrom2, record_list[0], record_list[1]))
                    continue
                
                for record_list in sv_dict1[sv_type][chrom][chrom2]:
                    flag = 0
                    temp_start = 0
                    temp_end = 0
                    for temp_list in sv_dict2[sv_type][chrom][chrom2]:
                        if check_same_variant(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1]) and temp_list[5] not in vis:
                            if temp_list[5] in vis:
                                flag = 1
                                file.write('another match on chr%s, (%d, %d), (%d, %d)\n'%(chrom, record_list[0], record_list[1], temp_list[0], temp_list[1]))
                            else:
                                file.write('same variant on chr%s (%d, %d), (%d, %d)\n'%(chrom, record_list[0], record_list[1], temp_list[0], temp_list[1]))
                                flag = 2
                                right += 1
                                vis.add(temp_list[5])
                                break
                        #if temp_list[0] - record_list[0] > 5000:
                        #    break
                if flag == 0:
                    file.write('wrong on chr%s %s, (%d, %d)\n'%(chrom, sv_type, record_list[0], record_list[1]))
                    wrong += 1 
    '''
    return right, right + wrong  # vcf1数量


if __name__ == '__main__':
    if len(sys.argv) > 1:
        vcf_standard = sys.argv[1]
        vcf_test = sys.argv[2]
        sv_dict1 = parse_chr_merge(vcf_standard)
        sv_dict2 = parse_vcf(vcf_test)
        output = sys.argv[3]
        sample_flag = int(sys.argv[4])
        right1, standard_cnt = solve_two_vcf(sv_dict1, sv_dict2, output, sample_flag)
        right2, test_cnt = solve_two_vcf(sv_dict2, sv_dict1, output, sample_flag)
        precision = right2 / test_cnt
        recall = right1 / standard_cnt
        F1 = 2 * precision * recall / (precision + recall)
        print('TP-base = %d'%(right1))
        print('TP-call = %d'%(right2))
        print('total-standard = %d'%(standard_cnt))
        print('total-input = %d'%(test_cnt))
        print('precision = %f'%(round(precision, 4)))
        print('recall = %f'%(round(recall, 4)))
        print('F1 = %f'%(round(F1, 4)))
    else:
        output = 'eval.tmp'
        sv_dict1 = parse_chr_merge('chr1_merge.vcf')
        for call_tool in ['cutesv', 'pbsv', 'sniffles', 'svim', 'cutesv0']:
            for merge_tool in ['merged1']:
                vcf_test = 'output/' + merge_tool + '_' + call_tool + '.vcf'
                sv_dict2 = parse_vcf(vcf_test)
                right1, standard_cnt = solve_two_vcf(sv_dict1, sv_dict2, output, 0)
                right2, test_cnt = solve_two_vcf(sv_dict2, sv_dict1, output, 0)
                precision = right2 / test_cnt
                recall = right1 / standard_cnt
                F1 = 2 * precision * recall / (precision + recall)
                print('======' + merge_tool + ', ' + call_tool + '======')
                print('TP-base = %d'%(right1))
                print('TP-call = %d'%(right2))
                print('base = %d'%(standard_cnt))
                print('call = %d'%(test_cnt))
                print('precision = %.4f'%(round(precision, 4)))
                print('recall = %.4f'%(round(recall, 4)))
                print('F1 = %.4f'%(round(F1, 4)))
    #'''

'''
def parse_vcf(filename):
    vcf_reader = VariantFile(filename, 'r')
    sv_dict = dict()
    for bx in ['INS', 'DEL', 'INV', 'DUP', 'TRA']:
        sv_dict[bx] = dict()
    for record in vcf_reader.fetch():
        sv_type, chrom1, chrom2, start, end, name = parse_record(record, 'NULL')
        if sv_type == 'TRA':
            if chrom1 in sv_dict[sv_type]:
                if chrom2 in sv_dict[sv_type][chrom1]:
                    sv_dict[sv_type][chrom1][chrom2].append([start, end, name])
                else:
                    sv_dict[sv_type][chrom1][chrom2] = [[start, end, name]]

            else:
                sv_dict[sv_type][chrom1] = dict()
                sv_dict[sv_type][chrom1][chrom2] = [[start, end, name]]
        else:
            if chrom1 in sv_dict[sv_type]:
                sv_dict[sv_type][chrom1].append([start, end, name])
            else:
                sv_dict[sv_type][chrom1] = [[start, end, name]]
    return sv_dict
'''
def parse_nstd(filename):
    sv_dict = dict()
    for bx in ['INS', 'DEL', 'INV', 'DUP', 'BND']:
        sv_dict[bx] = dict()
    with open(filename, 'r') as f:
        idx = 0
        pres = 0
        pree = 0
        pre_type = ''
        for line in f:
            if line[0] == '#':
                continue
            lineseq = line.strip().split('\t')
            chrom1 = lineseq[0]
            if 'nstd162' in filename and chrom1 != '1':
                break
            start = int(lineseq[1])
            sv_type = lineseq[7].split('SVTYPE=')[1].split(';')[0]
            if sv_type == 'TRA':
                sv_type = 'BND'
            if 'DUP' in sv_type:
                sv_type = 'DUP'
            idx += 1
            if sv_type == 'INS' or sv_type == 'DEL':
                end = abs(int(lineseq[7].split(';SVLEN=')[1].split(';')[0]))
            elif sv_type == 'BND':
                end = 0
            else:
                end = abs(int(lineseq[7].split(';END=')[1].split(';')[0]))
            name = 'variant' + str(idx)
            #sample = lineseq[7].split('SAMPLE=')[1].split(';')[0]
            if start == pres and pree == end and sv_type == pre_type and sv_type != 'BND':
                continue
            pres = start
            pree = end
            pre_type = sv_type
            if chrom1 in sv_dict[sv_type]:
                sv_dict[sv_type][chrom1].append([start, end, name])
            else:
                sv_dict[sv_type][chrom1] = [[start, end, name]]
    return sv_dict