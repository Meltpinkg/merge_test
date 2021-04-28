from pysam import VariantFile
import vcf
import sys
import os

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
    chrom2 = record.chrom
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
    return sv_type, chrom1, start, end, abs(svlen), supp_vec, chrom2


def check_same_variant(sv_type, start1, start2, end1, end2, len1, len2):
    #print('%s, %d, %d, %d, %d'%(sv_type, start1, start2, end1, end2))
    if sv_type == 'INS':
        if abs(start1 - start2) < 1000 and 0.7 < min(len1, len2) / max(len1, len2) <= 1:
            return True
    else:
        if abs(start1 - start2) < 1000 and 0.7 < min(len1, len2) / max(len1, len2) <= 1: # and max(start1, start2) <= min(end1, end2):
            return True
    return False

def check_same_bnd(sv_type, start1, start2, end1, end2):
    if abs(start1 - start2) < 1000 and abs(end1 - end2) < 1000:
        return True
    return False

def parse_vcf(filename):
    vcf_reader = VariantFile(filename, 'r')
    sv_dict = dict()
    for record in vcf_reader.fetch():
        sv_type, chrom1, start, end, svlen, supp_id, chrom2 = parse_record(record)
        if sv_type not in sv_dict:
            sv_dict[sv_type] = dict()
        if sv_type == 'BND':
            if chrom1 not in sv_dict[sv_type]:
                sv_dict[sv_type][chrom1] = dict()
            if chrom2 not in sv_dict[sv_type][chrom1]:
                sv_dict[sv_type][chrom1][chrom2] = list()
            sv_dict[sv_type][chrom1][chrom2].append([start, end, svlen, supp_id])
        else:
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
            if sv_type == 'INV':
                svlen = end - start + 1
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
    file = open(output, 'a')
    right = dict()
    wrong = dict()
    vis = set()
    idx = 0
    for sv_type in ['INS', 'DEL', 'INV', 'DUP', 'BND']:
        right[sv_type] = 0
        wrong[sv_type] = 0
    for sv_type in ['INS', 'DEL', 'INV', 'DUP']:
    #for sv_type in ['INV']:
        if sv_type not in sv_dict1:
            continue
        for chrom in sv_dict1[sv_type].keys():
            for record_list in sv_dict1[sv_type][chrom]:
                '''
                if sv_type == 'INV':
                    print('==')
                    print(record_list)
                    idx += 1
                '''
                flag = 0
                if chrom in sv_dict2[sv_type]:
                    for temp_list in sv_dict2[sv_type][chrom]:
                        '''
                        if record_list[0] == 64704785 and temp_list[0] == 64704786:
                            print(record_list[1])
                            print(temp_list[1])
                            print(check_same_variant(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1]))
                        '''
                        if check_same_variant(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1], record_list[2], temp_list[2]):
                            if sample_flag == 0:
                                file.write('same variant on chr%s (%d, %d, %d), (%d, %d, %d)\n'%(chrom, record_list[0], record_list[1], record_list[2], temp_list[0], temp_list[1], temp_list[2]))
                                flag = 1
                                right[sv_type] += 1
                            else:
                                for sample_id in record_list[2]:
                                    if sample_id in temp_list[2]:
                                        right[sv_type] += 1
                            break
                if flag == 0:
                    file.write('wrong on chr%s %s, (%d, %d, %d)\n'%(chrom, sv_type, record_list[0], record_list[1], record_list[2]))
                    if sample_flag == 0:
                        wrong[sv_type] += 1
                    else:
                        wrong[sv_type] += len(record_list[2])
                #if idx == 8:
                #    exit(0)
    if 'BND' in sv_dict1:
        sv_type = 'BND'
        for chrom in sv_dict1[sv_type].keys():
            for chrom2 in sv_dict1[sv_type][chrom].keys():
                for record_list in sv_dict1[sv_type][chrom][chrom2]:
                    flag = 0
                    if sv_type in sv_dict2 and chrom in sv_dict2[sv_type] and chrom2 in sv_dict2[sv_type][chrom]:
                        for temp_list in sv_dict2[sv_type][chrom][chrom2]:
                            if check_same_bnd(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1]):
                                file.write('same variant on chr(%s, %d, %d), (%s, %d, %d)\n'%(chrom, record_list[0], record_list[1], chrom2, temp_list[0], temp_list[1]))
                                flag = 1
                                right[sv_type] += 1
                                break
                    if flag == 0:
                        file.write('wrong on chr%s %s %s, (%d, %d)\n'%(chrom, chrom2, sv_type, record_list[0], record_list[1]))
                        wrong[sv_type] += 1
    return right, wrong

def solve_tp(sv_dict1, sv_dict2):  # sv_dict1是被比较的
    right = dict()
    wrong = dict()
    vis = set()
    idx = 0
    tp = list()
    fx = list()
    for sv_type in ['INS', 'DEL', 'INV', 'DUP', 'BND']:
        right[sv_type] = 0
        wrong[sv_type] = 0
    for sv_type in ['INS', 'DEL', 'INV', 'DUP']:
    #for sv_type in ['INV']:
        if sv_type not in sv_dict1:
            continue
        for chrom in sv_dict1[sv_type].keys():
            for record_list in sv_dict1[sv_type][chrom]:
                flag = 0
                if chrom in sv_dict2[sv_type]:
                    for temp_list in sv_dict2[sv_type][chrom]:
                        if check_same_variant(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1], record_list[2], temp_list[2]):
                            #file.write('same variant on chr%s (%d, %d, %d), (%d, %d, %d)\n'%(chrom, record_list[0], record_list[1], record_list[2], temp_list[0], temp_list[1], temp_list[2]))
                            tp.append([chrom, record_list[0], sv_type, record_list[1], record_list[2]]) # chrom, pos, type, end, len
                            flag = 1
                            right[sv_type] += 1
                            break
                if flag == 0:
                    #file.write('wrong on chr%s %s, (%d, %d, %d)\n'%(chrom, sv_type, record_list[0], record_list[1], record_list[2]))
                    fx.append([chrom, record_list[0], sv_type, record_list[1], record_list[2]])
                    wrong[sv_type] += 1
                #if idx == 8:
                #    exit(0)
    if 'BND' in sv_dict1:
        sv_type = 'BND'
        for chrom in sv_dict1[sv_type].keys():
            for chrom2 in sv_dict1[sv_type][chrom].keys():
                for record_list in sv_dict1[sv_type][chrom][chrom2]:
                    flag = 0
                    if sv_type in sv_dict2 and chrom in sv_dict2[sv_type] and chrom2 in sv_dict2[sv_type][chrom]:
                        for temp_list in sv_dict2[sv_type][chrom][chrom2]:
                            if check_same_bnd(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1]):
                                #file.write('same variant on chr(%s, %d, %d), (%s, %d, %d)\n'%(chrom, record_list[0], record_list[1], chrom2, temp_list[0], temp_list[1]))
                                tp.append([chrom, record_list[0], sv_type, record_list[1], record_list[2]])
                                flag = 1
                                right[sv_type] += 1
                                break
                    if flag == 0:
                        #file.write('wrong on chr%s %s %s, (%d, %d)\n'%(chrom, chrom2, sv_type, record_list[0], record_list[1]))
                        fx.append([chrom, record_list[0], sv_type, record_list[1], record_list[2]])
                        wrong[sv_type] += 1
    return right, wrong, tp, fx

if __name__ == '__main__':
    if len(sys.argv) == 4:
        vcf_standard = sys.argv[2]
        vcf_test = sys.argv[3]
        #sv_dict1 = parse_chr_merge(vcf_standard)
        sv_dict1 = parse_vcf(vcf_standard)
        sv_dict2 = parse_vcf(vcf_test)
        output = sys.argv[1]
        if not os.path.exists(output):
            os.mkdir(output)
        right1, wrong1, tpbase, fn = solve_tp(sv_dict1, sv_dict2) # standard_cnt
        right2, wrong2, tpcall, fp = solve_tp(sv_dict2, sv_dict1) # test_cnt
        #exit(0)
        right1_cnt = 0
        right2_cnt = 0
        standard_cnt = 0
        test_cnt = 0
        #print(right2)
        #print(wrong2)
        print('SVTYPE\tprecision\trecall\t\tF1\t\tTP-base\tTP-call\tbase\tcall\t')
        for sv_type in ['INS', 'DEL', 'INV', 'DUP', 'BND']:
        #for sv_type in ['DEL']:
            right1_cnt += right1[sv_type]
            right2_cnt += right2[sv_type]
            standard_cnt += right1[sv_type] + wrong1[sv_type]
            test_cnt += right2[sv_type] + wrong2[sv_type]
            precision = -1 if right2[sv_type] + wrong2[sv_type] == 0 else right2[sv_type] / (right2[sv_type] + wrong2[sv_type])
            recall = -1 if right1[sv_type] + wrong1[sv_type] == 0 else right1[sv_type] / (right1[sv_type] + wrong1[sv_type])
            F1 = 2 * precision * recall / (precision + recall)
            print('%s\t%f\t%f\t%f\t%d\t%d\t%d\t%d'%(sv_type, round(precision, 6), round(recall, 6), round(F1, 6), right1[sv_type], right2[sv_type], right1[sv_type] + wrong1[sv_type], right2[sv_type] + wrong2[sv_type]))
        precision = right2_cnt / test_cnt
        recall = right1_cnt / standard_cnt
        F1 = 2 * precision * recall / (precision + recall)
        print('%s\t%f\t%f\t%f\t%d\t%d\t%d\t%d'%('TOTAL', round(precision, 6), round(recall, 6), round(F1, 6), right1_cnt, right2_cnt, standard_cnt, test_cnt))
        # chrom, pos, type, end, len
        tpbase = sorted(tpbase, key = lambda x:(x[0],x[1]))
        tpcall = sorted(tpcall, key = lambda x:(x[0],x[1]))
        fp = sorted(fp, key = lambda x:(x[0],x[1]))
        fn = sorted(fn, key = lambda x:(x[0],x[1]))
        with open(output + '/tp-base.txt', 'w') as f:
            for x in tpbase:
                f.write('%s\t%d\tSVTYPE=%s;END=%d;SVLEN=%d\n'%(x[0], x[1], x[2], x[3], x[4]))
        with open(output + '/tp-call.txt', 'w') as f:
            for x in tpcall:
                f.write('%s\t%d\tSVTYPE=%s;END=%d;SVLEN=%d\n'%(x[0], x[1], x[2], x[3], x[4]))
        with open(output + '/fp.txt', 'w') as f:
            for x in fp:
                f.write('%s\t%d\tSVTYPE=%s;END=%d;SVLEN=%d\n'%(x[0], x[1], x[2], x[3], x[4]))
        with open(output + '/fn.txt', 'w') as f:
            for x in fn:
                f.write('%s\t%d\tSVTYPE=%s;END=%d;SVLEN=%d\n'%(x[0], x[1], x[2], x[3], x[4]))
        print('finish eval')
    elif len(sys.argv) == 5:
        vcf_standard = sys.argv[1]
        vcf_test = sys.argv[2]
        sv_dict1 = parse_chr_merge(vcf_standard)
        sv_dict2 = parse_vcf(vcf_test)
        output = sys.argv[3]
        file = open(output, 'w')
        file.close()
        sample_flag = int(sys.argv[4])
        right1, wrong1 = solve_two_vcf(sv_dict1, sv_dict2, output, sample_flag) # standard_cnt
        right2, wrong2 = solve_two_vcf(sv_dict2, sv_dict1, output, sample_flag) # test_cnt
        #exit(0)
        right1_cnt = 0
        right2_cnt = 0
        standard_cnt = 0
        test_cnt = 0
        #print(right2)
        #print(wrong2)
        print('SVTYPE\tprecision\trecall\t\tF1\t\tTP-base\tTP-call\tbase\tcall\t')
        for sv_type in ['INS', 'DEL', 'INV', 'DUP', 'BND']:
        #for sv_type in ['DEL']:
            right1_cnt += right1[sv_type]
            right2_cnt += right2[sv_type]
            standard_cnt += right1[sv_type] + wrong1[sv_type]
            test_cnt += right2[sv_type] + wrong2[sv_type]
            precision = -1 if right2[sv_type] + wrong2[sv_type] == 0 else right2[sv_type] / (right2[sv_type] + wrong2[sv_type])
            recall = -1 if right1[sv_type] + wrong1[sv_type] == 0 else right1[sv_type] / (right1[sv_type] + wrong1[sv_type])
            F1 = 2 * precision * recall / (precision + recall)
            print('%s\t%f\t%f\t%f\t%d\t%d\t%d\t%d'%(sv_type, round(precision, 6), round(recall, 6), round(F1, 6), right1[sv_type], right2[sv_type], right1[sv_type] + wrong1[sv_type], right2[sv_type] + wrong2[sv_type]))
        precision = right2_cnt / test_cnt
        recall = right1_cnt / standard_cnt
        F1 = 2 * precision * recall / (precision + recall)
        print('%s\t%f\t%f\t%f\t%d\t%d\t%d\t%d'%('TOTAL', round(precision, 6), round(recall, 6), round(F1, 6), right1_cnt, right2_cnt, standard_cnt, test_cnt))
    elif len(sys.argv) == 1:
        output = 'eval.tmp'
        file = open(output, 'w')
        file.close()
        sv_dict1 = parse_chr_merge('chr1_merge.vcf')
        #tool_list = ['output.bak311/merged', 'output322/merged', 'output.bak311/survivor']
        tool_list = ['output0/merged']
        #name_list = ['src', 'src1', 'sur']
        name_list = ['src', 'sur']
        for call_tool in ['cutesv', 'pbsv', 'sniffles', 'svim', 'cutesv0', 'svimS4']:
            print('===========================' + call_tool + '===========================')
            #print('SVTYPE\tprecision\trecall\t\tF1\t\tTP-base\tTP-call\tbase\tcall\t')
            for i in range(1):
                merge_tool = tool_list[i]
                vcf_test = merge_tool + '_' + call_tool + '.vcf'
                sv_dict2 = parse_vcf(vcf_test)
                right1, wrong1 = solve_two_vcf(sv_dict1, sv_dict2, output, 0)
                right2, wrong2 = solve_two_vcf(sv_dict2, sv_dict1, output, 0)
                #print(right1)
                #print(right2)
                right1_cnt = 0
                right2_cnt = 0
                standard_cnt = 0
                test_cnt = 0
                print('SVTYPE\tprecision\trecall\t\tF1\t\tTP-base\tTP-call\tbase\tcall\t')
                for sv_type in ['INS', 'DEL', 'INV', 'DUP']:
                #for sv_type in ['DEL']:
                    right1_cnt += right1[sv_type]
                    right2_cnt += right2[sv_type]
                    standard_cnt += right1[sv_type] + wrong1[sv_type]
                    test_cnt += right2[sv_type] + wrong2[sv_type]
                    precision = -1 if right2[sv_type] + wrong2[sv_type] == 0 else right2[sv_type] / (right2[sv_type] + wrong2[sv_type])
                    recall = -1 if right1[sv_type] + wrong1[sv_type] == 0 else right1[sv_type] / (right1[sv_type] + wrong1[sv_type])
                    F1 = 2 * precision * recall / (precision + recall)
                    print('%s\t%f\t%f\t%f\t%d\t%d\t%d\t%d'%(sv_type, round(precision, 6), round(recall, 6), round(F1, 6), right1[sv_type], right2[sv_type], right1[sv_type] + wrong1[sv_type], right2[sv_type] + wrong2[sv_type]))
                precision = right2_cnt / test_cnt
                recall = right1_cnt / standard_cnt
                F1 = 2 * precision * recall / (precision + recall)
                print('%s\t%f\t%f\t%f\t%d\t%d\t%d\t%d'%('TOTAL', round(precision, 6), round(recall, 6), round(F1, 6), right1_cnt, right2_cnt, standard_cnt, test_cnt))
    else:
        print('parameter error %d'%(len(sys.argv)))
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