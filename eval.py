from pysam import VariantFile
import vcf
import sys

def parse_record(record, sample_id):
    sv_type = record.info['SVTYPE']
    if sv_type == 'BND':
        sv_type = 'TRA'
    chrom1 = record.chrom
    name = record.id
    start = parse_to_int(record.pos)
    if record.info['SVTYPE'] == 'INS' and 'SVLEN' in record.info:
        end = parse_to_int(record.info['SVLEN'])
    else:
        try:
            end = parse_to_int(record.info['END'])
        except:
            try:
                end = parse_to_int(record.stop)
            except:
                pass   
    if record.info['SVTYPE'] == 'BND' or record.info['SVTYPE'] == 'TRA':
        tra_alt = str(record.alts[0])
        if tra_alt[0] == 'N':
            if tra_alt[1] == '[':
                tra_type = 'A'
            else:
                tra_type = 'B'
        elif tra_alt[0] == '[':
            tra_type = 'C'
        else:
            tra_type = 'D'
        if tra_alt[0] == 'N':
            tra_alt = tra_alt[2:-1]
        else:
            tra_alt = tra_alt[1:-2]
        chrom2 = tra_alt.split(':')[0]
        end = int(tra_alt.split(':')[1])
        strand = tra_type
    if record.info['SVTYPE'] != 'TRA' and record.info['SVTYPE'] != 'BND':
        chrom2 = record.chrom
        if 'STRAND' in record.info:
            strand = record.info['STRAND']
        elif 'STRANDS' in record.info:
            strand = record.info['STRANDS']
        else:
            strand = '.'
    return sv_type, chrom1, chrom2, start, end, name


def check_same_variant(sv_type1, start1, start2, end1, end2):
    if abs(start1 - start2) < 1000:
        if sv_type1 == 'INS':
            if 0.7 < min(end1, end2) / max(end1, end2) <= 1: 
                return True
        else:
            if abs(end1 - end2) < 1000:
                return True
    return False


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



def solve_two_vcf(vcf1, vcf2, output):  #vcf1是被比较的
    sv_dict1 = parse_nstd(vcf1)
    sv_dict2 = parse_nstd(vcf2)
    file = open(output, 'w')
    right = 0
    another = 0
    wrong = 0
    vis = set()
    for sv_type in ['INS', 'DEL', 'INV', 'DUP']:
        for chrom in sv_dict1[sv_type].keys():
            for record_list in sv_dict1[sv_type][chrom]:
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
                        if check_same_variant(sv_type, record_list[0], temp_list[0], record_list[1], temp_list[1]):
                            file.write('same variant on chr%s (%d, %d), (%d, %d)\n'%(chrom, record_list[0], record_list[1], temp_list[0], temp_list[1]))
                            flag = 1
                            right += 1
                            #vis.add(temp_list[2])
                            break
                if flag == 0:
                    file.write('wrong on chr%s %s, (%d, %d)\n'%(chrom, sv_type, record_list[0], record_list[1]))
                    wrong += 1 
                '''
                if flag == 2:
                    file.write('another match on chr%s %s, (%d, %d), (%d, %d)\n'%(chrom, sv_type, record_list[0], record_list[1], temp_start, temp_end))
                    another += 1 
                '''
    vis = set()
    if 'BND' in sv_dict1:
        for chrom in sv_dict1['BND']:
            for record_list in sv_dict1['BND'][chrom]:
                #print(record_list)
                wrong += 1
    
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
    '''
    print('right=%d, another=%d, wrong=%d'%(right, another, wrong))
    print('if another is right, precision=%f'%((right + another) / (right + another + wrong)))
    print('if another is wrong, pricision=%f'%(right  / (right + another + wrong)))
    '''
    return right, right + wrong  # vcf1数量


if __name__ == '__main__':
    vcf_standard = sys.argv[1]
    vcf_test = sys.argv[2]
    output = sys.argv[3]
    right1, standard_cnt = solve_two_vcf(vcf_standard, vcf_test, output)
    print(right1)
    print(standard_cnt)
    '''
    right2, test_cnt = solve_two_vcf(vcf_test, vcf_standard, output)
    precision = right2 / test_cnt
    recall = right1 / standard_cnt
    F1 = 2 * precision * recall / (precision + recall)
    print('precision = %f'%(round(precision, 4)))
    print('recall = %f'%(round(recall, 4)))
    print('F1 = %f'%(round(F1, 4)))
    '''
    '''
    vcf_standard = 'nstd162.GRCh37.variant_call.vcf'
    output = 'tmp'
    for call_tool in ['cutesv', 'pbsv', 'sniffles', 'svim']:
        for merge_tool in ['merged', 'survivor']:
            vcf_test = merge_tool + '_' + call_tool + '.vcf'
            right1, standard_cnt = solve_two_vcf(vcf_standard, vcf_test, output)
            right2, test_cnt = solve_two_vcf(vcf_test, vcf_standard, output)
            precision = right2 / test_cnt
            recall = right1 / standard_cnt
            F1 = 2 * precision * recall / (precision + recall)
            print('======' + merge_tool + ', ' + call_tool + '======')
            print('TP-base = %d'%(right1))
            print('TP-call = %d'%(right2))
            print('total-standard = %d'%(standard_cnt))
            print('total-input = %d'%(test_cnt))
            print('precision = %.4f'%(round(precision, 4)))
            print('recall = %.4f'%(round(recall, 4)))
            print('F1 = %.4f'%(round(F1, 4)))
    '''

            