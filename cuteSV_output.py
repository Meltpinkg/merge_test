import time

'''
    var_list: [[a, b, info], ...]
    return [info] for which a <= breakpoint < b
'''
def add_anotation_in(var_list, breakpoint):
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
            info.append(var_list[i][2])
        if breakpoint - var_list[i][0] > 3000000:
            break
    return info

'''
    var_list: [[a, b, info], ...]
    return [info] for which [start, end] and [a, b] overlap
'''
def add_anotation_overlap(var_list, start, end):
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
            info.append(var_list[i][2])
        if start - var_list[i][0] > 3000000:
            break
    for i in range(left + 1, len(var_list), 1):
        if var_list[i][0] < end:
            info.append(var_list[i][2])
        else:
            break
    return info


def parse_annotation_file(annotation_file):
    start_time = time.time()
    annotation_dict = dict()  # [chrom -> [[a, b, info], ...]]
    with open(annotation_file, 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            chr = seq[0]
            if chr[:3] == "chr":
                chr = chr[3:]
            if chr not in annotation_dict:
                annotation_dict[chr] = []
            annotation_dict[chr].append([int(seq[3]), int(seq[4]), seq[8]])
    print(time.time() - start_time)
    return annotation_dict
            


def output_result(semi_result, sample_ids, output_file, contiginfo, annotation_file):
    file = open(output_file, 'w')
    generate_header(file, contiginfo, sample_ids)
    #print('start writing to file ' + chrom)
    if annotation_file != None:
        annotation_dict = parse_annotation_file(annotation_file)
    for item in semi_result:  # [CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT]
        supp_vec = ''
        supp_id = []
        for i in range(len(sample_ids)):
            if i in item[5]:
                supp_vec += '1'
                supp_id.append(sample_ids[i])
            else:
                supp_vec += '0'
        can_record = item[2]
        if can_record.qual == "." or can_record.qual is None:
            filter_lable = "PASS"
        else:
            filter_lable = "PASS" if float(can_record.qual) >= 5.0 else "q5"
        annotation_list = []
        if annotation_file != None:
            if can_record.chrom1 in annotation_dict:
                if can_record.type == 'INS':
                    annotation_list = add_anotation_in(annotation_dict[can_record.chrom1], can_record.start)
                elif can_record.type == 'DEL':
                    annotation_list = add_anotation_overlap(annotation_dict[can_record.chrom1], can_record.start, can_record.start + can_record.end)
                elif can_record.type == 'INV' or can_record.type == 'DUP':
                    annotation_list = add_anotation_overlap(annotation_dict[can_record.chrom1], can_record.start, can_record.end)
                elif can_record.type == 'TRA':
                    annotation_list = add_anotation_in(annotation_dict[can_record.chrom1], can_record.start)
            if can_record.type == 'TRA' and can_record.chrom2 in annotation_dict:
                annotation_list += add_anotation_in(annotation_dict[can_record.chrom2], can_record.end)
        if can_record.type == 'INS':
            sv_len = can_record.end
            sv_end = can_record.start
            info_list = "SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN};SUPP={SUPP};SUPP_VEC={SUPP_VEC};SUPP_ID={SUPP_ID}".format(
                    SUPP = len(item[5]),
                    SUPP_ID = ','.join(supp_id),
                    SUPP_VEC = supp_vec,
                    SVTYPE = can_record.type, 
                    SVLEN = sv_len, 
                    END = sv_end, 
                    CIPOS = str(item[3]),
                    CILEN = str(item[4]))
            if can_record.strand != '.':
                info_list += ';STRAND=' + can_record.strand
            if annotation_list != []:
                info_list += ';ANNOTATION=' + '|'.join(annotation_list)
        elif can_record.type == 'DEL':
            sv_len = -can_record.end
            sv_end = can_record.start + can_record.end
            info_list = "SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CIEND={CIEND};SUPP={SUPP};SUPP_VEC={SUPP_VEC};SUPP_ID={SUPP_ID}".format(
                    SUPP = len(item[5]),
                    SUPP_ID = ','.join(supp_id),
                    SUPP_VEC = supp_vec,
                    SVTYPE = can_record.type, 
                    SVLEN = sv_len, 
                    END = sv_end, 
                    CIPOS = str(item[3]), 
                    CIEND = str(item[4]))
            if can_record.strand != '.':
                info_list += ';STRAND=' + can_record.strand
            if annotation_list != []:
                info_list += ';ANNOTATION=' + '|'.join(annotation_list)
        elif can_record.type == 'INV' or can_record.type == 'DUP':
            sv_len = can_record.end - can_record.start + 1
            sv_end = can_record.end
            info_list = "SVTYPE={SVTYPE};SVLEN={SVLEN};END={END}SUPP={SUPP};SUPP_VEC={SUPP_VEC};SUPP_ID={SUPP_ID}".format(
                    SUPP = len(item[5]),
                    SUPP_ID = ','.join(supp_id),
                    SUPP_VEC = supp_vec,
                    SVTYPE = can_record.type, 
                    SVLEN = sv_len, 
                    END = sv_end)
            if can_record.strand != '.':
                info_list += ';STRAND=' + can_record.strand
            if annotation_list != []:
                info_list += ';ANNOTATION=' + '|'.join(annotation_list)
        elif can_record.type == 'BND':
            info_list = "SVTYPE={SVTYPE};SUPP={SUPP};SUPP_VEC={SUPP_VEC};SUPP_ID={SUPP_ID}".format(
                    SUPP = len(item[5]),
                    SUPP_ID = ','.join(supp_id),
                    SUPP_VEC = supp_vec,
                    SVTYPE = can_record.type)
            if annotation_list != []:
                info_list += ';ANNOTATION=' + '|'.join(annotation_list)
            
        file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t".format(
            CHR = can_record.chrom1, 
            POS = can_record.start,
            ID = can_record.name,
            REF = can_record.ref,
            ALT = ','.join(can_record.alt), 
            QUAL = '.' if can_record.qual is None else can_record.qual,
            PASS = filter_lable,
            INFO = info_list, 
            FORMAT = "GT:ID:RAL:AAL:CO",
            ))
        for i in range(len(sample_ids)):
            if i in item[5]:
                file.write("%s:%s:%s:%s:"%(item[5][i].gt, item[5][i].name, item[5][i].ref, ','.join(item[5][i].alt)))
                file.write("%s_%d-%s_%d\t"%(item[5][i].chrom1, item[5][i].start, item[5][i].chrom2, item[5][i].end))
            else:
                file.write("./.:NAN:NAN:NAN:NAN\t")
        file.write('\n')
    file.close()


def generate_header(file, contiginfo, sample_ids):
    # General header
    file.write('##fileformat=VCFv4.2\n')
    file.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
    for i in contiginfo:
        file.write("##contig=<ID=%s,length=%d>\n"%(i, contiginfo[i]))

    # Specific header
    # ALT
    file.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n")
    file.write("##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n")
    file.write("##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n")
    file.write("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n")
    file.write("##ALT=<ID=BND,Description=\"Breakend of translocation\">\n")

    # INFO
    #file.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n")
    #file.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n")
    file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
    file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    file.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
    file.write("##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n")
    # file.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
    file.write("##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n")
    file.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
    #file.write("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n")
    file.write("##INFO=<ID=SUPP,Number=1,Type=String,Description=\"Number of samples supporting the variant\">\n")
    file.write("##INFO=<ID=SUPPID,Number=1,Type=String,Description=\"Samples supporting the variant\">\n")
    file.write("##FILTER=<ID=q5,Description=\"Quality below 5\">\n")
    
    # FORMAT
    file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    '''
    file.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
    file.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
    file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
    file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")
    '''
    file.write("##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Variant ID\">\n")
    file.write("##FORMAT=<ID=CO,Number=1,Type=String,Description=\"Sequence coordinates\">\n")
    file.write("##FORMAT=<ID=RAL,Number=1,Type=String,Description=\"Reference allele sequence\">\n")
    file.write("##FORMAT=<ID=AAL,Number=1,Type=String,Description=\"Alternative allele sequence\">\n")
    #file.write("##CommandLine=\"cuteSV %s\"\n"%(" ".join(argv)))
    file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for sample_id in sample_ids:
        file.write('\t' + sample_id)
    file.write('\n')