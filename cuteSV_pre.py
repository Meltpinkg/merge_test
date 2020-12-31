def output_result(semi_result, samples_id, output_file, contiginfo):
    file = open(output_file, 'w')
    generate_header(file, contiginfo)
    #print('start writing to file ' + chrom)
    for item in semi_result:  # [CHROM, POS, CANDIDATE_RECORD, CIPOS, CIEND, DICT]
        supp_id = ''
        for i in range(len(samples_id)):
            if i in item[5]:
                supp_id += '1'
            else:
                supp_id += '0'
        can_record = item[2]
        if can_record.qual == "." or can_record.qual is None:
            filter_lable = "PASS"
        else:
            filter_lable = "PASS" if float(can_record.qual) >= 5.0 else "q5"
        if can_record.type == 'INS':
            sv_len = can_record.end
            sv_end = can_record.start
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
            sv_len = -can_record.end
            sv_end = can_record.start + can_record.end
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
    file.close()


def generate_header(file, contiginfo):
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
    '''
    file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    file.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
    file.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
    file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
    file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")
    '''
    file.write("##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Variant ID from input.\">\n")
    file.write("##FORMAT=<ID=RAL,Number=1,Type=String,Description=\"Reference allele sequence reported from input.\">\n")
    file.write("##FORMAT=<ID=AAL,Number=1,Type=String,Description=\"Alternative allele sequence reported from input.\">\n")
    file.write("##FORMAT=<ID=CO,Number=1,Type=String,Description=\"Coordinates\">\n")
    #file.write("##CommandLine=\"cuteSV %s\"\n"%(" ".join(argv)))