from pysam import VariantFile
import sys

#gt: record.samples[i]['GT'] (tuple)
def parse_gt(gt):
    if gt == None or gt[0] == None:
        return 0
    if gt[0] + gt[1] == 0:
        return 0
    if gt[0] + gt[1] == 1:
        return 1
    if gt[0] + gt[1] == 2:
        return 2

# False -> inconsistency
def check(son, fa, ma):
    if fa == 0 and ma == 0:
        if son == 1 or son == 2:
            return False
    if fa == 0 and ma == 1 and son == 2:
        return False
    if fa == 1 and ma == 0 and son == 2:
        return False
    if fa == 0 and ma == 2:
        if son == 0 or son == 2:
            return False
    if fa == 2 and ma == 0:
        if son == 0 or son == 2:
            return False
    if fa == 1 and ma == 2 and son == 0:
        return False
    if fa == 2 and ma == 1 and son == 0:
        return False
    if fa == 2 and ma == 2:
        if son == 0 or son == 1:
            return False
    return True

def main(vcf_file):
    vcf_reader = VariantFile(vcf_file, 'r')
    '''
    sample1 = vcf_reader.header.samples[0]
    sample2 = vcf_reader.header.samples[1]
    sample3 = vcf_reader.header.samples[2]
    '''
    cnt = 0
    tot = 0
    for record in vcf_reader.fetch():
        tot += 1
        #supp_vec = record.info['SUPP_VEC']
        son = parse_gt(record.samples[0]['GT'])
        father = parse_gt(record.samples[1]['GT'])
        mother = parse_gt(record.samples[2]['GT'])
        res = check(son, father, mother)
        if res == False:
            cnt += 1
    print(cnt)
    print(tot)
    return cnt, tot

if __name__ == '__main__':
    main(sys.argv[1])