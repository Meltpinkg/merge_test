from pysam import VariantFile
import sys

def main(vcf_file):
    supp_dict = dict()
    vcf_reader = VariantFile(vcf_file, 'r')
    for record in vcf_reader.fetch():
        supp_vec = record.info['SUPP_VEC']
        if supp_vec in supp_dict:
            supp_dict[supp_vec] = supp_dict[supp_vec] + 1
        else:
            supp_dict[supp_vec] = 1
    print(supp_dict['001'])
    print(supp_dict['011'])
    print(supp_dict['101'])
    print(supp_dict['111'])
    print(supp_dict['110'])
    print(supp_dict['010'])
    print(supp_dict['100'])

if __name__ == '__main__':
    main(sys.argv[1])