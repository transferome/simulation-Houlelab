"""Finds the Minimum and Maximum Position of SNPs within a SNP.txt file
or a haplotypes.txt file"""


def region_min_max(snps_file, header=False):
    """ Finds the minimum and maximum coordinate position within the snp_txt file """
    # [1:] removes header info which is the chromosome
    if header:
        snp_positions = [s.split('\t')[0] for s in open(snps_file)][1:]
    else:
        snp_positions = [s.split('\t')[0] for s in open(snps_file)]
    position_integer = [int(s) for s in snp_positions]
    return [min(position_integer), max(position_integer)]


if __name__ == '__main__':
    test = region_min_max('snp_file.txt')
