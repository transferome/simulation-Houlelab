"""Commands to Prepare Forqs Simulations"""
import numpy as np
import sys
sys.path.append('/home/ltjones/forql/chromosome2L/pysrc')
import convert_recmap as cnrec
# import argparse


def organize_data(file, chromosome, region_min, region_max):
    """Takes the mimhap file and removes a few columns that get in the way"""
    data_filename = 'dgrp' + chromosome + '_subset.txt'
    newdata_filename = 'dgrp' + chromosome + '_rangesubset.txt'
    data = list()
    new_data = list()
    with open(file) as f:
        for line in f:
            pos = line.split('\t')[1]
            if region_min <= int(pos) <= region_max:
                # pos = line.split('\t')[1]
                data.append(line)
                hap_info = line.split('\t')[4:]
                hap_info.insert(0, pos)
                new_data.append(hap_info)
    with open(newdata_filename, 'w+') as f:
        for line in new_data:
            line = '\t'.join(line)
            f.write(line)
    with open(data_filename, 'w+') as f:
        for line in data:
            f.write(line)
    return newdata_filename


def data_matrix(chromosome, transpose_filename):
    """Takes the few column removed mimhap file, turns into matrix and transposes it.
    This makes it so that each row is a DGRP# and its haplotype sequence.  Writes it as a file"""
    # TODO: 107 is based off having 106 DGRP regionsubset_file
    subset_filename = 'dgrp' + chromosome + '_rangesubset.txt'
    subset_data = list()
    with open(subset_filename) as f:
        for line in f:
            line = line.rstrip('\n')
            subset_data.append(line.split('\t'))
    header_numbers = [str(x) for x in range(0, 106)]
    header_numbers.insert(0, 'position')
    subset_data.insert(0, header_numbers)
    hap_array = np.array(subset_data)
    hap_transpose = hap_array.transpose()
    haplist = hap_transpose.tolist()
    new_data_frame = ['\t'.join(lst) + '\n' for lst in haplist]
    with open(transpose_filename, 'w+') as f:
        for line in new_data_frame:
            f.write(line)


# def create_recomb_map(recombination_file, chromosome_range):
#     """Create recombination file based off how long chromosome will be"""
#     stem_name = recombination_file.split('rcc.txt')[0]
#     newfilename = stem_name + '_' + str(chromosome_range[0]) + '-' + str(chromosome_range[1]) + '.csv'
#     with open(recombination_file) as f, open(newfilename, 'w+') as f1:
#         for line in f:
#             if line.startswith('position'):
#                 f1.write(line)
#             else:
#                 if int(chromosome_range[0]) <= int(line.split(' ')[0]) <= int(chromosome_range[1]):
#                     f1.write(line)


def change_mix_name(mixed_haplotype_file):
    """Creates transpose filename"""
    new_filename = mixed_haplotype_file.split('_haplotypes')[0] + '_transpose.txt'
    return new_filename


def prepare_simulations(mixed_file, chr_arm, low_bound, high_bound, chromosome_length, recmap_file):
    """Main function for forqs begin"""
    transpose_name = change_mix_name(mixed_file)
    organize_data(mixed_file, chr_arm, int(low_bound), int(high_bound))
    data_matrix(chr_arm, transpose_name)
    cnrec.main_recmap(recmap_file, chr_arm, [low_bound, high_bound], chromosome_length)


if __name__ == '__main__':
    pass
