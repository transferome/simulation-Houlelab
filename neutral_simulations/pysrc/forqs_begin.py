"""Commands to Prepare Forqs Simulations"""
import numpy as np
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import convert_recmap as cnrec


def organize_data(file, chromosome, region_min, region_max):
    """Takes the mimhap file and removes a few columns that get in the way"""
    data_filename = 'dgrp{}_subset.txt'.format(chromosome)
    newdata_filename = 'dgrp{}_rangesubset.txt'.format(chromosome)
    data = list()
    new_data = list()
    with open(file) as f:
        filt_obj = list(filter(lambda x: int(region_min) <= int(x.split('\t')[1]) <= int(region_max), f))
        for line in filt_obj:
            data.append(line)
            linesp = line.split('\t')
            hap_info = linesp[4:]
            hap_info.insert(0, linesp[1])
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
    subset_filename = 'dgrp{}_rangesubset.txt'.format(chromosome)
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


def change_mix_name(mixed_haplotype_file):
    """Creates transpose filename"""
    new_filename = '{}_transpose.txt'.format(mixed_haplotype_file.split('_haplotypes')[0])
    return new_filename


def prepare_simulations(mixed_file, chr_arm, low_bound, high_bound, recmap_file):
    """Main function for forqs begin"""
    transpose_name = change_mix_name(mixed_file)
    organize_data(mixed_file, chr_arm, int(low_bound), int(high_bound))
    data_matrix(chr_arm, transpose_name)
    cnrec.main_recmap(recmap_file, chr_arm, [low_bound, high_bound])


if __name__ == '__main__':
    pass
