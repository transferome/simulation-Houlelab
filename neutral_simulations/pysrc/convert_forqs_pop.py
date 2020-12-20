"""Multiprocessing Module for Converting Multiple forqs Population Haplotype Files Into Snp tables for Simreads"""
import numpy as np
from ast import literal_eval as make_tuple
import itertools
from functools import partial
import os
from multiprocessing import Pool
import argparse
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def position_index(chromosome, lower_bound):
    """Gets list of positions, and indexes them with enumerate"""
    range_file = 'dgrp{}_rangesubset.txt'.format(chromosome)
    # pos_list = list()
    # subtractor is used to adjust the SNP positions from simulations that don't start at 0
    # for example 1,000,000 to 2,000,000 range has first SNP as 1005816
    # but in the forqs population file, 1,005,816 does not exist, from the start of the simulated chromosome
    # that actually exists at 5816 in forqs positions, so have to subtract 1,000,000
    subtractor = int(lower_bound)
    with open(range_file) as f:
        pos_list = [int(line.split('\t')[0]) - subtractor for line in f]
        # for line in f:
        #     pos = int(line.split('\t')[0]) - subtractor
        #     pos_list.append(pos)
    return list(enumerate(pos_list))


def hap_dictionary(chromosome):
    """Takes the transposed dataframe and creates a dictionary, where key is DGRP# and value
    is the haplotype sequence"""
    transpose_filename = 'dgrp{}_mixed_transpose.txt'.format(chromosome)
    data = [line.rstrip('\n') for line in open(transpose_filename)][1:]
    # TODO: for some reason forqs changed haplotype id's from 0 - 105, so all DGRP numbers need 1 number subtracted
    data_dict = {str(line.split('\t')[0]): ''.join(line.split('\t')[1:]) for line in data}
    return data_dict


def read_popfile(popfilename):
    """Reads in the popfile into usable data"""
    pop_file = [line.rstrip('\n') for line in open(popfilename)][2:]
    filtered_pop = list(filter(None, pop_file))
    new_list = list()
    for s in filtered_pop:
        if s.startswith('+ '):
            new_string = s.lstrip('+ { ')
            final_string = new_string.rstrip(' }')
            indi = final_string.split(' ')
            chunks = [make_tuple(s) for s in indi]
            new_list.append(chunks)
        if s.startswith('- '):
            new_string = s.lstrip('- { ')
            final_string = new_string.rstrip(' }')
            indi = final_string.split(' ')
            chunks = [make_tuple(s) for s in indi]
            new_list.append(chunks)
    return new_list


def index_getter(start_position, end_position, position_index_list):
    """Given an upper and lower bound for a position, this finds the index value of the maximum index value
    which pasts the equality testing"""
    index_min = [tup[0] for tup in position_index_list if tup[1] >= start_position]
    index_max = [tup[0] for tup in position_index_list if tup[1] <= end_position]
    ind_min = min(index_min)
    ind_max = max(index_max)
    return ind_min, ind_max


class ForqsHaplotype:
    """Forqs Haplotype Object
    Example [(0, 91), (1235052, 105), (1251805, 105)]
    or  [(0, 91)]
    """

    def __init__(self, forq_haplotype):
        """Forqs Haplotype Instance Creation"""
        self.info = forq_haplotype

        # tag object depending on it's structure
        self.length = len(self.info)
        self.founders = [info[1] for info in self.info]
        self.founders_set = list(set(self.founders))
        self.founder_number = len(self.founders_set)
        self.positions = [info[0] for info in self.info]
        self.struct = 'normal'
        if self.length == 1:
            self.struct = 'single'
        if self.length != self.founder_number:
            self.struct = 'contain duplicate'
        if self.length == 2 and self.founder_number == 1:
            self.struct = 'single recombined'
        self.position_map_offset = None
        self.temp_map = None
        self.map = list()

    def make_haplotype_map(self, position_index_list):
        """Create the haplotype guide for assembling the recombined forqs haplotype outputs"""
        position_indexes = [info[0] for info in position_index_list]
        position_values = [info[1] for info in position_index_list]
        min_index = min(position_indexes)
        max_index = max(position_indexes)
        if self.struct == 'single' or self.struct == 'single recombined':
            self.map.append([self.founders_set[0], min_index, max_index])
        if self.struct == 'normal' or self.struct == 'contain duplicate':
            if self.positions[1] < min(position_values):
                self.founders.pop(0)
                self.positions.pop(1)
            try:
                if self.positions[2] < min(position_values):
                    sys.exit("Two recombination points before the first SNP occurs")
            except IndexError:
                pass
            self.position_map_offset = self.positions[1:] + [max(position_values)]
            self.temp_map = list(zip(self.positions, self.position_map_offset))
            for hap in self.temp_map:
                min_id, max_id = index_getter(hap[0], hap[1], position_index_list)
                self.map.append([min_id, max_id])
            for i, hap_map in enumerate(self.map):
                if i < len(self.map) - 1:
                    current_end = hap_map[1]
                    next_start = self.map[i + 1][0]
                    # print(current_end, next_start)
                    if current_end == next_start:
                        self.map[i + 1][0] = int(next_start) + 1
            for line_id, hap_map in zip(self.founders, self.map):
                hap_map.insert(0, line_id)


def construct_haplotype(forq_haplotype_object, hap_dict):
    """Given the output from parse_individual, and the haplotype_dictionary, this creates
    an individuals haplotype"""
    if forq_haplotype_object.struct == 'single' or forq_haplotype_object.struct == 'single recombined':
        haplotype = hap_dict[str(forq_haplotype_object.founders_set[0])]
        sequence_chunk = haplotype[forq_haplotype_object.map[0][1]:forq_haplotype_object.map[0][2] + 1]
        return sequence_chunk
    if forq_haplotype_object.struct == 'normal' or forq_haplotype_object.struct == 'contain duplicate':
        sequence_list = list()
        for sub_map in forq_haplotype_object.map:
            haplotype = hap_dict[str(sub_map[0])]
            sequence_list.append(haplotype[sub_map[1]:sub_map[2] + 1])
        return ''.join(sequence_list)


def construct_individuals(poplist, position_index_list, hap_dict, forqs_stem):
    """Given the population list, this constructs haplotypes for all individuals"""
    # stem is meaningless but multiprocess needs to map different inputs to a function
    stem = forqs_stem
    out_list = list()
    for line in poplist:
        individual = ForqsHaplotype(line)
        individual.make_haplotype_map(position_index_list)
        out_list.append(construct_haplotype(individual, hap_dict))
    return out_list


def restructure_individuals(haplotypes):
    """Reorganizes so that haplotypes are columns instead of rows"""
    out_list = list()
    new_structure = list()
    for hap in haplotypes:
        new_hap = list(hap)
        new_structure.append(new_hap)
    hap_array = np.array(new_structure)
    hap_trans = hap_array.transpose()
    new_hap_list = hap_trans.tolist()
    for lst in new_hap_list:
        new_line = ','.join(lst)
        out_list.append(new_line)
    return out_list


def make_index_pair(haplotype_list):
    """Returns an Index List, so that haplotypes are paired maternal and paternal"""
    index_list = list()
    for x in range(0, len(haplotype_list)//2):
        index_list.extend([x] * 2)
    index_list = [str(x) for x in index_list]
    return index_list


def random_haplotypes(index_list):
    """Draws random sample of indexes (diploid individuals), given the sample number to sample"""
    hap_indexes = list(set(index_list))
    # rand_haps_indexes = random.sample(hap_indexes, k=sample_number)
    hap_indexes.sort(key=int)
    return hap_indexes


def haplotype_count_dict(random_haps):
    """Counts how many times each of the diploid individuals was sampled"""
    all_freq = dict()
    for s in random_haps:
        if s in all_freq:
            all_freq[s] += 1
        else:
            all_freq[s] = 1
    return all_freq


def haplotype_freq_list(hap_cnt_dict):
    """Gets frequency of each of the diploid individuals in the random sample.
    This will make the frequency list that simreads needs, this is the frequency of these in the pool, so haps not
    sampled do not get a freq of 0, because they were never poooled"""
    freq_dict = dict()
    for key in hap_cnt_dict.keys():
        count = hap_cnt_dict[key]
        freq = str(round(count / 200, 5))
        freq_dict[key] = freq
    freqs = list(freq_dict.values())
    freq_divide = [float(s) / 2 for s in freqs]
    final_freqs = list(itertools.chain.from_iterable(itertools.repeat(x, 2) for x in freq_divide))
    return final_freqs


def write_snp_table(hapfile, newfilename, final_haplotypes, chromosome):
    """Write the snp_table needed for simreads"""
    indi_nums = range(1, len(final_haplotypes[0].split(',')) + 1)
    indi_strs = [str(x) for x in indi_nums]
    indi_line = ','.join(indi_strs)
    header = ','.join([chromosome, 'Ref', indi_line])
    with open(hapfile) as f, open(newfilename, 'w+') as fo:
        fo.write(header)
        fo.write('\n')
        for org, new in zip(f, final_haplotypes):
            org_info = ','.join(org.split('\t')[1:3])
            new_line = ','.join([org_info, new])
            fo.write(new_line)
            fo.write('\n')


def write_simread_config(configfilename, frequencies, forqs_stem, chromosome):
    """Writes out the information for the simreads config file"""
    # TODO: not entirely sure how to match chromosome length and chromosome region.  Giving simreads region 0 - 1000000
    #  gives error
    snp_file = os.path.abspath('{}_haplotypes.txt'.format(forqs_stem))
    # chromosome_l = int(chromosome_length) + 4000
    str_freqs_list = [str(x) for x in frequencies]
    string_freqs = ' '.join(str_freqs_list)
    range_file = 'dgrp{}_rangesubset.txt'.format(chromosome)
    min_max = cfun.region_min_max(range_file)
    file_lines = ['#', '#Usage:harp sim_reads', '#', '',
                  'filename_refseq /home/dhoule/evoreseq/ref/dmel-majchr-norm-r6.24.fasta',
                  'filename_snps {}'.format(snp_file), 'filename_stem simreads_{}'.format(forqs_stem),
                  'region {}:{}-{}'.format(chromosome, str(min_max[0]), str(min_max[1])),
                  'haplotype_frequencies {}'.format(string_freqs), 'recombined haplotype frequencies',
                  'coverage 200', 'error_rate 0.2', 'read_length 150']
    with open(configfilename, 'w+') as f:
        for line in file_lines:
            f.write('{}\n'.format(line))

# def write_simread_config(configfilename, frequencies, forqs_stem, chromosome):
#     """Writes out the information for the simreads config file"""
#     # TODO: not entirely sure how to match chromosome length and chromosome region.  Giving simreads region 0 - 1000000
#     #  gives error
#     snp_file = os.path.abspath(forqs_stem + '_haplotypes.txt')
#     # chromosome_l = int(chromosome_length) + 4000
#     str_freqs_list = [str(x) for x in frequencies]
#     string_freqs = ' '.join(str_freqs_list)
#     filename_refseq = 'filename_refseq /home/dhoule/evoreseq/ref/dmel-majchr-norm-r6.24.fasta'
#     filename_snps = 'filename_snps {}'.format(snp_file)
#     filename_stem = 'filename_stem simreads_{}'.format(forqs_stem)
#     range_file = 'dgrp{}_rangesubset.txt'.format(chromosome)
#     min_max = cfun.region_min_max(range_file)
#     region = 'region {}:{}-{}'.format(chromosome, str(min_max[0]), str(min_max[1]))
#     haplotype_frequencies = 'haplotype_frequencies {}'.format(string_freqs)
#     recombined_haplotypes = 'recombined haplotype frequencies'
#     coverage = 'coverage 200'
#     error_rate = 'error_rate 0.2'
#     read_length = 'read_length 150'
#     file_lines = ['#', '#Usage:harp sim_reads', '#', '', filename_refseq, filename_snps, filename_stem, region,
#                   haplotype_frequencies, recombined_haplotypes, coverage, error_rate, read_length]
#     with open(configfilename, 'w+') as f:
#         for line in file_lines:
#             f.write('{}\n'.format(line))


def pop2_snp(chromosome, lower_bound, forqs_stem):
    """Main Function Which Using Code From Previously Defined Functions and Converts the file, and creates
    config file"""
    ref_file = 'dgrp{}_subset.txt'.format(chromosome)
    # chr_len = chromosome_length
    pos_id_list = position_index(chromosome, lower_bound)
    hap_dic = hap_dictionary(chromosome)
    pop_filename = '{}/population_final_pop1.txt'.format(forqs_stem)
    output_snptable = '{}_haplotypes.txt'.format(forqs_stem)
    simreads_config_filename = 'simreads_{}.config'.format(forqs_stem)
    forqs_pop = read_popfile(pop_filename)
    population_haplotypes = construct_individuals(forqs_pop, pos_id_list, hap_dic, forqs_stem)
    indxes = make_index_pair(population_haplotypes)
    rands = random_haplotypes(indxes)
    cnt_dict = haplotype_count_dict(rands)
    # currently setting sample size to 200
    freqies = haplotype_freq_list(cnt_dict)
    final_haps = restructure_individuals(population_haplotypes)
    write_snp_table(ref_file, output_snptable, final_haps, chromosome)
    write_simread_config(simreads_config_filename, freqies, forqs_stem, chromosome)


def multi_conversion(chromosome, lower_bound, target_func):
    """Convert population txt files to snp tables and create simreads config in parallel"""
    stems = cfun.list_stems(chromosome)
    pool = Pool(21, cfun.limit_cpu)
    target = partial(target_func, chromosome, lower_bound)
    pool.map(target, stems)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts forqs population outputs to snp tables')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    parser.add_argument('--low_bound', type=str, help='lower bound for position range')
    args = parser.parse_args()
    multi_conversion(args.chr, args.low_bound, pop2_snp)
