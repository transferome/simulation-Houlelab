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
    range_file = 'dgrp' + chromosome + '_rangesubset.txt'
    pos_list = list()
    # subtractor is used to adjust the SNP positions from simulations that don't start at 0
    # for example 1,000,000 to 2,000,000 range has first SNP as 1005816
    # but in the forqs population file, 1,005,816 does not exist, from the start of the simulated chromosome
    # that actually exists at 5816 in forqs positions, so have to subtract 1,000,000
    subtractor = int(lower_bound)
    with open(range_file) as f:
        for line in f:
            pos = int(line.split('\t')[0]) - subtractor
            pos_list.append(pos)
    return list(enumerate(pos_list))


def hap_dictionary(chromosome):
    """Takes the transposed dataframe and creates a dictionary, where key is DGRP# and value
    is the haplotype sequence"""
    transpose_filename = 'dgrp' + chromosome + '_mixed_transpose.txt'
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
            new_list.append(final_string)
        if s.startswith('- '):
            new_string = s.lstrip('- { ')
            final_string = new_string.rstrip(' }')
            new_list.append(final_string)
    return new_list


def index_getter(start_position, end_position, position_index_list):
    """Given an upper and lower bound for a position, this finds the index value of the maximum index value
    which pasts the equality testing"""
    index_min = [tup[0] for tup in position_index_list if tup[1] > start_position]
    index_max = [tup[0] for tup in position_index_list if tup[1] <= end_position]
    ind_min = min(index_min)
    ind_max = max(index_max)
    return str(ind_min), str(ind_max)


def parse_individual(popstring, position_index_list):
    """For an individual creates list of tuples, that include the position start and stop, and dgrp number"""
    indi = popstring.split(' ')
    tup_list = [make_tuple(s) for s in indi]
    positions = [tup[0] for tup in tup_list]
    dgrp_num = [tup[1] for tup in tup_list]
    if len(dgrp_num) > 1:
        positions_trim = positions[1:] + [positions[-1]]
        hap_guide = list(zip(dgrp_num, positions, positions_trim))
        out_list = list()
        for region in hap_guide[:-1]:
            try:
                minval, maxval = index_getter(region[1], region[2], position_index_list)
                new_tuple = (str(region[0]), int(minval), int(maxval))
                out_list.append(new_tuple)
            except ValueError:
                print('Value Error: Something wrong with ')
        max_possible_index_list = [tup[0] for tup in position_index_list]
        max_possible = max(max_possible_index_list)
        last_region = hap_guide[-1]
        fminval, fmaxval = index_getter(last_region[1], last_region[2], position_index_list)
        # Add 1 so there is not overlap with the previous region.
        ftuple = (str(last_region[0]), int(fminval), max_possible)
        out_list.append(ftuple)
        return out_list
    elif len(dgrp_num) == 1:
        out_list = list()
        max_possible_index_list = [tup[0] for tup in position_index_list]
        max_possible = max(max_possible_index_list)
        ftuple = (str(dgrp_num[0]), 0, max_possible)
        out_list.append(ftuple)
        return out_list
    else:
        print('Individual is of length less than 1')


def construct_haplotype(individual, hap_dict):
    """Given the output from parse_individual, and the haplotype_dictionary, this creates
    an individuals haplotype"""
    chunk = None
    if len(individual) > 1:
        chunk_list = list()
        for region in individual:
            try:
                haplotype = hap_dict[region[0]]
                # have to add one here and in the other similar code chunk, because splice not inclusive
                chunk = haplotype[region[1]:region[2] + 1]
                chunk_list.append(chunk)
            except KeyError:
                print('DGRP or Region Key error')
        return ''.join(chunk_list)
    if len(individual) == 1:
        for region in individual:
            try:
                haplotype = hap_dict[region[0]]
                chunk = haplotype[region[1]:region[2] + 1]
            except KeyError:
                print('DGRP or Region Key Error')
        return chunk
    if len(individual) is None:
        print('Error, Individual is of Length None')
        quit()


def construct_individuals(poplist, position_index_list, hap_dict, forqs_stem):
    """Given the population list, this constructs haplotypes for all individuals"""
    # stem is meaningless but multiprocess needs to map different inputs to a function
    stem = forqs_stem
    out_list = list()
    for line in poplist:
        individual = parse_individual(line, position_index_list)
        haplotype = construct_haplotype(individual, hap_dict)
        out_list.append(haplotype)
    return out_list


def make_index(haplotype_list):
    """Returns an Index List, so that haplotypes are paired maternal and paternal"""
    index_list = list()
    for x in range(0, len(haplotype_list)//2):
        rep_num = [x] * 2
        index_list.extend(rep_num)
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


def haplotype_freq_list(hap_cnt_dict, sample_number):
    """Gets frequency of each of the diploid individuals in the random sample.
    This will make the frequency list that simreads needs, this is the frequency of these in the pool, so haps not
    sampled do not get a freq of 0, because they were never poooled"""
    freq_dict = dict()
    for key in hap_cnt_dict.keys():
        count = hap_cnt_dict[key]
        freq = str(round(count / sample_number, 5))
        freq_dict[key] = freq
    freqs = list(freq_dict.values())
    freq_divide = [float(s) / 2 for s in freqs]
    final_freqs = list(itertools.chain.from_iterable(itertools.repeat(x, 2) for x in freq_divide))
    return final_freqs


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
    snp_file = os.path.abspath(forqs_stem + '_haplotypes.txt')
    # chromosome_l = int(chromosome_length) + 4000
    str_freqs_list = [str(x) for x in frequencies]
    string_freqs = ' '.join(str_freqs_list)
    filename_refseq = 'filename_refseq /home/dhoule/evoreseq/ref/dmel-majchr-norm-r6.24.fasta'
    filename_snps = 'filename_snps ' + snp_file
    filename_stem = 'filename_stem simreads_' + forqs_stem
    range_file = 'dgrp' + chromosome + '_rangesubset.txt'
    min_max = cfun.region_min_max(range_file)
    region = 'region 2L:' + str(min_max[0]) + '-' + str(min_max[1])
    haplotype_frequencies = 'haplotype_frequencies ' + string_freqs
    recombined_haplotypes = 'recombined haplotype frequencies'
    coverage = 'coverage 200'
    error_rate = 'error_rate 0.2'
    read_length = 'read_length 150'
    file_start = ['#', '#Usage:harp sim_reads', '#', '']
    file_addition = [filename_refseq, filename_snps, filename_stem, region, haplotype_frequencies,
                     recombined_haplotypes, coverage, error_rate, read_length]
    out_list = file_start + file_addition
    with open(configfilename, 'w+') as f:
        for line in out_list:
            f.write(line + '\n')


def pop2_snp(chromosome, lower_bound, forqs_stem):
    """Main Function Which Using Code From Previously Defined Functions and Converts the file, and creates
    config file"""
    ref_file = 'dgrp' + chromosome + '_subset.txt'
    # chr_len = chromosome_length
    pos_id_list = position_index(chromosome, lower_bound)
    hap_dic = hap_dictionary(chromosome)
    pop_filename = forqs_stem + '/population_final_pop1.txt'
    output_snptable = forqs_stem + '_haplotypes.txt'
    simreads_config_filename = 'simreads_' + forqs_stem + '.config'
    forqs_pop = read_popfile(pop_filename)
    population_haplotypes = construct_individuals(forqs_pop, pos_id_list, hap_dic, forqs_stem)
    indxes = make_index(population_haplotypes)
    rands = random_haplotypes(indxes)
    cnt_dict = haplotype_count_dict(rands)
    # currently setting sample size to 200
    freqies = haplotype_freq_list(cnt_dict, 200)
    final_haps = restructure_individuals(population_haplotypes)
    write_snp_table(ref_file, output_snptable, final_haps, chromosome)
    write_simread_config(simreads_config_filename, freqies, forqs_stem, chromosome)


def multi_conversion(chromosome, lower_bound, target_func):
    """Convert population txt files to snp tables and create simreads config in parallel"""
    stems = cfun.list_stems(chromosome)
    pool = Pool(18, cfun.limit_cpu)
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
