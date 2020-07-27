"""Main Forqs Simulation Script"""
import sys
import os
import subprocess
import shutil
import time
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import forqs_begin as fb
import common_functions as cfun


def harp_like_cleanup(chromosome):
    """Cleans up files created by forqs and harp_like"""
    # remove forqs.seed file
    try:
        os.remove('forqs.seed')
    except OSError:
        pass
    forqs_dirs = cfun.list_forqs_directories(chromosome)
    haps = cfun.list_forqs_haplotypes(chromosome)
    idxs = cfun.list_forqs_idx(chromosome)
    bams = cfun.list_sorted_bams(chromosome)
    bais = cfun.list_sorted_bams_idx(chromosome)
    haps_sorted = sorted(haps, key=lambda x: int(x.split('_haplotypes')[0].split('_')[-1]))
    idxs_sorted = sorted(idxs, key=lambda x: int(x.split('_haplotypes')[0].split('_')[-1]))
    bams_sorted = sorted(bams, key=lambda x: int(x.split('_sorted.bam')[0].split('_')[-1]))
    bai_sorted = sorted(bais, key=lambda x: int(x.split('_sorted.bam')[0].split('_')[-1]))
    len_check = [len(haps_sorted), len(idxs_sorted), len(bams_sorted), len(bai_sorted)]
    if all(len(forqs_dirs) == x for x in len_check):
        for hap, idx, bam, bai, forq_dir in zip(haps_sorted, idxs_sorted, bams_sorted, bai_sorted, forqs_dirs):
            shutil.move(hap, forq_dir)
            shutil.move(idx, forq_dir)
            shutil.move(bam, forq_dir)
            shutil.move(bai, forq_dir)
    else:
        print('File lists are of unequal length')


def harp_freq(chr_arm, step, width):
    """Runs harp freq in parallel"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/harp_freq_parallel.py', '--chr', chr_arm,
               '--step', step, '--width', width]
    return command


def remove_hlk(chr_arm):
    """Removes the hlk files, they are large"""
    hlks = cfun.list_hlk(chr_arm)
    for hlk in hlks:
        os.remove(hlk)


if __name__ == '__main__':
    start_time = time.time()
    # Set initial parameters and filenames
    chromo = '2L'
    mixed_haplotypes_file = 'dgrp{}_mixed_haplotypes.txt'.format(chromo)
    pos_range = ['8000000', '10000000']
    recombination_map = '{}rcc.txt'.format(chromo)
    recombination_map_subset = 'dmel_recRates_{}_{}-{}.csv'.format(chromo, pos_range[0], pos_range[1])
    simulation_number = 10
    chromosome_length = int(pos_range[1]) - int(pos_range[0])

    # prepare the recombination map and haplotype files for forqs simulation
    fb.prepare_simulations(mixed_haplotypes_file, chromo, pos_range[0], pos_range[1],
                           recombination_map)
    com = cfun.config_command(chromo, chromosome_length, recombination_map_subset, simulation_number)
    p1 = subprocess.call(com, shell=False)
    com2 = cfun.forq_command()
    p2 = subprocess.call(com2, shell=False)
    cfun.move_configs(chromo)
    com3 = cfun.convert_pops2snps(chromo, pos_range[0])
    p3 = subprocess.call(com3, shell=False)
    com4 = cfun.index_tables(chromo)
    p4 = subprocess.call(com4, shell=False)
    com5 = cfun.simreads_com()
    p5 = subprocess.call(com5, shell=False)
    cfun.simreads_cleanup(chromo)
    com6 = cfun.sam2bam(chromo)
    p6 = subprocess.call(com6, shell=False)
    com7 = cfun.bam_sort(chromo)
    p7 = subprocess.call(com7, shell=False)
    com8 = cfun.sambam_remove(chromo)
    p8 = subprocess.call(com8, shell=False)
    com9 = cfun.bam_index(chromo)
    p9 = subprocess.call(com9, shell=False)
    com10 = cfun.harp_like(chromo)
    p10 = subprocess.call(com10, shell=False)
    cfun.harp_like_cleanup(chromo)
    com11 = cfun.harp_freq(chromo, '50000', '100000')
    p11 = subprocess.call(com11, shell=False)
    com12 = cfun.harp_freq(chromo, '50000', '50000')
    p12 = subprocess.call(com12, shell=False)
    cfun.remove_hlk(chromo)
    print("--- %s seconds ---" % (time.time() - start_time))
