"""Module to test changes in haplotype frequency over the simulation process.
Step 1. Get Founding Haplotype frequencies in 100kb Windows at end of forqs simulation
Step 2. Check differences between actual.freqs and true.freqs in simreads
Step 3. Get HARP haplotype frequency estimates in 100kb Windows from the harp run on simreads output
ISSUE1 - Simreads is given frequency of constructed (recombined haplotype), not any DGRP haplotypes, thus simreads
true/actual freqs cannot be done in broken up windows"""
import os
import sys
import subprocess
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import forqs_begin as fb
import common_functions as cfun


def forqs_frequency_outputs(chromosome):
    """Gathers the forqs final frequency outputs and renames them appropriately"""
    os.mkdir('forqs_final_freqs')
    stems = cfun.list_stems(chromosome)
    for stem in stems:
        hapfreqs_filename = '{}/haplotype_frequencies_chr1_final_pop1.txt'.format(stem)
        newfilename = 'forqs_final_freqs/{}'.format(hapfreqs_filename.replace('chr1', stem).split('/')[1])
        os.rename(hapfreqs_filename, newfilename)


if __name__ == '__main__':
    chromo = '2L'
    mixed_haplotypes_file = 'dgrp{}_mixed_haplotypes.txt'.format(chromo)
    pos_range = ['4000000', '8000000']
    recombination_map = '{}rcc.txt'.format(chromo)
    recombination_map_subset = 'dmel_recRates_{}_{}-{}.csv'.format(chromo, pos_range[0], pos_range[1])
    simulation_number = 1000
    chromosome_length = int(pos_range[1]) - int(pos_range[0])

    # prepare the recombination map and haplotype files for forqs simulation
    fb.prepare_simulations(mixed_haplotypes_file, chromo, pos_range[0], pos_range[1],
                           recombination_map)
    com = cfun.config_command(chromo, chromosome_length, recombination_map_subset, simulation_number, configtype='HAPS')
    p1 = subprocess.call(com, shell=False)
    com2 = cfun.forq_command()
    p2 = subprocess.call(com2, shell=False)
    cfun.move_configs(chromo)
    forqs_frequency_outputs(chromo)
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
