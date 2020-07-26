"""Module to test changes in haplotype frequency over the simulation process.
Step 1. Get Founding Haplotype frequencies in 100kb Windows at end of forqs simulation
Step 2. Check differences between actual.freqs and true.freqs in simreads
Step 3. Get HARP haplotype frequency estimates in 100kb Windows from the harp run on simreads output
ISSUE1 - Simreads is given frequency of constructed (recombined haplotype), not any DGRP haplotypes, thus simreads
true/actual freqs cannot be done in broken up windows"""
import sys
import subprocess
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import forqs_begin as fb
import common_functions as cfun


# For step 1 will need to make a
if __name__ == '__main__':
    chromo = '2L'
    mixed_haplotypes_file = 'dgrp{}_mixed_haplotypes.txt'.format(chromo)
    pos_range = ['4000000', '8000000']
    recombination_map = '{}rcc.txt'.format(chromo)
    recombination_map_subset = 'dmel_recRates_{}_{}-{}.csv'.format(chromo, pos_range[0], pos_range[1])
    simulation_number = 10
    chromosome_length = int(pos_range[1]) - int(pos_range[0])

    # prepare the recombination map and haplotype files for forqs simulation
    fb.prepare_simulations(mixed_haplotypes_file, chromo, pos_range[0], pos_range[1],
                           recombination_map)
    com = cfun.config_command(chromo, chromosome_length, recombination_map_subset, simulation_number, configtype='HAPS')
    p1 = subprocess.call(com, shell=False)
    com2 = cfun.forq_command()
    p2 = subprocess.call(com2, shell=False)
    cfun.move_configs(chromo)

    with open('haplotype_frequencies_chr1_final_pop1.txt') as f:
        data = [line for line in f]