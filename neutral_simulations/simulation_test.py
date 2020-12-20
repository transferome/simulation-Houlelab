"""Module to test changes in haplotype frequency over the simulation process.
Step 1. Get Founding Haplotype frequencies in 100kb Windows at end of forqs simulation
Step 2. Check differences between actual.freqs and true.freqs in simreads
Step 3. Get HARP haplotype frequency estimates in 100kb Windows from the harp run on simreads output
ISSUE1 - Simreads is given frequency of constructed (recombined haplotype), not any DGRP haplotypes, thus simreads
true/actual freqs cannot be done in broken up windows"""
import os
import glob
import shutil
import sys
import subprocess
import time
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


def main(chromo, mixed_haplotypes_file, minimum_position, maximum_position, recombination_map,
         chromosome_length, recombination_map_subset, simulation_number):
    start_time = time.time()
    fb.prepare_simulations(mixed_haplotypes_file, chromo, minimum_position, maximum_position,
                           recombination_map)
    com = cfun.config_command(chromo, chromosome_length, recombination_map_subset, simulation_number, configtype='HAPS')
    subprocess.call(com, shell=False)
    com2 = cfun.forq_command()
    subprocess.call(com2, shell=False)
    cfun.move_configs(chromo)
    forqs_frequency_outputs(chromo)
    com3 = cfun.convert_pops2snps(chromo, minimum_position)
    subprocess.call(com3, shell=False)
    com4 = cfun.index_tables(chromo)
    subprocess.call(com4, shell=False)
    com5 = cfun.simreads_com()
    subprocess.call(com5, shell=False)
    cfun.simreads_cleanup(chromo)
    com6 = cfun.sam2bam(chromo)
    subprocess.call(com6, shell=False)
    com7 = cfun.bam_sort(chromo)
    subprocess.call(com7, shell=False)
    com8 = cfun.sambam_remove(chromo)
    subprocess.call(com8, shell=False)
    com9 = cfun.bam_index(chromo)
    subprocess.call(com9, shell=False)
    com10 = cfun.harp_like(chromo)
    subprocess.call(com10, shell=False)
    cfun.harp_like_cleanup(chromo)
    com11 = cfun.harp_freq(chromo, '25000', '50000')
    subprocess.call(com11, shell=False)
    com12 = cfun.harp_freq(chromo, '50000', '50000')
    subprocess.call(com12, shell=False)
    cfun.remove_hlk(chromo)
    total_time = "--- %s seconds ---\n" % (time.time() - start_time)
    with open('time_report.txt', 'a') as output_file:
        output_file.write(total_time)


def clean(chromosome):
    """Clean function keeps all .py files and files append to files_keep, deletes all other
    folders and files in the directory in which the script is within"""
    files_keep = glob.glob('*.py')
    files_keep.append('{}rcc.txt'.format(chromosome))
    files_keep.append('dgrp{}_mixed_haplotypes.txt'.format(chromosome))
    files_keep.append('time_report.txt')
    folder_contents = os.listdir()
    content_delete = [content for content in folder_contents if content not in files_keep]
    for content in content_delete:
        if os.path.isfile(content):
            os.remove(content)
        else:
            shutil.rmtree(content)
    print('Folder Cleaned')


if __name__ == '__main__':
    contig = '2L'
    mix_haplo_file = 'dgrp{}_mixed_haplotypes.txt'.format(contig)
    nuc_range = ['4000033', '8000000']
    recomb_map = '{}rcc.txt'.format(contig)
    recomb_map_subset = 'dmel_recRates_{}_{}-{}.csv'.format(contig, nuc_range[0], nuc_range[1])
    simu_number = 2
    contig_length = int(nuc_range[1]) - int(nuc_range[0])
    for simulator in range(0, 2):
        main(contig, mix_haplo_file, nuc_range[0], nuc_range[1], recomb_map, contig_length, recomb_map_subset, 2)
        clean(contig)
