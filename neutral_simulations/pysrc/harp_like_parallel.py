"""Multiprocessing of harp like command"""
import subprocess
from multiprocessing import Pool
from functools import partial
import argparse
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def like_process(chromosome, sorted_bam):
    """Execute harp like command and process"""
    # ref file is always the same
    forqs_stem = '_'.join(sorted_bam.split('.')[0].split('_')[1:3])
    ref_file = '/home/dhoule/evoreseq/ref/dmel-majchr-norm-r6.24.fasta'
    # snp_text will be based off chromosome arm
    snp_text = '/home/dhoule/evoreseq/ref/' + chromosome + '_snp_pool.txt'
    # use rangesubset to get min max SNP positions
    range_text = 'dgrp' + chromosome + '_rangesubset.txt'
    min_max = cfun.region_min_max(range_text)
    region = chromosome + ':' + str(min_max[0]) + '-' + str(min_max[1])
    harp_like_command = ['harp', 'like', '-I', '--bam', sorted_bam, '--region', region, '--refseq', ref_file, '--snps',
                         snp_text, '--stem', forqs_stem]
    subprocess.call(harp_like_command, shell=False)


def like_multi(chromosome, target_func):
    """Runs harp like function in parallel"""
    sorted_bams = cfun.list_sorted_bams(chromosome)
    pool = Pool(18, cfun.limit_cpu)
    adjusted_func = partial(target_func, chromosome)
    pool.map(adjusted_func, sorted_bams)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run harp like in parallel')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    args = parser.parse_args()
    like_multi(args.chr, like_process)
