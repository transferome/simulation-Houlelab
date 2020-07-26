"""Multiprocessing of samtools sort command on converted simulation bam files"""
import subprocess
from multiprocessing import Pool
import argparse
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def samtools_sort(bamfile):
    """Samtools sort command and process"""
    sorted_bam_file = bamfile.split('.bam')[0] + '_sorted.bam'
    command_list = ['samtools', 'sort', '-@4', bamfile, '-o', sorted_bam_file]
    command = ' '.join(command_list)
    subprocess.call(command, shell=True)


def samtools_sort_multi(chromosome, target_func):
    """runs samtools view in parallel"""
    bams = cfun.list_bams(chromosome)
    pool = Pool(21, cfun.limit_cpu)
    pool.map(target_func, bams)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sort Simulated Bams')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    args = parser.parse_args()
    samtools_sort_multi(args.chr, samtools_sort)
