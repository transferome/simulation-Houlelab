"""Multiprocessing of samtools index command"""
import subprocess
from multiprocessing import Pool
import argparse
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def samtools_index(bam_file):
    """Samtools index command and process"""
    command_list = ['samtools', 'index', bam_file]
    command = ' '.join(command_list)
    subprocess.call(command, shell=True)


def samtools_index_multi(chromosome, target_func):
    """runs samtools index in parallel"""
    bams = cfun.list_bams(chromosome)
    pool = Pool(18, cfun.limit_cpu)
    pool.map(target_func, bams)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Index Sorted Simulated Bam files')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    args = parser.parse_args()
    samtools_index_multi(args.chr, samtools_index)
