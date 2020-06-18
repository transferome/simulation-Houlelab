"""Multiprocessing of samtools view, sort, and index commands of the simulated reads"""
import subprocess
from multiprocessing import Pool
import argparse
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def samtools_view(sam_file):
    """Samtools view command and process"""
    bam_file = sam_file.split('.sam')[0] + '.bam'
    command_list = ['samtools', 'view', '-bhS', sam_file, '>', bam_file]
    command = ' '.join(command_list)
    subprocess.call(command, shell=True)


def samtools_view_multi(chromosome, target_func):
    """runs samtools view in parallel"""
    sams = cfun.list_sams(chromosome)
    pool = Pool(18, cfun.limit_cpu)
    pool.map(target_func, sams)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Simulated Bams to Sams')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    args = parser.parse_args()
    samtools_view_multi(args.chr, samtools_view)
