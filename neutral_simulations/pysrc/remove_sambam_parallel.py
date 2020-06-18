"""Removes old sams, and bams in parallel"""
import os
from multiprocessing import Pool
import argparse
import sys
sys.path.append('/home/ltjones/forql/chromosome2L/pysrc')
import common_functions as cfun


def remover(filename):
    """deletes filename from folder"""
    os.remove(filename)


def remove_sambam(chromosome, target_func):
    """removes old sam and bam files in parallel"""
    sams = cfun.list_sams(chromosome)
    bams = cfun.list_bams(chromosome)
    bams_old = [s for s in bams if not s.endswith('_sorted.bam')]
    remove_files = sams + bams_old
    pool = Pool(18, cfun.limit_cpu)
    pool.map(target_func, remove_files)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove old bam and sam files')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    args = parser.parse_args()
    remove_sambam(args.chr, remover)
