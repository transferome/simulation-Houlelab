"""Multiprocessing of Index Snp Table program from harp-simreads"""
import glob
import subprocess
from multiprocessing import Pool
import argparse
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def list_snp_files(chromosome):
    """Lists all the snp files in """
    return glob.glob('chr{}_*haplotypes.txt'.format(chromosome))


def index_snp_texts(snp_text):
    """index snp table process"""
    idx_com = '/usr/local/bin/index_snp_table'
    line_set = '10000'
    command = [idx_com, snp_text, line_set]
    subprocess.call(command, shell=False)


def index_snp_multi(chromosome, target_func):
    """Index multiple SNP files at same time """
    stems = list_snp_files(chromosome)
    pool = Pool(21, cfun.limit_cpu)
    pool.map(target_func, stems)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Index Snp Table in Parallel')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    args = parser.parse_args()
    index_snp_multi(args.chr, index_snp_texts)
