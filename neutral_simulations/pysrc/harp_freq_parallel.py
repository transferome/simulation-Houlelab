"""Multiprocessing of harp freq command"""
import subprocess
from multiprocessing import Pool
from functools import partial
import argparse
import sys
import os
import shutil
sys.path.append('/home/ltjones/forql/chromosome2L/pysrc')
import common_functions as cfun


def freq_process(chromosome, step, width, hlk_file):
    """Executes harps freq command"""
    range_text = 'dgrp' + chromosome + '_rangesubset.txt'
    min_max = cfun.region_min_max(range_text)
    region = chromosome + ':' + str(min_max[0]) + '-' + str(min_max[1])
    harp_freq_command = ['harp', 'freq', '--hlk', hlk_file, '--region', region, '--window_step', step,
                         '--window_width', width]
    subprocess.check_call(harp_freq_command, shell=False)


def freq_multi(chromosome, step, width, target_func):
    """runs harp freq in parallel multiprocess"""
    hlks = cfun.list_hlk(chromosome)
    pool = Pool(18, cfun.limit_cpu)
    adjusted_func = partial(target_func, chromosome, step, width)
    pool.map(adjusted_func, hlks)
    pool.close()
    pool.join()


def renamer(chromosome, step, width):
    """Renames files created from harp freq"""
    freqs = cfun.list_freqs(chromosome)
    for file in freqs:
        name_start = file.split('.freqs')[0]
        new_name = name_start + '_step' + step + '_width' + width + '.freqs'
        os.rename(file, new_name)


def mover(chromosome, width):
    """Moves all of the freq files into a folder"""
    current_dir = os.getcwd()
    new_dir = current_dir + '/freqs_' + chromosome + '_' + width + 'w'
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    freqs = cfun.list_freqs(chromosome)
    for file in freqs:
        shutil.move(file, new_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run harp like in parallel')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    parser.add_argument('--step', type=str, help='step size')
    parser.add_argument('--width', type=str, help='width size')
    args = parser.parse_args()
    freq_multi(args.chr, args.step, args.width, freq_process)
    renamer(args.chr, args.step, args.width)
    mover(args.chr, args.width)
