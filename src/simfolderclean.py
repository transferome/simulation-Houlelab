"""Removes files produced from a simulation, but not files needed to initiate the simulation"""
import os
import glob
import argparse
import shutil


def clean(chromosome):
    """Clean function keeps all .py files and files append to files_keep, deletes all other
    folders and files in the directory in which the script is within"""
    files_keep = glob.glob('*.py')
    files_keep.append('{}rcc.txt'.format(chromosome))
    files_keep.append('dgrp{}_mixed_haplotypes.txt'.format(chromosome))
    folder_contents = os.listdir()
    content_delete = [content for content in folder_contents if content not in files_keep]
    for content in content_delete:
        if os.path.isfile(content):
            os.remove(content)
        else:
            shutil.rmtree(content)
    print('Folder Cleaned')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Cleans Up Simulation Folders')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    args = parser.parse_args()
    clean(args.chr)
