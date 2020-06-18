"""Common Functions Used In Different Simulation Module scripts"""
import psutil
import os
import glob


def region_min_max(snp_txt):
    """ Finds the minimum and maximum coordinate position within the snp_txt file """
    # [1:] removes header info which is the chromosome
    first_column = [s.split('\t')[0] for s in open(snp_txt)]
    position_integer = [int(s) for s in first_column]
    return [min(position_integer), max(position_integer)]


def list_forqs_directories(chromosome):
    """List all the directories created by forqs"""
    current_dir = os.getcwd()
    dirs = glob.glob(current_dir + '/*/')
    forqs_dirs = [s for s in dirs if 'chr' + chromosome in s]
    forqs_dirs = [s for s in forqs_dirs if '.output' not in s]
    forqs_dirs_sorted = sorted(forqs_dirs, key=lambda x: int(x.split('/')[-2].split('_')[1]))
    return forqs_dirs_sorted


def list_stems(chromosome):
    """Converts forqs directories to list of forqs stems"""
    forqs_directories_list = list_forqs_directories(chromosome)
    stem_list = list()
    for s in forqs_directories_list:
        stem_list.append(s.split('/')[-2])
    return stem_list


def list_bams(chromosome):
    """Lists all the bam files in the folder"""
    return glob.glob('simreads_chr' + chromosome + '*.bam')


def list_sams(chromosome):
    """Lists all the sam files in the folder"""
    return glob.glob('simreads_chr' + chromosome + '*.sam')


def list_sorted_bams(chromosome):
    """Lists all sorted bam files"""
    return glob.glob('simreads_chr' + chromosome + '*_sorted.bam')


def list_sorted_bams_idx(chromosome):
    """Lists all sorted bam bai files"""
    return glob.glob('simreads_chr' + chromosome + '*_sorted.bam.bai')


def list_forqs_haplotypes(chromosome):
    """Lists all forqs haplotypes txt files"""
    return glob.glob('chr' + chromosome + '*_haplotypes.txt')


def list_forqs_idx(chromosome):
    """Lists all forqs haplotypes txt.idx files"""
    return glob.glob('chr' + chromosome + '*_haplotypes.txt.idx')


def list_hlk(chromosome):
    """Lists all hlk files in folder"""
    return glob.glob('chr' + chromosome + '_*.hlk')


def list_freqs(chromosome):
    """Lists all freq files created from harp freq"""
    return glob.glob('chr' + chromosome + '_*.freqs')


def limit_cpu():
    """Called At Every Process Started in a multiprcoess pool"""
    p = psutil.Process(os.getpid())
    p.nice(15)


if __name__ == '__main__':
    pass
