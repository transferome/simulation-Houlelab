"""Common Functions Used In Different Simulation Module scripts"""
import psutil
import os
import glob
import shutil


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
    p.nice(6)


def config_command(chr_arm, chr_length, recmap, n_size, configtype='Reg'):
    """Create the command to run the config multiprocess"""
    basic_command = ['python3', '/home/ltjones/neutral_sim/pysrc/forqs_configs.py', '--configtype', 'config_type_holder',
               '--chr', chr_arm, '--length', str(chr_length), '--recmap', recmap, '--size', str(n_size)]
    if configtype == 'Reg':
        basic_command[3] = 'Reg'
    if configtype == 'HAPS':
        basic_command[3] = 'HAPS'
    return basic_command


def forq_command():
    """Create the command to run the simulation multiprocess"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/forqs_parallel.py']
    return command


def move_configs(chromosome_arm):
    """Move all of the config files created for forqs, into the forqs directory the config file created"""
    forqs_dirs_sorted = list_forqs_directories(chromosome_arm)
    configs = glob.glob('*.config')
    configs_sorted = sorted(configs, key=lambda x: int(x.split('.')[0].split('_')[1]))
    if len(configs_sorted) == len(forqs_dirs_sorted):
        for config_file, forqs_dir in zip(configs_sorted, forqs_dirs_sorted):
            if config_file.split('.')[0].split('_')[1] == forqs_dir.split('/')[-2].split('_')[1]:
                shutil.move(config_file, forqs_dir)
            else:
                print('Forqs ID Numbers Dont Match')
    else:
        print('Different Number of Config Files and forqs Directories')


if __name__ == '__main__':
    pass
