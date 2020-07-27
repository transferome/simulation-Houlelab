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


def convert_pops2snps(chr_arm, lower_bound):
    """Converts all the population texts to snps in parallel"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/convert_forqs_pop.py', '--chr', chr_arm,
               '--low_bound', lower_bound]
    return command


def index_tables(chr_arm):
    """Indexes all snp tables in parallel prior to running simreads"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/index_snp_table_parallel.py', '--chr', chr_arm]
    return command


def simreads_com():
    """Runs simreads in parallel for all config files"""
    return ['python3', '/home/ltjones/neutral_sim/pysrc/simreads_parallel.py']


def simreads_cleanup(chromosome):
    """Cleans up files after simreads run, and puts them in the correct stem folder"""
    forqs_dirs = list_forqs_directories(chromosome)
    act_freq_list = glob.glob('simreads_*.actual.freqs')
    act_freq_sorted = sorted(act_freq_list, key=lambda x: int(x.split('.')[0].split('_')[-1]))
    true_freq_list = glob.glob('simreads_*.true.freqs')
    true_freq_sorted = sorted(true_freq_list, key=lambda x: int(x.split('.')[0].split('_')[-1]))
    config_list = glob.glob('simreads_*.config')
    config_sorted = sorted(config_list, key=lambda x: int(x.split('.')[0].split('_')[-1]))
    seed_list = glob.glob('simreads_*.seed')
    seed_sorted = sorted(seed_list, key=lambda x: int(x.split('.')[0].split('_')[-1]))
    len_check = [len(act_freq_sorted), len(true_freq_sorted), len(config_sorted), len(seed_sorted)]
    if all(len(forqs_dirs) == x for x in len_check):
        for act_freq, true_freq, config, seed, forq_dir in zip(act_freq_sorted, true_freq_sorted, config_sorted,
                                                               seed_sorted, forqs_dirs):
            shutil.move(act_freq, forq_dir)
            shutil.move(true_freq, forq_dir)
            shutil.move(config, forq_dir)
            shutil.move(seed, forq_dir)
    else:
        print('File lists are of unequal length')


def sam2bam(chr_arm):
    """Converts simulated sam files to bam files"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/samtools_view_parallel.py', '--chr', chr_arm]
    return command


def bam_sort(chr_arm):
    """Sorts the simulated bam files"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/samtools_sort_parallel.py', '--chr', chr_arm]
    return command


def sambam_remove(chr_arm):
    """Removes the old sam, and unsorted bam"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/remove_sambam_parallel.py', '--chr', chr_arm]
    return command


def bam_index(chr_arm):
    """Indexes all the sorted bam files"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/samtools_index_parallel.py', '--chr', chr_arm]
    return command


def harp_like(chr_arm):
    """Runs harp like in parallel"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/harp_like_parallel.py', '--chr', chr_arm]
    return command


def harp_like_cleanup(chromosome):
    """Cleans up files created by forqs and harp_like"""
    # remove forqs.seed file
    try:
        os.remove('forqs.seed')
    except OSError:
        pass
    forqs_dirs = list_forqs_directories(chromosome)
    haps = list_forqs_haplotypes(chromosome)
    idxs = list_forqs_idx(chromosome)
    bams = list_sorted_bams(chromosome)
    bais = list_sorted_bams_idx(chromosome)
    haps_sorted = sorted(haps, key=lambda x: int(x.split('_haplotypes')[0].split('_')[-1]))
    idxs_sorted = sorted(idxs, key=lambda x: int(x.split('_haplotypes')[0].split('_')[-1]))
    bams_sorted = sorted(bams, key=lambda x: int(x.split('_sorted.bam')[0].split('_')[-1]))
    bai_sorted = sorted(bais, key=lambda x: int(x.split('_sorted.bam')[0].split('_')[-1]))
    len_check = [len(haps_sorted), len(idxs_sorted), len(bams_sorted), len(bai_sorted)]
    if all(len(forqs_dirs) == x for x in len_check):
        for hap, idx, bam, bai, forq_dir in zip(haps_sorted, idxs_sorted, bams_sorted, bai_sorted, forqs_dirs):
            shutil.move(hap, forq_dir)
            shutil.move(idx, forq_dir)
            shutil.move(bam, forq_dir)
            shutil.move(bai, forq_dir)
    else:
        print('File lists are of unequal length')


def harp_freq(chr_arm, step, width):
    """Runs harp freq in parallel"""
    command = ['python3', '/home/ltjones/neutral_sim/pysrc/harp_freq_parallel.py', '--chr', chr_arm,
               '--step', step, '--width', width]
    return command


def remove_hlk(chr_arm):
    """Removes the hlk files, they are large"""
    hlks = list_hlk(chr_arm)
    for hlk in hlks:
        os.remove(hlk)


if __name__ == '__main__':
    pass
