"""Multiprocessing of Config File Creation"""
import os
import sys
from multiprocessing import Pool
from functools import partial
import argparse
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def make_stem(chromosome, sim_number):
    """Creates a stem for the simreads config file, based off directories already existing in folder
     program is executed in"""
    current_dir = os.getcwd()
    stem = 'chr{}_{}'.format(chromosome, str(sim_number))
    check_stem = current_dir + '/' + stem
    if os.path.exists(check_stem):
        print('Error With Stem Creation')
    else:
        return stem


# TODO: Using dummy_var because unsure what to pass map when all processes have same arguments,
# TODO: value of 200 for final population size is setting the sample size of the individuals sampled for sequencing
def forqs_config(chromosome, chromosome_length, recomb_map_file, sim_number):
    """Creates A Config File Given A  Chromosome"""
    stem = make_stem(chromosome, sim_number)
    output_dir_string = '\toutput_directory = ' + stem + '\n'
    chrom_length_string = '\tchromosome_lengths = ' + chromosome_length + '\n'
    recomb_map_string = '\tfilename = ' + recomb_map_file + '\n'
    config_file = ['Trajectory_Constant popsize_0\n', '    value = 53\n', '\n', 'Trajectory_Constant popsize_1\n',
                   '    value = 2000\n', '\n', 'Trajectory_Constant popsize_2\n', '    value = 30\n', '\n',
                   'Trajectory_Constant popsize_3\n', '    value = 200\n', '\n',
                   'Trajectory_GenerationComposite popsize\n', '    generation:trajectory = 0 popsize_0\n',
                   '    generation:trajectory = 1 popsize_1\n', '    generation:trajectory = 16 popsize_2\n',
                   '    generation:trajectory = 29 popsize_3\n', '\n',
                   'PopulationConfigGenerator_LinearSteppingStone pcg\n',
                   chrom_length_string, '    chromosome_pair_count = 1\n', '    generation_count = 30\n',
                   '    id_offset_step = 0\n', '    population_count = 1\n', '    population_size = popsize\n', '\n',
                   'RecombinationPositionGenerator_RecombinationMap rpg_map\n', recomb_map_string, '\n',
                   'Reporter_Population reporter_population\n', '\n', 'SimulatorConfig\n', output_dir_string,
                   '    population_config_generator = pcg\n', '    recombination_position_generator = rpg_map\n',
                   '    reporter = reporter_population\n']
    config_filename = stem + '.config'
    with open(config_filename, 'w+') as f:
        for s in config_file:
            f.write(s)
    # return stem


def multi_config(chromosome, chromosome_length, recomb_map_file, number_of_simulations):
    """Creats Config Files In Parallel"""
    sim_numbers = [x for x in range(number_of_simulations)]
    pool = Pool(18, cfun.limit_cpu)
    target_func = partial(forqs_config, chromosome, str(chromosome_length), recomb_map_file)
    pool.map(target_func, sim_numbers)
    pool.close()
    pool.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creat Config Files in Parallel')
    parser.add_argument('--chr', type=str, help='chromosome arm')
    parser.add_argument('--length', type=int, help='chromosome length')
    parser.add_argument('--recmap', type=str, help='recombination map file')
    parser.add_argument('--size', type=int, help='number of simulation iterations')
    args = parser.parse_args()
    multi_config(args.chr, args.length, args.recmap, args.size)
