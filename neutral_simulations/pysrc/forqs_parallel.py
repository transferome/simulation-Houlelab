"""Run multiple forqs simulations in parallel"""
import glob
import random
import subprocess
from multiprocessing import Pool
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def list_forq_configs():
    """Creates List of Forq Config Files"""
    return glob.glob('*.config')


def forqs_sim(forqs_config):
    """Executes forqs config file"""
    forqs_executable = '/usr/local/bin/forqs'
    # make a random seed, not really needed because forqs writes a new random seed at end of each simulation
    forq_seed = 'seed=' + str(random.randint(0, 32767))
    # command = ' '.join([forqs_executable, forqs_config, forq_seed])
    command = [forqs_executable, forqs_config, forq_seed]
    # print(command)
    subprocess.call(command, shell=False)


if __name__ == '__main__':
    # execute parallel processing of forqs simulations
    config_list = list_forq_configs()
    pool = Pool(18, cfun.limit_cpu)
    simulations = pool.map(forqs_sim, config_list)
    pool.close()
    pool.join()
