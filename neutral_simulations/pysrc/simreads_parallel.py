"""Multiprocessing of simreads config files, and execution"""
import glob
import subprocess
from multiprocessing import Pool
import sys
sys.path.append('/home/ltjones/neutral_sim/pysrc')
import common_functions as cfun


def list_config_files():
    """Lists all the simreads config files in folder"""
    return glob.glob('simreads_*.config')


def simreads(configfile):
    """The simreads command function"""
    sim_com = '/usr/local/bin/simreads'
    command = [sim_com, configfile]
    subprocess.call(command, shell=False)


def simreads_multi(target_func):
    """Multiprocess of the simreads commands"""
    configs = list_config_files()
    pool = Pool(21, cfun.limit_cpu)
    pool.map(target_func, configs)
    pool.close()
    pool.join()


if __name__ == '__main__':
    simreads_multi(simreads)
