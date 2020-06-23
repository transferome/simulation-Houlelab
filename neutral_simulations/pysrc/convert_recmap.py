"""Converts rcc map from webform to type used by forqs"""
import numpy as np


def read_rcc(filename):
    """Reads in rcc text"""
    with open(filename) as f:
        out_list = list()
        for line in f:
            line = line.replace(' ', '')
            line = line.rstrip('\n')
            out_list.append(line)
    return out_list


def subset_rcc(rcc_list, chr_arm, chromosome_range):
    """Subset the rcc text list"""
    out_list = list()
    newfile = '{}_{}-{}.csv'.format(chr_arm, chromosome_range[0], chromosome_range[1])
    with open(newfile, 'w+') as f:
        for line in rcc_list[1:]:
            if int(chromosome_range[0]) <= int(line.split('\t')[0].split('..')[0].split(':')[-1]) < int(chromosome_range[1]):
                out_list.append(line)
                f.write(line)
    return out_list


def get_midpoint(data, lower_bound):
    """Get midpoint value for each row"""
    out_list = list()
    for line in data[1:]:
        p2 = int(line.split('\t')[0].split('..')[-1])
        p1 = int(line.split('\t')[0].split('..')[0].split(':')[-1])
        mid_point = (p1 + p2)//2
        # TODO: is subtracting chromosome length from mid-point really what I want to fix range problem
        mid_point = mid_point - int(lower_bound)
        out_list.append(str(mid_point))
    return out_list


def get_midpoint_rate(data):
    """Gets the cameron midpoint rate"""
    out_list = list()
    for line in data:
        value = float(line.split('\t')[5])
        out_list.append(value)
    return out_list


def get_cumulative_rate(mid_rate_list):
    """Gets cumulate rate from midpoint_rate list"""
    mid_list = [x/10 for x in mid_rate_list]
    cum_array = np.cumsum(mid_list)
    out_list = cum_array.tolist()
    return [round(x, 5) for x in out_list]


def write_recmap(pos, mid, cum, chr_arm, chromosome_range):
    """Writes information to new file"""
    newfilename = 'dmel_recRates_{}_{}-{}.csv'.format(chr_arm, chromosome_range[0], chromosome_range[1])
    header = "position COMBINED_rate(cM / Mb) Genetic_Map(cM)\n"
    mid_str = [str(x) for x in mid]
    cum_str = [str(x) for x in cum]
    with open(newfilename, 'w+') as f:
        f.write(header)
        f.write("50000 0 0\n")
        for p, m, c in zip(pos, mid_str, cum_str):
            line = ' '.join([p, m, c]) + '\n'
            f.write(line)


def main_recmap(filename, chr_arm, chromosome_range):
    """Main functions combines previously defined functions into one function"""
    rcc_file = read_rcc(filename)
    dat = subset_rcc(rcc_file, chr_arm, chromosome_range)
    positions = get_midpoint(dat, chromosome_range[0])
    mid_rate = get_midpoint_rate(dat)
    cum_rate = get_cumulative_rate(mid_rate)
    write_recmap(positions, mid_rate, cum_rate, chr_arm, chromosome_range)


if __name__ == '__main__':
    pass
