"""Sum 10k forqs windows into 50k windows"""
import os
import glob
import pandas as pd


def list_final_csv(auto=True):
    """assuming working directory is directory containing final.csv's"""
    if auto:
        file_list = glob.glob('*final.csv')
        file_list.sort(key=lambda x: int(x.split('/')[-1].split('_')[1]))
        return file_list


def split_one_final(final_csv):
    """splits a final.csv into two with different start points to make comparison with overlapping harp freq
     files easier """
    file_start = final_csv.split('.csv')[0]
    outputfile1 = '{}_0start.csv'.format(file_start)
    outputfile2 = '{}_50kstart.csv'.format(file_start)
    with open(final_csv) as infile, open(outputfile1, 'w+') as outfile1, open(outputfile2, 'w+') as outfile2:
        data = [line for line in infile]
        for line in data:
            if int(line.split(',')[0]) < 3900000:
                outfile1.write(line)
        for line in data:
            if 50000 <= int(line.split(',')[0]) <= 3940000:
                outfile2.write(line)


def split_finals(auto=True):
    """splits all the final.csv files from forqs into files with two different start points"""
    if auto:
        csv_files = list_final_csv()
        for csv in csv_files:
            split_one_final(csv)


def list_zero_starts(auto=True):
    """Lists the zero start csv's"""
    if auto:
        file_list = glob.glob('*0start.csv')
        file_list.sort(key=lambda x: int(x.split('/')[-1].split('_')[1]))
        return file_list


def list_fifty_starts(auto=True):
    """Lists the 50k start csv's"""
    if auto:
        file_list = glob.glob('*50kstart.csv')
        file_list.sort(key=lambda x: int(x.split('/')[-1].split('_')[1]))
        return file_list


def final_average_csv(pandas_csv, region_start, region_stop):
    """Average the 10k windwos to 50k winoows"""
    new_lines = list()
    row_groups = list()
    for x in range(0, 385, 10):
        group = (x, x + 10)
        row_groups.append(group)
    average_rows = list()
    for tup in row_groups:
        rows_sum_temp = pandas_csv[tup[0]:tup[1]].sum(axis=0).tolist()
        # remove first element it is sum of region value
        rows_sum = rows_sum_temp[1:]
        rows_avg = [x/10 for x in rows_sum]
        average_rows.append(rows_avg)
    regions = range(region_start, region_stop, 100000)
    for region, row in zip(regions, average_rows):
        line = [str(region)] + [str(x) for x in row]
        new_lines.append(','.join(line))
    return new_lines


def average_csvs(auto=True):
    """Writing the average csv's"""
    if auto:
        zeroes = list_zero_starts()
        fifties = list_fifty_starts()
        for z, f in zip(zeroes, fifties):
            dat_zero = pd.read_csv(z, header=None)
            dat_50k = pd.read_csv(f, header=None)
            dat_zero_avg = final_average_csv(dat_zero, 0, 3900000)
            dat_50k_avg = final_average_csv(dat_50k, 50000, 3950000)
            new_filename_zero = z.split('.csv')[0] + '_avg.csv'
            new_filename_fifty = f.split('.csv')[0] + '_avg.csv'
            with open(new_filename_zero, 'w+') as zout:
                for line in dat_zero_avg:
                    zout.write('{}\n'.format(line))
            with open(new_filename_fifty, 'w+') as fout:
                for line in dat_50k_avg:
                    fout.write('{}\n'.format(line))


if __name__ == '__main__':
    os.chdir('/home/solid-snake/Simulations/2017_Example/forqs_final_freqs')
    split_finals()
    average_csvs()
