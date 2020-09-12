"""Module For Commonly Used Functions To Speed Up Coding Efficiency"""


def reader(filename, newline=True):
    """Read Text File Into List"""
    if newline:
        with open(filename) as text_file:
            text_lines = [line for line in text_file]
    else:
        with open(filename) as text_file:
            text_lines = [line.rstrip('\n') for line in text_file]
    return text_lines


def writer(filename, text_list, newline=True):
    """Write Text List To Text File"""
    if newline:
        with open(filename, 'w+') as text_file:
            for line in text_list:
                text_file.write(line)
    else:
        with open(filename, 'w+') as text_file:
            for line in text_list:
                text_file.write('{}\n'.format(line))


if __name__ == '__main__':
    test = reader('/home/solid-snake/Programs/ssdTrim')
    print(test)
