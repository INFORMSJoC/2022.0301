#!/usr/bin/env python

from random import randint
import numpy as np
import argparse
import os


class OptionParser(argparse.ArgumentParser):
    """OptionParser"""

    def __init__(self):
        super(OptionParser, self).__init__()
        self.add_argument('-c', '--class_inst', type=int, nargs='+',
                          help='class instance',
                          default=[1, 2, 3, 4, 5, 6])
        self.add_argument('-n', '--number_jobs', type=int, nargs='+',
                          help='number of jobs', default=[20, 50, 100, 150])
        self.add_argument('-m', '--number_machines', type=int, nargs='+',
                          help='number of machines', default=[3, 5, 8, 10, 12])


class ValidationError(Exception):
    """Vaildation Error"""

    def __init__(self, msg, err):
        super(ValidationError, self).__init__()
        self.msg = msg
        self.err = err


def give_weight_duration(n):
    if n > 6:
        raise ValidationError('Error occured number to large', n)
    elif n == 1:
        return randint(10, 100), randint(1, 10)
    elif n == 2:
        return randint(1, 100), randint(1, 100)
    elif n == 3:
        return randint(10, 20), randint(10, 20)
    elif n == 4:
        return randint(90, 100), randint(90, 100)
    elif n == 5:
        p = randint(90, 100)
        w = randint(p - 5, p + 5)
        return w, p
    elif n == 6:
        p = randint(10, 100)
        w = randint(p - 5, p + 5)
        return w, p


def main():
    parser = OptionParser()
    args = parser.parse_args()

    for n in args.number_jobs:
        directory = './wt%03d' % (n)
        os.mkdir(directory)
        i = 1
        for rdd in [0.2 + 0.2 * i for i in range(0, 5)]:
            for tdf in [0.2 + 0.2 * i for i in range(0, 5)]:
                for j in range(0, 5):
                    p = np.random.randint(low=1, high=100, size=n)
                    w = np.random.randint(low=1, high=10, size=n)
                    sum_p = p.sum()
                    a = np.max([0, sum_p * (1 - tdf - rdd / 2)])
                    b = np.max([0, sum_p * (1 - tdf + rdd / 2)])
                    d = np.random.randint(
                        low=a, high=b, size=n)
                    f = open(directory + '/wt%03d_%03d.txt' %
                             (n, i), 'w')
                    f.write("%d\n" % (n))
                    for it in zip(p, d, w):
                        line = "%d %d %d\n" % (it[0], it[1], it[2])
                        f.write(line)
                    f.close()
                    i = i + 1


if __name__ == '__main__':
    main()
