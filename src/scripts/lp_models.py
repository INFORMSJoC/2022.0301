import argparse
import re
import time
import os
import sys
import datetime
from gurobipy import *


class DataParser(argparse.ArgumentParser):

    def __init__(self):
        super(DataParser, self).__init__()
        self.add_argument('-m', '--machines', nargs='+', type=int,
                          help='store a list of  number of jobs that' +
                          ' you want to calculate',
                          default=[3, 5, 8, 10, 12])
        self.add_argument('-d', '--directory', type=str,
                          help='Give the directory of the instances',
                          default='./instances_and_or/wt040/0.2')
        self.add_argument('-f', '--file', type=str,
                          help='Give the instance',
                          default='/home/daniel/Dropbox/PapersOR/PM/implementation/resources/instances/wt40/wt040/wt040_001.dat')
        self.add_argument('-M', '--model', nargs='+', type=str,
                          help='TI Model',
                          default='Laura')


def read_instance_rcp(file):
    with open(file) as f:
        lines = f.readlines()
        n = int(lines[1].split()[0]) - 2
        m = int(lines[1].split()[1])
        and_constr = list()
        for idx, line in enumerate(lines[5:n + 5]):
            a = map(lambda x: x - 1, map(int, line.split())[m + 2:])
            for i in a:
                if i is not n + 1:
                    and_constr.append((idx + 1, i))

    f.close()
    return n, and_constr


def read_instance_file_type(file):
    with open(file) as f:
        lines = f.readlines()
        line = lines[0]
        matchObj = re.match(r'^\#\sn=(\d+)', line)
        n = int()
        and_constr = list()
        if matchObj is None:
            print("Error")
        else:
            n = int(matchObj.group(1))
        for line in lines[1:]:
            a = map(int, line.split())
            and_constr.append((a[0] + 1, a[1] + 1))

    f.close()
    return n, and_constr


def read_instance_salim(file):
    with open(file) as f:
        lines = f.readlines()
        and_constr = list()
        n = int(lines[0].split()[0])
        for idx, line in enumerate(lines[1:]):
            a = map(int, line.split()[5:])
            for i in a:
                and_constr.append((idx + 1, i + 1))

    f.close()
    return n, and_constr


def read_instance_TW(file, m):
    with open(file) as f:
        lines = f.readlines()
        n = int(lines[0].split()[0])
        p = [0]
        w = [0]
        d = [0]
        p_max = 0
        T = 0
        for idx, line in enumerate(lines[1:]):
            a = map(int, line.split())
            print(a)
            pp = a[0]
            dd = a[1] / m
            ww = a[2]
            if pp > dd:
                d[0] += ww * (pp - dd)
                dd = pp
            if p_max < pp:
                p_max = pp
            p.append(pp)
            d.append(dd)
            w.append(ww)
            p[0] += pp

        T = p[0] - p_max
        T /= m
        T += p_max

    f.close()

    return n, T, p, w, d


class TI_model(object):
    """docstring for TI_model"""

    def __init__(self, file, file_type, m):
        super(TI_model, self).__init__()
        self.m = Model()
        self.x = dict()
        self.p = None
        self.w = None
        self.d = None
        self.u = None
        self.T = None
        self.cputime = None
        self.file = file
        self.nb_mach = float(m)
        self.statistics_file = str()
        if file_type is 0:
            self.n, self.and_constr = read_instance_rcp(file)
        elif file_type is 1:
            self.n, self.T, self.p, self.w, self.d = read_instance_TW(file, m)
        elif file_type is 2:
            self.n, self.and_constr = read_instance_salim(file)

        self.iterator = range(1, self.n + 1)

    def print_statistics(self):
        print(os.getcwd())
        cur_path = os.path.join(os.getcwd(), os.curdir)
        os.chdir("/home/daniel/")
        f = open(self.statistics_file, "a")
        if os.stat(self.statistics_file).st_size == 0.0:
            f.write("name, status, cputime, runtime, gap,\
                    objBound, objVal, data, time\n")

        name = os.path.basename(self.file)
        status = self.m.getAttr(GRB.Attr.Status)
        cputime = str(self.cputime)
        runtime = self.m.getAttr(GRB.Attr.Runtime)
        gap = self.m.getAttr(GRB.Attr.MIPGap)
        objBound = self.m.getAttr(GRB.Attr.ObjBound)
        objValue = self.m.getAttr(GRB.Attr.ObjVal)
        now = datetime.datetime.now()
        date = "%s/%s/%s" % (now.day, now.month, now.year)
        time = "%s:%s:%s" % (now.hour, now.minute, now.second)

        s = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            name, status, cputime, runtime, gap,
            objBound, objValue, date, time)

        f.write(s)
        f.close()
        os.chdir(cur_path)

    def solve_model(self):
        self.m.setParam(GRB.Param.Threads, 1.0)
        # self.m.setParam(GRB.Param.OutputFlag, 0)
        self.cputime = time.clock()
        self.m.optimize()
        self.cputime = time.clock() - self.cputime


class TI_model_Laura(TI_model):
    """docstring for TI_model_Laura"""

    def __init__(self, file, file_type, m):
        TI_model.__init__(self, file, file_type, m)
        self.statistics_file = "%d_%d_TI_Laura.csv" % (self.n, m)

    def build_model(self):
        # construct all the variables
        for i in self.iterator:
            for t in self.iterator:
                self.x[i, t] = self.m.addVar(
                    obj=0, vtype="B", name="x[%d,%d]" % (i, t))

        self.u = self.m.addVar(obj=1.0, vtype="C", lb=0.0)
        self.m.update()

        # assign every job to a time period
        for i in self.iterator:
            coef = [1.0 for t in self.iterator]
            vars = [self.x[i, t] for t in self.iterator]
            linexpr = LinExpr(coef, vars)
            self.m.addConstr(linexpr, "=", 1.0)
        self.m.update()

        # No more than m jobs at a time
        for t in self.iterator:
            coef = [1.0 for i in self.iterator]
            vars = [self.x[i, t] for i in self.iterator]
            linexpr = LinExpr(coef, vars)
            self.m.addConstr(linexpr, "<", self.nb_mach)
        self.m.update()

        # add the and-constraints
        for t in self.iterator:
            for e in self.and_constr:
                coef = [1.0 for i in range(1, t)]
                vars = [self.x[e[0], i] for i in range(1, t)]
                if len(coef) is 0:
                    self.m.addConstr(0.0, ">", self.x[e[1], t])
                else:
                    linexpr = LinExpr(coef, vars)
                    self.m.addConstr(linexpr, ">", self.x[e[1], t])
        self.m.update()

        # bound the objective value
        for t in self.iterator:
            for i in self.iterator:
                self.m.addConstr(float(t) * self.x[i, t], "<", self.u)

        self.m.update()

    def print_solution(self):
        if self.m.status is GRB.Status.OPTIMAL:
            sol = self.m.getAttr(GRB.Attr.X, self.x)
            obj = int(self.m.getAttr(GRB.Attr.ObjVal))
            for t in range(1, obj + 1):
                a = str()
                for j in self.iterator:
                    if sol[j, t] == 1.0:
                        a += str(j) + ' '

                print("at time %d: %s" % (t, a))


class TI_model_General(TI_model):
    """docstring for TI_model_Peter"""

    def __init__(self, file, file_type, m):
        TI_model.__init__(self, file, file_type, m)
        self.statistics_file = "%d_%d_TW.csv" % (self.n, m)

    def build_model(self):
        # construct all the variables
        for i in self.iterator:
            for t in range(0, self.T):
                self.x[i, t] = self.m.addVar(
                    obj=self.w[i] * max(0, t + self.p[i] - self.d[i]),
                    vtype="C", name="x[%d,%d]" % (i, t))

        for i in self.iterator:
            coef = [1.0 for t in range(0, self.T - self.p[i] + 1)]
            vars = [self.x[i, t] for t in range(0, self.T)]
            linexpr = LinExpr(coef, vars)
            self.m.addConstr(linexpr, "=", 1.0)
        self.m.update()

        for t in range(0, self.T):
            coef = [1.0 for i in self.iterator for j in range(
                max(0, t - self.p[i] + 1), t + 1)]
            vars = [self.x[i, j]
                    for i in self.iterator for j in range(
                        max(0, t - self.p[i] + 1), t + 1)]
            linexpr = LinExpr(coef, vars)
            self.m.addConstr(linexpr, "<", self.nb_mach)
        self.m.update()

    def print_solution(self):
        if self.m.status is GRB.Status.OPTIMAL:
            sol = self.m.getAttr(GRB.Attr.X, self.x)
            obj = int(self.m.getAttr(GRB.Attr.ObjVal))
            for t in range(0, self.T):
                a = str()
                for j in self.iterator:
                    if sol[j, t] > 0:
                        a += str(j) + ' ' + str(sol[j, t]) + ' '

                if len(a) > 0:
                    print("at time %d: %s" % (t, a))


def main():
    try:
        parser = DataParser()
    except argparse.ArgumentError:
        print("Error")
        sys.exit(2)
    else:
        args = parser.parse_args()

    for m in args.machines:
        name = args.file
        print("Solving instance %s" % name)
        file_type = None
        model = None
        if name.endswith("rcp"):
            file_type = 0
        elif name.endswith("dat"):
            file_type = 1
        elif name.endswith("txt"):
            file_type = 2

        model = TI_model_General(name, file_type, m)

        model.build_model()
        model.solve_model()
        model.print_solution()
        # model.print_statistics()


if __name__ == '__main__':
    main()
