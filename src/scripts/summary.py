
# %%
import re
import pandas as pd
import numpy as np
import os
import sys
from shutil import copy, copyfile
from pathlib import Path
from pandas.core import groupby

from tikzplotlib import save as tikz_save
import matplotlib.pyplot as plt

workdir = Path.cwd().parent
results = workdir.joinpath(Path("./results"))

# %% Load the data of the new results
file_name = "~/results_2021_12_10.csv"
file_path = Path(file_name)
data = pd.read_csv(file_name, index_col=False)
match = re.search(r'(CG_overall|results)\_(20\d{2})\_(\d{2})\_(\d{2})\.csv', file_name)
year = match.group(2)
month = match.group(3)
day = match.group(4)

# %% Calculate some extra columns
data['gap'] = (data['global_upper_bound'] - data['global_lower_bound']
               )/(data['global_lower_bound'] + 0.00001)
data['opt'] = data['global_lower_bound'] == data['global_upper_bound']
data['reduction'] = (data['first_size_graph'] -
                     data['size_after_reduced_cost'])/(data['first_size_graph'] + 0.000001)
data['Inst'] = data.NameInstance.apply(
    lambda x:  int(re.search(r'.*\_(\d+)', x).group(1)))

for it in data.columns[3:12]:
    data[it] = 0.7*data[it]

# %% create result directory and copy results to that directory

results_path = results.joinpath("./results_{}_{}_{}".format(year, month, day))

if results_path.exists() == False:
    os.mkdir(results_path)

# copy(file_path, results_path.joinpath(match.group(0)))
tex_file = str()

# %% Create tex files for Column generation results
template_dir_path = results.joinpath("./template_dir")
for lst in template_dir_path.iterdir():
    if lst.name == "CG_tables_template.tex":
        copy(lst, results_path.joinpath(
            "CG_tables_{}_{}_{}.tex".format(year, month, day)))
        tex_file = str(results_path.joinpath(
            "CG_tables_{}_{}_{}.tex".format(year, month, day)))
    else:
        copy(lst, results_path.joinpath(lst.name))

os.popen("sd  \"CG_summary_20191024.csv\" \"CG_summary_{}_{}_{}.csv\" ".format(
    year, month, day)+tex_file)
os.popen("sd  \"CG_allinstances_20191024.csv\" \"CG_allinstances_{}_{}_{}.csv\" ".format(
    year, month, day)+tex_file)


# %% Compute summary results for CG over all solvers
summary_grouped = data.groupby(['pricing_solver', 'n', 'm'])
aggregation = {"tot_lb": {np.max, np.mean},
               "gap": {np.max, np.mean},
               "first_size_graph": {np.max, np.mean},
               "opt": np.sum,
               "reduction": {np.max, np.mean},
               "tot_cputime": {np.max, np.mean}}
summary_write = summary_grouped.agg(aggregation).pivot_table(index=['n', 'm'], values=[
    'tot_lb', 'gap', 'first_size_graph', 'reduction', 'opt'], columns=['pricing_solver'])
print(summary_write.columns)
summary_write.columns.set_levels(
    ['AFBC', 'TI_BDD'], level=2, inplace=True)
summary_write.columns = ["_".join(x) for x in summary_write.columns.ravel()]
summary_write.to_csv(results_path.joinpath(
    "/CG_summary_{}_{}_{}.csv".format(year, month, day)))

# %% pivot results for all pricing_solvers
all_instances = data.pivot_table(values=['tot_lb', 'gap', 'first_size_graph', 'reduction', 'opt', 'rel_error', 'nb_generated_col',
                                         'global_lower_bound', 'global_upper_bound', 'tot_cputime', 'tot_bb'], index=['n', 'm', 'Inst'], columns=['pricing_solver'])
all_instances.columns.set_levels(
    ['AFBC'], level=1, inplace=True)
all_instances.columns = ["_".join(x) for x in all_instances.columns.ravel()]
all_instances.to_csv(
    results_path.joinpath(
        "CG_allinstances_{}_{}_{}.csv".format(year, month, day)))

# %% Load results of Pessoa et al. and Oliveira qnd Pessoa
df_pessoa = pd.read_csv(results.joinpath("all_pessoa.csv"))
df_pessoa.Opt = df_pessoa.Opt.apply(str)
df_pessoa['best'] = df_pessoa.apply(lambda x: re.search(
    r'[^0-9]?(\d+)', x['Opt']).group(1), axis=1)
df_pessoa.best = df_pessoa.best.apply(pd.to_numeric)

df_oliveira = pd.read_csv(results.joinpath("oliveira_overall.csv"))

# %% Merge our results with results of Oliveira
df_all = pd.merge(data, df_oliveira, on=['Inst', 'n', 'm'])


# %% Compute overall performance profile curve
df_all['best_solver'] = df_all[['tot_bb', 'TimeOliveira']].min(axis=1)
df_all['ratio_tot_bb_best'] = df_all['tot_bb'] / df_all['best_solver']

df_all['ratio_TimeOliveira_best'] = df_all['TimeOliveira'] / \
    df_all['best_solver']

sorted_ratio_tot_bb = df_all[['ratio_tot_bb_best']
                             ].sort_values(by='ratio_tot_bb_best')
yvals = np.arange(len(sorted_ratio_tot_bb)) / \
    float(len(sorted_ratio_tot_bb) - 1.0)

sorted_ratio_TimeOliveira = df_all[['ratio_TimeOliveira_best']].sort_values(
    by='ratio_TimeOliveira_best')

yvalues = np.arange(len(sorted_ratio_TimeOliveira)) / \
    float(len(sorted_ratio_TimeOliveira) - 1.0)

width, height = plt.figaspect(1.68)
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
ax.step(sorted_ratio_tot_bb, yvals, label='BDD')
ax.step(sorted_ratio_TimeOliveira, yvalues, label='ATIF')
ax.set_xlim([10**0, 100])
# ax.set_title(
#     r"Performance profile for instances with $m = %d$ and $n = %d$"
#     % (i, j))
ax.set_xlabel(r"$\tau$")
ax.set_ylabel(r"$P(r_{p,s} \leq \tau)$")
ax.legend(loc='lower right')
# plt.savefig('profile_curve%d_%d.pdf' % (i, j), dpi=200)
name_file = 'profile_curve_overall_{}_{}_{}.tex'.format(year, month, day)
tikz_save(results_path.joinpath(name_file))
plt.savefig(results_path.joinpath('profile_curve_overall_{}_{}_{}.pdf'.format(
    year, month, day)), dpi=200)
# %% Compute performance profile curves per instance class
for n in [40, 50]:
    for m in [2, 4]:
        sorted_ratio_tot_bb = df_all.loc[(df_all['n'] == n) & (
            df_all["m"] == m), "ratio_tot_bb_best"].sort_values()
        yvals = np.arange(len(sorted_ratio_tot_bb)) / \
            float(len(sorted_ratio_tot_bb) - 1.0)

        sorted_ratio_TimeOliveira = df_all.loc[(df_all['n'] == n) & (
            df_all["m"] == m), "ratio_TimeOliveira_best"].sort_values()

        yvalues = np.arange(len(sorted_ratio_TimeOliveira)) / \
            float(len(sorted_ratio_TimeOliveira) - 1.0)

        width, height = plt.figaspect(1.68)
        fig, ax = plt.subplots(figsize=(width, height), dpi=200)
        ax.step(sorted_ratio_tot_bb, yvals, label='BDD')
        ax.step(sorted_ratio_TimeOliveira, yvalues, label='ATIF')
        ax.set_xlim([10**0, 100])
        ax.set_title(
            "Performance profile for instances with $m = {}$ and $n = $".format(m, n))
        ax.set_xlabel(r"$\tau$")
        ax.set_ylabel(r"$P(r_{p,s} \leq \tau)$")
        ax.legend(loc='lower right')
        name_file = 'profile_curve_overall_{}_{}_{}_{}_{}.tex'.format(
            n, m, year, month, day)
        tikz_save(results_path.joinpath(name_file))
        plt.savefig(results_path.joinpath('profile_curve_{}_{}_{}_{}_{}.pdf'.format(
            n, m, year, month, day)), dpi=200)
