# %%
import re
import pandas as pd
import numpy as np
import os
from shutil import copy, copyfile
from pathlib import Path
from pandas.core import groupby
import summary_methods

from tikzplotlib import save as tikz_save
import matplotlib.pyplot as plt

plt.style.use('ggplot')

# %%
data_file_info = summary_methods.DataFileInfo("./results_2022_08_25.csv")

# %%
data_file_info.write_cg_results()

# %%
data_file_info.pivot_cg_table()

# %%
# data_file_info.get_list_of_instances_solved_by_one_of_the_solvers('ArcTimeIndexed', 'TimeIndexed')
# data_file_info.get_list_of_instances_solved_by_one_of_the_solvers('ArcTimeIndexed', 'BddBackward')
# data_file_info.get_list_of_instances_solved_by_one_of_the_solvers('ArcTimeIndexed', 'BddBackwardCycle')
# data_file_info.get_list_of_instances_solved_by_one_of_the_solvers('TimeIndexed', 'BddBackwardCycle')
# data_file_info.get_list_of_instances_solved_by_one_of_the_solvers('TimeIndexed', 'BddBackward')
# data_file_info.get_list_of_instances_solved_by_one_of_the_solvers('BddBackward', 'BddBackwardCycle')

# %% 
data_file_info_exact = summary_methods.DataFileInfo("./CG_overall_2021_07_12.csv")
data_file_info_exact.create_templates()

# %%
data_file_info_exact.draw_overall_profile_curve()

# %%
data_file_info_exact.draw_profile_curve_per_instance_class()

# %%
data_file_info_exact.write_summary_exact_opt()

# %%
data_file_info_exact.write_all_instances()


# %%
data_file_info_exact.list_of_instance_only_solved_by_bdd()

# %%
