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

# %% This are the results from the CG algorithm
data_file_info = summary_methods.DataFileInfo("./results_2022_08_25.csv")
data_file_info.create_templates()

# %%
data_file_info.write_cg_results()

# %%
data_file_info.pivot_cg_table()

# %% This are the results from the exact algorithm
data_file_info_exact = summary_methods.DataFileInfo("./results_2021_07_12.csv")
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
data_file_info = summary_methods.DataFileInfo("./results_2023_07_01.csv")

# %%
data_file_info.create_templates()
# %%
data_file_info.summary_variable_heuristic()
# %%
