import os
from pathlib import Path
import re
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from shutil import copy
from tikzplotlib import save as tikz_save
plt.style.use('ggplot')

def mean_opt_func(x):
	s = x.sum()
	return s/(x > 0).sum()

class DataFileInfo:
    names_pricing_solvers = names_pricing_solvers = {"ArcTimeIndexed", "BddBackward", "TimeIndexed", "BddBackwardCycle" }
    number_instances_solved_by_oliveira = {(40, 2): 25, (40, 4): 25, (50, 2): 25, (50, 4): 25, (100, 2): 24, (100, 4): 22}
    
    def __init__(self, file_path):
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")
        self.file_path = file_path
        self.year, self.month, self.day = self.extract_date_info(file_path)
        self.data = pd.read_csv(self.file_path, index_col=False)
        self.calculate_extra_cols()
        self.best = self.data.sort_values(by=['tot_bb']).drop_duplicates(['Inst','n','m'],keep='first')
        self.workdir = Path.cwd()
        self.results = self.workdir.joinpath(Path("./results"))
        self.results_path = self.create_results_path(self.results, self.year, self.month, self.day)
        self.df_oliveira = pd.read_csv(self.workdir.joinpath("./oliveira_overall.csv"), index_col=False)
        self.df_all_opt = None
        self.all_df = None
        self.summary_exact_opt = None
        self.all_instances_cg = None
        self.list_all_solvers = self.data.pricing_solver.unique()
    
    def create_results_path(self, results, year, month, day):
        if results.exists() == False:
            os.mkdir(results)
        results_path = results.joinpath("./results_{}_{}_{}".format(year, month, day))
        print(results_path)
        if results_path.exists() == False:
            os.mkdir(results_path)
        return results_path

    def extract_date_info(self, file_name):
        match = re.search(r'(results|CG_overall)\_(20\d{2})\_(\d{2})\_(\d{2})\.csv', file_name)
        if match is None:
            raise ValueError(f"The file name {file_name} does not match the pattern '(results|CG_overall)_(20\d{2})\_(\d{2})\_(\d{2})\.csv'")
        
        year = match.group(2)
        month = match.group(3)
        day = match.group(4)
        return year, month, day

    def calculate_extra_cols(self):
        self.data['gap'] = (self.data['global_upper_bound'] - self.data['global_lower_bound']
               )/(self.data['global_lower_bound'] + 0.00001)
        self.data['opt'] = self.data['global_lower_bound'] == self.data['global_upper_bound']
        self.data['reduction'] = (self.data['first_size_graph'] -
                     self.data['size_after_reduced_cost'])/(self.data['first_size_graph'] + 0.000001)
        self.data['Inst'] = self.data.NameInstance.apply(
            lambda x:  int(re.search(r'.*\_(\d+)', x).group(1)))
    
    def write_cg_results(self):
        summary_grouped = self.data.groupby(['pricing_solver', 'n', 'm'])
        list_all_solvers = self.data.pricing_solver.unique()
        if not set(list_all_solvers).issubset(set(self.names_pricing_solvers)):
            raise ValueError(f"Not all solvers are in the data. The following solvers are not in the data: {set(list_all_solvers).difference(self.names_pricing_solvers)}")
        aggregation = {"tot_lb": {np.max, np.mean},
                    "tot_lb_root": {np.max, np.mean},
                    "gap": {np.max, np.mean},
                    "first_size_graph": {np.max, np.mean},
                    "size_after_reduced_cost": {np.max, np.mean},
                    "opt": np.sum,
                    "reduction": {np.max, np.mean},
                    "tot_cputime": {np.max, np.mean}}
        summary_write = summary_grouped.agg(aggregation).pivot_table(index=['n', 'm'], values=[
            'tot_lb', 'tot_lb_root', 'size_after_reduced_cost', 'gap', 'first_size_graph', 'reduction', 'opt'], columns=['pricing_solver'])
        summary_write.columns =  summary_write.columns.set_levels(
            list_all_solvers, level=2)
        summary_write.columns = ["_".join(x) for x in summary_write.columns.ravel()]
        summary_write.to_csv(self.results_path.joinpath(
            "CG_summary.csv"))
    
    def pivot_cg_table(self):
        self.all_instances_cg = self.data.pivot_table(values=['tot_lb', 'gap', 'first_size_graph', 'reduction', 'opt', 'rel_error', 'nb_generated_col',
                                         'global_lower_bound', 'global_upper_bound', 'tot_cputime', 'tot_bb','first_rel_error','nb_generated_col_root','nb_nodes_explored','size_after_reduced_cost','global_lowerbound_root','global_upper_bound_root','tot_lb_root'], index=['n', 'm', 'Inst'], columns=['pricing_solver'])
        list_all_solvers = self.data.pricing_solver.unique()
        self.all_instances_cg.columns = self.all_instances_cg.columns.set_levels(
            list_all_solvers, level=1)
        self.all_instances_cg.columns = ["_".join(x) for x in self.all_instances_cg.columns.ravel()]
        self.all_instances_cg.to_csv(
        self.results_path.joinpath("TUE_server_results_CG.csv"))
    
    def get_list_of_instances_solved_by_one_of_the_solvers(self, solver1, solver2):
        lst = self.all_instances_cg[(self.all_instances_cg[ 'opt_' + solver1] == 1) & (self.all_instances_cg['opt_' + solver2] == 0)]
        lst.to_csv(self.results_path.joinpath("list_of_instances_solved_by_{}_not_{}.csv".format(solver1, solver2)))
        lst = self.all_instances_cg[(self.all_instances_cg['opt_' + solver1] == 0) & (self.all_instances_cg[ 'opt_' + solver2] == 1)]
        lst.to_csv(self.results_path.joinpath("list_of_instances_solved_by_{}_not_{}.csv".format(solver2, solver1)))

        
    
    def create_templates(self):
        template_dir_path = self.workdir.joinpath("./template_dir")
        for lst in template_dir_path.iterdir():
            # if lst.name == "template_table.tex":
            #     copy(lst, self.results_path.joinpath(
            #         "CG_tables_{}_{}_{}.tex".format(self.year, self.month, self.day)))
            # else:
            copy(lst, self.results_path.joinpath(lst.name))
    
    
    def calculate_df_all_opt(self):
        for it in ['tot_real_time', 'tot_cputime', 'tot_bb', 'tot_lb', 'tot_lb_root', 'tot_heuristic', 'tot_build_dd', 'tot_pricing', 'tot_reduce_cost_fixing']:
            self.best["adjusted_{}".format(it)] = 0.7*self.best[it]
        df_oliveira_opt = self.df_oliveira[(self.df_oliveira['OptFound'] == 1)]
        best_opt = self.best[self.best['opt']]
        self.df_all_opt = pd.merge(best_opt, df_oliveira_opt, on=['Inst', 'n', 'm'])
    
    def calculate_all_df(self):
        self.all_df = pd.merge(self.best, self.df_oliveira,on=['Inst','n','m'],how='outer')
        weight_oliveira = self.all_df['OptFound'].fillna(0.0)
        self.all_df["opt_TimeOliveira"] = weight_oliveira * self.all_df["TimeOliveira"]
        self.all_df["opt_tot_bb"] = self.all_df["opt"]*self.all_df["adjusted_tot_bb"]
        self.all_df["opt_found_tot_bb"] = weight_oliveira*self.all_df["opt"]*self.all_df["adjusted_tot_bb"]
        self.all_df["opt_found"] = weight_oliveira*self.all_df["opt"]
    
    def write_all_instances(self):
        all_instances = self.best.pivot_table(values=['tot_lb', 'gap', 'first_size_graph', 'reduction', 'opt', 'rel_error', 'nb_generated_col',
                                                'global_lower_bound', 'global_upper_bound', 'tot_cputime', 'tot_bb','first_rel_error','nb_generated_col_root','nb_nodes_explored','size_after_reduced_cost','global_lowerbound_root','global_upper_bound_root','tot_lb_root'], index=['n', 'm', 'Inst'], columns=['pricing_solver'])
        all_instances.columns = all_instances.columns.set_levels(
            ['AFBC'], level=1)
        all_instances.columns = ["_".join(x) for x in all_instances.columns.ravel()]
        all_instances.to_csv(
            self.results_path.joinpath(
                "TUE_server_results.csv"))
    
    def write_summary_exact_opt(self):
        if self.all_df is None:
            self.calculate_all_df()
        agg = {"opt_tot_bb": {mean_opt_func}, "opt_TimeOliveira": {mean_opt_func}, "opt_found_tot_bb": {mean_opt_func}, "opt_found":{np.sum},"opt": {np.sum}}
        grouped = self.all_df.groupby(['n', 'm'])
        if self.summary_exact_opt is None:
            self.summary_exact_opt = grouped.agg(agg)
            self.summary_exact_opt.columns = ["_".join(x) for x in self.summary_exact_opt.columns.ravel()]
            self.summary_exact_opt["OptFound_sum"] = self.number_instances_solved_by_oliveira.values()
            self.summary_exact_opt.to_csv(self.results_path.joinpath(
                "TUEserver_summary_opt.csv"))
    
    def draw_overall_profile_curve(self):
        if self.df_all_opt is None:
            self.calculate_df_all_opt()
        
        self.df_all_opt['best_solver'] = self.df_all_opt[['adjusted_tot_bb', 'TimeOliveira']].min(axis=1)
        self.df_all_opt['ratio_tot_bb_best'] = self.df_all_opt['adjusted_tot_bb'] / \
            self.df_all_opt['best_solver']

        self.df_all_opt['ratio_TimeOliveira_best'] = self.df_all_opt['TimeOliveira'] / \
            self.df_all_opt['best_solver']

        sorted_ratio_tot_bb = self.df_all_opt[['ratio_tot_bb_best']
                                        ].sort_values(by='ratio_tot_bb_best').to_numpy()
        yvals = np.arange(len(sorted_ratio_tot_bb)) / \
            float(len(sorted_ratio_tot_bb) - 1.0)

        sorted_ratio_TimeOliveira = self.df_all_opt[['ratio_TimeOliveira_best']].sort_values(
            by='ratio_TimeOliveira_best').to_numpy()

        yvalues = np.arange(len(sorted_ratio_TimeOliveira)) / \
            float(len(sorted_ratio_TimeOliveira) - 1.0)

        width, height = plt.figaspect(1.68)
        fig, ax = plt.subplots(figsize=(width, height), dpi=200)
        ax.step(sorted_ratio_tot_bb, yvals, label='BDD')
        ax.step(sorted_ratio_TimeOliveira, yvalues, label='ATIF')
        ax.set_xlim([10.0**0, 20])
        ax.set_xlabel(r"$\tau$")
        ax.set_ylabel(r"$P(r_{p,s} \leq \tau)$")
        ax.legend(loc='lower right')
        name_file = 'profile_curve_overall.tex'
        tikz_save(self.results_path.joinpath(name_file))
        plt.savefig(self.results_path.joinpath('profile_curve_overall.pdf'), dpi=200)
    
    def draw_profile_curve_per_instance_class(self):
        if self.df_all_opt is None:
            self.calculate_df_all_opt()
        
        for n in [40, 50, 100]:
            for m in [2, 4]:
                sorted_ratio_tot_bb = self.df_all_opt.loc[(self.df_all_opt['n'] == n) & (
                    self.df_all_opt["m"] == m), "ratio_tot_bb_best"].sort_values().to_numpy()
                yvals = np.arange(len(sorted_ratio_tot_bb)) / \
                    float(len(sorted_ratio_tot_bb) - 1.0)

                sorted_ratio_TimeOliveira = self.df_all_opt.loc[(self.df_all_opt['n'] == n) & (
                    self.df_all_opt["m"] == m), "ratio_TimeOliveira_best"].sort_values().to_numpy()

                yvalues = np.arange(len(sorted_ratio_TimeOliveira)) / \
                    float(len(sorted_ratio_TimeOliveira) - 1.0)

                width, height = plt.figaspect(1.68)
                fig, ax = plt.subplots(figsize=(width, height), dpi=200)
                ax.step(sorted_ratio_tot_bb, yvals, label='BDD')
                ax.step(sorted_ratio_TimeOliveira, yvalues, label='ATIF')
                ax.set_xlim([10.0**0, 20])
                ax.set_title(
                    "Performance profile for instances with $m = {}$ and $n = $".format(m, n))
                ax.set_xlabel(r"$\tau$")
                ax.set_ylabel(r"$P(r_{p,s} \leq \tau)$")
                ax.legend(loc='lower right')
                name_file = 'profile_curve_{}_{}.tex'.format(
                    n, m)
                tikz_save(self.results_path.joinpath(name_file))
                plt.savefig(self.results_path.joinpath('profile_curve_{}_{}.pdf'.format(
                    n, m )), dpi=200)
    
    def list_of_instance_only_solved_by_bdd(self):
        if self.all_df is None:
            self.calculate_all_df()
        lst = self.all_df[(self.all_df['opt'] == 1) & (self.all_df['OptFound'] == 0)]
        # output = lst[['Inst', 'n', 'm', 'tot_bb', 'TimeOliveira', 'opt', 'OptFound']]
        output = lst[['Inst', 'n', 'm', 'tot_bb', 'TimeOliveira', 'opt', 'OptFound']]
        output.to_csv(self.results_path.joinpath("list_of_instance_only_solved_by_bdd.csv"), index=False)
    
    def summary_variable_heuristic(self):
        grouped = self.data.groupby(['ijoc_review','n','m'])
        lst = self.data['ijoc_review'].unique()
        aggregation = {"tot_lb": {np.max, np.mean},
                    "tot_lb_root": {np.max, np.mean},
                    "gap": {np.max, np.mean},
                    "first_size_graph": {np.max, np.mean},
                    "size_after_reduced_cost": {np.max, np.mean},
                    "opt": np.sum,
                    "tot_bb": {np.max, np.mean},
                    "reduction": {np.max, np.mean},
                    "tot_cputime": {np.max, np.mean}}
        summary_write = grouped.agg(aggregation).pivot_table(index=['n', 'm'], values=[
            'tot_lb', 'tot_lb_root', 'size_after_reduced_cost', 'gap', 'first_size_graph', 'reduction', 'opt', 'tot_bb'], columns=['ijoc_review'])
        summary_write.columns =  summary_write.columns.set_levels(
            lst, level=2)
        summary_write.columns = ["_".join(x) for x in summary_write.columns.ravel()]
        summary_write.to_csv(self.results_path.joinpath(
            "heuristic_variable_summary.csv"))
        
