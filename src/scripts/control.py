# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import re


# %%
df = pd.read_csv('../CG_overall_20200421.csv')


# %%
df.info()


# %%
df['Inst'] = df.NameInstance.apply(lambda x:  int(re.search(r'.*\_(\d+)', x).group(1)))


# %%
df_pessoa = pd.read_csv('../results/all_pessoa.csv')


# %%
df_pessoa.head()


# %%
df_pessoa.Opt= df_pessoa.Opt.apply(str)


# %%
df_pessoa['best'] = df_pessoa.apply(lambda x: re.search(r'[^0-9]?(\d+)', x['Opt']).group(1),axis=1)


# %%
df_pessoa.best = df_pessoa.best.apply(pd.to_numeric)


# %%
df_pessoa.info()


# %%
df_all = pd.merge(df, df_pessoa, on=['Inst', 'n','m'])


# %%
df_all[df_all['global_lower_bound'] > df_all['best']]



# %%
