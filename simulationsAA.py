#%%
import pandas as pd # type: ignore
import numpy as np
import seaborn as sns # type: ignore
import matplotlib.pyplot as plt
import os
import scipy.stats as ss
import scipy as sp

## ---
# Choose the best parameter
path = '/data/simulations/'
files = os.listdir(path)

# First create the file with the community growth rate
community = pd.DataFrame(columns=['Parameter','sample','Growth_rate'])
computeFluxes = 'data/simulations/compute.fluxes.040225_0.005.out'
with open(computeFluxes)as file:
    enter = False
    enter2 = False
    for line in file:   
        if 'Running (' in line:
            l = line.split('(')[1]
            parameter = l.split(')')[0].replace(', ','_')
        if '>' in line:
            sample = line.replace('>','').replace('\n','').strip()
            enter = True
            enter2 = False
        if '--- CoCo ---\n' == line:
            enter2 = True
        if enter2 == True and '----------------------------------\n' == line:
            gr = line0.replace('\n','')
            try:
                gr = float(gr)
                community.loc[len(community)] = (parameter,sample,gr)
                enter = False
                enter2 = False
            except:
                print (parameter, sample, gr)
        elif enter2 == True and 'Finished iter' in line:
            gr = line0.replace('\n','')
            try:
                gr = float(gr)
                community.loc[len(community)] = (parameter,sample,gr)
                enter = False
                enter2 = False
            except:
                print (parameter,sample,gr)
        elif enter2 == True:
            line0 = line

# plot the community and relate them with OD
lnOD = pd.read_csv(f'data/tables/ln_od.tsv',sep='\t', index_col=0)
cols = ['Control', 'Control.1', 'Control.2', 'AA', 'AA.1', 'AA.2','AA+Ab', 'AA+Ab.1','AA+Ab.2']
lnOD.drop([0,2],inplace=True)
lnlOD = lnOD[cols]
# Calculate growth rate for each combination
index = list(lnOD.index)

gr_dict = {}
grSimulation = {}
for col in cols:
    gr_table = []
    dftmp = lnOD[col]
    i = 0
    while i < len(index)-1:
        j = 0
        while j < len(index):
            if index[j] != index[i]:
                if j > i:
                    gr = (float(dftmp.loc[index[j]]) - float(dftmp.loc[index[i]]))/((index[j]-index[i])*24)
                    if index[j] == 21 and index[i] == 19:
                        grSimulation[col] = gr
                    gr_table.append(gr)
            j += 1
        i += 1
    gr_dict[col] = gr_table
            
gr_average = {}
for col in gr_dict.keys():
    
    gr = gr_dict[col]
    gr = gr[1:]
    gr = gr[1:]
    gr_i = np.mean(gr)
    gr_average[col] = gr_i

grPD = pd.DataFrame.from_dict(gr_average, orient='index')
samplesSynonims = {'Control_1_B':'Control', 'Control_2_B':'Control.1', 'Control_3_B':'Control.2', 'AA_1_B':'AA', 'AA_2_B':'AA.1','AA_3_B':'AA.2',
            'AAAb_1_B':'AA+Ab','AAAb_2_B':'AA+Ab.1','AAAb_3_B':'AA+Ab.2'}
listas = []
for par in community.Parameter.unique():
    dftmp = community[community['Parameter'] == par]
    dftmp2 = dftmp.copy()
    for i in dftmp.itertuples():
        sample = i[2]
        sample2 = samplesSynonims[sample]
        dftmp2.loc[i[0],'GR_OD'] = gr_average[sample2]
    listas.append(dftmp2)
# %%
dftmp = pd.DataFrame()
for l in listas:
    dftmp = pd.concat([dftmp,l])

def annotate(data, **kws):
    r, p = ss.pearsonr(data['Growth_rate'], data['GR_OD'])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
            transform=ax.transAxes)
# %%
# Create Correlation analysis between parameters
# Add the correlation of the Archaea, consider also the growth rate 
# should be proportional to the community since they have more than 85% of the total community
import scipy.stats as ss
corrParameters = pd.DataFrame(columns=['Parameter','CorrGRvsOD','pval1','CorrGRvsMt','pval1','CorrODvsMt','pval3','min-max','#samples'])
for parameter in dftmp.Parameter.unique():
    dfP = dftmp[dftmp['Parameter'] == parameter]
    d,g,ct,ct2,std = parameter.split("_")
	
    try: 
        dfGR = pd.read_csv(f'{path}coco_growth_rate_members_gamma.{g}.delta.{d}.ct_{ct}_aa_{std}std_blocked.none.csv')
        dfGR = dfGR.rename(columns={dfGR.columns[2]:'growth_rate'}) # Just in case the name is not the same
        idx = set(dfGR['Unnamed: 0'].unique()) - set(dfP['sample'].values)
        dfGR = dfGR[~dfGR['Unnamed: 0'].isin(idx)]
        z = dfGR[dfGR['compartments'] == 'maxbin2_spades_001'].growth_rate.values	
        x = dfP.Growth_rate.values
        y = dfP.GR_OD.values
        try:
            corr,pval = ss.pearsonr(x,y)
            corr2,pval2 = ss.pearsonr(x,z)
            corr3,pval3 = ss.pearsonr(y,z)
            a = (z*100/x)
            min_ = np.round(np.min(a),2)
            max_ = np.round(np.max(a),2)
            corrParameters.loc[len(corrParameters)] = (parameter,corr,pval,corr2,pval2,corr3,pval3,(min_,max_),len(x))
        except:
            pass
    except:
        print (f'coco_growth_rate_members_gamma.{g}.delta.{d}.ct_{ct}_aa_{std}std_blocked.none.csv')
# %%
corrParameters.to_clipboard()
# %%
best = '2.0_1.0_0.7_0.7_1.0'
# %%
# r, p = sp.stats.pearsonr(dftmp['Growth_rate'], dftmp['GR_OD'])
to_check = []
for parameter in dftmp.Parameter.unique():
    dftmp2 = dftmp[dftmp['Parameter'] == parameter]
    r,p = sp.stats.pearsonr(dftmp2['Growth_rate'], dftmp2['GR_OD'])
    if r > 0.5 and p < 0.1:
        to_check.append([parameter,np.round(r,2),np.round(p,4)])
# %%
# Graph the best
savepath='figures/'
df_best = dftmp[dftmp['Parameter'] == best]
g = sns.lmplot(x='Growth_rate', y='GR_OD', data=df_best, row='Parameter',
               height=3, aspect=1)

def annotate(data, **kws):
    r, p = sp.stats.pearsonr(data['Growth_rate'], data['GR_OD'])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
            transform=ax.transAxes)

g.map_dataframe(annotate)
g.savefig(f'{savepath}compute.fluxes.040225_0.005.best.pdf')
# %%
def CheckFluxes(cpd,eof='m',path=path,sp=None):
    files = os.listdir(path)
    indices = ["Control_1_B","Control_2_B","Control_3_B","AA_1_B","AA_2_B","AA_3_B","AAAb_1_B","AAAb_2_B","AAAb_3_B",]
    df = pd.DataFrame(index=indices)
    columnas = []
    if eof == 'm':
        check = "exchanges"
    elif eof == "f":
        check = "fluxes"
    for f in files:
        if f.endswith(".csv") and check in f:
            # print (f, check)
            if check == "exchanges":
                df2 = pd.read_csv(f"{path}{f}", sep=",", index_col=0)
            elif check == "fluxes":
                df2 = pd.read_csv(f"{path}{f}", sep=",", index_col=[0,1])
            else:
                break
            if f.startswith(f"coco_{check}_"):
                f2 = f.replace(f"coco_{check}_gamma.","").replace(".0.delta","").replace("0.ct_","").replace("_aa_","_").replace("std_blocked.none.csv","_c")
                
            elif f.startswith(f"base_{check}_"):
                f2 = f.replace(f"base_{check}_gamma.","").replace(".0.delta","").replace("0.ct_","").replace("_aa_","_").replace("std_blocked.none.csv","_b")

            columnas.append(f2)
            for i in indices:
                try:
                    if check == "exchanges":
                        c = df2.loc[cpd][i]
                    elif check == "fluxes":
                        c = df2.loc[i,sp][cpd]
                except:
                    c = None
                # print (i,f2,c)
                df.loc[i,f2] = c

        columnas.sort()
    return df[columnas]
# %%
##########################################################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import PowerNorm
import pandas as pd
import seaborn as sns

##########################################################################
d,g,ct,ct2,std = best.split("_")
df = pd.read_csv(f"{path}coco_fluxes_gamma.{g}.delta.{d}.ct_{ct}_aa_{std}std_blocked.none.csv", index_col=[0,1])
coverage = pd.read_csv('data/tables/coverage.mags.aminoacids.tsv',sep='\t',index_col=0)
aa_names_p = 'data/tables/AAModelSeed2.tsv'
mags_names_f = 'data/tables/names.list.tsv'

# Select only AA_1_B and AA_2_B
aa_samples = ['AA_1_B','AA_2_B','AA_3_B']
control_samples = ['Control_1_B','Control_2_B','Control_3_B']
aaab_samples = ['AAAb_1_B','AAAb_2_B']
samples = control_samples + aa_samples + aaab_samples
# Create a dictionary for the AA names
aa_seed = {}
aminoacids = []
with open(aa_names_p)as f:
    for line in f:
        seed = 'EX_'+line.split("\t")[0]+'_e'
        name = line.split("\t")[1].replace("\n","")
        aa_seed[seed] = name
        aminoacids.append(seed)

# Create a dictionary to replace the names of the MAGs
mags_names = {}
with open(mags_names_f)as f:
    for line in f:
        mag = line.split("\t")[0].replace(".","_")
        name = line.split("\t")[1].replace("\n","")
        mags_names[mag] = name

mags = set()
for samples in [aa_samples,control_samples,aaab_samples]:
    ms = list(coverage[samples].mean(axis=1).sort_values().tail(n=10).index)
    mags = mags | set(ms)

# Create functionts to plot the AAs    
def CreateTable(samples,MAGs,cols=''):
    # Select the most abundant mags in those samples
    sure_to_keep = set(['L-Leucine','L-Aspartate','L-Valine'])
    mags = dict(coverage.loc[list(MAGs)][samples].mean(axis=1))
    amino = list(set(aminoacids) & set(df.columns))
    df2 = pd.DataFrame(columns=amino)
    # Create the table for the plot
    for m in mags.keys():
        v = mags[m]
        m = m.replace(".","_")
        ss = []
        for s in aa_samples:
            if (s,m) in df.index:
                ss.append((s,m))
        dftmp = df.loc[ss][amino]
        df1 = pd.DataFrame(dftmp.mean(),columns=[m])*(v/100)
        df2 = pd.concat([df2,df1.T])
    df2 = df2.rename(columns=aa_seed,index=mags_names)
    if cols == '':
        filter = np.abs(df2).sum() > 0.001
        cols = df2.loc[:,filter].columns.tolist()
        cols = list(set(cols) | sure_to_keep)
    df3 = df2[cols]
    return df3

# Plot several plots
def CreatePlot(dataframes,names,number):
    fig, axs = plt.subplots(nrows=number)
    i = 0
    for df in dataframes:
        normalized_df = (df - df.mean()) / df.std()
        sns.clustermap(normalized_df, row_cluster=False,col_cluster=False,cmap="RdBu_r",ax=axs[i])
        i += 1

    return fig

# Print samples names
samples = []
for i in df.index:
    if i[0] not in samples:
        samples.append(i[0])

aa_samples = list(set(samples) & set(aa_samples))
control_samples = list(set(samples) & set(control_samples))
aaab_samples = list(set(samples) & set(aaab_samples))

print (samples)
# %%
# Who is consuming CO2 and H2 and producing acetate?
df_sorted = df.sort_index()
for sample in samples:
    dftmp = df_sorted.loc[sample,][['EX_cpd00011_e','EX_cpd11640_e','EX_cpd00029_e']]
    for i in dftmp.itertuples():
        co2 = np.round(i[1],2)
        h2 = np.round(i[2],2)
        acetate = np.round(i[3],2)
        if co2 < 0 and acetate > 0:
            print (i[0],mags_names[i[0]],co2,h2,acetate)
#%%
# Who is consuming Leucine? --> cpd00107
# Who is consuming Valine? --> cpd00156
# Who is consuming Cysteine? --> cpd00084
for sample in samples:
    dftmp = df_sorted.loc[sample,][['EX_cpd00107_e','EX_cpd00156_e','EX_cpd00084_e']]
    for i in dftmp.itertuples():
        leucine = np.round(i[1],2)
        valine = np.round(i[2],2)
        cysteine = np.round(i[3],2)
        if leucine < 0 or valine < 0 or cysteine < 0:
            print(i[0],mags_names[i[0]],leucine,valine,cysteine)
# %%
for sample in samples:
    print(df_sorted.loc[sample,'metabat2_spades_18'][['EX_cpd00011_e','EX_cpd00029_e']])

# %%
def cpd_modelSeed_to_name(cpd_list):
    seed_metabolites_df = pd.read_csv("data/db/seed_metabolites_edited.tsv", sep="\t", index_col=0)
    cpd_modelSeed_to_name_dict = dict(seed_metabolites_df.loc[cpd_list, 'name'])
    return cpd_modelSeed_to_name_dict

def rxn_modelSeed_to_name(rxn_list):
    seed_reactions_f = 'data/db/seed_reactions_corrected.tsv'
    seed_reactions_df = pd.read_csv(seed_reactions_f, sep='\t', index_col=0)
    return dict(seed_reactions_df.loc[rxn_list, 'name'])
# %%
# Load and normalize the fluxes of bacteria based on relative abundance
# Load exchanges and fluxes to medium
d,g,ct,ct2,std = best.split("_")
exchanges_sim_f = f'{path}coco_fluxes_gamma.{g}.delta.{d}.ct_{ct}_aa_{std}std_blocked.none.csv'
import_sim_f = f'{path}coco_exchanges_gamma.{g}.delta.{d}.ct_{ct}_aa_{std}std_blocked.none.csv'
exchanges_sim_df = pd.read_csv(exchanges_sim_f, index_col=[0,1])

samples = []
for idx in exchanges_sim_df.index:
    if idx[0] not in samples:
        samples.append(idx)

import_sim_df = pd.read_csv(import_sim_f, index_col=[0])
# load MAG relative abundances 
abundances_f = f"data/tables/coverage.mags.aminoacids.tsv"
abundances_df = pd.read_csv(abundances_f, sep="\t", index_col=0)
abundances_df = abundances_df/100
abundance_threshold = 0.00001
# get the abundance of each species in each sample in order 
abundances = []
considered_samples = abundances_df.columns[abundances_df.columns.str.contains("_B")]


for idx in exchanges_sim_df.index:
    sample,mag = idx
    mag2 = mag.replace("_",".")
    if 'concoct' in mag2:
        mag2 = mag
    exchanges_sim_df.loc[idx] = exchanges_sim_df.loc[idx] * abundances_df.loc[mag2][sample]

# sanity checks on normalization
cpd = "cpd01024" 
ex = exchanges_sim_df.columns[exchanges_sim_df.columns.str.contains("_e")]
print("""
      This is a sanity check, if the numbers are closed the normalization was correct
      The total import from the medium is {}
      The total exchanges are {}
      """.format(
                round(import_sim_df.loc["EX_"+cpd+"_m",considered_samples[0]], 2), 
                round(sum(exchanges_sim_df.loc[considered_samples[0],"EX_"+cpd+"_e"]), 2)))

# set to zero fluxes below the limit of sensitivity of the solver 
exchanges_sim_df[exchanges_sim_df.abs() < 1e-6] = 0
exchanges_sim_df = exchanges_sim_df.drop(columns = exchanges_sim_df.columns[exchanges_sim_df.isna().all()]) # remove rxn with not predicted fluxes
exchanges_sim_df = exchanges_sim_df.drop(columns = exchanges_sim_df.columns[(abs(exchanges_sim_df) < 1e-6).all()]) # remove rxn with fluxes under the sensitivity limit
# %%
exchanges_sim_e_df = exchanges_sim_df.loc[:,exchanges_sim_df.columns[exchanges_sim_df.columns.str.contains('_e')]]
dt_plot = exchanges_sim_e_df.loc[(slice(None), 'maxbin2_spades_001'), :].copy()
idx = dt_plot.abs().max(axis=0) > 0.01
cpd_to_drop = ['EX_cpd00001_e']

dt_plot.drop(columns=cpd_to_drop, inplace=True)
col_name = cpd_modelSeed_to_name([cpd.split('_')[1] for cpd in idx.index])
new_keys = ['EX_'+cpd+'_e' for cpd in col_name.keys()]
col_name = dict(zip(new_keys, col_name.values()))
dt_plot = dt_plot.loc[:, idx]
dt_plot.rename(columns=col_name, inplace=True)

dt_plot_control = dt_plot.loc[dt_plot.index[dt_plot.index.get_level_values(0).str.contains('Control_')], :]
dt_plot_aa = dt_plot.loc[dt_plot.index[dt_plot.index.get_level_values(0).str.contains('AA_')], :]
dt_plot_aaab = dt_plot.loc[dt_plot.index[dt_plot.index.get_level_values(0).str.contains('AAAb_')], :]

fig, axes = plt.subplots(1, 3, figsize=(10, 12.), tight_layout=True)
sns.barplot(data=dt_plot_control, orient='h', color='lightblue', ax=axes[0])
axes[0].set_title('Control')
axes[0].set_xlim(left=-10, right=7.5) 
sns.barplot(data=dt_plot_aa, orient='h', color='orange', ax=axes[1])
axes[1].set_title('AA')
axes[1].set_xlim(left=-10, right=7.5) 
sns.barplot(data=dt_plot_aaab, orient='h', color='grey', ax=axes[2])
axes[2].set_title('AA+Ab')
axes[2].set_xlim(left=-10, right=7.5)

plt.savefig(f'{savepath}mt.fluxes.pdf')
# %%
# Who es consuming most the acetate?
samples = import_sim_df.columns
acetate_df = pd.DataFrame(index=samples,columns=['TotalConsumed','PrimaryConsumer','Consumed','perc'])
for sample in samples:
	# acetate_S = import_sim_df.loc['EX_cpd00029_m'][sample]
	dftmp = exchanges_sim_df.loc[sample,]['EX_cpd00029_e'].dropna()
	sum_negative_acetate = dftmp[dftmp < 0].sum()
	mag = exchanges_sim_df.loc[sample,]['EX_cpd00029_e'].dropna().sort_values().index[0]
	val = exchanges_sim_df.loc[sample,]['EX_cpd00029_e'].dropna().sort_values()[0]
	perc = (val * 100)/ sum_negative_acetate
	acetate_df.loc[sample] = (sum_negative_acetate,mag,val,perc)
# %%
acetate_df # Here I corroborate acetate is consumed by mt and is the organism who most consumes this compounds
# %%
def sortDF(df):
    source = df.columns[0]
    target = df.columns[1]
    
    sources = list(df[source].unique())
    sources.sort()
    
    targets = list(df[target].unique())
    targets.sort()

    idx_left = []
    idx_right = []
    mets = []
    
    for s in sources:
        if s.startswith('EX'):
            if s not in mets:
                mets.append(s)
        else:
            idx_left = idx_left + list(df[df[source] == s].index)
    
    for t in targets:
        if t.startswith('EX'):
            if t not in mets:
                mets.append(t)
            
        else:
            idx_right = idx_right + list(df[df[target] == t].index)
    
    idx = idx_left + idx_right
    return idx
		
def pl_sankey(df, label_color, categories, value, colores, colorX, title='Sankey Diagram', fname=None, width=3000, height=1600, scale=2):
    from IPython.display import Image
    import plotly.graph_objects as go
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib

    df = df.copy()
    metabolites = df['Metabolite'].unique()
    labels = []
    colors = []

    # associate labels to colors
    for k, v in label_color.items():
        labels += [k]
        colors += [v]

    # transform df into a source-target pair
    st_df = None

    for i in range(len(categories)-1):
        _st_df = df[[categories[i],categories[i+1],value,colorX]]
        _st_df.columns = ['source', 'target', 'count','Color']
        st_df = pd.concat([st_df, _st_df])

    st_df = st_df.reset_index().drop('index', axis=1)

    for i in st_df.itertuples():
        idx = i[0]
        s = i[1]
        t = i[2]
        if s in metabolites:
            t = t.replace("_o", "")
            c = colores[t]
            st_df.loc[idx, 'Color'] = c
        else:
            c = colores[s]
            st_df.loc[idx, 'Color'] = colores[s]

    # sort
    st_df = st_df.sort_values(by=['source', 'target'])
    idxs = sortDF(st_df)
    print(len(idxs),len(st_df))
    st_df = st_df.loc[idxs]

    labels = list(st_df['source'].unique()) + list(st_df['target'].unique())
    labels = list(set(labels))  # Eliminar duplicados y preservar el orden

    # add index for source-target pair
    st_df['sourceID'] = st_df['source'].apply(lambda x: labels.index(str(x)))
    st_df['targetID'] = st_df['target'].apply(lambda x: labels.index(str(x)))

    # Creating the sankey diagram
    data = dict(
        type='sankey',
        node=dict(
            pad=15, thickness=20, line=dict(color='black', width=0.5), label=labels, color=colors,
        ),
        link=dict(source=st_df['sourceID'], target=st_df['targetID'], value=st_df['count'], color=st_df['Color']),
    )
    
    layout = dict(title=title, font=dict(size=16, family='Arial'))

    # Creating figure
    fig = go.Figure(dict(data=[data], layout=layout))
    
    if fname:
        fig.write_image(f'{fname}.svg', format='svg', width=width, height=height, scale=scale)
    
    return Image(fig.to_image(format='png', width=width, height=height, scale=scale))
# %%
import matplotlib.pyplot as plt
import seaborn as sns
import random
import pandas as pd
import numpy as np

it = []
control = samples[0:3]
AAc = samples[3:6]
AAAb = samples[6:]
ss = [control,AAc,AAAb]
col_name = cpd_modelSeed_to_name([cpd.split('_')[1] for cpd in idx.index])

cpd_to_drop = ['EX_cpd00001_e','EX_cpd00002_e','EX_cpd00003_e','EX_cpd00004_e','EX_cpd00005_e','EX_cpd00006_e',\
	'EX_cpd00007_e','EX_cpd00008_e','EX_cpd00009_e','EX_cpd00010_e']

interactionsSankey = pd.DataFrame(columns=['MAG','CPD','Name','Flux','Sample'])
sampleNames = ['control','AAc','AAAb']
i = 0
for s in  ss:
    sample = sampleNames[i]
    muestras = list(s)
    dfS = exchanges_sim_df.loc[muestras].groupby(level=1).mean()
    cols = []
    for c in dfS.columns:
        if '_e' in c:
            cols.append(c)
    dfSE = dfS[cols].dropna(how='all',axis=1)
    dfSE[np.isnan(dfSE)] = 0
    dfSE[np.abs(dfSE) < 1e-6] = 0
    dfSE = dfSE[dfSE != 0].dropna(how='all',axis=1)
    dfSE = dfSE.loc[:, ~dfSE.columns.isin(cpd_to_drop)]

    cpdsToKeep = set()
    for mag in dfSE.index:
        consumed = set(dfSE.loc[mag].dropna().sort_values().head(n=5).index)
        produced = set(dfSE.loc[mag].dropna().sort_values().tail(n=5).index)
        cpdsToKeep = consumed | produced | cpdsToKeep
	
    dfSE = dfSE[list(cpdsToKeep)]

    for mag in dfSE.index:
        for cpd in dfSE.columns:
            name = col_name[cpd.replace('EX_','').replace('_e','')]
            val = dfSE.loc[mag,cpd]
            if ~np.isnan(val):
                interactionsSankey.loc[(len(interactionsSankey))] = (mag,cpd,name,val,sample)

    print (f'Total cpds for sample {sample}: {len(cpdsToKeep)}')
    i += 1
# %%
totalMetabolites = set()
cpd_to_drop2 = []
for x in cpd_to_drop:
    cpd_to_drop2.append(x.replace("_e","_m"))

for sample in sampleNames:
    mediumES = import_sim_df[samples[samples.str.contains(sample)]].mean(axis=1).dropna()
    mediumES = mediumES.loc[~mediumES.index.isin(cpd_to_drop2)]
    mostConsumed = set(mediumES[np.abs(mediumES) > 0].sort_values().tail(n=10).index)
    mostProduced = set(mediumES[np.abs(mediumES) > 0].sort_values().head(n=10).index)
    totalMetabolites = mostConsumed | mostProduced | totalMetabolites

# %%
# Change the name of the cpds
cpds = []
for cpd in list(totalMetabolites):
    cpds.append(cpd.replace("_m","_e"))
interactions = interactionsSankey[interactionsSankey['CPD'].isin(cpds)]
# %%
# Add the names of the MAGs
for mag in interactions.MAG.unique():
    idxs = list(interactions[interactions['MAG'] == mag].index)
    interactions.loc[idxs,'MAGName'] = mags_names[mag]
#%%
# Create Gephi Table
gephi = pd.DataFrame(columns=['Source','Target','Flux','Direction'])

for i in interactions.itertuples():
    idx,mag,cpd,name,flux,sample,magname = i
    if flux > 0: # Source is the MAG and target the medium
        gephi.loc[len(gephi)] = (magname,name,flux,'Production')
    elif flux < 0:
        gephi.loc[len(gephi)] = (name,magname,flux,'Consumption')
# Create a Weight for the flux
gephi['Weight1'] = np.abs(gephi['Flux']) / np.max(np.abs(gephi['Flux']))
# %%
gephi.to_csv('data/tables/interactions_syntrophy_AA.tsv',sep='\t',index=False)
# %%
# ADD Nodes tables
# identify the consumers and producers of co2
mags = []
co2Consumer = []
cpds = []
for i in interactions.itertuples():
    idx,mag,cpd,name,flux,sample,magname = i
    if magname not in mags:
        mags.append(magname)
    if name not in cpds:
        cpds.append(name)
    if name == 'CO2':
        if flux < 0:
            co2Consumer.append(magname)

nodes = pd.DataFrame(columns=['ID','Latitude','Longitude'])
longitudA = 0
latitudA = -40

longitudB = 0
latitudB = 40

longitudC = 30
latitudC = 0
done = []
for node in list(set(gephi.Source.unique()) | set(gephi.Target.unique())):
    if node not in done:
        done.append(node)
        if node in mags:
            if node in co2Consumer:
                nodes.loc[len(nodes)] = (node,latitudA,longitudA)
                longitudA += 20
                if latitudA == -40:
                    latitudA = -45
                else:
                    latitudA = -40

            else:
                nodes.loc[len(nodes)] = (node,latitudB,longitudB)
                longitudB += 20
                if latitudB == 40:
                    latitudB = 45
                else:
                    latitudB = 40

        else:
            nodes.loc[len(nodes)] = (node,latitudC,longitudC)
            longitudC += 20
            if latitudC == 0:
                latitudC = -5
            else:
                latitudC = 0
nodes.to_csv('data/tables/nodes.tsv',sep='\t',index=False)
#%%
AA_ids = ['cpd00035',
          'cpd00051',
          'cpd00132',
          'cpd00041',
          'cpd00084',
          'cpd00023',
          'cpd00053',
          'cpd00033',
          'cpd00119',
          'cpd00322',
          'cpd00107',
          'cpd00039',
          'cpd00060',
          'cpd00066',
          'cpd00129',
          'cpd00054',
          'cpd00161',
          'cpd00065',
          'cpd00069',
          'cpd00156',
          'cpd00281']

# Which are the most consumed - produced compounds?
sumImport = import_sim_df.mean(axis=1)
sumImport = sumImport[~sumImport.index.isin(cpd_to_drop2)]
sumImport = sumImport[np.abs(sumImport) > 0]

def cpd_modelSeed_to_name(cpd_list,change='medium'):
    seed_metabolites_df = pd.read_csv("data/db/seed_metabolites_edited.tsv", sep="\t", index_col=0)
    # Remove EX_ _m or _e
    cpd_list2 = []
    for cpd in cpd_list:
        cpd_list2.append(cpd.replace("EX_","").replace("_m","").replace("_e",""))
    cpd_modelSeed_to_name_dict = dict(seed_metabolites_df.loc[cpd_list2, 'name'])
    cpd_modelSeed_to_name_dict2 = {}
    for cpd in cpd_modelSeed_to_name_dict.keys():
        name = cpd_modelSeed_to_name_dict[cpd]
        if change == 'medium':
            cpd = 'EX_'+cpd+'_m'
        elif change == 'exchange':
            cpd = 'EX_'+cpd+'_e'
        cpd_modelSeed_to_name_dict2[cpd] = name
    return cpd_modelSeed_to_name_dict2
cpdNames = cpd_modelSeed_to_name(list(sumImport.index))
sumImport.rename(index=cpdNames,inplace=True)

# %%
exchanges_mean_df = pd.DataFrame()
samples_dict = {'control':['Control_1_B','Control_2_B','Control_3_B'],
                'AA':['AA_1_B', 'AA_2_B','AA_3_B',],
                 'AAAb':['AAAb_1_B', 'AAAb_2_B']}
mostAbundant = {}
for sample in samples_dict.keys():
    samplesNames = samples_dict[sample]
    dftmp = exchanges_sim_e_df.loc[samplesNames,]#.dropna(axis=1,how='all')
    dftmp = dftmp.iloc[:,~dftmp.columns.isin(cpd_to_drop)]
    dftmp.dropna(axis=1,how='all',inplace=True)
    dftmp = dftmp.reset_index()
    dftmp = dftmp.fillna(0)
    dftmpMean = dftmp.drop('level_0',axis=1).groupby(['compartment']).mean()
    # most abundantant  
    mostAbundant[sample] = set(list(np.abs(dftmpMean).sum().sort_values().tail(n=115).index))
    dftmpMean.reset_index(inplace=True)
    dftmpMean['Sample'] = sample
    exchanges_mean_df = pd.concat([exchanges_mean_df,dftmpMean],ignore_index=True)
#%%
mets = list(mostAbundant['control'] | mostAbundant['AA'] | mostAbundant['AAAb'])
names = cpd_modelSeed_to_name(mets,'exchange')
met_names = []
met_names_dict = {}
for met in mets:
    met_names.append(names[met])
    met_names_dict[met] = names[met]
plot_mets = pd.Series(mets,index=met_names)

mets.append('compartment')
mets.append('Sample')

exchanges_mean_subset_df = exchanges_mean_df[mets]
exchanges_mean_subset_df = exchanges_mean_subset_df.set_index(['Sample','compartment'])
exchanges_mean_subset_df[exchanges_mean_subset_df.abs() < 1e-6] = 0

# set MAG order from most to least abundant
plot_MAGs = pd.Series(exchanges_mean_subset_df.reset_index().compartment.unique())
coverage2 = coverage.copy()
for mag in coverage.index:
    mag2 = mag.replace(".","_")
    coverage2.rename(index={mag:mag2},inplace=True)
plot_abundances = pd.Series(coverage2.mean(axis=1).loc[plot_MAGs])
plot_MAGs = plot_MAGs.iloc[plot_abundances.argsort()[::-1]]
colorscale = [[0, '#b3003b'], [0.5, '#994e83'], [1, '#4775d1']]

plot_df_A = exchanges_mean_subset_df.loc['control',].dropna(axis=1, how='all')
plot_df_B = exchanges_mean_subset_df.loc['AA',].dropna(axis=1, how='all')
plot_df_C = exchanges_mean_subset_df.loc['AAAb',].dropna(axis=1, how='all')

representation_df = exchanges_mean_subset_df.rename(columns=met_names_dict)
representation_df
#%%
# Normalizar los flujos
df = representation_df.copy()
normalized_df = pd.DataFrame(index=df.index,columns=df.columns)
for sample,row in df.iterrows():
    for cpd, flux in row.items():
        if flux > 0:
            norm_flux = float(np.log(flux+1))
        elif flux < 0:
            norm_flux = float(-np.log(-flux+1))
        else:
            norm_flux = np.nan
        normalized_df.loc[sample,cpd] = norm_flux
# Create heatmap
cols = []
for cpd in normalized_df.columns:
    cols.append(cpd+'_1')
    cols.append(cpd+'_2')
    cols.append(cpd+'_3')

samples_names = list(exchanges_mean_df['Sample'].unique())
mags = list(exchanges_mean_df['compartment'].unique())
normalized_df_heatmap = pd.DataFrame(columns=cols,index=list(mags_names.values()))

for sampleMAG,row in normalized_df.iterrows():
    sample, bin_ = sampleMAG
    mag = mags_names[bin_]
    for cpd,flux in row.items():
        if sample == 'control':
            cpd_to_plot = cpd + '_1'
        elif sample == 'AAAb':
            cpd_to_plot = cpd + '_3'
        else:
            cpd_to_plot = cpd+ '_2'
        normalized_df_heatmap.loc[mag,cpd_to_plot] = flux
normalized_df_heatmap         
#%%
normalized_df_heatmap = normalized_df_heatmap.apply(pd.to_numeric, errors='coerce')
most_abundant_mags = list(abundances_df.mean(axis=1).sort_values().tail(n=32).index)
most_abundant_mags2 = []
for mag in most_abundant_mags:
    mag2 = mags_names[mag.replace(".","_")]
    most_abundant_mags2.append(mag2)
normalized_df_heatmap2 = normalized_df_heatmap.loc[most_abundant_mags2]
#%%
normalized_df_heatmap2 = normalized_df_heatmap2.apply(pd.to_numeric, errors='coerce')
sns.heatmap(normalized_df_heatmap2, cmap="RdYlGn", mask=normalized_df_heatmap2.isnull())
#%%
aminoacids = []
for aa_ in cpd_modelSeed_to_name(AA_ids).values():
    aminoacids.append(aa_+'_1')
    aminoacids.append(aa_+'_2')
    aminoacids.append(aa_+'_3')
normalized_df_heatmap_aa = normalized_df_heatmap2[aminoacids]
normalized_df_heatmap3 = normalized_df_heatmap2[normalized_df_heatmap2.columns[~normalized_df_heatmap2.columns.isin(aminoacids)]]

#%%
cpds_maximum = {}
columns = normalized_df_heatmap3.columns
for i,col in enumerate(columns):
    if i == 0 or i%3 == 0:
        value = np.abs(normalized_df_heatmap3[columns[i:i+3]]).mean(axis=1).max()
        if np.isnan(value) == False:
            cpds_maximum[col.split("_")[0]] = value
cpds_maximum = {k: v for k, v in sorted(cpds_maximum.items(), key=lambda item: item[1])}

# get the same order for the most abundant 
for_heatmap = []
for c in list(cpds_maximum.keys())[-18:]:
    for_heatmap.append(c+'_1')
    for_heatmap.append(c+'_2')
    for_heatmap.append(c+'_3')

normalized_df_heatmap4 = normalized_df_heatmap3[for_heatmap]
filtered_df = normalized_df_heatmap4.copy()
#%%

############################################################
############################################################
############################################################
#################### HEATMAP GRAPH #########################
############################################################
############################################################
############################################################
# Plot the aminoacids
mags_to_plot = ['T. xylanilyticum PD19',
'Tepidanaerobacteraceae sp. PD32',
'Firmicutes sp. PD7',
'Firmicutes sp. PD31',
'Tepidimicrobiaceae sp. PD23',
'Tepidanaerobacteraceae sp. PD3',
'Limnochordia sp. PD22',
'Carbobacillus sp. PD11',
'Lutisporaceae sp. PD17',
'Tepidanaerobacteraceae sp. PD27',
'Peptococcales sp. PD24',
'Moorellia sp. PD5',
'Tepidimicrobiaceae sp. PD26',
'Symbiobacterium sp. PD15',
'Limnochordia sp. PD14',
'Tissierellaceae sp. PD12',
'C. subterraneus PD13',
'M. thermautotrophicus PD4']
normalized_df_heatmap_aa_plot = normalized_df_heatmap_aa.loc[mags_to_plot[::-1]]
# Now i transpose the table in order to have it like the rest of the figure
magsTotalFluxes = []
for mag in normalized_df_heatmap_aa_plot.index:
    magsTotalFluxes.append(mag+'_1')
    magsTotalFluxes.append(mag+'_2')
    magsTotalFluxes.append(mag+'_3')

columnsToPlot = []
for col in normalized_df_heatmap_aa_plot.columns:
    col = col.split("_")[0]
    if col not in columnsToPlot:
        columnsToPlot.append(col)

normalized_df_heatmap_aa_plot2 = pd.DataFrame(index=columnsToPlot,columns=magsTotalFluxes)
normalized_df_heatmap_aa_plot2.fillna(0,inplace=True)
for mag in normalized_df_heatmap_aa_plot.index:
    for col in normalized_df_heatmap_aa_plot.columns:
        value = normalized_df_heatmap_aa_plot.loc[mag][col].astype('object')
        mag2 = mag+'_'+col.split("_")[1]
        normalized_df_heatmap_aa_plot2.loc[col.split("_")[0],mag2] = value

fig, ax = plt.subplots(figsize=(12, 5))

# Create the heatmap
sns.heatmap(
    normalized_df_heatmap_aa_plot2,
    cmap="RdBu_r",
    linewidths=0.5,
    linecolor="gray",
    center=0, vmin=-0.75,vmax=0.75,
    cbar_kws={"shrink": 0.8},
    ax=ax
)

# Add vertical lines every three columns to separate conditions
for i in range(3, normalized_df_heatmap_aa_plot2.shape[1], 3):
    ax.vlines(i, 0, normalized_df_heatmap_aa_plot2.shape[0], color="gray", linestyle="--", linewidth=2)

# Adjust x-axis ticks
xtick_positions = np.arange(1, normalized_df_heatmap_aa_plot2.shape[1], 3)  # Center of each triplet
xtick_labels = [col.split('_')[0] for col in normalized_df_heatmap_aa_plot2.columns[::3]]  # Compound base names
ax.set_xticks(xtick_positions + 0.5)  # Center ticks
ax.set_xticklabels(xtick_labels, rotation=90, fontsize=10)

# Adjust y-axis ticks
ytick_positions = np.arange(normalized_df_heatmap_aa_plot2.shape[0])
ax.set_yticks(ytick_positions + 0.5)  # Center ticks

# Save the figure
plt.savefig(f"{savepath}flux_exchanges_transposed_aa.svg")
plt.savefig(f"{savepath}flux_exchanges_transposed_aa.pdf")


#%%
############################################################
############################################################
############################################################
#################### HEATMAP GRAPH 2 #######################
############################################################
############################################################
############################################################
normalized_df_heatmap_all_plot = filtered_df.loc[mags_to_plot[::-1]]
fig, ax = plt.subplots(figsize=(12, 5))

# Create the heatmap
sns.heatmap(
    normalized_df_heatmap_all_plot,
    cmap="RdBu_r",
    linewidths=0.5,
    linecolor="gray",
    center=0, vmin=-1.15,vmax=1.15,
    cbar_kws={"shrink": 0.8},
    ax=ax
)

# Add vertical lines every three columns to separate conditions
for i in range(3, normalized_df_heatmap_all_plot.shape[1], 3):
    ax.vlines(i, 0, normalized_df_heatmap_all_plot.shape[0], color="gray", linestyle="--", linewidth=2)

# Adjust x-axis ticks
xtick_positions = np.arange(1, normalized_df_heatmap_all_plot.shape[1], 3)  # Center of each triplet
xtick_labels = [col.split('_')[0] for col in normalized_df_heatmap_all_plot.columns[::3]]  # Compound base names
ax.set_xticks(xtick_positions + 0.5)  # Center ticks
ax.set_xticklabels(xtick_labels, rotation=45, fontsize=10)

# Adjust y-axis ticks
ytick_positions = np.arange(normalized_df_heatmap_all_plot.shape[0])
ax.set_yticks(ytick_positions + 0.5)  # Center ticks
# Save the figure
plt.savefig(f"{savepath}flux_exchanges_transposed_all.svg")
plt.savefig(f"{savepath}flux_exchanges_transposed_all.pdf")
#%%
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Plot heatmap
fig, ax = plt.subplots(figsize=(12, 5))
sns.heatmap(
    normalized_df_heatmap_all_plot,
    cmap="RdBu_r",
    vmin=-1.5, vmax=1.5,
    linewidths=0.5,
    linecolor='gray',
    center=0,
    cbar_kws={"shrink": 0.8},  # Shrink colorbar
    ax=ax
)

# Add thick gridlines every 3 squares
for i in range(3, normalized_df_heatmap_all_plot.shape[1],3):
    ax.vlines(i, 0,normalized_df_heatmap_all_plot.shape[0], color='grey', linestyle='--', linewidth=2)

xtick_positions = np.arange(1, normalized_df_heatmap_all_plot.shape[1], 3)
ax.set_xticks(xtick_positions)  # Center ticks
ax.set_xticklabels([f'{normalized_df_heatmap_all_plot.columns[i].split("_")[0]}' for i in xtick_positions], 
                   rotation=45, fontsize=10)

# Custom y-ticks
ytick_positions = np.arange(normalized_df_heatmap_all_plot.shape[0])
ax.set_yticks(ytick_positions + 0.5)  # Center ticks

plt.savefig(f'{savepath}flux_exchanges.svg')
plt.savefig(f'{savepath}flux_exchanges.pdf')

#%%
producers = ['Environment']*len(mets)
metabolites = []
counts = []
for cpd in exchanges_mean_subset_df.columns:
    counts.append(np.abs(plot_df_A[cpd]).sum())

colors = [0, 0.5, 1,  0, 0.5, 1] # A -> 0, C -> 0.5, B -> 1
for m1 in plot_MAGs:
    for m2 in plot_mets:
        # A
        if m1 in plot_df_A.index:
            flux = plot_df_A.loc[m1, m2]
            if flux > 0:
                producers.append(m1)
                metabolites.append(plot_mets.index[plot_mets==m2][0])
                counts.append(flux)
                colors.append(0)
        # C
        if m1 in plot_df_C.index:
            flux = plot_df_C.loc[m1, m2]
            if flux > 0:
                producers.append(m1)
                metabolites.append(plot_mets.index[plot_mets==m2][0])
                counts.append(flux)
                colors.append(0.5)
        # B
        if m1 in plot_df_B.index:
            flux = plot_df_B.loc[m1, m2]
            if flux > 0:
                producers.append(m1)
                metabolites.append(plot_mets.index[plot_mets==m2][0])
                counts.append(flux)
                colors.append(1)

fig = go.Figure(go.Parcats(
    dimensions=[
        {'label': 'Producer',
         'values': producers},
        {'label': 'Metabolite',
         'values': metabolites}],
    counts=counts,
    line={'color': colors, 'colorscale': colorscale},
    labelfont={'size': 8, 'family': 'Arial'},
    tickfont={'size': 8, 'family': 'Arial'}),
    go.Layout(margin=dict(t=1, b=3)))
left_height = len(np.unique(producers))

fig.update_layout(width=380, height=313)
fig.show()
fig.write_image("data/tables/sankey_producers.pdf")
#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = exchanges_mean_subset_df.copy()
df['Sample_compartment'] = df.index.get_level_values('Sample') + '-' + df.index.get_level_values('compartment')
df.reset_index(drop=True, inplace=True)

# Set the new combined column as index.
df.set_index('Sample_compartment', inplace=True)

# Prepare the data for the heatmap
heatmap_data = df  # Rows: Sample-compartment, Columns: EX_cpd...

# Plot the heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(
    heatmap_data,
    cmap='coolwarm',  # Choose a color map
    center=0,         # Center around 0
    annot=False,      # Set True if annotations are needed
    cbar_kws={'label': 'Flux Value'}  # Add a label to the color bar
)

# Customize labels and title
plt.title("Heatmap of Exchange Reactions by Sample and Compartment")
plt.xlabel("Exchange Reactions")
plt.ylabel("Sample-Compartment")
plt.tight_layout()
plt.show()
# %%
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

# Change the plot according to each sample
interactions = interactions.reset_index()
interactions = interactions.drop('index',axis=1)

# Create medium for import and export
samples = ['control','AAc','AAAb']
for cpd in interactions.CPD.unique():
    cpdFlux = interactions[interactions['CPD'] == cpd]
    name = cpdFlux.Name.values[0]

    importFlux = cpdFlux[cpdFlux['Flux'] < 0]
    exportFlux = cpdFlux[cpdFlux['Flux'] > 0]
    for sample in samples:
        importFlux_sum = importFlux[importFlux['Sample'] == sample].Flux.sum()
        exportFlux_sum = exportFlux[exportFlux['Sample'] == sample].Flux.sum()
        importFlux_sum *= -1
        exportFlux_sum *= -1
        
        interactions.loc[len(interactions)] = ('medium',cpd,name,importFlux_sum,sample,'Medium')
        interactions.loc[len(interactions)] = ('medium',cpd,name,exportFlux_sum,sample,'Medium')

importI_df = interactions[interactions['Flux'] < 0]
exportI_df = interactions[interactions['Flux'] > 0]

#%%
importI_df['Flux'] = importI_df.Flux * -5
importI_df = importI_df.loc[importI_df.index.repeat(importI_df['Flux'].astype(int))]

exportI_df['Flux'] = exportI_df.Flux * 5
exportI_df = exportI_df.loc[exportI_df.index.repeat(exportI_df['Flux'].astype(int))]

# Crear las dimensiones
mag_dim = go.parcats.Dimension(
    values=exportI_df['MAGName'], label="Source"
)

compound_dim = go.parcats.Dimension(
    values=exportI_df['Name'], label="Compound"
)

sample_dim = go.parcats.Dimension(
    values=exportI_df['Sample'], label="Sample"
)
line_thickness = exportI_df['Flux'] / exportI_df['Flux'].max()  # Normalize by flux 
color_mapping = {'control': 0, 'AAc': 1, 'AAAb': 2}
exportI_df['Color'] = exportI_df['Sample'].map(color_mapping)

colorscale = [[0, 'lightblue'], [0.5, 'orange'], [1, 'purple']]  # colors for the samples

fig = go.Figure(data=[go.Parcats(
    dimensions=[mag_dim, compound_dim, sample_dim],
    line={'color': exportI_df['Color'], 'colorscale': colorscale},#, 'width': line_thickness * 10},  
    hoveron='color',
    hoverinfo='count+probability',
    labelfont={'size': 18, 'family': 'Times'},
    tickfont={'size': 16, 'family': 'Times'},
    arrangement='freeform'
)])

fig.update_layout(title_text="Flux Export")
fig.show()

#%%
# Dimensions
mag_dim = go.parcats.Dimension(
    values=importI_df['MAGName'], label="Source"
)

compound_dim = go.parcats.Dimension(
    values=importI_df['Name'], label="Compound"
)

sample_dim = go.parcats.Dimension(
    values=importI_df['Sample'], label="Sample"
)

# Definir el color basado en las muestras
color_mapping = {'control': 0, 'AAc': 1, 'AAAb': 2}
importI_df['Color'] = importI_df['Sample'].map(color_mapping)

# Crear el gráfico parallel categories
colorscale = [[0, 'lightblue'], [0.5, 'orange'], [1, 'purple']]  # Colores para las muestras

fig = go.Figure(data=[go.Parcats(
    dimensions=[sample_dim,compound_dim, mag_dim],
    line={'color': importI_df['Color'], 'colorscale': colorscale},
    hoveron='color',
    hoverinfo='count+probability',
    labelfont={'size': 18, 'family': 'Times'},
    tickfont={'size': 16, 'family': 'Times'},
    arrangement='freeform'
)])

fig.update_layout(title_text="Flux Import")
fig.show()
#%%
# Crear las dimensiones para el gráfico
mag_import_dim = go.parcats.Dimension(
    values=importI_df['MAGName'], label="Import"
)

compound_dim = go.parcats.Dimension(
    values=interactions['Name'], label="Compound"
)

mag_export_dim = go.parcats.Dimension(
    values=exportI_df['MAGName'], label="Export"
)

# Color for each sample
color_mapping = {'control': 0, 'AAc': 1, 'AAAb': 2}
idxs = list(set(importI_df.index) | set(exportI_df.index))
interactions2 = interactions.loc[idxs]
interactions2['Flux'] = np.abs(interactions2['Flux'])
interactions2['Color'] = interactions2['Sample'].map(color_mapping)

# Line
line_thickness = interactions2['Flux'] / interactions2['Flux'].max()  # Normalizando Flux para grosor

# Crear el gráfico Parallel Categories
colorscale = [[0, 'lightblue'], [0.5, 'orange'], [1, 'purple']]  # Colores para las muestras

fig = go.Figure(data=[go.Parcats(
    dimensions=[mag_export_dim, compound_dim, mag_import_dim],
    line={'color': interactions2['Color'], 'colorscale': colorscale},#, 'width': line_thickness * 10},
    hoveron='color',
    hoverinfo='count+probability',
    labelfont={'size': 18, 'family': 'Times'},
    tickfont={'size': 16, 'family': 'Times'},
    arrangement='freeform'
)])

fig.update_layout(title_text="Main Fluxes")
fig.show()

#%%
interaction2 = interactions.copy()
interaction2['Source'] = ''
for i in interactions.itertuples():
    idx = i[0]
    flux = i[4]
    if flux > 0:
        interaction2.loc[idx,'Source'] = 'Export'
    elif flux < 0:
        interaction2.loc[idx,'Source'] = 'Import'
interaction2['Flux'] = np.abs(interaction2.Flux) * 5
interaction2 = interaction2.loc[interaction2.index.repeat(interaction2['Flux'].astype(int))]

# Crear las dimensiones
export_dim = go.parcats.Dimension(
    values=interaction2[interaction2['Source'] =='Export'].MAGName, label="Export"
)

compound_dim = go.parcats.Dimension(
    values=interactions2['Name'], label="Compound"
)

import_dim = go.parcats.Dimension(
    values=interaction2[interaction2['Source'] =='Import'].MAGName, label="Import"
)
line_thickness = interaction2['Flux'] / interaction2['Flux'].max()  # Normalizando Flux para grosor
# Definir el color basado en las muestras
color_mapping = {'control': 0, 'AAc': 1, 'AAAb': 2}
interaction2['Color'] = interaction2['Sample'].map(color_mapping)

# Crear el gráfico parallel categories
colorscale = [[0, 'lightblue'], [0.5, 'orange'], [1, 'purple']]  # Colores para las muestras
fig = go.Figure(data=[go.Parcats(
    dimensions=[export_dim, compound_dim, import_dim],
    line={'color': interaction2['Color'], 'colorscale': colorscale},#, 'width': line_thickness * 10},  # Ajuste del grosor
    hoveron='color',
    hoverinfo='count+probability',
    labelfont={'size': 18, 'family': 'Times'},
    tickfont={'size': 16, 'family': 'Times'},
    arrangement='freeform'
)])

fig.update_layout(title_text="Flux Export")
fig.show()
#%%
# Normalize Flux
interaction2['Flux'] = np.abs(interaction2['Flux']).astype(int)
interaction2 = interaction2.loc[interaction2.index.repeat(interaction2['Flux'])]

# Dimensions
export_values = interaction2.loc[interaction2['Source'] == 'Export', 'MAGName']
import_values = interaction2.loc[interaction2['Source'] == 'Import', 'MAGName']
compound_values = interaction2['Name']

min_length = min(len(export_values), len(import_values), len(compound_values))

export_dim = go.parcats.Dimension(values=export_values[:min_length], label="Export")
import_dim = go.parcats.Dimension(values=import_values[:min_length], label="Import")
compound_dim = go.parcats.Dimension(values=compound_values[:min_length], label="Compound")

# Color normalization
interaction2['Color'] = interaction2['Color'] / interaction2['Color'].max()

# Graph
colorscale = [[0, 'lightblue'], [0.5, 'orange'], [1, 'purple']]
fig = go.Figure(data=[go.Parcats(
    dimensions=[export_dim, compound_dim, import_dim],
    line={'color': interaction2['Color'], 'colorscale': colorscale},
    hoveron='color',
    hoverinfo='count+probability',
    labelfont={'size': 18, 'family': 'Times'},
    tickfont={'size': 16, 'family': 'Times'},
    arrangement='perpendicular'
)])

fig.update_layout(title_text="Flux Export")
fig.show()        
#%%
# Ploting final graph
# Plot M. thermautotrophicus barplots
dt_plot = exchanges_sim_df.loc[(slice(None), 'maxbin2_spades_001'), :].copy()
idx = dt_plot.abs().max(axis=0) > 0.01
rxn_names = pd.read_csv('data/db/seed_reactions_corrected.tsv',sep='\t').set_index('id')['name'].to_dict()
names = pd.read_csv('data/db/seed_metabolites_edited.tsv',sep='\t').set_index('id')['name'].to_dict()

dt_plot.reset_index(inplace=True)
dt_plot.drop('compartment',axis=1,inplace=True)
dt_plot['Sample'] = dt_plot['level_0'].str.split("_").str[0]

exchanges_to_plot = {'EX_cpd01024_e':'Methane','EX_cpd00011_e':'CO2',
                'EX_cpd11640_e':'H2','EX_cpd00029_e':'Acetate',
                'EX_cpd00141_e':'Propionate','EX_cpd00137_e':'Citrate',
                'EX_cpd00130_e':'L-Malate','EX_cpd00106_e':'Fumarate',
                'EX_cpd00047_e':'Formate','EX_cpd00035_e':'L-Alanine',
                'EX_cpd00041_e':'L-Aspartate','EX_cpd00023_e':'L-Glutamate',
                'EX_cpd00053_e':'L-Glutamine','EX_cpd00033_e':'Glycine',
                'EX_cpd00060_e':'L-Methionine','EX_cpd00054_e':'L-Serine',
                'EX_cpd00161_e':'L-Threonine','EX_cpd00065_e':'L-Tryptophan'}

reactions_to_plot = ['rxn00278_c','rxn05938_c','rxn00175_c','rxn90072_c', 'rxn00250_c','rxn00260_c',
                       'rxn15962_c','rxn07189_c','rxn00260_c','rxn03127_c','rxn03020_c']
custom_palette = {
    'Control': '#56B4E9',
    'AA': '#D55E00',
    'AAAb': '#009E73',
}

import matplotlib.pyplot as plt
import seaborn as sns # type: ignore

sns.set_theme(style="whitegrid")

cpd = list(exchanges_to_plot.keys())[0]
g = sns.catplot(
    data=dt_plot, kind="bar",
    x="Sample",y=cpd, hue="Sample",
    errorbar="sd", palette=custom_palette, alpha=1, height=6
)
g.despine(left=True)
g.set_axis_labels("", "Flux")

g.fig.suptitle(exchanges_to_plot[cpd],fontsize=10, fontdict={"weight": "bold"})
# %%
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set_theme(style="whitegrid")
savepath = 'figures/'
# Convert exchanges_to_plot dictionary into a list of tuples (id, name)
exchanges_list = list(exchanges_to_plot.items())

# Reshape the data for plotting
dt_melted = dt_plot.melt(id_vars=["Sample"], value_vars=list(exchanges_to_plot.keys()), var_name="Exchange", value_name="Flux")

# Map Exchange IDs to their real names
dt_melted["Exchange Name"] = dt_melted["Exchange"].map(exchanges_to_plot)

# Plot using FacetGrid
g = sns.FacetGrid(dt_melted, col="Exchange Name", col_wrap=4, sharey=False, height=2, aspect=0.7)  
g.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

# Adjust aesthetics
g.set_titles(col_template="{col_name}")  
g.set_axis_labels("", "Flux")
g.add_legend()

# Show the plot
plt.savefig(f'{savepath}mt.cpds.pdf')
plt.show()
#%%
# Melt the dataframe for reactions
dt_melted_rxns = dt_plot.melt(id_vars=["Sample"], value_vars=reactions_to_plot, var_name="Reaction", value_name="Flux")

# # Plot using FacetGrid
g_rxn = sns.FacetGrid(dt_melted_rxns, col="Reaction", col_wrap=4, sharey=False, height=2, aspect=0.7)
g_rxn.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

# # Adjust aesthetics
g_rxn.set_titles(col_template="{col_name}")
g_rxn.set_axis_labels("", "Flux")
g_rxn.add_legend()

# # Show the plot
savepath = 'figures/'
plt.savefig(f'{savepath}mt.rxns.pdf')
plt.show()
#%%
##############################################################
##############################################################
##############################################################
############### Symbiobacterium sp. PD15 #####################
##############################################################
##############################################################
##############################################################
def run_pd15():
    dt_plot = exchanges_sim_df.loc[(slice(None), 'metabat2_spades_25'), :].copy()
    idx = dt_plot.abs().max(axis=0) > 0.000001

    dt_plot.reset_index(inplace=True)
    dt_plot.drop('compartment',axis=1,inplace=True)
    dt_plot['Sample'] = dt_plot['level_0'].str.split("_").str[0]

    exchanges_to_plot = {'EX_cpd00047_e':'Formate','EX_cpd00011_e':'CO2',
                    'EX_cpd11640_e':'H2','EX_cpd00029_e':'Acetate',
                    'EX_cpd00020_e':'Pyruvate','EX_cpd00161_e':'Threonine',
                    'EX_cpd00054_e':'Serine','EX_cpd00023_e':'Glutamate',
                    'EX_cpd00013_e':'NH3','EX_cpd00141_e':'Propionate',
                    'EX_cpd00048_e':'Sulfate','EX_cpd00239_e':'H2S',
                    'EX_cpd00036_e':'Succinate','EX_cpd00137_e':'Citrate',
                    'EX_cpd00033_e':'Glycine','EX_cpd00221_e':'D-Lactate',
                    'EX_cpd00084_e':'Cysteine'}

    reactions_to_plot = ['rxn08094_c','rxn00260_c','rxn00347_c','rxn00799_c', 'rxn00248_c','rxn00161_c',
                        'rxn00566_c','rxn05902_c','rxn05256_c','rxn00379_c']
    custom_palette = {
        'Control': '#56B4E9',
        'AA': '#D55E00',
        'AAAb': '#009E73',
    }

    import matplotlib.pyplot as plt
    import seaborn as sns # type: ignore

    sns.set_theme(style="whitegrid")

    cpd = list(exchanges_to_plot.keys())[0]
    g = sns.catplot(
        data=dt_plot, kind="bar",
        x="Sample",y=cpd, hue="Sample",
        errorbar="sd", palette=custom_palette, alpha=1, height=6
    )
    g.despine(left=True)
    g.set_axis_labels("", "Flux")

    g.fig.suptitle(exchanges_to_plot[cpd],fontsize=10, fontdict={"weight": "bold"})
    
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    sns.set_theme(style="whitegrid")
    savepath = 'figures/'
    # Convert exchanges_to_plot dictionary into a list of tuples (id, name)
    exchanges_list = list(exchanges_to_plot.items())

    # Reshape the data for plotting
    dt_melted = dt_plot.melt(id_vars=["Sample"], value_vars=list(exchanges_to_plot.keys()), var_name="Exchange", value_name="Flux")

    # Map Exchange IDs to their real names
    dt_melted["Exchange Name"] = dt_melted["Exchange"].map(exchanges_to_plot)

    # Plot using FacetGrid
    g = sns.FacetGrid(dt_melted, col="Exchange Name", col_wrap=4, sharey=False, height=2, aspect=0.7)  
    g.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

    # Adjust aesthetics
    g.set_titles(col_template="{col_name}")  
    g.set_axis_labels("", "Flux")
    g.add_legend()

    # Show the plot
    plt.savefig(f'{savepath}pd15.cpds.pdf')
    plt.show()
    
    # Melt the dataframe for reactions
    dt_melted_rxns = dt_plot.melt(id_vars=["Sample"], value_vars=reactions_to_plot, var_name="Reaction", value_name="Flux")

    # # Plot using FacetGrid
    g_rxn = sns.FacetGrid(dt_melted_rxns, col="Reaction", col_wrap=4, sharey=False, height=2, aspect=0.7)
    g_rxn.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

    # # Adjust aesthetics
    g_rxn.set_titles(col_template="{col_name}")
    g_rxn.set_axis_labels("", "Flux")
    g_rxn.add_legend()

    # # Show the plot
    plt.savefig(f'{savepath}pd15.rxns.pdf')
    plt.show()
run_pd15()
#%%
def run_pd14():
    dt_plot = exchanges_sim_df.loc[(slice(None), 'metabat2_spades_18'), :].copy()
    idx = dt_plot.abs().max(axis=0) > 0.000001

    dt_plot.reset_index(inplace=True)
    dt_plot.drop('compartment',axis=1,inplace=True)
    dt_plot['Sample'] = dt_plot['level_0'].str.split("_").str[0]

    exchanges_to_plot = {'EX_cpd00047_e':'Formate','EX_cpd00011_e':'CO2',
                    'EX_cpd11640_e':'H2','EX_cpd00029_e':'Acetate',
                    'EX_cpd00054_e':'Serine','EX_cpd00023_e':'Glutamate',
                    'EX_cpd00141_e':'Propionate','EX_cpd00035_e':'Alanine',
                    'EX_cpd00033_e':'Glycine','EX_cpd00041_e':'Aspartate',
                    'EX_cpd00084_e':'Cysteine'}

    reactions_to_plot = ['rxn00371_c','rxn00225_c','rxn00175_c']
    custom_palette = {
        'Control': '#56B4E9',
        'AA': '#D55E00',
        'AAAb': '#009E73',
    }

    import matplotlib.pyplot as plt
    import seaborn as sns # type: ignore

    sns.set_theme(style="whitegrid")

    cpd = list(exchanges_to_plot.keys())[0]
    g = sns.catplot(
        data=dt_plot, kind="bar",
        x="Sample",y=cpd, hue="Sample",
        errorbar="sd", palette=custom_palette, alpha=1, height=6
    )
    g.despine(left=True)
    g.set_axis_labels("", "Flux")

    g.fig.suptitle(exchanges_to_plot[cpd],fontsize=10, fontdict={"weight": "bold"})
    
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    sns.set_theme(style="whitegrid")
    savepath = 'figures/'
    # Convert exchanges_to_plot dictionary into a list of tuples (id, name)
    exchanges_list = list(exchanges_to_plot.items())

    # Reshape the data for plotting
    dt_melted = dt_plot.melt(id_vars=["Sample"], value_vars=list(exchanges_to_plot.keys()), var_name="Exchange", value_name="Flux")

    # Map Exchange IDs to their real names
    dt_melted["Exchange Name"] = dt_melted["Exchange"].map(exchanges_to_plot)

    # Plot using FacetGrid
    g = sns.FacetGrid(dt_melted, col="Exchange Name", col_wrap=4, sharey=False, height=2, aspect=0.7)  
    g.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

    # Adjust aesthetics
    g.set_titles(col_template="{col_name}")  
    g.set_axis_labels("", "Flux")
    g.add_legend()

    # Show the plot
    plt.savefig(f'{savepath}pd15.cpds.pdf')
    plt.show()
    
    # Melt the dataframe for reactions
    dt_melted_rxns = dt_plot.melt(id_vars=["Sample"], value_vars=reactions_to_plot, var_name="Reaction", value_name="Flux")

    # # Plot using FacetGrid
    g_rxn = sns.FacetGrid(dt_melted_rxns, col="Reaction", col_wrap=4, sharey=False, height=2, aspect=0.7)
    g_rxn.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

    # # Adjust aesthetics
    g_rxn.set_titles(col_template="{col_name}")
    g_rxn.set_axis_labels("", "Flux")
    g_rxn.add_legend()

    # # Show the plot

    plt.savefig(f'{savepath}pd15.rxns.pdf')
    plt.show()
run_pd14()

# %%
def run_pd13():
    dt_plot = exchanges_sim_df.loc[(slice(None), 'metabat2_spades_15'), :].copy()
    idx = dt_plot.abs().max(axis=0) > 0.000001

    dt_plot.reset_index(inplace=True)
    dt_plot.drop('compartment',axis=1,inplace=True)
    dt_plot['Sample'] = dt_plot['level_0'].str.split("_").str[0]

    exchanges_to_plot = {'EX_cpd00047_e':'Formate','EX_cpd00011_e':'CO2',
                    'EX_cpd11640_e':'H2','EX_cpd00029_e':'Acetate',
                    'EX_cpd00054_e':'Serine','EX_cpd00023_e':'Glutamate',
                    'EX_cpd00141_e':'Propionate','EX_cpd00035_e':'Alanine',
                    'EX_cpd00033_e':'Glycine','EX_cpd00041_e':'Aspartate',
                    'EX_cpd00084_e':'Cysteine','EX_cpd00156_e':'Valine'}

    reactions_to_plot = ['rxn00165_c']
    custom_palette = {
        'Control': '#56B4E9',
        'AA': '#D55E00',
        'AAAb': '#009E73',
    }

    import matplotlib.pyplot as plt
    import seaborn as sns # type: ignore

    sns.set_theme(style="whitegrid")

    cpd = list(exchanges_to_plot.keys())[0]
    g = sns.catplot(
        data=dt_plot, kind="bar",
        x="Sample",y=cpd, hue="Sample",
        errorbar="sd", palette=custom_palette, alpha=1, height=6
    )
    g.despine(left=True)
    g.set_axis_labels("", "Flux")

    g.fig.suptitle(exchanges_to_plot[cpd],fontsize=10, fontdict={"weight": "bold"})
    
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    sns.set_theme(style="whitegrid")
    savepath = 'figures/'
    # Convert exchanges_to_plot dictionary into a list of tuples (id, name)
    exchanges_list = list(exchanges_to_plot.items())

    # Reshape the data for plotting
    dt_melted = dt_plot.melt(id_vars=["Sample"], value_vars=list(exchanges_to_plot.keys()), var_name="Exchange", value_name="Flux")

    # Map Exchange IDs to their real names
    dt_melted["Exchange Name"] = dt_melted["Exchange"].map(exchanges_to_plot)

    # Plot using FacetGrid
    g = sns.FacetGrid(dt_melted, col="Exchange Name", col_wrap=4, sharey=False, height=2, aspect=0.7)  
    g.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

    # Adjust aesthetics
    g.set_titles(col_template="{col_name}")  
    g.set_axis_labels("", "Flux")
    g.add_legend()

    # Show the plot
    
    plt.savefig(f'{savepath}pd13.cpds.pdf')
    plt.show()
    
    # Melt the dataframe for reactions
    dt_melted_rxns = dt_plot.melt(id_vars=["Sample"], value_vars=reactions_to_plot, var_name="Reaction", value_name="Flux")

    # Plot using FacetGrid
    g_rxn = sns.FacetGrid(dt_melted_rxns, col="Reaction", col_wrap=4, sharey=False, height=2, aspect=0.7)
    g_rxn.map_dataframe(sns.barplot, x="Sample", y="Flux", hue="Sample", palette=custom_palette, errorbar="sd")

    # # Adjust aesthetics
    g_rxn.set_titles(col_template="{col_name}")
    g_rxn.set_axis_labels("", "Flux")
    g_rxn.add_legend()

    # # Show the plot
    
    plt.savefig(f'{savepath}pd13.rxns.pdf')
    plt.show()
run_pd13()