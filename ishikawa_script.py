# %%
#!/usr/bin/env python3

## "# %%" allows running cells independently of each other

# Import modules
import os
import sys
import csv
import functools as fn
import matplotlib
from matplotlib.pyplot import yticks
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib_venn import venn2
from numpy.core.fromnumeric import transpose
import numpy as np
import pandas as pd
import seaborn as sns
import random
import cptac
import cptac.utils as ut
import statistics
import gtfparse as gtfp
import scipy.stats as stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler, StandardScaler
import random
import math
import subprocess
# %%
## Custom functions

def getFile(fileIN):
    data = pd.read_csv(fileIN)
    print(data)
    return data

def fileListing(dir):
    files=[]
    filenames=[]
    for filename in sorted(os.listdir(dir)):
        if filename.endswith('counts'):
            file=os.path.join(dir+'/'+filename)
            if file not in files:
                files.append(file)
                filenames.append(filename)
    return files,filenames


def gtfParser(file):
    '''Get protein coding genes from gtf file'''
    df = gtfp.read_gtf(file)
    df_protCod = df[(df["feature"] == "gene") & (df["gene_type"] == "protein_coding")]
    prot_genes=df_protCod[["gene_id","gene_name"]]
    prot_genes[['ens_id','version']] = prot_genes.gene_id.str.split(".",expand=True)
    prot_genes.drop(columns =["gene_id","version"], inplace = True)
    return prot_genes


def parseClin(clin_data,stages_collapsed='yes'):
    df = pd.read_csv(clin_data, delimiter = "\t")
    df[['name','ext']]=df.file.str.split(".",n=1,expand=True)
    df.drop(columns=['file','ext'],inplace=True)
    df['name'] = 's_' + df['name'].astype(str) # add prefix for compatibility with R
    df['name'] = df['name'].str.replace('-','.')
    df.sort_values("name", inplace=True)
    if stages_collapsed == 'yes':
        df['stage']=df['stage'].replace({"Stage I":"Stage_I",
                    "Stage IA":"Stage_I","Stage IB":"Stage_I", 
                    "Stage IC":"Stage_I","Stage II":"Stage_II",
                    "Stage IIA":"Stage_II","Stage IIB":"Stage_II",
                    "Stage III":"Stage_III","Stage IIIA":"Stage_III",
                    "Stage IIIB":"Stage_III","Stage IIIC":"Stage_III",
                    "Stage IIIC1":"Stage_III","Stage IIIC2":"Stage_III",
                    "Stage IV":"Stage_IV","Stage IVA":"Stage_IV",
                    "Stage IVB":"Stage_IV"})
        
    return df


def seekName(row,dict_names):
    for name,info in dict_names.items():
        if row['name'].startswith(name):
            metadata=[]
            for key in info:
                metadata.append(info[key])
            return metadata


def prepVSD(vsd_csv):
    df=pd.read_csv(vsd_csv)
    df.set_index('gene_name',inplace=True)
    df_t = df.transpose()
    df_t.reset_index(inplace=True)
    df_t.rename(columns={'index':'name'},inplace=True)
    return df_t


def plot_PCA(data_csv_dir,labels_csv_dir):
    '''Plot PCA from normalize count data. Need separate files
    for counts and labels'''
    datafile=pd.read_csv(data_csv_dir)
    datafile.set_index('name',inplace=True)
    x=datafile.values

    labels_file=pd.read_csv(labels_csv_dir)
    y=labels_file.loc[:,['stage']].values

    x=StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                , columns = ['principal component 1', 'principal component 2'])

    finalDF=pd.concat([principalDf,labels_file[['stage']]],axis=1)
    pca.explained_variance_ratio_

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)
    targets = ['Stage_I', 'Stage_II', 'Stage_III','Stage_IV']
    colors = ['r', 'g', 'b','y']
    for target, color in zip(targets,colors):
        indicesToKeep = finalDF['stage'] == target
        ax.scatter(finalDF.loc[indicesToKeep, 'principal component 1']
                , finalDF.loc[indicesToKeep, 'principal component 2']
                , c = color
                , s = 50)
    ax.legend(targets)
    ax.grid()


def getKmeans(data_csv,labels_csv,upper_limit=19419):
    '''Use normalized data and the label column in your
        label_csv file should be the first (0). Upper limit
        is defined by number of variables (genes) included in 
        your count matrix, default is 19419 after filtering 
        by pritein coding genes.
        '''
    data = np.genfromtxt(
        data_csv,
        delimiter=",",
        usecols=range(1, upper_limit),
        skip_header=1)


    true_label_names = np.genfromtxt(
        labels_csv,
        delimiter=",",
        usecols=(0,),
        skip_header=1,
        dtype="str"
    )

    label_encoder = LabelEncoder()

    true_labels = label_encoder.fit_transform(true_label_names)

    label_encoder.classes_
    n_clusters = len(label_encoder.classes_)


    preprocessor = Pipeline(
        [
            ("scaler", MinMaxScaler()),
            ("pca", PCA(n_components=2, random_state=None)),
        ]
    )

    clusterer = Pipeline(
    [
        (
            "kmeans",
            KMeans(
                n_clusters=n_clusters,
                init="k-means++",
                n_init=50,
                max_iter=500,
                random_state=None,
            ),
        ),
    ]
    )

    pipe = Pipeline(
        [
            ("preprocessor", preprocessor),
            ("clusterer", clusterer)
        ]
    )

    pipe.fit(data)

    preprocessed_data = pipe["preprocessor"].transform(data)

    predicted_labels = pipe["clusterer"]["kmeans"].labels_

    silhouette_score(preprocessed_data, predicted_labels)

    adjusted_rand_score(true_labels, predicted_labels)

    pcadf = pd.DataFrame(
        pipe["preprocessor"].transform(data),
        columns=["component_1", "component_2"],
    )

    pcadf["predicted_cluster"] = pipe["clusterer"]["kmeans"].labels_
    pcadf["true_label"] = label_encoder.inverse_transform(true_labels)

    plt.style.use("fivethirtyeight")
    plt.figure(figsize=(8, 8))

    scat = sns.scatterplot(
        "component_1",
        "component_2",
        s=50,
        data=pcadf,
        hue="predicted_cluster",
        style="true_label",
        palette="Set2",
    )

    scat.set_title(
        "Clustering results from TCGA-UCEC\nGene Expression Data"
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

    plt.show()



def interGene(list1,list2):
    inter=[]
    for element in list1:
        if element in list2:
            inter.append(element)
    return inter


def my_percentile(data, percentile):
    n = len(data)
    p = n * percentile / 100
    if p.is_integer():
        return sorted(data)[int(p)]
    else:
        return sorted(data)[int(math.ceil(p)) - 1]

    
def bootsModel(list1,list2,size,iter,ref):
    '''Takes two lists and the number of different random variables     
        to prepare distribution (without replacement).
        list1 contains protein coding genes in TADs with PgCRs, while 
        list2 is used to prepare random samples of size=size. Iter is the
        number of random samples to produce for distribution and ref is
        the shared number of genes between the original datasets.'''
    n = 0
    shared_nlist=[]
    while n <= iter:
        r_list=random.sample(list2,size)
        shared=interGene(r_list,list1)
        l_shared=len(shared)
        shared_nlist.append(l_shared)
        n+=1
    
    # calculate confidence interval (2.5 and 97.5th)
    index=[2.5, 97.5]
    perc_func=[my_percentile(shared_nlist,i) for i in index]
    perc2_5=perc_func[0]
    print('percentile 2.5 --> '+str(perc2_5))
    perc97_5=perc_func[1]
    print('percentile 97.5 --> '+str(perc97_5))
    print('Mean: '+ str(statistics.mean(shared_nlist)))

    # plot distribution of shared genes after bootstraping ## intersection with actual group is 75
    plt.figure(figsize=(7,7))
    p=sns.displot(shared_nlist, kind="kde")
    (p.map(plt.axvline, x=ref,ls='--',linewidth=2, alpha=.7,color='black')
        .map(plt.axvline, x=perc2_5,ls='-',linewidth=1, alpha=.9,color='red')
        .map(plt.axvline, x=perc97_5,ls='-',linewidth=1, alpha=.9,color='red')
        .tight_layout(w_pad=0))
    p.fig.text(0.55, 0.9,'95% CI', fontsize=12)
    plt.savefig('bootstrapping_distribution_withCI.pdf')

    #sns.boxplot(shared_nlist)
    return shared_nlist

def fixedWidthClusterMap(dataFrame, cellSizePixels=50):
    # Calulate the figure size, this gets us close, but not quite to the right place
    dpi = matplotlib.rcParams['figure.dpi']
    marginWidth = matplotlib.rcParams['figure.subplot.right']-matplotlib.rcParams['figure.subplot.left']
    marginHeight = matplotlib.rcParams['figure.subplot.top']-matplotlib.rcParams['figure.subplot.bottom']
    Ny,Nx = dataFrame.shape
    figWidth = (Nx*cellSizePixels/dpi)/0.8/marginWidth
    figHeigh = (Ny*cellSizePixels/dpi)/0.8/marginHeight

    # do the actual plot
    grid = sns.clustermap(dataFrame, figsize=(figWidth, figHeigh),z_score=0,cmap='vlag',
        row_cluster=True,col_cluster=False, method='median',metric='euclidean')

    # calculate the size of the heatmap axes
    axWidth = (Nx*cellSizePixels)/(figWidth*dpi)
    axHeight = (Ny*cellSizePixels)/(figHeigh*dpi)

    # resize heatmap
    ax_heatmap_orig_pos = grid.ax_heatmap.get_position()
    grid.ax_heatmap.set_position([ax_heatmap_orig_pos.x0, ax_heatmap_orig_pos.y0, 
                                  axWidth, axHeight])

    # resize dendrograms to match
    ax_row_orig_pos = grid.ax_row_dendrogram.get_position()
    grid.ax_row_dendrogram.set_position([ax_row_orig_pos.x0, ax_row_orig_pos.y0, 
                                         ax_row_orig_pos.width, axHeight])
    ax_col_orig_pos = grid.ax_col_dendrogram.get_position()
    grid.ax_col_dendrogram.set_position([ax_col_orig_pos.x0, ax_heatmap_orig_pos.y0+axHeight,
                                         axWidth, ax_col_orig_pos.height])
    return grid # return ClusterGrid object


def get_expression_perGO(expr_df,clusters_dirPath,go_term,cell_size):
    clust_dic={}
    for clus in clusters_dirPath:
        clus_name=clus.split("/")[-1]
        df_clust=pd.read_csv(clus,sep="\t",
            usecols=[1,5,7,8])
        clust_dic.update({clus_name:df_clust})

    id_list=[]
    for key,value in clust_dic.items():
        value=value.loc[value['Description'] == go_term]
        value_melt = value.assign(geneID=value.geneID.str.split("/"))
        geneids=value_melt['geneID'].tolist()
        id_list.append(geneids)
        
    genes=[]
    for id in id_list:
        for gene in id:
            for gn in gene:
                genes.append(gn)

    expr_df_copy=expr_df.copy()
    expr_df_copy.reset_index(inplace=True)
    go_values=expr_df_copy[expr_df_copy['gene_name'].isin(genes)]
    go_values.set_index('gene_name',inplace=True)
    
    sns.set(font_scale=0.4)
    plt.figure()
    cm=fixedWidthClusterMap(go_values,cell_size)
    plt.savefig('/Users/alejandrolagreca/Google Drive/Figuras final Dai/daiana_proteomics/python_analysis/heatmap_'+go_term+'.pdf')

    return go_values


# %% [markdown]
# # Evaluate PAXs expression in Ishikawa data
# %%
ishi=pd.read_csv('ishi_expression_dataset.txt', sep='\t') # expression dataset named "Figure 1 - Source data 1"
ishi.rename(columns={'NAME':'ens'},inplace=True)
ishi.drop(columns=['DESCRIPTION'],inplace=True)
#ishi[['ens_id','genename']] = ishi.ens.str.split("_",expand=True)
ishi_melt=ishi.melt(id_vars=['ens'])
ishi_melt=ishi_melt[ishi_melt['ens'].str.contains(r'PAX\d{1}')]
#ishi_melt.set_index('ens',inplace=True)
ishi_melt['value'] = np.log2(ishi_melt['value']+1)

plt.figure()
g=sns.catplot(data=ishi_melt,x='variable',y='value',col='ens',kind='bar',sharey=False)
g.savefig('pax_expression_ishi_data.pdf')
plt.show()


# %% [markdown]
# # TCGA-UCEC
# %%
###Preparing data

cts='tcga_data/' # directory with raw count files from tcga-ucec project
gtf='gencode.v38.annotation.gtf' # gene annotation from gencode
clin='stage_an_type_EC.txt' # tab delimited file with FIGO stage\tHistologic type\tfilename

## Parsing and merging count files
fileList,fileNames=fileListing(cts)

dataframes = [ pd.read_csv( f, sep='\t',header=None,names=("A","B") ) for f in fileList ] # add arguments as necessary to the read_csv method
merged = fn.reduce(lambda left,right: pd.merge(left,right,on='A', how='outer'), dataframes)
merged[['ens_id','version']] = merged.A.str.split(".",expand=True)
merged.drop(columns =["version","A"], inplace = True)
merged.set_index("ens_id", inplace=True)
merged.set_axis(fileNames, inplace=True, axis=1)
merged.columns = merged.columns.str.replace(r'.htseq.counts$', '')

##Parsing and filtering gtf file
gtf_prot=gtfParser(gtf)

##Merge count matrix and gtf file to place gene_names (keep protein coding)
dfs=[gtf_prot,merged]
merged_gtf=fn.reduce( lambda left,right: pd.merge(left,right,on="ens_id",how="outer"),dfs)
merged_gtf.dropna(subset=['gene_name'],inplace=True)
merged_gtf.drop(columns=['ens_id'],inplace=True)
merged_gtf.fillna(0,inplace=True)
merged_gtf.set_index("gene_name", inplace=True)
merged_gtf=merged_gtf.round(0).astype(int)

merged_gtf_f = merged_gtf[(merged_gtf.T != 0).any()]
merged_gtf_f=merged_gtf_f.add_prefix('s_') # add prefix for compatibility with R
merged_gtf_f.columns = merged_gtf_f.columns.str.replace('-', '.')
merged_gtf_f.reset_index(inplace=True)
merged_gtf_f.sort_values("gene_name", inplace=True)
merged_gtf_f.drop_duplicates(subset="gene_name", keep='first', inplace=True)

#Save merged count matrix as csv file (go to deseq2 in R)
merged_gtf_f.to_csv('ucec_count_matrix.csv',index=False)

# %%
#Parse and clinical data
clinical=parseClin(clin)

#Save clinical data as csv file (go to deseq2 in R)
clinical.to_csv('ucec_clinical_data.csv',index=False)

# %%
### Using R ###
path_to_counts='ucec_count_matrix.csv'
path_to_metadata='ucec_clinical_data.csv'
job_title='TCGA'
feature='stage'
condition1='Stage_IV'
condition2='Stage_I'

subprocess.call(['deseq_analysis/run_deseq2.R',path_to_counts,path_to_metadata,job_title,
                feature,condition1,condition2])


# %%
###Load normalized count (deseq2) and clinical data
norm_counts=pd.read_csv('count_matrix_vstransformation.csv')
clinical=pd.read_csv('ucec_clinical_data.csv')

#Melt count matrix by gene name
norm_counts_melt=norm_counts.melt(id_vars=['gene_name'])
norm_counts_melt.rename(columns={'variable':'name'},inplace=True)

merged_clin=pd.merge(norm_counts_melt, clinical, on='name')

#Select PR and ER genes and transpose
merged_clin.set_index("gene_name", inplace=True)
genes=merged_clin.loc[['PGR','ESR1']]
genes.reset_index(inplace=True)
genes_org=genes.pivot(index='name',columns='gene_name')
genes_org.columns = ['_'.join(col) for col in genes_org.columns]
genes_org.drop(columns=['type_ESR1','stage_ESR1'],inplace=True)
genes_org.rename(columns={'type_PGR':'type','stage_PGR':'stage'},inplace=True)


# %%
#Plot scatter all stages
pr_cut=10
er_cut=10
plt.figure(figsize=(7,7))
p=sns.relplot(data=genes_org, x="value_PGR", y="value_ESR1", 
                hue="stage",style='stage',kind="scatter")

(p.map(plt.axhline, y=er_cut,ls='--',linewidth=2, alpha=.7,color='black')
    .map(plt.axvline, x=pr_cut,ls='--',linewidth=2, alpha=.7,color='black')
    .set_axis_labels("value_PGR", "value_ESR1")
    .tight_layout(w_pad=0))
p.fig.text(0.25, 0.95,'High ER - Low PR', fontsize=12)
p.fig.text(0.52, 0.95,'High ER - High PR', fontsize=12)
p.fig.text(0.52, 0.15,'Low ER - High PR', fontsize=12)
p.fig.text(0.25, 0.15,'Low ER - Low PR', fontsize=12)
plt.show()

#Plot scatter stage I
genes_org.reset_index(inplace=True)
stageI=genes_org.set_index("stage")
stageI=stageI.loc[['Stage_I']]
stageI.reset_index(inplace=True)

plt.figure(figsize=(7,7))
sI=sns.relplot(data=stageI, x="value_PGR", y="value_ESR1", 
                hue="stage",style='stage',kind="scatter")

(sI.map(plt.axhline, y=er_cut,ls='--',linewidth=2, alpha=.7,color='black')
    .map(plt.axvline, x=pr_cut,ls='--',linewidth=2, alpha=.7,color='black')
    .set_axis_labels("value_PGR", "value_ESR1")
    .tight_layout(w_pad=0))
p.fig.text(0.25, 0.95,'High ER - Low PR', fontsize=12)
p.fig.text(0.52, 0.95,'High ER - High PR', fontsize=12)
p.fig.text(0.52, 0.15,'Low ER - High PR', fontsize=12)
p.fig.text(0.25, 0.15,'Low ER - Low PR', fontsize=12)
plt.show()

# %%
#####--##### (use above data to filter outlier stage I samples)
#Filter count matrix by High ER - High PR
rejects=stageI.loc[(stageI['value_ESR1'] < er_cut) | (stageI['value_PGR'] < pr_cut)]
rejects_list = rejects['name'].tolist()

merged_gtf_f.set_index('gene_name', inplace=True)
theones = merged_gtf_f[merged_gtf_f.columns.difference(rejects_list)]

#Save stage I-filtered count matrix as csv file
theones.to_csv('ucec_stageI-filt_count_matrix.csv',index=True)

#Prepare and save clinical data for deseq2
clinical_filt=clinical[~clinical['name'].isin(rejects_list)]
clinical_filt.to_csv('ucec_clinical_data_filtered.csv',index=False)

# %% [markdown]
### Using R ###
path_to_counts='ucec_stageI-filt_count_matrix.csv'
path_to_metadata='ucec_clinical_data_filtered.csv'
job_title='TCGA-filt'
feature='stage'
condition1='Stage_IV'
condition2='Stage_I' # Stage_II and Stage_III

subprocess.call(['deseq_analysis/run_deseq2.R',path_to_counts,path_to_metadata,job_title,
                feature,condition1,condition2])


################

# %%
vsd=prepVSD('count_matrix_vstransformation.csv')
clin=pd.read_csv('ucec_clinical_data_filtered.csv')

# %%
merged=pd.merge(vsd,clin,on='name',how='inner')
merged

# %%
graphingSite = 'PAX2'

sns.catplot(x='stage', y=graphingSite, kind='box', data=merged, showfliers=False,
        order=['Stage_I','Stage_II','Stage_III','Stage_IV'])
sns.stripplot(x='stage', y=graphingSite, data=merged, color='.3',
        order=['Stage_I','Stage_II','Stage_III','Stage_IV'])

plt.show()

# %%
graphingSite = 'PAX2'

sns.catplot(x='type', y=graphingSite, kind='box', data=merged, showfliers=False)
sns.stripplot(x='type', y=graphingSite, data=merged, color='.3')
plt.xticks(rotation=45)

plt.show()

# %%
# Bootstrapping using genes in TADs with PgCRs as fixed and randomizing deseq results (s4 vs all the rest, pooled)

df_s4_1=pd.read_csv('deseq2_results_S4vSI.csv',usecols=['row'])
df_s4_2=pd.read_csv('deseq2_results_S4vSII.csv',usecols=['row'])
df_s4_3=pd.read_csv('deseq2_results_S4vSIII.csv',usecols=['row'])
merged=[df_s4_1,df_s4_2,df_s4_3]
merged_dfs=pd.concat(merged)

list_merged=merged_dfs['row'].tolist()
list_merged_uniq=list(set(list_merged))
l_list_merged_uniq=len(list_merged_uniq) 
print('lenght of list of merged DE stages --> '+str(l_list_merged_uniq)) 

with open('tad-pgcr_genes', 'r') as f: # provided as source data for figure 6
        content = f.read()
        tadlist=content.split('\n')

l_tadlist=len(tadlist)


common=interGene(tadlist,list_merged_uniq)
l_common=len(common)
print('lenght of common genes --> '+str(l_common))
#with open('shared_genes.csv', 'w') as f:
 #   write=csv.writer(f)
  #  write.writerow(common)

list_prot=gtf_prot['gene_name'].tolist()

iterations=10000

len_list_boot=bootsModel(tadlist,list_prot,l_list_merged_uniq,iterations,l_common)

# %%
statistics.mean(len_list_boot)

plt.figure(figsize=(5,5))
venn2(subsets = (l_list_merged_uniq, l_tadlist, l_common), set_labels = ('Endometrial Cancer DE genes', 'Genes in TADs with PgCRs'),
        set_colors=('r', 'g'), alpha = 0.7)
plt.show()
###################################
# %%
shared = pd.read_csv('shared_genes') # provided as source data for figure 6
shared_list=shared['gene'].tolist()

# %% [markdown]
#Testing Dou endometrial data

# %%
cptac.download(dataset='endometrial')

# %%
endo=cptac.Endometrial()
endo.list_data()

# %%
transc=endo.get_transcriptomics()
samples=transc.index
genes=transc.columns

# %%
transc

# %%
#Set desired attribute to variable 'clinical_attribute'
clinical_attribute = "Histologic_type"

#Join attribute with acetylation dataframe
clinical_and_transc = endo.join_metadata_to_omics(metadata_df_name="clinical", omics_df_name="transcriptomics",
                                                     metadata_cols=clinical_attribute)

# Use the cptac.utils.reduce_multiindex function to combine the
# multiple column levels, so it's easier to graph our data
clinical_and_transc = ut.reduce_multiindex(df=clinical_and_transc, flatten=True)

clinical_and_transc.head()

# %%
clinical_and_transc

# %%
#Make dataframes with only endometrioid and only serous data in order to compare 
endom = clinical_and_transc.loc[clinical_and_transc[clinical_attribute] == "Endometrioid"]
serous = clinical_and_transc.loc[clinical_and_transc[clinical_attribute] == "Serous"]
clinical_and_transc[[clinical_attribute]] = clinical_and_transc[[clinical_attribute]].fillna(
    value="Non_Tumor")

nont = clinical_and_transc.loc[clinical_and_transc[clinical_attribute] == "Non_Tumor"]

# %%
serous

# %%
graphingSite = 'PAX2_transcriptomics'
print(stats.ttest_ind(endom[graphingSite], serous[graphingSite]))

plt.figure()
sns.catplot(x=clinical_attribute, y=graphingSite, kind='box', data=clinical_and_transc, showfliers=False, 
            order=["Non_Tumor", "Endometrioid", "Serous"])
sns.stripplot(x=clinical_attribute, y=graphingSite, data=clinical_and_transc, color='.3', 
            order=["Non_Tumor", "Endometrioid", "Serous"])
plt.show()

# %%
clinical_and_transc.reset_index(inplace=True)
transc_melt = clinical_and_transc.melt(id_vars=['Patient_ID',"Histologic_type"])
pax_genes=transc_melt[transc_melt['Name'].str.contains(r'^PAX\d{1}_')]


# %%
pax_genes

# %%
plt.figure()
sns.catplot(x='Histologic_type', y='value', kind='box', data=pax_genes, #showfliers=False, 
            order=["Non_Tumor", "Endometrioid", "Serous"],col='Name',col_wrap=3)
#sns.stripplot(x='Histologic_type', y='value', data=pax_genes, color='.3', 
 #          order=["Non_Tumor", "Endometrioid", "Serous"])
plt.show()

# %%
box = sns.catplot(x='Histologic_type',y='value', kind='box', data=pax_genes,linewidth=0.2,
fliersize=0.3)
# #fig = box.get_figure()
box.set(xticklabels=[])
#box.savefig('boxplot_BBK.png')
plt.show()

# %%
# pax2 expression in histologic types (+ ANOVA)
pax2=pax_genes[pax_genes['Name'].str.contains('PAX2')]
pax2=pax2[['Histologic_type','value']]

# %%
pax2

# %%

# Test Normality
unique_majors = pax2['Histologic_type'].unique()
for major in unique_majors:
    stats.probplot(pax2[pax2['Histologic_type'] == major]['value'], dist="norm", plot=plt)
    plt.title("Probability Plot - " +  major)
    plt.show()

# %%
# Test homogeneity (value must be under 2)
ratio = pax2.groupby('Histologic_type').std().max() / pax2.groupby('Histologic_type').std().min()
ratio

# %%
pax2[['value']] = pax2[['value']].fillna(
    value=0)

# %%
stats.f_oneway(pax2['value'][pax2['Histologic_type'] == 'Endometrioid'],
                pax2['value'][pax2['Histologic_type'] == 'Serous'],
                pax2['value'][pax2['Histologic_type'] == 'Non_Tumor'])

# %%
import statsmodels.stats.multicomp as mc

comp = mc.MultiComparison(pax2['value'], pax2['Histologic_type'])
post_hoc_res = comp.tukeyhsd()
post_hoc_res.summary()

# %%
post_hoc_res.plot_simultaneous(ylabel= "Histologic type", xlabel= "Score Difference")