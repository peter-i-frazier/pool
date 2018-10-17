import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#
# Functions for producing distances between recommended point and the closest training data,
#

def closest_data(x,y,data):
    # pass either data_sfp or data_acps for data
    # looks for the closest point in data using the euclidean distance
    dist_sqr = (data['x']-x)**2 + (data['y']-y)**2
    min_dist = np.sqrt(dist_sqr.min())
    idx = dist_sqr.idxmin()
    return min_dist

def add_distances(method_df,data_df):
    method_df['closest_data_dist'] = pd.Series()
    for i,row in method_df.iterrows():
        method_df['closest_data_dist'][i]=closest_data(row['x'],row['y'],data_df)

#
# Functions for plot with line between recommended point and the closest training data,
#

def closest_idx(x,y,data):
    # pass either data_sfp or data_acps for data
    # looks for the closest point in data using the euclidean distance
    dist_sqr = (data['x']-x)**2 + (data['y']-y)**2
    idx = dist_sqr.idxmin()
    return idx

def plot_lines(method_df,data_df,col):
    for i,row in method_df.iterrows():
        idx=closest_idx(row['x'],row['y'],data_df)
        x = data_df.loc[idx]['x']
        y = data_df.loc[idx]['y']
        plt.plot([row['x'],x],[row['y'],y],col,alpha=1,linewidth=1)

#
# Get data and add distances
#

sfp_coord_data = pd.read_csv('figures/sfp_specific_simulation_coord.csv')
pool_sfp = sfp_coord_data.iloc[:100,].copy()
mutation_sfp = sfp_coord_data.iloc[100:200,].copy()
predict_sfp = sfp_coord_data.iloc[200:300,].copy()
data_sfp = sfp_coord_data.iloc[300:,].copy()

acps_coord_data = pd.read_csv('figures/acps_specific_simulation_coord.csv')
pool_acps = acps_coord_data.iloc[:100,].copy()
mutation_acps = acps_coord_data.iloc[100:200,].copy()
predict_acps = acps_coord_data.iloc[200:300,].copy()
data_acps = acps_coord_data.iloc[300:,].copy()

add_distances(pool_sfp,data_sfp)
add_distances(mutation_sfp,data_sfp)
add_distances(pool_acps,data_acps)
add_distances(mutation_acps,data_acps)

# The data is a number between 0 and 1, where 1 indicates a 100% success rate
benchmark_sfp = pd.read_csv('figures/sfp_specific_simulation_quality.csv')
benchmark_acps = pd.read_csv('figures/acps_specific_simulation_quality.csv')



#
# Produce distance statistics
#
stat_file=open('figures/distance_statistics.txt','w')

print('Statistics for distance to closest datapoint in training data', file=stat_file)
print('POOL Sfp', file=stat_file)
pool_sfp['closest_data_dist'].describe().to_csv(stat_file)
print('\n',file=stat_file)

print('Mutation Sfp', file=stat_file)
mutation_sfp['closest_data_dist'].describe().to_csv(stat_file)
print('\n',file=stat_file)

print('POOL AcpS',file=stat_file)
pool_acps['closest_data_dist'].describe().to_csv(stat_file)
print('\n',file=stat_file)

print('Mutation AcpS',file=stat_file)
mutation_acps['closest_data_dist'].describe().to_csv(stat_file)
print('\n',file=stat_file)

stat_file.close()


#
# Formatting details for plots
# colors from http://colorbrewer2.org/#type=qualitative&scheme=Accent&n=3
# 

pool_col = '#ca0020'
predict_col = '#f4a582'
mutation_col='#0571b0'

pool_mrk='o'
predict_mrk='d'
mutation_mrk='^'
data_mrk='s'
mrk_sz=4



# Fig 2C and 2D, ROC curves
sfp_roc_data = pd.read_csv('figures/sfp_specific_roc_data.csv')
plt.figure()
plt.plot(sfp_roc_data['x'], sfp_roc_data['y'])
plt.plot([0,1], [0,1], 'k--')
plt.legend(['Sfp','50:50 probability'])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Sfp Specific')
plt.savefig('figures/fig2c.eps')
plt.close()


acps_roc_data = pd.read_csv('figures/acps_specific_roc_data.csv')
plt.figure()
plt.plot(acps_roc_data['x'], acps_roc_data['y'])
plt.plot([0,1], [0,1], 'k--')
plt.legend(['AcpS','50:50 probability'])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('AcpS Specific')
plt.savefig('figures/fig2d.eps')
plt.close()

#
# Fig 3A and 3B, benchmark comparisons
#

subset_to_plot = list(range(0,100,5))
plt.plot(benchmark_sfp['pool'][subset_to_plot]*100,'o',color=pool_col,markeredgecolor='k')
plt.plot(benchmark_sfp['mutate'][subset_to_plot]*100,'^',color=mutation_col,markeredgecolor='k')
plt.plot(benchmark_sfp['predict_optimize'][subset_to_plot]*100,'d',color=predict_col,markeredgecolor='k')

plt.plot(benchmark_sfp['pool']*100,'-',color=pool_col)
plt.plot(benchmark_sfp['mutate']*100,'-',color=mutation_col)
plt.plot(benchmark_sfp['predict_optimize']*100,'-',color=predict_col)

plt.plot(benchmark_sfp['pool'][subset_to_plot]*100,'o',color=pool_col,markeredgecolor='k')
plt.plot(benchmark_sfp['mutate'][subset_to_plot]*100,'^',color=mutation_col,markeredgecolor='k')
plt.plot(benchmark_sfp['predict_optimize'][subset_to_plot]*100,'d',color=predict_col,markeredgecolor='k')

plt.xlabel('# of peptides recommended')
plt.ylabel('Probability of selectivity (%)')
#plt.text(0,6,'Sfp',fontsize=14)
plt.title('Sfp')
plt.ylim(0,60)
plt.legend(['POOL','Mutation','Predict-then-Optimize'],frameon=False,loc='upper center',ncol=3)
plt.savefig('figures/fig3a.eps')
plt.close()

subset_to_plot = list(range(0,100,5))
plt.plot(benchmark_acps['pool'][subset_to_plot]*100,'o',color=pool_col,markeredgecolor='k')
plt.plot(benchmark_acps['mutate'][subset_to_plot]*100,'^',color=mutation_col,markeredgecolor='k')
plt.plot(benchmark_acps['predict_optimize'][subset_to_plot]*100,'d',color=predict_col,markeredgecolor='k')

plt.plot(benchmark_acps['pool']*100,'-',color=pool_col)
plt.plot(benchmark_acps['mutate']*100,'-',color=mutation_col)
plt.plot(benchmark_acps['predict_optimize']*100,'-',color=predict_col)

plt.plot(benchmark_acps['pool'][subset_to_plot]*100,'o',color=pool_col,markeredgecolor='k')
plt.plot(benchmark_acps['mutate'][subset_to_plot]*100,'^',color=mutation_col,markeredgecolor='k')
plt.plot(benchmark_acps['predict_optimize'][subset_to_plot]*100,'d',color=predict_col,markeredgecolor='k')

plt.xlabel('# of peptides recommended')
plt.ylabel('Probability of selectivity (%)')
#plt.text(0,15,'AcpS',fontsize=14)
plt.title('AcpS')
plt.ylim(0,2.5)
plt.legend(['POOL','Mutation','Predict-then-Optimize'],frameon=False,loc='upper center',ncol=3)
plt.savefig('figures/fig3b.eps')
plt.close()



#
# FIG 3C (2-dimensional embedding of training data and recommended points for Sfp)
#

plt.plot(pool_sfp['x'],pool_sfp['y'],pool_mrk,color=pool_col,markersize=mrk_sz)
plt.plot(mutation_sfp['x'],mutation_sfp['y'],mutation_mrk,color=mutation_col,markersize=mrk_sz)
plt.plot(data_sfp['x'],data_sfp['y'],data_mrk,color='0.75',markersize=mrk_sz)
plt.plot(predict_sfp['x'],predict_sfp['y'],predict_mrk,color=predict_col,markersize=mrk_sz)

plot_lines(pool_sfp,data_sfp,pool_col)
plot_lines(mutation_sfp,data_sfp,mutation_col)

# Overplot with training data first so that everything goes on top of the lines,
# and pool and mutation last so they are the top-most layer
plt.plot(data_sfp['x'],data_sfp['y'],data_mrk,color='0.75',markersize=mrk_sz,markeredgecolor='k')
plt.plot(predict_sfp['x'],predict_sfp['y'],predict_mrk,color=predict_col,markersize=mrk_sz,markeredgecolor='k')
plt.plot(mutation_sfp['x'],mutation_sfp['y'],mutation_mrk,color=mutation_col,markersize=mrk_sz,markeredgecolor='k')
plt.plot(pool_sfp['x'],pool_sfp['y'],pool_mrk,color=pool_col,markersize=mrk_sz,markeredgecolor='k')

plt.legend(['POOL','Mutation','Training Data','Predict-then-Optimize'],frameon=False,markerfirst=False)

#plt.xlim(-60,35)
#plt.ylim(-35,40)
plt.xlabel('x coordinate')
plt.ylabel('y coordinate')
#plt.text(-50,15,'Sfp',fontsize=14)
plt.title('Sfp')
plt.savefig('figures/fig3c.eps',transparent=True)
plt.close()


#
# FIG 3D (2-dimensional embedding of training data and recommended points for AcpS)
#

plt.plot(pool_acps['x'],pool_acps['y'],pool_mrk,color=pool_col,markersize=mrk_sz)
plt.plot(mutation_acps['x'],mutation_acps['y'],mutation_mrk,color=mutation_col,markersize=mrk_sz)
plt.plot(data_acps['x'],data_acps['y'],data_mrk,color='0.75',markersize=mrk_sz)
plt.plot(predict_acps['x'],predict_acps['y'],predict_mrk,color=predict_col,markersize=mrk_sz)

plot_lines(pool_acps,data_acps,pool_col)
plot_lines(mutation_acps,data_acps,mutation_col)

# Overplot with training data first so that everything goes on top of the lines,
# and pool and mutation last so they are the top-most layer
plt.plot(data_acps['x'],data_acps['y'],data_mrk,color='0.75',markersize=mrk_sz,markeredgecolor='k')
plt.plot(predict_acps['x'],predict_acps['y'],predict_mrk,color=predict_col,markersize=mrk_sz,markeredgecolor='k')
plt.plot(mutation_acps['x'],mutation_acps['y'],mutation_mrk,color=mutation_col,markersize=mrk_sz,markeredgecolor='k')
plt.plot(pool_acps['x'],pool_acps['y'],pool_mrk,color=pool_col,markersize=mrk_sz,markeredgecolor='k')

plt.legend(['POOL','Mutation','Training Data','Predict-then-Optimize'],frameon=False,markerfirst=False)
#plt.xlim(-60,35)
#plt.ylim(-40,40)
plt.xlabel('x coordinate')
plt.ylabel('y coordinate')
#plt.text(-40,15,'AcpS',fontsize=14)
plt.title('AcpS')
plt.savefig('figures/fig3d.eps',transparent=True)
plt.close()


#
# FIG S3A (Sfp, Histogram of distance to closest point in training data, used to be called Fig 3E)
#

max_val = 40 # Largest value to include in the histogram
bins = np.arange(0,max_val,step=1)

plt.hist(pool_sfp['closest_data_dist'],bins,density=None,alpha=0.5,color=pool_col,histtype='bar', ec='black')
plt.hist(mutation_sfp['closest_data_dist'],bins,density=None,alpha=0.5,color=mutation_col,histtype='bar', ec='black')
plt.legend(['POOL','Mutation'])
plt.xlabel('distance to closest point in original data')
plt.ylabel('frequency')
plt.title('Sfp')
plt.savefig('figures/figS3A.pdf',transparent=True)
plt.close()

#
# FIG S3B (AcpS, Histogram of distance to closest point in training data, used to be called Fig 3F)
#

plt.hist(pool_acps['closest_data_dist'], bins, density=None,alpha=0.5,color=pool_col,histtype='bar', ec='black')
plt.hist(mutation_acps['closest_data_dist'],bins, density=None,alpha=0.5,color=mutation_col,histtype='bar', ec='black')
plt.legend(['POOL','Mutation'])
plt.xlabel('distance to closest point in original data')
plt.ylabel('frequency')
plt.title('AcpS')
plt.savefig('figures/figS3B.pdf',transparent=True)
plt.close()

