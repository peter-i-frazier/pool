import matplotlib.pyplot as plt
import numpy as np

def plot_visualization(data, name, fig_path=None):
    n_orig = data['n_orig_{0}'.format(name)]
    n_benchmark = data['n_{0}'.format(name)]
    coords = data[name]
    plt.figure()
    colors = np.concatenate((np.repeat(['#3333cc'], n_orig), np.repeat(['#33ccff', '#ffff00', '#cc3300'], n_benchmark)))
    plt.scatter(coords[:, 0], coords[:, 1], s=40, c=colors, alpha=.7, edgecolors='none')
    if fig_path:
        plt.savefig(fig_path)

def plot_coord(data, name, fig_path=None):
    hit_idx = np.array(data["orig_hit_"+name])
    plt.figure(figsize=(10,10))
    plt.scatter(data["orig_"+name][hit_idx==0,0], data["orig_"+name][hit_idx==0,1], c='#E1E5E4', label="Data miss")
    plt.scatter(data["orig_"+name][hit_idx==1,0], data["orig_"+name][hit_idx==1,1], s=40, c='#E5000A', label="Data hit")
    plt.scatter(data["pool_"+name][:,0], data["pool_"+name][:,1], s=40, c='#338DFF', label="POOL")
    plt.scatter(data["mutate_"+name][:,0], data["mutate_"+name][:,1], c='#E9FD17', label="Mutation")
    plt.scatter(data["predict_optimize_"+name][:,0], data["predict_optimize_"+name][:,1], c='#00990C', label="Predict-then-optimize")
    ax = plt.gca()
    ax.legend(loc='upper left')
    if fig_path:
        plt.savefig(fig_path+".pdf")

def plot_benchmark(table, fig_path=None):
    plt.figure()
    plt.plot(range(1, len(table.index)+1), table['pool'] * 100, label="POOL", color=(0, 112./255, 192./255))
    plt.plot(range(1, len(table.index)+1), table['mutate'] * 100, label="Mutate", color="darkgreen")
    plt.plot(range(1, len(table.index)+1), table['predict_optimize'] * 100, label="Predict-then-optimize", color="red")
    plt.xlabel("number of peptides recommended", fontsize=26)
    plt.ylabel("Probability %", fontsize=26)
    plt.tick_params(labelsize=12)
    plt.legend(loc="upper left", fontsize=12)
    plt.subplots_adjust(bottom=0.115, top=0.9)
    if fig_path:
        plt.savefig(fig_path+".pdf")

def plot_roc(table, fig_path=None):
    plt.figure()
    plt.plot(table.x, table.y)
    plt.plot([0,1], [0,1], linestyle='--', color='grey')
    plt.xlabel("False Positive Rate (FPR)", fontsize=26)
    plt.ylabel("True Positive Rate (TPR)", fontsize=26)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tick_params(labelsize=12)
    plt.subplots_adjust(bottom=0.115, top=0.9)
    if fig_path:
        plt.savefig(fig_path)
