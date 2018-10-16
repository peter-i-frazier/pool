import sys
from sklearn.manifold import TSNE
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas
import pickle

"""
Functions computing and visualizing distance between peptides
"""
np.random.seed(0)
################### Global variables ##########################
# DO NOT modify if you do not know what you are doing
class_dict = {'A': 5,
              'R': 4,
              'N': 2,
              'D': 1,
              'C': 8,
              'Q': 2,
              'E': 1,
              'G': 6,
              'H': 4,
              'I': 5,
              'L': 5,
              'K': 4,
              'M': 5,
              'F': 3,
              'P': 6,
              'S': 7,
              'T': 7,
              'W': 3,
              'Y': 3,
              'V': 5}

########################## Methods #############################
def compute_dist_one_side(str1, str2, penalties):
    # alpha, beta, gamma = penalties
    pen = 0.  # total penalty
    if len(str1) < len(str2):
        short = str1
        long = str2
    else:
        short = str2
        long = str1
    pen += penalties[2] * (len(long) - len(short))
    for j in range(len(short)):
        s = short[j]
        l = long[j]
        if s != l:
            if class_dict[s] != class_dict[l]:
                pen += penalties[0]
            else:
                pen += penalties[1]
    return pen

def compute_dist(strlist1, strlist2, penalties):
    #print '******Computing dist*******'
    # Substrings before Serine must be reversed
    # before being parsed to function
    pen1 = compute_dist_one_side(strlist1[0][::-1], strlist2[0][::-1], penalties)
    #print 'left pen', pen1
    pen2 = compute_dist_one_side(strlist1[1], strlist2[1], penalties)
    #print 'right pen', pen2
    return pen1 + pen2

def compute_dist_mat(X, penalties):
    """
    X: numpy.ndarray of str, data

    return: numpy.ndarray, distance between strings
            a square and symmetric matrix with dim len(X)
    """
    print('Computing distance matrix...')
    mat = np.zeros((len(X), len(X)))
    for i in range(len(mat)):
        for j in range(i):
            #d = compute_dist(X[i], X[j])
            d = compute_dist(X[i, :], X[j, :], penalties)
            mat[i, j] = d
            mat[j, i] = d
    print('Done')
    return mat

def compute_rel_pos(X, penalties):
    dist_mat = compute_dist_mat(X, penalties)
    print('Computing coordinates for plotting...')
    tsne = TSNE(metric='precomputed', random_state=0)
    rel_pos = tsne.fit_transform(dist_mat)
    print('Done')
    return rel_pos

def plot(X, n_orig, penalties):
    assert (len(X) - n_orig) % 3 == 0
    plt.style.use('ggplot')
    rel_pos = compute_rel_pos(X, penalties)
    N = (len(X) - n_orig) / 3  # size of other groups
    # colors = np.concatenate((np.zeros(n_orig), np.repeat([1, 2, 3], N)))
    colors = np.concatenate((np.repeat(['b'], n_orig), np.repeat(['r', 'y', 'c'], N)))
    end_idx_list = range(n_orig, len(X)+1)
    with PdfPages('plots.pdf') as pdf:
        for end_idx in end_idx_list:
            plt.figure()
            plt.scatter(rel_pos[:end_idx, 0], rel_pos[:end_idx, 1],
                        s=40, c=colors[:end_idx], alpha=.7, edgecolors='none')
            # plt.title('alpha={0}, beta={1}, gamma={2}'.format(penalties[0], penalties[1], penalties[2]))
            pdf.savefig()
            plt.close()
            # plt.savefig('plot.pdf')

def plot_animation(X, n_orig, penalties):
    assert (len(X) - n_orig) % 3 == 0
    plt.style.use('ggplot')
    rel_pos = compute_rel_pos(X, penalties)
    N = (len(X) - n_orig) / 3  # size of other groups
    colors = np.concatenate((np.repeat('#3C3C3C', n_orig), np.repeat(['#0070C0', 'darkgreen', '#FF0000'], N)))
    pdf_names = ['POOL', 'mutation', 'predict-then-optimize']
    for method_idx in range(3):
        with PdfPages('{0}_animation.pdf'.format(pdf_names[method_idx])) as pdf:
            for itr in range(N):
                include_indices = np.concatenate((np.arange(n_orig), np.arange(n_orig+method_idx*N, n_orig+method_idx*N+itr+1)))
                plt.figure()
                plt.scatter(rel_pos[include_indices, 0], rel_pos[include_indices, 1],
                            s=40, c=colors[include_indices], alpha=.7, edgecolors='none')
                pdf.savefig()
                plt.close()

## No need to call this when data is generated from get_data.py
def format_data(X):
    # Compute number of original peptides
    orig = X[:, 0]
    n_orig = len(orig[orig != ''])
    # Reformat data into two columns
    X_clean = np.vstack(np.split(X, [2, 4, 6], axis=1))
    return X_clean[~np.any(X_clean == '', axis=1)], n_orig

if __name__ == '__main__':
    # Penalties, must have length of 3 and that
    # penalties[0] > penalties[1] > penalties[2] > 0
    # Modify to make plot perfect
    penalties = [10., 1.5, 1.]
    enzyme = sys.argv[1]
    total_X = [[], []]
    seq_table = pandas.read_csv('figures/{0}_specific_simulation_seq.csv'.format(enzyme))
    total_X[0].extend(seq_table['pool_nterm'].values)
    total_X[0].extend(seq_table['mutation_nterm'].values)
    total_X[0].extend(seq_table['naive_nterm'].values)
    total_X[1].extend(seq_table['pool_cterm'].values)
    total_X[1].extend(seq_table['mutation_cterm'].values)
    total_X[1].extend(seq_table['naive_cterm'].values)
    train_data = pandas.DataFrame()
    for i in range(3):
        subdata = pandas.read_csv('data/TS{0}.csv'.format(i))
        subdata = subdata.loc[subdata['nterm'].notnull(), :]
        train_data = train_data.append(subdata, ignore_index=True)
    total_X[0].extend(train_data['nterm'].values)
    total_X[1].extend(train_data['cterm'].values)
    total_X = np.array(total_X).T
    rel_pos = compute_rel_pos(total_X, penalties)
    table = pandas.DataFrame(rel_pos, columns=['x', 'y'])
    table.to_csv('figures/{0}_specific_simulation_coord.csv'.format(enzyme), index=False)
