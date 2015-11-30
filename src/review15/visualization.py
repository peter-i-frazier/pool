from sklearn.manifold import TSNE
import numpy as np
import matplotlib.pyplot as plt

"""
Functions computing and visualizing distance between peptides
"""

################### Global variables ##########################
# Penalties, must have length of 3 and that
# pens[0] > pens[1] > pens[2] > 0
# Modify to make plot perfect
pens = [10., 1.5, 1.]
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
def compute_dist_one_side(str1, str2):
    # alpha, beta, gamma = pens
    pen = 0.  # total penalty
    if len(str1) < len(str2):
        short = str1
        long = str2
    else:
        short = str2
        long = str1
    pen += pens[2] * (len(long) - len(short))
    for j in xrange(len(short)):
        s = short[j]
        l = long[j]
        if s != l:
            if class_dict[s] != class_dict[l]:
                pen += pens[0]
            else:
                pen += pens[1]
    return pen

# This one is obsolete due to change in input format
#def compute_dist(str1, str2):
#    print '******Computing dist*******'
#    print str1, str2
#    # Split string on the first occurrence of S
#    strlist1 = str1.split('S', 1)
#    strlist2 = str2.split('S', 1)
#    # Substrings before Serine must be reversed
#    # before being parsed to function
#    pen1 = compute_dist_one_side(strlist1[0][::-1],
#                                 strlist2[0][::-1])
#    print 'left pen', pen1
#    pen2 = compute_dist_one_side(strlist1[1],
#                                 strlist2[1])
#    print 'right pen', pen2
#    return pen1 + pen2

def compute_dist(strlist1, strlist2):
    #print '******Computing dist*******'
    # Substrings before Serine must be reversed
    # before being parsed to function
    pen1 = compute_dist_one_side(strlist1[0][::-1],
                                 strlist2[0][::-1])
    #print 'left pen', pen1
    pen2 = compute_dist_one_side(strlist1[1],
                                 strlist2[1])
    #print 'right pen', pen2
    return pen1 + pen2

def compute_dist_mat(X):
    """
    X: numpy.ndarray of str, data

    return: numpy.ndarray, distance between strings
            a square and symmetric matrix with dim len(X)
    """
    print 'Computing distance matrix...'
    mat = np.zeros((len(X), len(X)))
    for i in xrange(len(mat)):
        for j in xrange(i):
            #d = compute_dist(X[i], X[j])
            d = compute_dist(X[i, :], X[j, :])
            mat[i, j] = d
            mat[j, i] = d
    print 'Done'
    return mat

def compute_rel_pos(X):
    dist_mat = compute_dist_mat(X)
    print 'Computing coordinates for plotting...'
    tsne = TSNE(metric='precomputed', random_state=0)
    rel_pos = tsne.fit_transform(dist_mat)
    print 'Done'
    np.savetxt('rel_pos.csv', rel_pos, delimiter=',')
    return rel_pos

def plot(X, n_orig):
    assert (len(X) - n_orig) % 3 == 0
    plt.style.use('ggplot')
    rel_pos = compute_rel_pos(X)
    N = (len(X) - n_orig) / 3  # size of other groups
    colors = np.concatenate((np.zeros(n_orig), np.repeat([1, 2, 3], N)))
    plt.figure()
    plt.scatter(rel_pos[:, 0], rel_pos[:, 1],
                s=40, c=colors, alpha=.7, edgecolors='none')
    plt.title('alpha={0}, beta={1}, gamma={2}'.format(pens[0], pens[1], pens[2]))
    plt.savefig('plot.pdf')

## No need to call this when data is generated from get_data.py
def format_data(X):
    # Compute number of original peptides
    orig = X[:, 0]
    n_orig = len(orig[orig != ''])
    # Reformat data into two columns
    X_clean = np.vstack(np.split(X, [2, 4, 6], axis=1))
    return X_clean[~np.any(X_clean == '', axis=1)], n_orig

if __name__ == '__main__':
    n_orig = 835
    # numpy array of numpy strings
    #X = np.genfromtxt('test_input.csv', delimiter=',', dtype=str)
    X = np.genfromtxt('data.csv', delimiter=',', dtype=str)
    #X = np.genfromtxt('toy_visualization_data.csv',
    #                  delimiter=',', dtype=str, skip_header=1)
    #X, n_orig = format_data(X)
    #print X
    #print compute_dist_mat(X)
    plot(X, n_orig)
    print 'Penalties', pens
    print 'n_orig =', n_orig
