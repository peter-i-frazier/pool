import sqlalchemy
import pandas
import numpy as np

USER = "jialeiwang"
PASSWORD = "wangjialei123"
HOST = "work.cxcjqzn7ydtp.us-east-1.rds.amazonaws.com"
DB = "sfp_AcpS"

engine = sqlalchemy.create_engine("mysql+pymysql://{0}:{1}@{2}/{3}".format(USER, PASSWORD, HOST, DB))
benchmark_seq_table = pandas.read_sql_table('benchmark_seq', engine)
original_seq_table = pandas.read_sql_table('original_seq', engine)

# below each object is a list of strings, where each string is nterm / cterm of the corresponding sequence
#pool_nterm = benchmark_seq_table['pool_nterm'].values
#pool_cterm = benchmark_seq_table['pool_cterm'].values
#
#mutation_nterm = benchmark_seq_table['mutation_nterm'].values
#mutation_cterm = benchmark_seq_table['mutation_cterm'].values
#
naive_nterm = benchmark_seq_table['naive_nterm'].values
#naive_cterm = benchmark_seq_table['naive_cterm'].values

original_nterm = original_seq_table['original_nterm'].values
#original_cterm = original_seq_table['original_cterm'].values

## Reformat data
X_ben = np.asarray(benchmark_seq_table)
X_ben = np.vstack(np.split(X_ben, [2, 4], axis=1))
X_orig = np.asarray(original_seq_table)
X = np.vstack((X_orig, X_ben))
#print X.shape
#print original_nterm[:5]
#print X[:5, :]
#print X
#print naive_nterm

### Need to use this number in visualization code
print 'n_orig =', len(original_nterm)
np.savetxt('data.csv', X, delimiter=',', fmt='%s')


