__author__ = 'jialeiwang'
import sqlalchemy
import pandas
import matplotlib.pyplot as plt

USER = "jialeiwang"
PASSWORD = "wangjialei123"
HOST = "work.cxcjqzn7ydtp.us-east-1.rds.amazonaws.com"
DB = "sfp_AcpS"

engine = sqlalchemy.create_engine("mysql+pymysql://{0}:{1}@{2}/{3}".format(USER, PASSWORD, HOST, DB))
plot_path = "./"

# plot ROC
def plot_roc(table_name):
    df = pandas.read_sql_table(table_name, engine)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    df.plot(x='x', y='y', legend=False)
    # ax.legend_.remove()
    plt.xlabel("False Positive Rate (FPR)", fontsize=30)
    plt.ylabel("True Positive Rate (TPR)", fontsize=30)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tick_params(labelsize=20)
    plt.subplots_adjust(bottom=0.15, top=0.9)
    plt.savefig(plot_path + table_name + ".pdf")


# plot benchmark
def plot_benchmark(table_name):
    df = pandas.read_sql_table(table_name, engine)
    fig = plt.figure()
    plt.plot(range(1, len(df.index)+1), df['pool'] * 100, label="POOL", color=(0, 112./255, 192./255))
    plt.plot(range(1, len(df.index)+1), df['mutate'] * 100, label="Mutate", color="darkgreen")
    plt.plot(range(1, len(df.index)+1), df['predict_optimize'] * 100, label="Predict-then-optimize", color="red")
    plt.xlabel("number of peptides recommended", fontsize=26)
    plt.ylabel("P(at least one peptide is a hit) %", fontsize=26)
    plt.tick_params(labelsize=12)
    plt.legend(loc="upper left", fontsize=12)
    plt.subplots_adjust(bottom=0.115, top=0.9)
    plt.savefig(plot_path + table_name + ".pdf")
    # plt.show()

# plot local search
def plot_local_search():
    df = pandas.read_sql_table("local_search", engine)
    fig = plt.figure()
    plt.plot(range(1, len(df.index)+1), df['pool_start'] * 100, label="start from POOL recommended peptides", color=(0, 112./255, 192./255))
    plt.plot(range(1, len(df.index)+1), df['random_start'] * 100, label="start from random set of peptides", color="red")
    plt.xlabel("Iteration", fontsize=26)
    plt.ylabel("P(at least one peptide is a hit) %", fontsize=26)
    plt.tick_params(labelsize=12)
    plt.legend(loc="lower right", fontsize=12)
    plt.subplots_adjust(bottom=0.115, top=0.9)
    plt.savefig(plot_path + "local_search.pdf")

if __name__ == "__main__":
    # plot_roc("ROC_sfp")
    # plot_roc("ROC_AcpS")
    # plot_roc("ROC_PfAcpH")
    #
    # plot_benchmark("benchmark_type1")
    # plot_benchmark("benchmark_type2")

    # plot_local_search()

    df_random_71 = pandas.read_sql_table("local_search", engine)
    df_pool_71 = pandas.read_sql_table("local_search_pool_71", engine)
    df_random_500 = pandas.read_sql_table("local_search_random_500", engine)
    df_pool_500 = pandas.read_sql_table("local_search_pool_500", engine)
    fig = plt.figure()
    plt.plot(range(len(df_pool_71.index)), df_pool_71['PI'] * 100, label="start from POOL recommended peptides", color=(0, 112./255, 192./255))
    plt.plot(range(len(df_random_71.index)), df_random_71['random_start'] * 100, label="start from random set of peptides", color="red")
    # plt.plot(range(len(df_pool_500.index)), df_pool_500['PI'] * 100, label="start from POOL recommended peptides", color=(0, 112./255, 192./255))
    # plt.plot(range(len(df_random_500.index)), df_random_500['PI'] * 100, label="start from random set of peptides", color="red")
    plt.xlabel("Iteration", fontsize=26)
    plt.ylabel("P(at least one peptide is a hit) %", fontsize=26)
    plt.tick_params(labelsize=12)
    plt.legend(loc="lower right", fontsize=12)
    plt.subplots_adjust(bottom=0.115, top=0.9)
    plt.savefig(plot_path + "local_search_71_peptides.pdf")
