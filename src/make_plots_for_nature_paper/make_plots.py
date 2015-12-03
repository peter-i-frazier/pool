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

plot_roc("ROC_sfp")
plot_roc("ROC_AcpS")
plot_roc("ROC_PfAcpH")

# plot benchmark
def plot_benchmark(table_name):
    df = pandas.read_sql_table(table_name, engine)
    fig = plt.figure()
    plt.plot(range(len(df.index)), df['pool'] * 100, label="POOL")
    plt.plot(range(len(df.index)), df['mutate'] * 100, label="Mutate")
    plt.plot(range(len(df.index)), df['predict_optimize'] * 100, label="Predict-then-optimize")
    plt.xlabel("number of peptides recommended", fontsize=26)
    plt.ylabel("P(at least one peptide is a hit) %", fontsize=26)
    plt.tick_params(labelsize=12)
    plt.legend(loc="upper left", fontsize=12)
    plt.subplots_adjust(bottom=0.115, top=0.9)
    plt.savefig(plot_path + table_name + ".pdf")
plot_benchmark("benchmark_type1")
plot_benchmark("benchmark_type2")
