import matplotlib.pyplot as plt


def plt_style():
    plt.style.use('seaborn')
    plt.style.use("my538_2020_03_09_cnf.mplstyle")
    plt.style.use("tableau-colorblind10")


def savefig_check(name):
    plt.savefig("Figures/checks/%s.png" % name)
