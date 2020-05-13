import os
import glob
import matplotlib.pyplot as plt


def plt_style():
    """
    Set personal matplotlib settings
    """
    plt.style.use('seaborn')
    lst_styles = glob.glob('*/*mplstyle')
    plt.style.use(lst_styles[-1])
    plt.style.use("tableau-colorblind10")


def savefig_check(name, path='None'):
    """
    Save a matplotlib figure as sanity check
    :param name:name of figure
    """
    if not path:
        path = os.path.join('Figure', 'checks')

    if not os.path.exists(path):
        os.makedirs(path)

    plt.savefig(os.path.join(path, "%s.png" % name))


def savefig_solution(dir, name):
    """
    Save a matplotlib figure as a png ain this path
    :param dir: path
    :param name: figure name
    """
    plt.savefig("%s/%s.png" % (dir, name), dpi=100)


def clean_directories(params_IO):
    """
    deletes files before new run:
    File types deleted:
      - png
      - hdf5
      - xmf
    """
    print(params_IO)
    out_dir = params_IO["directory"]

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    png_files = glob.glob(os.path.join(out_dir, "*.png"))
    hdf_files = glob.glob(os.path.join(out_dir, "*.h5"))
    xmf_files = glob.glob(os.path.join(out_dir, "*.xmf"))
    csv_files = glob.glob(os.path.join(out_dir, "*.csv"))

    for file_list in [png_files, hdf_files, xmf_files, csv_files]:
        for my_file in file_list:
            try:
                os.remove(my_file)
            except:
                print("Error while deleting file : ", my_file)
