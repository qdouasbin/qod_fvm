import os
import glob
import matplotlib.pyplot as plt


# import input_1d_solver as inp


def plt_style():
    plt.style.use('seaborn')
    plt.style.use("my538_2020_03_09_cnf.mplstyle")
    plt.style.use("tableau-colorblind10")


def savefig_check(name):
    plt.savefig("Figures/checks/%s.png" % name)


def savefig_solution(dir, name):
    plt.savefig("%s/%s.png" % (dir, name), dpi=100)


def clean_directories(inp):
    """
    deletes files before new run:
    File types deleted:
      - png
      - hdf5
      - xmf
    """
    png_files = glob.glob(os.path.join(inp.output_dir, "*.png"))
    hdf_files = glob.glob(os.path.join(inp.output_dir, "*.h5"))
    xmf_files = glob.glob(os.path.join(inp.output_dir, "*.xmf"))
    csv_files = glob.glob(os.path.join(inp.output_dir, "*.csv"))

    for file_list in [png_files, hdf_files, xmf_files, csv_files]:
        for my_file in file_list:
            try:
                os.remove(my_file)
            except:
                print("Error while deleting file : ", my_file)
