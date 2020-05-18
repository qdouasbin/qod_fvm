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


def savefig_solution(name):
    """
    Save a matplotlib figure as a png ain this path
    :param dir: path
    :param name: figure name
    """
    print("\t--> Saving %s" % name)
    plt.savefig(name, bbox_inches='tight', pad_inches=0.01, dpi=200)



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


def plot_prim_var_field(dom_1D, fig_name):
    """
    Plot the field of primitive variable with matplotlib
    :param dom_1D: domain (Field object)
    :param fig_name: file name
    """
    xi = dom_1D.get_xi()
    area = dom_1D.get_area()
    area_max = area.max()
    temp = dom_1D.get_T()
    temp_max = temp.max()
    pres = dom_1D.get_P()
    pres_max = pres.max()
    rho = dom_1D.get_rho()
    rho_max = rho.max()

    u = dom_1D.get_u()
    u_max = np.amax(u)

    plt.figure()
    plt.plot(xi, area / area_max, ':.', label='A/%2.2e' % area_max)
    plt.plot(xi, temp / temp_max, '--+', label='T/%2.2e' % temp_max)
    plt.plot(xi, pres / pres_max, '--x', label='P/%2.2e' % pres_max)
    plt.plot(xi, rho / rho_max, '--', label=r'$\rho$/%2.2e' % rho_max)

    if (u_max):
        plt.plot(xi, u / u_max, '--', label=r'u/%2.2e' % u_max)
    else:
        plt.plot(xi, u, '--', label=r'u [m/s]')

    plt.legend()
    plt.xlabel(r'x [m]')
    plt.savefig("Figures/checks/%s.png" % fig_name)
    plt.show()
