#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
# Import personal libraries
# from finite_difference_library import first_derivative, fd_non_uniform_grid
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
#=====================================================
# Helper functions
#=====================================================
def plot_transient(filenames_,labels_,which_lines_dashed_):
    expected_mean_value_for_skin_friction_coefficient = 6.25e-3 # from Lodato's source term paper

    # data store
    time_store=[]
    skin_friction_coefficient_store=[]
    wall_shear_stress_store=[]
    skin_friction_coefficient_store=[]
    bulk_mass_flow_store=[]
    # plot function inputs
    labels_store=[]
    which_lines_dashed_store=[]

    for i,filename in enumerate(filenames_):
        # load data
        time,wall_shear_stress,skin_friction_coefficient,bulk_density,bulk_velocity = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        # compute the bulk mass flow
        bulk_mass_flow = bulk_density*bulk_velocity
        # compute the skin friction coefficient
        skin_friction_coefficient /= expected_mean_value_for_skin_friction_coefficient

        # store the data
        time_store.append(time)
        labels_store.append(labels_[i])
        wall_shear_stress_store.append(wall_shear_stress)
        skin_friction_coefficient_store.append(skin_friction_coefficient)
        bulk_mass_flow_store.append(bulk_mass_flow)

    # plot the quantities
    figure_filename="skin_friction_coefficient"
    qp.plotfxn(time_store,skin_friction_coefficient_store,
        figure_filename=figure_filename,
        figure_size=(6,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        title_label="Transient Flow Convergence",
        xlabel="$t$",
        ylabel="$C_{f}(t)$/$C_{f}^{expected}$",
        which_lines_dashed=which_lines_dashed_store)
    figure_filename="wall_shear_stress"
    qp.plotfxn(time_store,wall_shear_stress_store,
        figure_filename=figure_filename,
        figure_size=(6,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        title_label="Transient Flow Convergence",
        xlabel="$t$",
        ylabel="$\\tau_{w}$",
        which_lines_dashed=which_lines_dashed_store)
    figure_filename="bulk_mass_flow"
    qp.plotfxn(time_store,bulk_mass_flow_store,
        figure_filename=figure_filename,
        figure_size=(6,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        title_label="Transient Flow Convergence",
        xlabel="$t$",
        ylabel="$\\rho_{b}U_{b}$",
        which_lines_dashed=which_lines_dashed_store)
    return
#=====================================================
# Plot the transient quantities
#=====================================================
filenames=[\
"turbulent_quantities-1.txt",\
"/home/julien/Codes/dummy_dir_for_testing/turbulent_quantities.txt",\
]
labels=[\
"failed",\
"local",\
]
which_lines_dashed=[\
False,\
True,\
]
plot_transient(filenames,labels,which_lines_dashed)

