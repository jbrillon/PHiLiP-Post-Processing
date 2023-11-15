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
#-----------------------------------------------------
def plotfxn(x_store,y_store,x_label,y_label,
    figure_filename,labels_store,which_lines_dashed_store,
    log_axes=None):
    qp.plotfxn(x_store,y_store,
        figure_filename=figure_filename,
        figure_size=(6,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        title_label="Transient Flow Convergence",
        xlabel=x_label,
        ylabel=y_label,
        which_lines_dashed=which_lines_dashed_store,
        transparent_legend=True,
        legend_border_on=False,
        grid_lines_on=True,
        log_axes=log_axes)
    return
#-----------------------------------------------------
def plot_transient(filenames_,labels_,which_lines_dashed_,
    plot_skin_friction_coefficient=True,
    plot_wall_shear_stress=True,
    plot_bulk_mass_flow=True,
    starting_data_index_for_plot=8
    ):
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
        time_store.append(time[starting_data_index_for_plot:])
        labels_store.append(labels_[i])
        wall_shear_stress_store.append(wall_shear_stress[starting_data_index_for_plot:])
        skin_friction_coefficient_store.append(skin_friction_coefficient[starting_data_index_for_plot:])
        bulk_mass_flow_store.append(bulk_mass_flow[starting_data_index_for_plot:])

    # plot the quantities
    if(plot_skin_friction_coefficient):
        plotfxn(time_store,skin_friction_coefficient_store,\
            "$t$","$C_{f}(t)$/$C_{f}^{expected}$","skin_friction_coefficient",\
            labels_store,which_lines_dashed_store)
    if(plot_wall_shear_stress):
        plotfxn(time_store,wall_shear_stress_store,\
            "$t$","$\\tau_{w}$","wall_shear_stress",\
            labels_store,which_lines_dashed_store)
    if(plot_bulk_mass_flow):
        plotfxn(time_store,bulk_mass_flow_store,\
            "$t$","$\\rho_{b}U_{b}$","bulk_mass_flow",\
            labels_store,which_lines_dashed_store,log_axes="y")
    return
#-----------------------------------------------------
#=====================================================
# Plot the transient quantities
#=====================================================
filenames=[\
"/home/julien/Codes/dummy_dir_for_testing/turbulent_quantities.txt",\
"turbulent_quantities-22892028.txt",\
"turbulent_quantities-22907920.txt",\
"turbulent_quantities-22909869.txt",\
"turbulent_quantities-22909869.txt",\
]
labels=[\
"running: $\\Delta t=6.8\\times10^{-5}$, $\\alpha=0.3$",\
"failed: $\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.3$",\
"running: $\\Delta t=1.75\\times10^{-4}$, $\\alpha=0.3$",\
"running: $\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.2$",\
"failed: $\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.05$",\
]#,
which_lines_dashed=[\
False,\
False,\
False,\
False,\
]
plot_transient(filenames,labels,which_lines_dashed)
