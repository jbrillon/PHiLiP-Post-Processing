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
    log_axes=None,xlimits=[],ylimits=[]):
    qp.plotfxn(x_store,y_store,
        figure_filename=figure_filename,
        figure_size=(8,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        title_label="WMLES: TCF $c_{DG}$ NSFR.IR.GLL 20x10x10 p4, $CFL\\approx0.2$",
        xlabel=x_label,
        ylabel=y_label,
        xlimits=xlimits,
        ylimits=ylimits,
        which_lines_dashed=which_lines_dashed_store,
        transparent_legend=False,
        legend_border_on=True,
        grid_lines_on=True,
        log_axes=log_axes,
        legend_location="best")
    return
#-----------------------------------------------------
def plot_transient(filenames_,labels_,which_lines_dashed_=[],
    plot_skin_friction_coefficient=True,
    plot_wall_shear_stress=True,
    plot_bulk_mass_flow=True,
    starting_data_index_for_plot=0,
    friction_velocity_based_reynolds_number=[]
    ):
    # expected_mean_value_for_skin_friction_coefficient = 6.25e-3 # from Lodato's source term paper
    # expected_mean_value_for_skin_friction_coefficient = 0.00347754 # for Re=5200
    # expected_mean_value_for_wall_shear_stress = 0.00161876 # for Re=5200
    # update above using Dean's expression
    # data store
    time_store=[]
    skin_friction_coefficient_store=[]
    wall_shear_stress_store=[]
    skin_friction_coefficient_store=[]
    bulk_mass_flow_store=[]
    # plot function inputs
    labels_store=[]

    for i,filename in enumerate(filenames_):
        # load data
        time,wall_shear_stress,skin_friction_coefficient,bulk_density,bulk_velocity = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        # compute the bulk mass flow
        bulk_mass_flow = bulk_density*bulk_velocity

        expected_mean_value_for_wall_shear_stress = 0.00161876 # for Re=5200
        expected_mean_value_for_skin_friction_coefficient = 0.00347754 # for Re=5200
        if(friction_velocity_based_reynolds_number[i]==5200):
            print("here")
            expected_mean_value_for_wall_shear_stress = 0.00161876 # for Re=5200
            # expected_mean_value_for_skin_friction_coefficient = 0.00347754 # for Re=5200
            expected_mean_value_for_skin_friction_coefficient = 2.0*expected_mean_value_for_wall_shear_stress
        elif(friction_velocity_based_reynolds_number[i]==395):
            print("here2")
            expected_mean_value_for_wall_shear_stress = 0.003380747 # for Re=395
            expected_mean_value_for_skin_friction_coefficient = 2.0*expected_mean_value_for_wall_shear_stress

        # store the data
        time_store.append(time[starting_data_index_for_plot:])
        labels_store.append(labels_[i])
        wall_shear_stress_store.append(wall_shear_stress[starting_data_index_for_plot:]/expected_mean_value_for_wall_shear_stress)
        skin_friction_coefficient_store.append(skin_friction_coefficient[starting_data_index_for_plot:]/expected_mean_value_for_skin_friction_coefficient)
        bulk_mass_flow_store.append(bulk_mass_flow[starting_data_index_for_plot:]/bulk_mass_flow[0])

    # plot the quantities
    if(plot_skin_friction_coefficient):
        plotfxn(time_store,skin_friction_coefficient_store,\
            "$t^{*}$","Normalized Skin Friction Coefficient, $C_{f}/C^{expected}_{f}$","skin_friction_coefficient",\
            labels_store,which_lines_dashed_)
            #,xlimits=[30,70],ylimits=[3e-4,4e-4])
        # plotfxn(time_store,skin_friction_coefficient_store/expected_mean_value_for_skin_friction_coefficient,\
        #     "$t$","$C_{f}(t)$/$C_{f}^{expected}$","skin_friction_coefficient",\
        #     labels_store,which_lines_dashed_)
    if(plot_wall_shear_stress):
        plotfxn(time_store,wall_shear_stress_store,\
            "$t^{*}$","Normalized Nondim. Wall Shear Stress, $\\tau_{w}/\\tau^{expected}_{w}$","wall_shear_stress",\
            labels_store,which_lines_dashed_)
    if(plot_bulk_mass_flow):
        plotfxn(time_store,bulk_mass_flow_store,\
            "$t^{*}$","Normalized Nondim. Bulk Mass Flow Rate, $\\rho_{b}U_{b}/(\\rho_{b}U_{b})_{initial}$","bulk_mass_flow",\
            labels_store,which_lines_dashed_,log_axes="y")
    return
#-----------------------------------------------------
#-----------------------------------------------------
def plot_boundary_layer_profile(filenames_,labels_,which_lines_dashed_=[],
    ):
    # expected_mean_value_for_skin_friction_coefficient = 6.25e-3 # from Lodato's source term paper
    # expected_mean_value_for_skin_friction_coefficient = 0.00347754 # for Re=5200
    # expected_mean_value_for_wall_shear_stress = 0.00161876 # for Re=5200
    # update above using Dean's expression
    # data store
    y_coordinate_store=[]
    average_x_velocity_store=[]
    # plot function inputs
    labels_store=[]

    for i,filename in enumerate(filenames_):
        # load data
        y_coordinate, average_x_velocity = np.loadtxt(filename,dtype=np.float64,unpack=True)

        # store the data
        y_coordinate_store.append(y_coordinate)
        labels_store.append(labels_[i])
        average_x_velocity_store.append(average_x_velocity)
        print(average_x_velocity)

    # swap the x and y data here
    plotfxn(y_coordinate_store,average_x_velocity_store,\
        "$y^{*}$","Average x-velocity, $u_{avg.}$","boundary_layer_profile",\
        labels_store,which_lines_dashed_,log_axes=None)
    return
#-----------------------------------------------------
#=====================================================
# Plot the transient quantities
#=====================================================
filenames=[\
# "turbulent_quantities-22907920.txt",\
# "turbulent_quantities-22909869.txt",\
# "turbulent_quantities-22909869.txt",\
# "turbulent_quantities-22911528.txt",\
# "turbulent_quantities-22916636.txt",\
# "turbulent_quantities-22916869.txt",\
# "turbulent_quantities-22931922.txt",\
# only below this is good to plot
# "turbulent_quantities-22892028.txt",\
# "turbulent_quantities-23080330.txt",\
# "turbulent_quantities-23077970.txt",\
# "turbulent_quantities-23117286.txt",\
# "turbulent_quantities-23117307.txt",\
# "turbulent_quantities-23134131.txt",\
# "turbulent_quantities-23135133.txt",\
# "turbulent_quantities-23117307.txt",\
# "turbulent_quantities-nov25-local-roe.txt",\
# "turbulent_quantities-23226879.txt",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GL_OI-0_Re395_p3/turbulent_quantities.txt",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_Re395_p3/turbulent_quantities.txt",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_Re395_p3/turbulent_quantities.txt",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cPlus_IR_2PF-Roe_GL_OI-0_Re395_p3/turbulent_quantities.txt",\
"turbulent_quantities.txt",\
]
labels=[\
# "running: $\\Delta t=1.75\\times10^{-4}$, $\\alpha=0.3$",\
# "running: $\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.2$",\
# "failed: $\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.05$",\
# "running: $\\Delta t=6.8\\times10^{-5}$, $\\alpha=0.3$",\
# "running: $\\Delta t=3.5\\times10^{-4}$, $\\alpha=0.1$",\
# "running: $\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.05$",\
# "running: $\\Delta t=1.0\\times10^{-4}$, $\\alpha=0.3$",\
# only below this is good to plot
# "$\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.3$, before fixes, Turb. IC",\
# "$\\Delta t=5\\times10^{-5}$, $\\alpha=0.3$",\
# "$\\Delta t=5\\times10^{-5}$, $\\alpha=0.3$, $|\\tau_{w}|$",\
# "$\\Delta t=1\\times10^{-5}$, $\\alpha=0.3$",\
# "$\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.3$, Chao fix, Turb. IC",\
# "$\\Delta t=1.5\\times10^{-4}$, $\\alpha=0.3$, Chao fix, Turb. IC",\
# "$\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.3$, Chao fix, Lam. IC",\
# "$\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.3$, Chao",\
# "$\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.3$, Chao, Laminar",\
# "$\\Delta t=3.44\\times10^{-4}$, $\\alpha=0.0$, Chao, Brian, Roe, Laminar",\
# "$\\Delta t=6.8\\times10^{-5}$, $\\alpha=0.0$, Chao, Brian, 2PF, Laminar",\
"$c_{DG}$ NSFR.IR",\
"$c_{+}$ NSFR.IR",\
"$c_{DG}$ NSFR.IR.Roe",\
"$c_{+}$ NSFR.IR.Roe",\
"constant source term"
]
# uncomment for the old results
# plot_transient(filenames,labels,which_lines_dashed=[2,3])

'''
filenames=[\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/turbulent_quantities.txt",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re395_p4_20x10x10_turbulent_initialization/turbulent_quantities.txt",\
]
labels=[\
"$Re_{\\tau}\\approx5200$",\
"$Re_{\\tau}\\approx395$",\
]
which_lines_dashed=[]
friction_velocity_based_reynolds_number=[5200,395]
plot_transient(filenames,labels,starting_data_index_for_plot=0,friction_velocity_based_reynolds_number=friction_velocity_based_reynolds_number)
'''

# plot boundary layer profile
filenames=[\
"/Users/Julien/NarvalFiles/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile.dat",\
]
labels=[\
"$Re_{\\tau}\\approx5200$",\
]
which_lines_dashed=[]
plot_boundary_layer_profile(filenames,labels)