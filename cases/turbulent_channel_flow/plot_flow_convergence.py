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
    log_axes=None,xlimits=[],ylimits=[],which_lines_black=[]):
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
        which_lines_black=which_lines_black,
        which_lines_dashed=which_lines_dashed_store,
        transparent_legend=True,
        legend_border_on=False,
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
            expected_mean_value_for_wall_shear_stress = 0.00161876 # for Re=5200
            # expected_mean_value_for_skin_friction_coefficient = 0.00347754 # for Re=5200
            expected_mean_value_for_skin_friction_coefficient = 2.0*expected_mean_value_for_wall_shear_stress
        elif(friction_velocity_based_reynolds_number[i]==395):
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
def plot_boundary_layer_profile(filenames_,labels_,which_lines_dashed_=[]):
    # expected_mean_value_for_skin_friction_coefficient = 6.25e-3 # from Lodato's source term paper
    # expected_mean_value_for_skin_friction_coefficient = 0.00347754 # for Re=5200
    # expected_mean_value_for_wall_shear_stress = 0.00161876 # for Re=5200
    # update above using Dean's expression
    # data store
    average_y_plus_store=[]
    average_u_plus_store=[]
    average_velocity_fluctuation_rms_store=[]
    # plot function inputs
    labels_store=[]

    # load reference data
    y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/lee_moser_2015_dns.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    average_y_plus_store.append(y_plus_reference)
    average_u_plus_store.append(average_u_plus_reference)
    labels_store.append("DNS [Lee \\& Moser, 2015]")

    for i,filename in enumerate(filenames_):
        # load data
        y_coordinate, average_x_velocity, average_kinematic_viscosity, average_x_velocity_fluctuation_rms, average_y_velocity_fluctuation_rms, average_z_velocity_fluctuation_rms = np.loadtxt(filename,dtype=np.float64,unpack=True)

        number_of_points = np.size(y_coordinate)
        y_plus = []
        u_plus = []
        u_plus_fluctuation = []
        v_plus_fluctuation = []
        w_plus_fluctuation = []
        for j in range(0,number_of_points):
            y_value = y_coordinate[j] # nondimensional
            distance_from_wall = 0.0 # nondimensional
            if(y_value < 0.0):
                distance_from_wall = 1.0 + y_value
            else:
                distance_from_wall = 1.0 - y_value
            Re_tau = 5200.0
            Re_inf = 129245.0
            nondimensional_density = 1.0
            nondimensional_viscosity = 1.0
            nondimensional_friction_velocity = Re_tau/Re_inf
            nondimensional_kinematic_viscosity = nondimensional_viscosity/nondimensional_density
            y_plus_local = Re_inf*distance_from_wall*nondimensional_friction_velocity/nondimensional_kinematic_viscosity
            # y_plus_local = Re_inf*distance_from_wall*nondimensional_friction_velocity/average_kinematic_viscosity[j]
            y_plus.append(y_plus_local)
            u_plus_local = average_x_velocity[j]/nondimensional_friction_velocity
            u_plus.append(u_plus_local)
            # fluctuations
            u_plus_fluctuation_local = average_x_velocity_fluctuation_rms[j]/nondimensional_friction_velocity
            v_plus_fluctuation_local = average_y_velocity_fluctuation_rms[j]/nondimensional_friction_velocity
            w_plus_fluctuation_local = average_z_velocity_fluctuation_rms[j]/nondimensional_friction_velocity
            u_plus_fluctuation.append(u_plus_fluctuation_local)
            v_plus_fluctuation.append(v_plus_fluctuation_local)
            w_plus_fluctuation.append(w_plus_fluctuation_local)

        y_plus = np.array(y_plus)
        u_plus = np.array(u_plus)
        u_plus_fluctuation = np.array(u_plus_fluctuation)
        v_plus_fluctuation = np.array(v_plus_fluctuation)
        w_plus_fluctuation = np.array(w_plus_fluctuation)

        number_of_unique_BL_points = int(0.5*np.size(y_plus))+1
        index_centerline = number_of_unique_BL_points - 1
        unique_y_plus=np.zeros(number_of_unique_BL_points)
        unique_u_plus=np.zeros(number_of_unique_BL_points)
        unique_u_plus_fluctuation=np.zeros(number_of_unique_BL_points)
        unique_v_plus_fluctuation=np.zeros(number_of_unique_BL_points)
        unique_w_plus_fluctuation=np.zeros(number_of_unique_BL_points)
        unique_y_plus[index_centerline] = 1.0*y_plus[index_centerline]
        unique_u_plus[index_centerline] = 1.0*u_plus[index_centerline]
        unique_u_plus_fluctuation[index_centerline] = 1.0*u_plus_fluctuation[index_centerline]
        unique_v_plus_fluctuation[index_centerline] = 1.0*v_plus_fluctuation[index_centerline]
        unique_w_plus_fluctuation[index_centerline] = 1.0*w_plus_fluctuation[index_centerline]
        for j in range(0,index_centerline):
            unique_y_plus[j] = y_plus[j]
            index_for_second_BL = -(j+1)
            if(np.abs(y_plus[j]-y_plus[index_for_second_BL])>1.0e-4):
                print("ERROR in boundary layer profile averaging! Aborting...")
                print(y_plus[j])
                print(y_plus[index_for_second_BL])
                print(j)
                print(index_for_second_BL)
                exit()
            # average the two boundary layer profiles (because channel has two walls)
            unique_u_plus[j] = 0.5*(u_plus[j]+u_plus[index_for_second_BL])
            unique_u_plus_fluctuation[j] = 0.5*(u_plus_fluctuation[j]+u_plus_fluctuation[index_for_second_BL])
            unique_v_plus_fluctuation[j] = 0.5*(v_plus_fluctuation[j]+v_plus_fluctuation[index_for_second_BL])
            unique_w_plus_fluctuation[j] = 0.5*(w_plus_fluctuation[j]+w_plus_fluctuation[index_for_second_BL])

        # marker for the input
        y_plus_wall_model_input = 0.2*Re_tau # 0.2 is the uniform delta y from the grid
        u_plus_wall_model_input = np.average(u_plus[np.where(np.abs(y_plus - y_plus_wall_model_input)<1e-4)])

        # store the data
        y_plus_after_wall_model_input_point = unique_y_plus[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        u_plus_after_wall_model_input_point = unique_u_plus[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        u_plus_fluctuation_after_wall_model_input_point = unique_u_plus_fluctuation[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        v_plus_fluctuation_after_wall_model_input_point = unique_v_plus_fluctuation[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        w_plus_fluctuation_after_wall_model_input_point = unique_w_plus_fluctuation[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        
        # add to plot
        average_y_plus_store.append(y_plus_after_wall_model_input_point)
        labels_store.append(labels_[i])
        average_u_plus_store.append(u_plus_after_wall_model_input_point)

        average_velocity_fluctuation_rms_store.append(u_plus_fluctuation_after_wall_model_input_point)
        average_velocity_fluctuation_rms_store.append(v_plus_fluctuation_after_wall_model_input_point)
        average_velocity_fluctuation_rms_store.append(w_plus_fluctuation_after_wall_model_input_point)

        # add the marker to the plot
        average_y_plus_store.append(np.array(y_plus_wall_model_input))
        labels_store.append("Wall Model Input")
        average_u_plus_store.append(np.array(u_plus_wall_model_input))

    # load reference data
    y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/frere_p3_fig7a.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    average_y_plus_store.append(y_plus_reference)
    average_u_plus_store.append(average_u_plus_reference)
    labels_store.append("p3 DG-WMLES [Fr\\`ere, 2017]")

    # load reference data
    y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/frere_p4_fig7a.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    average_y_plus_store.append(y_plus_reference)
    average_u_plus_store.append(average_u_plus_reference)
    labels_store.append("p4 DG-WMLES [Fr\\`ere, 2017]")

    average_y_plus_store.append(unique_y_plus)
    labels_store.append("averaged")
    average_u_plus_store.append(unique_u_plus)

    qp.plotfxn(average_y_plus_store,average_u_plus_store,
        figure_filename="boundary_layer_profile",
        figure_size=(7,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx5200$",
        xlabel="$\\left\\langle y^{+}\\right\\rangle$",
        ylabel="$\\left\\langle u^{+}\\right\\rangle$",
        xlimits=[1.0e0,5200.0],
        ylimits=[0,30],
        which_lines_black=[0],
        which_lines_only_markers=[2,3,4],
        transparent_legend=False,
        legend_border_on=True,
        grid_lines_on=True,
        log_axes="x",
        legend_location="best",
        vertical_lines=[y_plus_wall_model_input])

    qp.plotfxn(average_y_plus_store,average_u_plus_store,
        figure_filename="boundary_layer_profile_zoom",
        figure_size=(7,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx5200$",
        xlabel="$\\left\\langle y^{+}\\right\\rangle$",
        ylabel="$\\left\\langle u^{+}\\right\\rangle$",
        xlimits=[1.0e2,5200.0],
        ylimits=[18,28],
        which_lines_black=[0],
        which_lines_only_markers=[2,3,4],
        which_lines_dashed=[5],
        transparent_legend=False,
        legend_border_on=True,
        grid_lines_on=True,
        log_axes="x",
        legend_location="best",
        vertical_lines=[y_plus_wall_model_input])

    # fluctuations
    xdata_for_fluctuations=[y_plus_after_wall_model_input_point,y_plus_after_wall_model_input_point,y_plus_after_wall_model_input_point]
    labels_store = ["$\\left\\langle u'^{+}\\right\\rangle _{rms}$",\
                    "$\\left\\langle v'^{+}\\right\\rangle _{rms}$",\
                    "$\\left\\langle w'^{+}\\right\\rangle _{rms}$"]
    qp.plotfxn(xdata_for_fluctuations,average_velocity_fluctuation_rms_store,
        figure_filename="boundary_layer_profile_of_fluctuations",
        figure_size=(7,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx5200$",
        xlabel="$\\left\\langle y^{+}\\right\\rangle$",
        ylabel="$\\left\\langle u'^{+}_{i}\\right\\rangle _{rms}$",
        xlimits=[1.0e0,5200.0],
        ylimits=[],
        which_lines_black=[],
        which_lines_only_markers=[],
        transparent_legend=False,
        legend_border_on=True,
        grid_lines_on=True,
        log_axes=None,
        legend_location="best",
        vertical_lines=[y_plus_wall_model_input])
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
# plot_transient(filenames,labels,starting_data_index_for_plot=0,friction_velocity_based_reynolds_number=friction_velocity_based_reynolds_number)

# plot boundary layer profile
filenames=[\
# filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0200.dat",\
# filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0240.dat",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0330.dat",\
]
labels=[\
"$c_{DG}$ NSFR.IR.GLL 20x10x10 p4"\
]
which_lines_dashed=[]
plot_boundary_layer_profile(filenames,labels)