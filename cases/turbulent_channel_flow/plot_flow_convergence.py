#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
# Import personal libraries
# from finite_difference_library import first_derivative, fd_non_uniform_grid
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
# from sys import platform
# if platform == "linux" or platform == "linux2":
#     # linux
#     sys.path.append("/home/julien/Codes/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes
# elif platform == "darwin":
#     # OS X
#     sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes
#-----------------------------------------------------
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
import matplotlib;from matplotlib.lines import Line2D
#=====================================================
# Helper functions
#=====================================================
#-----------------------------------------------------
def plotfxn(x_store,y_store,x_label,y_label,
    figure_filename,labels_store,which_lines_dashed_store,
    log_axes=None,xlimits=[],ylimits=[],which_lines_black=[]):
    qp.plotfxn(x_store,y_store,
        figure_filename=figure_filename,
        figure_size=(6,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        # title_label="WMLES: TCF $c_{DG}$ NSFR.IR.GLL 20x10x10 p4, $CFL\\approx0.2$",
        xlabel=x_label,
        ylabel=y_label,
        xlimits=xlimits,
        ylimits=ylimits,
        which_lines_black=which_lines_black,
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
            "$t^{*}$","Normalized Skin Friction Coefficient, $C_{f}(t^{*})/C^{expected}_{f}$","skin_friction_coefficient",\
            labels_store,which_lines_dashed_,
            # xlimits=[0,460]
            )
        # plotfxn(time_store,skin_friction_coefficient_store/expected_mean_value_for_skin_friction_coefficient,\
        #     "$t$","$C_{f}(t)$/$C_{f}^{expected}$","skin_friction_coefficient",\
        #     labels_store,which_lines_dashed_)
    if(plot_wall_shear_stress):
        plotfxn(time_store,wall_shear_stress_store,\
            "$t^{*}$","Normalized Nondim. Wall Shear Stress, $\\tau_{w}/\\tau^{expected}_{w}$","wall_shear_stress",\
            labels_store,which_lines_dashed_)
    if(plot_bulk_mass_flow):
        plotfxn(time_store,bulk_mass_flow_store,\
            "$t^{*}$","Normalized Bulk Mass Flow Rate, $\\rho_{b}U_{b}/(\\rho_{b}U_{b})_{0}$","bulk_mass_flow",\
            labels_store,which_lines_dashed_,log_axes="y")
    return
#-----------------------------------------------------
#-----------------------------------------------------
def plot_boundary_layer_profile(filenames_,labels_,friction_velocity_based_reynolds_number,which_lines_dashed_=[]):
    # expected_mean_value_for_skin_friction_coefficient = 6.25e-3 # from Lodato's source term paper
    # expected_mean_value_for_skin_friction_coefficient = 0.00347754 # for Re=5200
    # expected_mean_value_for_wall_shear_stress = 0.00161876 # for Re=5200
    # update above using Dean's expression
    # data store
    average_y_plus_store=[]
    average_u_plus_store=[]
    average_x_velocity_fluctuation_rms_store=[]
    average_y_velocity_fluctuation_rms_store=[]
    average_z_velocity_fluctuation_rms_store=[]
    average_reynolds_stress_uv_store=[]
    nondim_y_store = []
    average_u_by_ubulk_store = []
    # plot function inputs
    labels_store=[]

    if(friction_velocity_based_reynolds_number==5200):
        Re_tau = 5200.0
        Re_inf = 129245.0
    elif(friction_velocity_based_reynolds_number==395):
        Re_tau = 395.0
        Re_inf = 6793.46
    else:
        print("Error: Invalid friction_velocity_based_reynolds_number. Aborting...")
        exit()

    if(friction_velocity_based_reynolds_number==5200):
        # load reference data
        # y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/lee_moser_2015_dns.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",") # digitized
        # taken from https://turbulence.oden.utexas.edu/
        y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/LM_Channel_5200_mean_prof.dat",skiprows=72,usecols=(1,2),dtype=np.float64,unpack=True)
        average_y_plus_store.append(y_plus_reference)
        average_u_plus_store.append(average_u_plus_reference)
        nondim_y, uuplus, vvplus, wwplus, uvplus = np.loadtxt("./data/reference/LM_Channel_5200_vel_fluc_prof.dat",skiprows=75,usecols=(0,2,3,4,5),dtype=np.float64,unpack=True)
        nondim_y_store.append(nondim_y)
        average_x_velocity_fluctuation_rms_store.append(np.sqrt(uuplus))
        average_y_velocity_fluctuation_rms_store.append(np.sqrt(vvplus))
        average_z_velocity_fluctuation_rms_store.append(np.sqrt(np.abs(wwplus)))
        average_reynolds_stress_uv_store.append(-uvplus)
        labels_store.append("DNS [Lee \\& Moser, 2015]\n$(\\Delta x^{+},\\Delta y_{c}^{+},\\Delta y_{w}^{+},\\Delta z^{+})=(12.7,10.3,0.498,6.4)$")
        if(False):
            # load reference data
            y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/frere_p3_fig7a.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            average_y_plus_store.append(y_plus_reference)
            average_u_plus_store.append(average_u_plus_reference)
            labels_store.append("p3 DG-WMLES\n [Fr\\`ere et al., 2017]")

            # load reference data
            y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/frere_p4_fig7a.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            average_y_plus_store.append(y_plus_reference)
            average_u_plus_store.append(average_u_plus_reference)
            labels_store.append("p4 DG-WMLES\n [Fr\\`ere et al., 2017]")

    elif(friction_velocity_based_reynolds_number==395):
        # # load reference data
        # y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/moser_et_al_1999_dns_Re395.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        # average_y_plus_store.append(y_plus_reference)
        # average_u_plus_store.append(average_u_plus_reference)
        # labels_store.append("DNS [Moser et al., 1999]")
        # load Kim data; available at http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case032
        if(False):
            nondim_y, y_plus_reference, average_u_plus_reference, uuplus, vvplus, wwplus, uvplus = np.loadtxt("./data/reference/kim_1989_data_Re395.txt",usecols=(1,2,3,4,5,6,7),skiprows=11,max_rows=97,dtype=np.float64,unpack=True)
            average_y_plus_store.append(y_plus_reference)
            average_u_plus_store.append(average_u_plus_reference)
            nondim_y_store.append(nondim_y)
            average_x_velocity_fluctuation_rms_store.append(np.sqrt(uuplus))
            average_y_velocity_fluctuation_rms_store.append(np.sqrt(vvplus))
            average_z_velocity_fluctuation_rms_store.append(np.sqrt(wwplus))
            average_reynolds_stress_uv_store.append(-uvplus)
            labels_store.append("DNS [Lee, 1989]")
        # MOSER KIM MANSOUR DATA
        nondim_y, y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/moser_1999_chan395_means.txt",usecols=(0,1,2),skiprows=25,dtype=np.float64,unpack=True)
        uuplus, vvplus, wwplus, uvplus = np.loadtxt("./data/reference/moser_1999_chan395_reystress.txt",usecols=(2,3,4,5),skiprows=25,dtype=np.float64,unpack=True)
        average_y_plus_store.append(y_plus_reference)
        average_u_plus_store.append(average_u_plus_reference)
        nondim_y_store.append(nondim_y)
        average_x_velocity_fluctuation_rms_store.append(np.sqrt(uuplus))
        average_y_velocity_fluctuation_rms_store.append(np.sqrt(vvplus))
        average_z_velocity_fluctuation_rms_store.append(np.sqrt(wwplus))
        average_reynolds_stress_uv_store.append(-uvplus)
        labels_store.append("DNS [Moser et al., 1999]")
        average_u_by_ubulk_store.append(average_u_plus_reference*(Re_tau/Re_inf))

    for i,filename in enumerate(filenames_):
        # load data
        y_coordinate, average_x_velocity, average_kinematic_viscosity, average_x_velocity_fluctuation_rms, average_y_velocity_fluctuation_rms, average_z_velocity_fluctuation_rms, average_reynolds_stress_uv = np.loadtxt(filename,dtype=np.float64,unpack=True)

        number_of_points = np.size(y_coordinate)
        y_plus = []
        u_plus = []
        u_plus_fluctuation = []
        v_plus_fluctuation = []
        w_plus_fluctuation = []
        reynolds_stress_uv = []
        for j in range(0,number_of_points):
            y_value = y_coordinate[j] # nondimensional
            distance_from_wall = 0.0 # nondimensional
            if(y_value < 0.0):
                distance_from_wall = 1.0 + y_value
            else:
                distance_from_wall = 1.0 - y_value
            if(friction_velocity_based_reynolds_number==5200):
                Re_tau = 5200.0
                Re_inf = 129245.0
            elif(friction_velocity_based_reynolds_number==395):
                Re_tau = 395.0
                Re_inf = 6793.46
            else:
                print("Error: Invalid friction_velocity_based_reynolds_number. Aborting...")
                exit()
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
            # Reynolds stresses
            reynolds_stress_uv_local = np.abs(average_reynolds_stress_uv[j])/nondimensional_friction_velocity/nondimensional_friction_velocity
            reynolds_stress_uv.append(reynolds_stress_uv_local)

        y_plus = np.array(y_plus)
        u_plus = np.array(u_plus)
        u_plus_fluctuation = np.array(u_plus_fluctuation)
        v_plus_fluctuation = np.array(v_plus_fluctuation)
        w_plus_fluctuation = np.array(w_plus_fluctuation)
        reynolds_stress_uv = np.array(reynolds_stress_uv)
        number_of_unique_BL_points = int(0.5*np.size(y_plus))+1
        index_centerline = number_of_unique_BL_points - 1
        unique_y_plus=np.zeros(number_of_unique_BL_points)
        unique_u_plus=np.zeros(number_of_unique_BL_points)
        unique_u_plus_fluctuation=np.zeros(number_of_unique_BL_points)
        unique_v_plus_fluctuation=np.zeros(number_of_unique_BL_points)
        unique_w_plus_fluctuation=np.zeros(number_of_unique_BL_points)
        unique_reynolds_stress_uv=np.zeros(number_of_unique_BL_points)
        unique_y_plus[index_centerline] = 1.0*y_plus[index_centerline]
        unique_u_plus[index_centerline] = 1.0*u_plus[index_centerline]
        unique_u_plus_fluctuation[index_centerline] = 1.0*u_plus_fluctuation[index_centerline]
        unique_v_plus_fluctuation[index_centerline] = 1.0*v_plus_fluctuation[index_centerline]
        unique_w_plus_fluctuation[index_centerline] = 1.0*w_plus_fluctuation[index_centerline]
        unique_reynolds_stress_uv[index_centerline] = 1.0*reynolds_stress_uv[index_centerline]
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
            unique_reynolds_stress_uv[j] = 0.5*(reynolds_stress_uv[j]+reynolds_stress_uv[index_for_second_BL])

        # marker for the input
        y_plus_wall_model_input = 0.2*Re_tau # 0.2 is the uniform delta y from the grid # WARNING -- SHOULD BE AN INPUT -- TO DO
        u_plus_wall_model_input = np.average(u_plus[np.where(np.abs(y_plus - y_plus_wall_model_input)<1e-4)])
        index_wall_model_input = np.abs(unique_y_plus - y_plus_wall_model_input).argmin()

        # store the data
        y_plus_after_wall_model_input_point = unique_y_plus[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        u_plus_after_wall_model_input_point = unique_u_plus[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        u_plus_fluctuation_after_wall_model_input_point = unique_u_plus_fluctuation[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        v_plus_fluctuation_after_wall_model_input_point = unique_v_plus_fluctuation[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        w_plus_fluctuation_after_wall_model_input_point = unique_w_plus_fluctuation[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]
        reynolds_stress_uv_after_wall_model_input_point = unique_reynolds_stress_uv[np.where(unique_y_plus >= (y_plus_wall_model_input-1e-4))]

        # vbulk_store=[9.8272033399675807e-01,1.0242263253573591e+00] # to account for imperfect bulk velocity due to source term
        # vbulk_store=[1.0,1.0]
        instantaneous_nondim_bulk_velocity = 1.0 #vbulk_store[i]
        instantaneous_nondim_bulk_velocity_sqr = instantaneous_nondim_bulk_velocity*instantaneous_nondim_bulk_velocity

        # add to plot
        average_y_plus_store.append(y_plus_after_wall_model_input_point)
        labels_store.append(labels_[i])
        average_u_plus_store.append(u_plus_after_wall_model_input_point/instantaneous_nondim_bulk_velocity)

        average_x_velocity_fluctuation_rms_store.append(u_plus_fluctuation_after_wall_model_input_point/instantaneous_nondim_bulk_velocity)
        average_y_velocity_fluctuation_rms_store.append(v_plus_fluctuation_after_wall_model_input_point/instantaneous_nondim_bulk_velocity)
        average_z_velocity_fluctuation_rms_store.append(w_plus_fluctuation_after_wall_model_input_point/instantaneous_nondim_bulk_velocity)

        average_reynolds_stress_uv_store.append(reynolds_stress_uv_after_wall_model_input_point/instantaneous_nondim_bulk_velocity_sqr)
        nondim_y_store.append(y_plus_after_wall_model_input_point/np.float64(friction_velocity_based_reynolds_number))
        # add the marker to the plot for wall model input location -- can uncomment to add the marker
        # average_y_plus_store.append(np.array(y_plus_wall_model_input))
        # labels_store.append("Wall Model Input")
        # average_u_plus_store.append(np.array(u_plus_wall_model_input))

        if(friction_velocity_based_reynolds_number==5200):
            Re_tau = 5200.0
            Re_inf = 129245.0
        elif(friction_velocity_based_reynolds_number==395):
            Re_tau = 395.0
            Re_inf = 6793.46
        else:
            print("Error: Invalid friction_velocity_based_reynolds_number. Aborting...")
            exit()
        nondimensional_friction_velocity = Re_tau/Re_inf
        average_u_by_ubulk_store.append(u_plus_after_wall_model_input_point*np.float64(nondimensional_friction_velocity)/instantaneous_nondim_bulk_velocity)

    if(friction_velocity_based_reynolds_number==5200):
        xlimits_=[1.0e0,5200.0]
        ylimits_=[0,30]
        xlimits_zoom=[1e2,5200.0]
        ylimits_zoom=[18,28]
        which_lines_black_=[0]
        # which_lines_only_markers_=[1,2] # uncomment if plotting Frere's results too
        which_lines_only_markers_=[]
        which_lines_dashed_=[5]
        y_plus_val_for_reynolds_stress_textboxes=1100.0
    elif(friction_velocity_based_reynolds_number==395):
        xlimits_=[]
        ylimits_=[]
        xlimits_zoom=[5e1,5e2]
        ylimits_zoom=[12.5,22.5]
        which_lines_black_=[0]#[0,1]
        which_lines_only_markers_=[]
        which_lines_dashed_=[]#[1]
        y_plus_val_for_reynolds_stress_textboxes=180.0

    if(friction_velocity_based_reynolds_number==395):
        qp.plotfxn(nondim_y_store,average_u_by_ubulk_store,
            figure_filename="boundary_layer_profile_scaled_by_bulk_velocity_Re"+str(friction_velocity_based_reynolds_number),
            figure_size=(7,6),
            legend_labels_tex=labels_store,
            figure_filetype="pdf",
            # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
            # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx5200$",
            # xlabel="$\\left\\langle y^{+}\\right\\rangle$",
            xlabel="$y/\\delta$",
            ylabel="$\\left\\langle u/U_{b}\\right\\rangle$",
            xlimits=xlimits_,
            ylimits=ylimits_,
            # xlimits=xlimits_zoom,
            # ylimits=ylimits_zoom,
            which_lines_black=which_lines_black_,
            which_lines_only_markers=which_lines_only_markers_,
            which_lines_dashed=which_lines_dashed_,
            transparent_legend=True,
            legend_border_on=False,
            grid_lines_on=False,
            # log_axes="x",
            legend_location="best",
            vertical_lines=[y_plus_wall_model_input/np.float64(friction_velocity_based_reynolds_number)])

    if(friction_velocity_based_reynolds_number==395):
        y_plus_reference, average_u_plus_reference = np.loadtxt("./data/reference/vermeire2016cpr/p4_fine_uplus_vs_yplus.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        average_y_plus_store.insert(1,y_plus_reference)
        average_u_plus_store.insert(1,average_u_plus_reference)
        # labels_store.insert(1,"p4 fine CPR\n[Vermeire et al., 2016]")
        labels_store.insert(1,"p$4$ WR-ILES via CPR\n[Vermeire et al., 2016]")
        which_lines_black_=[0,1]
        which_lines_only_markers_=[1]


    qp.plotfxn(average_y_plus_store,average_u_plus_store,
        figure_filename="boundary_layer_profile_Re"+str(friction_velocity_based_reynolds_number),
        figure_size=(7,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx5200$",
        # xlabel="$\\left\\langle y^{+}\\right\\rangle$",
        xlabel="$y^{+}$",
        ylabel="$\\left\\langle u^{+}\\right\\rangle$",
        xlimits=xlimits_,
        ylimits=ylimits_,
        which_lines_black=which_lines_black_,
        which_lines_only_markers=which_lines_only_markers_,
        which_lines_dashed=which_lines_dashed_,
        transparent_legend=True,
        legend_border_on=False,
        grid_lines_on=False,
        log_axes="x",
        legend_location="upper left",
        vertical_lines=[y_plus_wall_model_input])

    # uncomment for unresolved section
    # average_y_plus_store.append(unique_y_plus[:(index_wall_model_input+1)])
    # labels_store.append("Unresolved")
    # average_u_plus_store.append(unique_u_plus[:(index_wall_model_input+1)])

    qp.plotfxn(average_y_plus_store,average_u_plus_store,
        figure_filename="boundary_layer_profile_zoom_Re"+str(friction_velocity_based_reynolds_number),
        figure_size=(7,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx5200$",
        # xlabel="$\\left\\langle y^{+}\\right\\rangle$",
        xlabel="$y^{+}$",
        ylabel="$\\left\\langle u^{+}\\right\\rangle$",
        xlimits=xlimits_zoom,
        ylimits=ylimits_zoom,
        which_lines_black=which_lines_black_,
        which_lines_only_markers=which_lines_only_markers_,
        which_lines_dashed=which_lines_dashed_,
        transparent_legend=True,
        legend_border_on=False,
        grid_lines_on=False,
        log_axes="x",
        legend_location="upper left",
        vertical_lines=[y_plus_wall_model_input])

    if(friction_velocity_based_reynolds_number==395):
        # remove vermeire data
        average_y_plus_store.pop(1)
        average_u_plus_store.pop(1)
        labels_store.pop(1)
        which_lines_black_=[0]
        which_lines_only_markers_=[]

    if(True):
        plot_vs_nondim_y=True # FLAG

        xlabel_="$\\left\\langle y^{+}\\right\\rangle$"
        xdata_=average_y_plus_store
        if(plot_vs_nondim_y==True):
            xlabel_="$y/\\delta$"
            xdata_=nondim_y_store

        if(friction_velocity_based_reynolds_number==395 and plot_vs_nondim_y==True):
            nondim_y, value = np.loadtxt("./data/reference/vermeire2016cpr/p4_fine_uvrms_vs_y.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            nondim_y_store.insert(1,nondim_y)
            average_reynolds_stress_uv_store.insert(1,value)
            # labels_store.insert(1,"p$4$ WR-ILES CPR\n[Vermeire et al., 2016]")
            labels_store.insert(1,"p$4$ WR-ILES via CPR\n[Vermeire et al., 2016]")
            which_lines_black_=[0,1]
            which_lines_only_markers_=[1]
            xdata_=nondim_y_store # as a precaution; to avoid relying on deep copy
        qp.plotfxn(
            xdata_,
            average_reynolds_stress_uv_store,
            figure_filename="reynolds_stress_uv_profile_Re"+str(friction_velocity_based_reynolds_number),
            figure_size=(7,6),
            legend_labels_tex=labels_store,
            figure_filetype="pdf",
            # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
            # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx5200$",
            xlabel=xlabel_,
            ylabel="$-u'v'/u_{\\tau}^{2}$",
            # xlimits=[],
            # ylimits=ylimits_,
            which_lines_black=which_lines_black_,
            which_lines_only_markers=which_lines_only_markers_,
            which_lines_dashed=which_lines_dashed_,
            transparent_legend=True,
            legend_border_on=False,
            grid_lines_on=False,
            # log_axes="x",
            legend_location="best",
            vertical_lines=[y_plus_wall_model_input/np.float64(friction_velocity_based_reynolds_number)],
            # secondary_vertical_lines=[0.3]
            )
        
        if(friction_velocity_based_reynolds_number==395 and plot_vs_nondim_y==True):
            nondim_y, value = np.loadtxt("./data/reference/vermeire2016cpr/p4_fine_urms_vs_y.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            nondim_y_store[1]=nondim_y
            average_x_velocity_fluctuation_rms_store.insert(1,value)
            xdata_=nondim_y_store # as a precaution; to avoid relying on deep copy

        # fluctuations
        # xdata_for_fluctuations=[y_plus_after_wall_model_input_point,y_plus_after_wall_model_input_point,y_plus_after_wall_model_input_point]
        # labels_store = ["$\\left\\langle u'^{+}\\right\\rangle _{rms}$",\
        #                 "$\\left\\langle v'^{+}\\right\\rangle _{rms}$",\
        #                 "$\\left\\langle w'^{+}\\right\\rangle _{rms}$"]
        qp.plotfxn(
            xdata_,
            average_x_velocity_fluctuation_rms_store,
            figure_filename="boundary_layer_profile_of_x_velocity_fluctuations_Re"+str(friction_velocity_based_reynolds_number),
            figure_size=(7,6),
            legend_labels_tex=labels_store,
            figure_filetype="pdf",
            # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
            # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx395$",
            # xlabel="$\\left\\langle y^{+}\\right\\rangle$",
            xlabel=xlabel_,
            ylabel="$\\left\\langle u'^{+}\\right\\rangle _{rms}$",
            # xlimits=xlimits_,
            # ylimits=[],
            which_lines_black=which_lines_black_,
            which_lines_dashed=which_lines_dashed_,
            which_lines_only_markers=which_lines_only_markers_,
            transparent_legend=True,
            legend_border_on=False,
            grid_lines_on=False,
            # log_axes="x",
            legend_location="best",
            vertical_lines=[y_plus_wall_model_input/np.float64(friction_velocity_based_reynolds_number)],
            # secondary_vertical_lines=[0.3]
            )
        
        if(friction_velocity_based_reynolds_number==395 and plot_vs_nondim_y==True):
            nondim_y, value = np.loadtxt("./data/reference/vermeire2016cpr/p4_fine_vrms_vs_y.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            nondim_y_store[1]=nondim_y
            average_y_velocity_fluctuation_rms_store.insert(1,value)
            xdata_=nondim_y_store # as a precaution; to avoid relying on deep copy
        qp.plotfxn(
            xdata_,
            average_y_velocity_fluctuation_rms_store,
            figure_filename="boundary_layer_profile_of_y_velocity_fluctuations_Re"+str(friction_velocity_based_reynolds_number),
            figure_size=(7,6),
            legend_labels_tex=labels_store,
            figure_filetype="pdf",
            # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
            # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx395$",
            xlabel=xlabel_,
            ylabel="$\\left\\langle v'^{+}\\right\\rangle _{rms}$",
            # xlimits=xlimits_,
            # ylimits=[],
            which_lines_black=which_lines_black_,
            which_lines_dashed=which_lines_dashed_,
            which_lines_only_markers=which_lines_only_markers_,
            transparent_legend=True,
            legend_border_on=False,
            grid_lines_on=False,
            # log_axes="x",
            legend_location="best",
            vertical_lines=[y_plus_wall_model_input/np.float64(friction_velocity_based_reynolds_number)],
            # secondary_vertical_lines=[0.3]
            )
        
        if(friction_velocity_based_reynolds_number==395 and plot_vs_nondim_y==True):
            nondim_y, value = np.loadtxt("./data/reference/vermeire2016cpr/p4_fine_wrms_vs_y.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            nondim_y_store[1]=nondim_y
            average_z_velocity_fluctuation_rms_store.insert(1,value)
            xdata_=nondim_y_store # as a precaution; to avoid relying on deep copy
        qp.plotfxn(
            xdata_,
            average_z_velocity_fluctuation_rms_store,
            figure_filename="boundary_layer_profile_of_z_velocity_fluctuations_Re"+str(friction_velocity_based_reynolds_number),
            figure_size=(7,6),
            legend_labels_tex=labels_store,
            figure_filetype="pdf",
            # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
            # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx395$",
            xlabel=xlabel_,
            ylabel="$\\left\\langle w'^{+}\\right\\rangle _{rms}$",
            # xlimits=xlimits_,
            # ylimits=[],
            which_lines_black=which_lines_black_,
            which_lines_only_markers=which_lines_only_markers_,
            which_lines_dashed=which_lines_dashed_,
            transparent_legend=True,
            legend_border_on=False,
            grid_lines_on=False,
            # log_axes="x",
            legend_location="best",
            vertical_lines=[y_plus_wall_model_input/np.float64(friction_velocity_based_reynolds_number)],
            # secondary_vertical_lines=[0.3]
            )

        if(friction_velocity_based_reynolds_number==395 and plot_vs_nondim_y==True):
            # remove vermeire data
            nondim_y_store.pop(1)
            average_reynolds_stress_uv_store.pop(1)
            average_x_velocity_fluctuation_rms_store.pop(1)
            average_y_velocity_fluctuation_rms_store.pop(1)
            average_z_velocity_fluctuation_rms_store.pop(1)
            labels_store.pop(1)
            which_lines_black_=[0]
            which_lines_only_markers_=[]

    # plot with all profiles in one
    number_of_profiles = 4
    number_of_runs = len(filenames_)
    xdata_for_plot = []
    ydata_for_plot = []
    clr_reference = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    mrkr_reference=['None','None','None','None','None','None','None','None','None']
    lnstl_reference=['solid','solid','solid','dashed','dashdot','solid','dashed']
    clr_input_store = []
    mrkr_input_store = []
    lnstl_input_store = []
    leg_elements_input = []
    # labels_store_all_profiles = []
    # load in the reference results
    leg_elements_input.append(Line2D([0],[0], label=labels_store[0], color='k', marker='None', markersize=6, mfc='None', linestyle='solid'))
    for i in range(0,number_of_profiles):
        clr_input_store.append('k')
        mrkr_input_store.append('None')
        lnstl_input_store.append('solid')
        xdata_for_plot.append(average_y_plus_store[0])
    # ydata_for_plot.append(average_u_plus_store[0])
    ydata_for_plot.append(average_reynolds_stress_uv_store[0])
    ydata_for_plot.append(average_x_velocity_fluctuation_rms_store[0])
    ydata_for_plot.append(average_y_velocity_fluctuation_rms_store[0])
    ydata_for_plot.append(average_z_velocity_fluctuation_rms_store[0])

    for j in range(number_of_runs):
        leg_elements_input.append(Line2D([0],[0], label=labels_store[j+1], color=clr_reference[j], marker=mrkr_reference[j], markersize=6, mfc='None', linestyle=lnstl_reference[j]))
        for i in range(0,number_of_profiles):
            clr_input_store.append(clr_reference[j])
            mrkr_input_store.append(mrkr_reference[j])
            lnstl_input_store.append(lnstl_reference[j])
            xdata_for_plot.append(average_y_plus_store[j+1])
        # ydata_for_plot.append(average_u_plus_store[j+1])
        ydata_for_plot.append(average_reynolds_stress_uv_store[j+1])
        ydata_for_plot.append(average_x_velocity_fluctuation_rms_store[j+1])
        ydata_for_plot.append(average_y_velocity_fluctuation_rms_store[j+1])
        ydata_for_plot.append(average_z_velocity_fluctuation_rms_store[j+1])


    # index_loc_for_text = np.argmin(np.abs(average_y_plus_store[0]-y_plus_wall_model_input))

    index_loc_for_text = np.argmin(np.abs(average_y_plus_store[0]-y_plus_val_for_reynolds_stress_textboxes))
    x_loc_for_text = []
    y_loc_for_text = []
    strings_for_textbox = []
    for i in range(0,number_of_profiles):
        x_loc_for_text.append(average_y_plus_store[0][index_loc_for_text])
    y_loc_for_text.append(average_reynolds_stress_uv_store[0][index_loc_for_text])
    y_loc_for_text.append(average_x_velocity_fluctuation_rms_store[0][index_loc_for_text])
    y_loc_for_text.append(average_y_velocity_fluctuation_rms_store[0][index_loc_for_text])
    y_loc_for_text.append(average_z_velocity_fluctuation_rms_store[0][index_loc_for_text])
    strings_for_textbox.append("$-\\left\\langle u'v'\\right\\rangle^{+}$")
    strings_for_textbox.append("$\\left\\langle u'\\right\\rangle^{+}_{rms}$")
    strings_for_textbox.append("$\\left\\langle v'\\right\\rangle^{+}_{rms}$")
    strings_for_textbox.append("$\\left\\langle w'\\right\\rangle^{+}_{rms}$")
    horizontal_alignment_for_textbox_input=['right','left','right','left']
    vertical_alignment_for_textbox_input=['top','bottom','top','bottom']

    qp.plotfxn(xdata_for_plot,ydata_for_plot,
        figure_filename="reynolds_stress_profiles_Re"+str(friction_velocity_based_reynolds_number),
        figure_size=(7,6),
        # legend_labels_tex=labels_store,
        figure_filetype="pdf",
        # title_label="Turbulent Channel Flow $Re_{\\tau}\\approx395$, $CFL\\approx0.2$, $\\alpha=0.0$",
        # title_label="WMLES Approach to Turbulent Channel Flow at $Re_{\\tau}\\approx395$",
        xlabel="$\\left\\langle y^{+}\\right\\rangle$",
        ylabel="Reynolds Stress Profiles",
        xlimits=xlimits_,
        # ylimits=[0.0,3.0],
        clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
        leg_elements_input=leg_elements_input,
        transparent_legend=True,
        legend_border_on=False,
        grid_lines_on=False,
        log_axes=None,
        legend_location="best",
        # vertical_lines=[y_plus_wall_model_input/np.float64(friction_velocity_based_reynolds_number)],
        vertical_lines=[y_plus_wall_model_input],
        # secondary_vertical_lines=[0.3],
        strings_for_textbox=strings_for_textbox,
        x_positions_for_textbox=x_loc_for_text,
        y_positions_for_textbox=y_loc_for_text,
        horizontal_alignment_for_textbox=horizontal_alignment_for_textbox_input,
        vertical_alignment_for_textbox=vertical_alignment_for_textbox_input)


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
# filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/turbulent_quantities.txt",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re395_p4_20x10x10_turbulent_initialization/turbulent_quantities.txt",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re395_p4_20x10x10_turbulent_initialization_second_element_input/turbulent_quantities.txt",\
]
labels=[\
# "p4 $c_{DG}$ NSFR.IR.GLL-WMLES\n $(\\Delta x^{+},\\Delta y^{+},\\Delta z^{+})=(400,250,400)$",\
# "$Re_{\\tau}\\approx5200$",\
"$Re_{\\tau}\\approx395$",\
"$Re_{\\tau}\\approx395$ 2nd el. input",\
]
which_lines_dashed=[1]
# friction_velocity_based_reynolds_number=[5200,395]
friction_velocity_based_reynolds_number=[395,395]
# plot_transient(filenames,labels,which_lines_dashed_=which_lines_dashed,starting_data_index_for_plot=0,friction_velocity_based_reynolds_number=friction_velocity_based_reynolds_number)
# exit()

filenames=[\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re395_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0300.dat",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re395_p4_20x10x10_turbulent_initialization_second_element_input/flow_field_files/velocity_vorticity-0_boundary_layer_profile.dat",\
# filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re395_p4_20x10x10_turbulent_initialization_second_element_input/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t670.dat",\
]
labels=[\
# "P$4$ NSFR-WMLES\n $n_{x}\\times n_{y}\\times n_{z} = 20\\times 10\\times 10$\n$(t^{*}=300)$",\
# "P$4$ NSFR-WMLES 2nd el. input\n $n_{x}\\times n_{y}\\times n_{z} = 20\\times 10\\times 10$\n$(t^{*}=400)$",\
# "P$4$ NSFR-WMLES 2nd el. input\n $n_{x}\\times n_{y}\\times n_{z} = 20\\times 10\\times 10$\n$(t^{*}=670)$",\
"p$4$ WM-ILES-1st-element",\
"p$4$ WM-ILES-2nd-element",\
# "3",\
]
which_lines_dashed=[2]
plot_boundary_layer_profile(filenames,labels,395,which_lines_dashed_=which_lines_dashed)
exit()
# plot boundary layer profile
filenames=[\
# filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0200.dat",\
# filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0240.dat",\
# filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0330.dat",\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0450.dat",\
]
labels=[\
# "p4 $c_{DG}$ NSFR.IR.GLL-WMLES\n $(\\Delta x^{+},\\Delta y^{+},\\Delta z^{+})=(400,250,400)$",\
# "$t^{*}=200$",\
# "$t^{*}=240$",\
# "$t^{*}=330$",\
"$t^{*}=450$",\
]
which_lines_dashed=[2]
plot_boundary_layer_profile(filenames,labels,5200,which_lines_dashed_=which_lines_dashed)
exit()
filenames=[\
filesystem+"NarvalFiles/2024_AIAA/turbulent_channel_flow/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re395_p4_20x10x10_turbulent_initialization/flow_field_files/velocity_vorticity-0_boundary_layer_profile_t0300.dat",\
]
labels=[\
"$t^{*}=300$",\
]
which_lines_dashed=[]
plot_boundary_layer_profile(filenames,labels,395)