import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../src");
import numpy as np
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
#-----------------------------------------------------
# define functions
#-----------------------------------------------------
from scipy.interpolate import splrep, splev
def get_dissipation_discrete(time,kinetic_energy,smoothing=False):
    if(smoothing):
        # smoothing if noisy
        dissipation_rate_val = splrep(time,kinetic_energy,k=5,s=0.0000001)
        dissipation_rate_val = -splev(time,dissipation_rate_val,der=1)
        return dissipation_rate_val
    else:
        dissipation_rate_val = -np.gradient(kinetic_energy,time)
        return dissipation_rate_val
def smooth_data_np_cumsum_my_average(arr, span):
    cumsum_vec = np.cumsum(arr)
    moving_average = (cumsum_vec[2 * span:] - cumsum_vec[:-2 * span]) / (2 * span)

    # The "my_average" part again. Slightly different to before, because the
    # moving average from cumsum is shorter than the input and needs to be padded
    front, back = [np.average(arr[:span])], []
    for i in range(1, span):
        front.append(np.average(arr[:i + span]))
        back.insert(0, np.average(arr[-i - span:]))
    back.insert(0, np.average(arr[-2 * span:]))
    return np.concatenate((front, moving_average, back))
#=====================================================
def get_smoothing_parameters_from_subdirectories(subdirectories_for_plot):
    smoothing_parameters_store=[]
    number_of_result_curves=len(subdirectories_for_plot)
    for i in range(0,number_of_result_curves):
        name=subdirectories_for_plot[i]
        if(name=="supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512"):
            smoothing_parameter=450
        elif(name=="supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p3_procs512"):
            smoothing_parameter=450 # copied
        elif(name=="supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512"):
            smoothing_parameter=60
        elif(name=="supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512"):
            smoothing_parameter=60
        elif(name=="supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128"):
            smoothing_parameter=15 # copied
        elif(name=="supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128"):
            smoothing_parameter=15 # copied
        elif(name=="supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128"):
            smoothing_parameter=500 # copied
        elif("dofs0256" in name):
            smoothing_parameter=450 # copied
        elif("dofs0128" in name):
            smoothing_parameter=60 # copied
        elif("dofs0064" in name):
            smoothing_parameter=15 # copied
        else:
            smoothing_parameter=1 #dummy
        smoothing_parameters_store.append(smoothing_parameter)
    return smoothing_parameters_store
#=====================================================
# Input variables for plotting
#=====================================================
global subdirectories, filenames, labels, black_line_flag, \
dashed_line_flag, figure_filename_postfix, figure_title, \
ylimits_kinetic_energy_input, ylimits_dissipation_input, \
log_axes_input, legend_on_input, legend_inside_input, \
plot_reference_result, nlegendcols_input, \
figure_subdirectory, data_directory_base, figure_directory_base
#=====================================================
def plot_for_presentation(
    subdirectories_for_plot,
    labels_for_plot,
    black_line_flag_for_plot,
    dashed_line_flag_for_plot,
    plotting_subsonic_result=False,
    which_was_ran_with_corrected_quantites=[],
    number_of_degrees_of_freedom=[],
    compare_with_reference_result_at_same_degrees_of_freedom=False,
    lnstl_input_store_=[],
    smooth_dilatational_dissipation_rate=[],
    smoothing_parameters_input=[]):
    
    global subdirectories, filenames, labels, black_line_flag, \
    dashed_line_flag, figure_filename_postfix, figure_title, \
    ylimits_kinetic_energy_input, ylimits_dissipation_input, \
    log_axes_input, legend_on_input, legend_inside_input, \
    plot_reference_result, nlegendcols_input, \
    figure_subdirectory, data_directory_base, figure_directory_base
    #-----------------------------------------------------
    # data store
    time_store = []
    kinetic_energy_store = []
    solenoidal_dissipation_store = []
    dilatational_dissipation_store = []
    pressure_dissipation_store = []
    dissipation_store = []
    #-----------------------------------------------------
    if(plotting_subsonic_result):
        time, kinetic_energy = np.loadtxt("./data/chapelier2024/subsonic/kinetic_energy.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    else:
        time, kinetic_energy = np.loadtxt("./data/chapelier2024/kinetic_energy.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    time_store.append(time)
    kinetic_energy_store.append(kinetic_energy)
    dissipation = get_dissipation_discrete(time,kinetic_energy)
    dissipation_store.append(dissipation)
    if(plotting_subsonic_result):
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/subsonic/solenoidal_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    else:
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    solenoidal_dissipation_store.append(solenoidal_dissipation)
    if(plotting_subsonic_result):
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/subsonic/dilatational_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    else:
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    dilatational_dissipation_store.append(dilatational_dissipation)
    pressure_dissipation_store.append(np.nan*dilatational_dissipation)
    if(plotting_subsonic_result):
        labels.append("SPADE $512^3$ DOF\n[Chapelier et al.]")
    else:
        labels.append("$2048^3$ DOF\n[Chapelier et al.]")
    black_line_flag.append(True)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    clr_input_store = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    mrkr_input_store = ['o','None','None','None','None','None','None','None']
    lnstl_input_store = ['None','solid','solid','solid','solid','dashed','solid','dashed','solid']
    if(plotting_subsonic_result):
        lnstl_input_store = ['None','solid','dashed','solid','solid','dashed','solid','dashed','solid']
    if(256 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        clr_input_store = ['k','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        # clr_input_store = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        lnstl_input_store = ['None','solid','solid','solid','dashed','solid','solid','solid']
    if(128 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        # clr_input_store = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        clr_input_store = ['k','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        lnstl_input_store = ['None','solid','solid','solid','solid','solid','solid','dashed','solid']
    if(64 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        # clr_input_store = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        clr_input_store = ['k','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        lnstl_input_store = ['None','solid','solid','solid','solid','solid','solid','dashed','solid']
    #-----------------------------------------------------
    number_of_result_curves=len(subdirectories_for_plot)
    for i in range(0,number_of_result_curves):
        figure_filename_postfix_input=figure_filename_postfix
        #-----------------------------------------------------
        subdirectories.append(subdirectories_for_plot[i])
        labels.append(labels_for_plot[i])
        black_line_flag.append(black_line_flag_for_plot[i])
        dashed_line_flag.append(dashed_line_flag_for_plot[i])
        filenames.append("turbulent_quantities.txt")
        #-----------------------------------------------------
        # load file
        if("/Volumes/KAUST/" in subdirectories[i]):
            filename = subdirectories[i]+"/"+filenames[i]
        else:
            filename = data_directory_base+"/"+subdirectories[i]+"/"+filenames[i]
        if(which_was_ran_with_corrected_quantites!=[] and i in which_was_ran_with_corrected_quantites):
            # time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation, corrected_pressure_dilatation_based_dissipation, corrected_dilatational_dissipation, uncorrected_pressure_dilatation_based_dissipation, uncorrected_dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
            time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,usecols=(0,1,2,3,4,5,6,7,8),dtype=np.float64,unpack=True)
            pressure_dissipation_store.append(pressure_dilatation_based_dissipation)
            if(smooth_dilatational_dissipation_rate!=[] and smooth_dilatational_dissipation_rate[i]==True):
                dilatational_dissipation = smooth_data_np_cumsum_my_average(dilatational_dissipation,smoothing_parameters_input[i])
            dilatational_dissipation_store.append(dilatational_dissipation)
            # pressure_dissipation_store.append(corrected_pressure_dilatation_based_dissipation)
            # dilatational_dissipation_store.append(uncorrected_dilatational_dissipation)
            # print(dilatational_dissipation[-1])
            # print(uncorrected_dilatational_dissipation[-1])
            # print(corrected_dilatational_dissipation[-1])            
        else:
            time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,usecols=(0,1,2,3,4,5,6,7,8),dtype=np.float64,unpack=True)
            pressure_dissipation_store.append(pressure_dilatation_based_dissipation)
            if(smooth_dilatational_dissipation_rate!=[] and smooth_dilatational_dissipation_rate[i]==True):
                dilatational_dissipation = smooth_data_np_cumsum_my_average(dilatational_dissipation,smoothing_parameters_input[i])
            dilatational_dissipation_store.append(dilatational_dissipation)
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        dissipation = get_dissipation_discrete(time,kinetic_energy)
        dissipation_store.append(dissipation)
        solenoidal_dissipation_store.append(solenoidal_dissipation)
        # dilatational_dissipation_calc = 2.0*(strain_rate_based_dissipation - deviatoric_strain_rate_based_dissipation)
        # print(np.linalg.norm(dilatational_dissipation_calc-dilatational_dissipation))

    final_time_for_plot = 20.0
    if(plotting_subsonic_result):
        final_time_for_plot = 10.0

    ylimits_for_plot = [0.02,0.14]
    if(plotting_subsonic_result):
        ylimits_for_plot = [0.04,0.13]

    if(lnstl_input_store_!=[]):
        lnstl_input_store = lnstl_input_store_

    qp.plotfxn(xdata=time_store,#[time,time],
            ydata=kinetic_energy_store,#[kinetic_energy,kolmogorov_slope],
            ylabel='Nondimensional Kinetic Energy, $K^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'kinetic_energy_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0,final_time_for_plot],
            ylimits=ylimits_for_plot,
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")
    
    #-----------------------------------------------------
    ylimits_for_plot = []
    qp.plotfxn(xdata=time_store,#[time,time],
            ydata=dissipation_store,#[kinetic_energy,kolmogorov_slope],
            ylabel='Nondimensional Dissipation Rate, $\\epsilon^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'dissipation_rate_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0,final_time_for_plot],
            ylimits=ylimits_for_plot,
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")

    #-----------------------------------------------------
    ylimits_for_plot = [0.0,0.012]
    if(plotting_subsonic_result):
        ylimits_for_plot = [0.0,0.012]
    time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    if(plotting_subsonic_result):
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/subsonic/solenoidal_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    time_store[0] = time # replace it -- this is a hack
    if(256 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation_256_flexi.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store.insert(1,time)
        solenoidal_dissipation_store.insert(1,solenoidal_dissipation)
        # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels.insert(1,"FLEXI\n(p$3$ DGSEM, sub-cell FV)\n[Chapelier et al.]")
        dashed_line_flag.insert(1,True)
        # clr_input_store.insert(1,"k")
        # lnstl_input_store.insert(1,"dashed")
        clr_input_store.insert(1,'tab:blue')
        lnstl_input_store.insert(1,"solid")
        # NS3D
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation_256_ns3d.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store.insert(2,time)
        solenoidal_dissipation_store.insert(2,solenoidal_dissipation)
        # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels.insert(2,"NS3D\n(FD-6, HO filter)\n[Chapelier et al.]")
        dashed_line_flag.insert(2,True)
        # clr_input_store.insert(2,"k")
        # lnstl_input_store.insert(2,"dotted")
        clr_input_store.insert(2,'tab:red')
        lnstl_input_store.insert(2,"solid")
    if(128 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation_128_flexi.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store.insert(1,time)
        solenoidal_dissipation_store.insert(1,solenoidal_dissipation)
        # labels.append("FLEXI $128^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels.insert(1,"FLEXI\n(p$3$ DGSEM, sub-cell FV)\n[Chapelier et al.]")
        dashed_line_flag.insert(1,True)
        # clr_input_store.insert(1,"k")
        # lnstl_input_store.insert(1,"dashed")
        clr_input_store.insert(1,'tab:blue')
        lnstl_input_store.insert(1,"solid")
        # NS3D
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation_128_ns3d.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store.insert(2,time)
        solenoidal_dissipation_store.insert(2,solenoidal_dissipation)
        # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels.insert(2,"NS3D\n(FD-6, HO filter)\n[Chapelier et al.]")
        dashed_line_flag.insert(2,True)
        # clr_input_store.insert(2,"k")
        # lnstl_input_store.insert(2,"dotted")
        clr_input_store.insert(2,'tab:red')
        lnstl_input_store.insert(2,"solid")
    if(64 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation_64_flexi.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store.insert(1,time)
        solenoidal_dissipation_store.insert(1,solenoidal_dissipation)
        # labels.append("FLEXI $128^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels.insert(1,"FLEXI\n(p$3$ DGSEM, sub-cell FV)\n[Chapelier et al.]")
        dashed_line_flag.insert(1,True)
        # clr_input_store.insert(1,"k")
        # lnstl_input_store.insert(1,"dashed")
        clr_input_store.insert(1,'tab:blue')
        lnstl_input_store.insert(1,"solid")
        # NS3D
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation_64_ns3d.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store.insert(2,time)
        solenoidal_dissipation_store.insert(2,solenoidal_dissipation)
        # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels.insert(2,"NS3D\n(FD-6, HO filter)\n[Chapelier et al.]")
        dashed_line_flag.insert(2,True)
        # clr_input_store.insert(2,"k")
        # lnstl_input_store.insert(2,"dotted")
        clr_input_store.insert(2,'tab:red')
        lnstl_input_store.insert(2,"solid")
    if(lnstl_input_store_!=[]):
        lnstl_input_store = lnstl_input_store_
    if(plotting_subsonic_result==False):
        labels[0]= "$2048^3$ DOF [Chapelier et al.]"
    qp.plotfxn(xdata=time_store,
            ydata=solenoidal_dissipation_store,
            ylabel='Nondimensional Solenoidal Dissipation, $\\varepsilon_{s}^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'solenoidal_dissipation_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0,final_time_for_plot],
            ylimits=ylimits_for_plot,
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="upper left")
    #-----------------------------------------------------
    ylimits_for_plot = [0.0,0.0012]
    if("supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128" in subdirectories_for_plot):
        ylimits_for_plot = [0.0,0.0013]
    if(plotting_subsonic_result):
        # ylimits_for_plot = []
        ylimits_for_plot = [0.0,5.0e-6]
    time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    if(plotting_subsonic_result):
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/subsonic/dilatational_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")    
    time_store[0] = time # replace it -- this is a hack

    if(plotting_subsonic_result==False):
        labels[0]= "$2048^3$ DOF\n[Chapelier et al.]"
    if(256 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        # FLEXI
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation_256_flexi.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store[1] = time
        dilatational_dissipation_store.insert(1,dilatational_dissipation)
        pressure_dissipation_store.insert(1,np.nan*dilatational_dissipation)
        # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels[1] = "FLEXI\n(p$3$ DGSEM, sub-cell FV)\n[Chapelier et al.]"
        # NS3D
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation_256_ns3d.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store[2] = time
        dilatational_dissipation_store.insert(2,dilatational_dissipation)
        pressure_dissipation_store.insert(2,np.nan*dilatational_dissipation)
        # labels.append("NS3D $256^3$\n($6^{th}$-order FD with HO Filter)\n[Chapelier et al.]")
        labels[2] = "NS3D\n(FD-6, HO filter)\n[Chapelier et al.]"
        dashed_line_flag[2]=True
        # dashed_line_flag.insert(1,True)
        # clr_input_store.insert(1,"k")
        # lnstl_input_store.insert(1,"dotted")
    if(128 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        # FLEXI
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation_128_flexi.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store[1] = time
        dilatational_dissipation_store.insert(1,dilatational_dissipation)
        pressure_dissipation_store.insert(1,np.nan*dilatational_dissipation)
        # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels[1] = "FLEXI\n(p$3$ DGSEM, sub-cell FV)\n[Chapelier et al.]"
        # NS3D
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation_128_ns3d.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store[2] = time
        dilatational_dissipation_store.insert(2,dilatational_dissipation)
        pressure_dissipation_store.insert(2,np.nan*dilatational_dissipation)
        # labels.append("NS3D $256^3$\n($6^{th}$-order FD with HO Filter)\n[Chapelier et al.]")
        labels[2] = "NS3D\n(FD-6, HO filter)\n[Chapelier et al.]"
        dashed_line_flag[2]=True
        # dashed_line_flag.insert(1,True)
        # clr_input_store.insert(1,"k")
        # lnstl_input_store.insert(1,"dotted")
    if(64 in number_of_degrees_of_freedom and (compare_with_reference_result_at_same_degrees_of_freedom==True)):
        # FLEXI
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation_64_flexi.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store[1] = time
        dilatational_dissipation_store.insert(1,dilatational_dissipation)
        pressure_dissipation_store.insert(1,np.nan*dilatational_dissipation)
        # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
        labels[1] = "FLEXI\n(p$3$ DGSEM, sub-cell FV)\n[Chapelier et al.]"
        # NS3D
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation_64_ns3d.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store[2] = time
        dilatational_dissipation_store.insert(2,dilatational_dissipation)
        pressure_dissipation_store.insert(2,np.nan*dilatational_dissipation)
        # labels.append("NS3D $256^3$\n($6^{th}$-order FD with HO Filter)\n[Chapelier et al.]")
        labels[2] = "NS3D\n(FD-6, HO filter)\n[Chapelier et al.]"
        dashed_line_flag[2]=True
        # dashed_line_flag.insert(1,True)
        # clr_input_store.insert(1,"k")
        # lnstl_input_store.insert(1,"dotted")

    qp.plotfxn(xdata=time_store,
            ydata=dilatational_dissipation_store,
            ylabel='Nondimensional Dilatational Dissipation, $\\varepsilon_{d}^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'dilatational_dissipation_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0,final_time_for_plot],
            ylimits=ylimits_for_plot,
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="upper right")

    #-----------------------------------------------------
    ylimits_for_plot = [0.0,0.002]
    
    # time_store.pop(0);
    qp.plotfxn(xdata=time_store,
            ydata=pressure_dissipation_store,
            ylabel='Nondimensional Pressure Dissipation, $\\varepsilon_{p}^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'pressure_dissipation_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0,final_time_for_plot],
            # ylimits=ylimits_for_plot,
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=True,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="upper right")
    if(64 in number_of_degrees_of_freedom):
        time, kinetic_energy = np.loadtxt("./data/chapelier2024/kinetic_energy.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
        time_store[0] = time # replace it -- this is a hack
        ylimits_for_plot = [0.02,0.14]
        if(compare_with_reference_result_at_same_degrees_of_freedom==True):
            # FLEXI
            time, kinetic_energy = np.loadtxt("./data/chapelier2024/kinetic_energy_64_flexi.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            time_store[1] = time
            kinetic_energy_store.insert(1,kinetic_energy)
            # kinetic_energy_store.insert(1,time)
            # labels.append("FLEXI $256^3$\n($4^{th}$-order DGSEM with subgrid FV)\n[Chapelier et al.]")
            labels[1] = "FLEXI\n(p$3$ DGSEM, sub-cell FV)\n[Chapelier et al.]"
            # NS3D
            time, kinetic_energy = np.loadtxt("./data/chapelier2024/kinetic_energy_64_ns3d.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
            time_store[2] = time
            kinetic_energy_store.insert(2,kinetic_energy)
            # kinetic_energy_store.insert(2,time)
            labels[2] = "NS3D\n(FD-6, HO filter)\n[Chapelier et al.]"
            dashed_line_flag[2]=True
        qp.plotfxn(xdata=time_store,#[time,time],
                ydata=kinetic_energy_store,#[kinetic_energy,kolmogorov_slope],
                ylabel='Nondimensional Kinetic Energy, $K^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'kinetic_energy_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels,
                black_lines=False,
                xlimits=[0,final_time_for_plot],
                ylimits=ylimits_for_plot,
                log_axes=log_axes_input,
                which_lines_black=black_line_flag,
                which_lines_dashed=dashed_line_flag,
                which_lines_only_markers=[0],
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=True,#transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=12,#14
                legend_location="best")
#-----------------------------------------------------
#=====================================================
def reinit_inputs():
    global subdirectories, filenames, labels, black_line_flag, \
    dashed_line_flag, figure_filename_postfix, figure_title, \
    ylimits_kinetic_energy_input, ylimits_dissipation_input, \
    log_axes_input, legend_on_input, legend_inside_input, \
    plot_reference_result, nlegendcols_input, \
    figure_subdirectory, data_directory_base, figure_directory_base, \
    smoothing_input, plot_PHiLiP_DNS_result_as_reference_input

    subdirectories = []
    filenames = []
    labels = []
    black_line_flag = []
    dashed_line_flag = []
    figure_filename_postfix = "" # default
    figure_title = "" # default
    ylimits_kinetic_energy_input = [] # default
    ylimits_dissipation_input = [] # default
    log_axes_input=None # default
    legend_on_input=True # default
    legend_inside_input=False # default
    plot_reference_result=False # default
    nlegendcols_input=1
    figure_subdirectory="" # default
    data_directory_base = "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex"
    # figure_directory_base = "/Users/Julien/julien_phd/post_processing/figures/taylor_green_vortex"
    figure_directory_base = "figures"
    plot_PHiLiP_DNS_result_as_reference_input=True # default
#=====================================================
#-----------------------------------------------------
smooth_dilatational_dissipation_rate_input=[True,True,True,True,True,True,True,True]
# smooth_dilatational_dissipation_rate_input=[False,False,False,False,False,False,False,False]
#=====================================================
# All results
#-----------------------------------------------------
#-----------------------------------------------------
# clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
reinit_inputs()
data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
date_for_runs="."
figure_subdirectory="./"
# figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
figure_filename_postfix = "_128_p7_time_step_advantage"
legend_inside_input=True
plot_reference_result=True
plot_PHiLiP_DNS_result_as_reference_input=False
#-----------------------------------------------------
subdirectories_for_plot=[\
"supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
"supersonic_viscous_TGV_ILES_NSFR_cHU_Ra_2PF_GLL_OI-0_dofs0128_p7_procs512",\
"/Volumes/KAUST/NarvalFiles/2024_JCP/time_step_advantage/supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p7_procs512_CFL-0point25",\
"/Volumes/KAUST/NarvalFiles/2024_JCP/time_step_advantage/supersonic_viscous_TGV_ILES_NSFR_cHU_Ra_2PF_GLL_OI-0_dofs0128_p7_procs512_CFL-0point3",\
]
# labels
labels_for_plot=[\
"$128^{3}$, $c_{DG}$, CFL$=0.1$",\
"$128^{3}$, $c_{HU}$, CFL$=0.1$",\
"$128^{3}$, $c_{DG}$, CFL$=0.25$",\
"$128^{3}$, $c_{HU}$, CFL$=0.3$",\
]
black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
dashed_line_flag_for_plot=[False,False,False,True,False,True,True]
which_was_ran_with_corrected_quantites=[1]
number_of_degrees_of_freedom_input=[]
compare_with_ref_result_at_same_dof=False
plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
    which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
    number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
    compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
    smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
    smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))

exit()
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_128_high_poly_degree"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128",\
    ]
    labels_for_plot=[\
    "p$3$",\
    "p$7$",\
    "p$15$",\
    # "p$15$ smoothed",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[0]
    number_of_degrees_of_freedom_input=[128]
    compare_with_ref_result_at_same_dof=False
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        # smooth_dilatational_dissipation_rate=[True,True,False,True],
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_p7_convergence"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0512_p7_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$64^{3}$",\
    # "$64^{3}$ smoothed",\
    "$128^{3}$",\
    # "$128^{3}$ smoothed",\
    "$256^{3}$",\
    # "$256^{3}$ smoothed",\
    # "$512^{3}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,True,False,True,True]
    which_was_ran_with_corrected_quantites=[1]
    number_of_degrees_of_freedom_input=[]
    compare_with_ref_result_at_same_dof=False
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_p3_convergence"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p3_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p3_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0512_p7_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$64^{3}$",\
    # "$64^{3}$ smoothed",\
    "$128^{3}$",\
    # "$128^{3}$ smoothed",\
    "$256^{3}$",\
    # "$256^{3}$ smoothed",\
    # "$512^{3}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[0]
    number_of_degrees_of_freedom_input=[]
    compare_with_ref_result_at_same_dof=False
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_64_two_point_flux"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_KG_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_CH_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    ]
    # labels
    # labels_for_plot=[\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $16$p$7$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL\n $32$p$3$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $8$p$15$ ($128^3$ DOF) CFL=0.01",\
    # ]
    labels_for_plot=[\
    "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{DG}$ NSFR.IR",\
    "p$3$ $c_{DG}$ NSFR.KG",\
    "p$3$ $c_{DG}$ NSFR.CH",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[2]
    number_of_degrees_of_freedom_input=[64]
    compare_with_ref_result_at_same_dof=False
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        lnstl_input_store_ = ['None','solid','dashed','solid','dashed','dashed','solid','dashed','solid'],
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_64_correction_parameter"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cHU_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    ]
    # labels
    # labels_for_plot=[\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $16$p$7$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL\n $32$p$3$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $8$p$15$ ($128^3$ DOF) CFL=0.01",\
    # ]
    labels_for_plot=[\
    "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$7$ $c_{HU}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{+}$ NSFR.CH$_{\\mathrm{RA}}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[2]
    number_of_degrees_of_freedom_input=[64]
    compare_with_ref_result_at_same_dof=False
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        lnstl_input_store_ = ['None','solid','dashed','solid','dashed','solid','dashed','solid','dashed','solid'],
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    # figure_filename_postfix = "_064_correction_parameter"
    figure_filename_postfix = "_64"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128",\
    # "supersonic_viscous_TGV_ILES_NSFR_cHU_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    # "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_KG_2PF_GLL_OI-0_dofs0064_p3_procs128",\
    ]
    # labels
    # labels_for_plot=[\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $16$p$7$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL\n $32$p$3$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $8$p$15$ ($128^3$ DOF) CFL=0.01",\
    # ]
    labels_for_plot=[\
    "p$7$ $c_{DG}$ NSFR",\
    # "p$7$ $c_{HU}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{DG}$ NSFR",\
    # "p$3$ $c_{+}$ NSFR.CH$_{\\mathrm{RA}}$",\
    # "p$3$ $c_{DG}$ NSFR.KG",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[2]
    number_of_degrees_of_freedom_input=[64]
    compare_with_ref_result_at_same_dof=True
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    # exit()
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_128_correction_parameter"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cHU_Ra_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    ]
    # labels
    # labels_for_plot=[\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $16$p$7$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL\n $32$p$3$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $8$p$15$ ($128^3$ DOF) CFL=0.01",\
    # ]
    labels_for_plot=[\
    "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$7$ $c_{HU}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{+}$ NSFR.CH$_{\\mathrm{RA}}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,True,False,True,False,True,True]
    which_was_ran_with_corrected_quantites=[2]
    number_of_degrees_of_freedom_input=[128]
    compare_with_ref_result_at_same_dof=False
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        lnstl_input_store_ = ['None','solid','dashed','solid','dashed','solid','dashed','solid','dashed','solid'],
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_256"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p3_procs512",\
    ]
    # labels
    labels_for_plot=[\
    # "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GLL-Roe-PPL",\
    # "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-Roe\n(with PPL)",\
    "p$7$ $c_{DG}$ NSFR",\
    "p$3$ $c_{DG}$ NSFR",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[0]
    number_of_degrees_of_freedom_input=[256]
    compare_with_ref_result_at_same_dof=True
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_128"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    ]
    # labels
    # labels_for_plot=[\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $16$p$7$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL\n $32$p$3$ ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL\n $8$p$15$ ($128^3$ DOF) CFL=0.01",\
    # ]
    labels_for_plot=[\
    # "p$15$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GLL-Roe-PPL",\
    # "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GLL-Roe-PPL",\
    # "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GLL-Roe-PPL",\
    # "p$3$ $c_{+}$ NSFR.CH$_{\\mathrm{RA}}$-GLL-Roe-PPL",\
    # "p$15$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$7$ $c_{DG}$ NSFR",\
    "p$3$ $c_{DG}$ NSFR",\
    # "p$3$ $c_{+}$ NSFR.CH$_{\\mathrm{RA}}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[1]
    number_of_degrees_of_freedom_input=[128]
    compare_with_ref_result_at_same_dof=True
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites,
        number_of_degrees_of_freedom=number_of_degrees_of_freedom_input,
        compare_with_reference_result_at_same_degrees_of_freedom=compare_with_ref_result_at_same_dof,
        smooth_dilatational_dissipation_rate=smooth_dilatational_dissipation_rate_input,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))

#=====================================================
# DOFs: 256^3 | Subsonic case
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_subsonic"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "subsonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512",\
    # "subsonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512_no_limiter",\
    # "subsonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512_corrected_quantities",\
    ]
    # labels
    labels_for_plot=[\
    "PHiLiP $256^3$ DOF\n(p$7$ $c_{DG}$ NSFR, PPL)",\
    #"$c_{DG}$ NSFR.CH$_{RA}$-GLL-Roe-PPL $32^{3}$p$7$\n ($256^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe $32^{3}$p$7$\n ($256^3$ DOF) CFL=0.1",\
    # "with correction",\
    # "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$3$\n ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $8$p$15$\n ($128^3$ DOF) CFL=0.01",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$7$\n ($256^3$ DOF) CFL=0.1",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,True,False,False,False,True,True]
    which_was_ran_with_corrected_quantites=[2]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,plotting_subsonic_result=True,which_was_ran_with_corrected_quantites=which_was_ran_with_corrected_quantites)

exit()
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

filename = filesystem+"NarvalFiles/2024_JCP/"+"supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512/turbulent_quantities.txt"
time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
ax.semilogx(time, solenoidal_dissipation)

filename = filesystem+"NarvalFiles/2024_JCP/"+"supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128/turbulent_quantities.txt"
time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)

ax.semilogx(time, solenoidal_dissipation)
# ax.plot(time, kinetic_energy)

filename = filesystem+"NarvalFiles/2024_JCP/"+"supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128_cfl_1e-2/turbulent_quantities.txt"
time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)

ax.semilogx(time, solenoidal_dissipation)

filename = filesystem+"NarvalFiles/2024_JCP/"+"supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128_cfl_5e-2/turbulent_quantities.txt"
time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)

ax.semilogx(time, solenoidal_dissipation)

filename = filesystem+"NarvalFiles/2024_JCP/"+"supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128_cfl_25e-3/turbulent_quantities.txt"
time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)

ax.semilogx(time, solenoidal_dissipation)
# ax.plot(time, kinetic_energy)

ax.set(xlabel='t', ylabel='$\\varepsilon_{s}$',
       title='check')
# ax.grid()

# fig.savefig("test.png")
plt.show()

#=====================================================
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(False):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_128_overint10"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512_overint-10-for-unsteady-quantities",\
    "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512_overint-10-for-unsteady-quantities",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $16$p$7$\n ($128^2$ DOF) CFL=0.1",\
    "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$3$\n ($128^2$ DOF) CFL=0.1",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot)

#=====================================================
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(False):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_128_check"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512_overint-10-for-unsteady-quantities",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $16$p$7$\n ($128^2$ DOF) CFL=0.1",\
    "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $16$p$7$\n ($128^2$ DOF) CFL=0.1, OI-10",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,True,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot)
