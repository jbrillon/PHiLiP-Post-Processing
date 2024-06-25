#-----------------------------------------------------
# Python code for plotting unsteady periodic turbulent flow quantities
# - Post-processing code
# - Written by Julien Brillon
# - Email: julien.brillon@mail.mcgill.ca
# - McGill University
#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
import scipy # SciPy: contains additional numerical routines to numpy
#-----------------------------------------------------
# Import personal libraries
# from finite_difference_library import first_derivative, fd_non_uniform_grid
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../submodules/quickplotlib/lib"); import quickplotlib as qp
# from sys import platform
# if platform == "linux" or platform == "linux2":
#     # linux
#     sys.path.append("/home/julien/Codes/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes
# elif platform == "darwin":
#     # OS X
#     sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes

import matplotlib;from matplotlib.lines import Line2D
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
#-----------------------------------------------------
def plot_periodic_turbulence(
    figure_subdirectory,
    subdirectories,
    filenames,
    labels,
    black_line_flag,
    dashed_line_flag,
    figure_directory_base,
    data_directory_base,
    plot_reference_result,
    figure_filename_postfix,
    figure_title,
    log_axes_input,
    legend_on_input,
    legend_inside_input,
    nlegendcols_input,
    transparent_legend_input=True,
    plot_kinetic_energy=True,
    plot_dissipation_rate=True,
    plot_enstrophy=True,
    plot_dissipation_components=False,
    plot_numerical_viscosity=False,
    clr_input=[],mrkr_input=[],lnstl_input=[],
    legend_fontSize_input=16,
    tmax=20.0,
    solid_and_dashed_lines=False,
    dashed_and_solid_lines=False,
    reference_result_author="Dairay et al.",
    plot_numerical_dissipation=False,
    plot_PHiLiP_DNS_result_as_reference=False,
    dissipation_rate_smoothing=[],
    plot_filtered_dns=False,
    plot_zoomed_section_dissipation_rate=False,
    plot_zoomed_section_numerical_dissipation_components=False,
    plot_zoomed_section_enstrophy=False,
    dofs_for_zoomed_section=256,
    ):
    # plotting parameters store
    labels_store = []
    which_lines_black_input=[]
    which_lines_dashed_input=[]
    # data store
    time_store = []
    kinetic_energy_store = []
    dissipation_store = []
    enstrophy_store = []
    vorticity_based_dissipation_store = []
    strain_rate_based_dissipation_store = []
    deviatoric_strain_rate_based_dissipation_store = []
    pressure_dilatation_based_dissipation_store = []
    eps_K_minus_eps_S_plus_eps_p_store = []
    eps_K_minus_eps_Sd_plus_eps_p_store = []
    eps_K_minus_eps_S_store = []
    eps_K_minus_eps_Sd_store = []
    eps_S_minus_eps_Sd_store = []
    eps_S_plus_eps_p_store = []
    eps_p_store = []
    eps_S_store = []
    numerical_viscosity_store = []
    clr_input_store=[]
    mrkr_input_store=[]
    lnstl_input_store=[]
    i_curve = 0
    figure_subdirectory = figure_subdirectory+"/"
    # post-fix
    if(figure_filename_postfix!=""):
        figure_filename_postfix = "_" + figure_filename_postfix
    # reference result
    if(plot_reference_result):
        if(reference_result_author=="Vermeire"):
            # DNS - Kinetic Energy
            path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/vermiere"
            filename=path_to_reference_result+"/"+"kinetic_energy"+".txt"
            time, kinetic_energy = np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
            time_store.append(time)
            kinetic_energy_store.append(kinetic_energy)
            # labels_store.append("Spectral DNS (P4, $260^{3}$ DOFs)")
            # labels_store.append("DNS ($260^{3}$ DOFs)\n [Vermeire, 2014]")
            labels_store.append("DNS [Vermeire]")
            # DNS - dissipation
            filename=path_to_reference_result+"/"+"dissipation"+".txt"
            time, dissipation = np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
            dissipation_store.append(dissipation)
            # DNS - enstrophy
            filename=path_to_reference_result+"/"+"enstrophy"+".txt"
            time, enstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
            enstrophy_store.append(enstrophy)
        elif(reference_result_author=="Dairay et al."):
            labels_store.append("DNS [Dairay et al.]")
            filename=CURRENT_PATH+"../cases/taylor_green_vortex/data/TGV_Re1600.dat"
            time, kinetic_energy, dissipation, enstrophy = np.loadtxt(filename,skiprows=43,dtype=np.float64,unpack=True,usecols=(0,1,2,4))
            time_store.append(time)
            kinetic_energy_store.append(kinetic_energy)
            dissipation_store.append(dissipation)
            enstrophy_store.append(enstrophy)
            vorticity_based_dissipation_store.append((2.0/1600.0)*enstrophy)
            numerical_viscosity_store.append(1600.0*dissipation/(2.0*enstrophy))
        else:
            print("ERROR: Invalid value passed for reference_result_author. Aborting...")
            exit()
        # black_line_flag.append(True) # inputs
        # dashed_line_flag.append(True) # inputs
        # if(black_line_flag[i_curve]):
        #   which_lines_black_input.append(i_curve)
        # if(dashed_line_flag[i_curve]):
        #   which_lines_dashed_input.append(i_curve)
        which_lines_black_input.append(i_curve)
        # which_lines_dashed_input.append(i_curve) # uncomment for dashed DNS result
        i_curve += 1
        if(clr_input!=[]):
            clr_input_store.append('k')
        if(mrkr_input!=[]):
            mrkr_input_store.append('None')
        if(lnstl_input!=[]):
            lnstl_input_store.append('solid') # supported values are '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
    elif(plot_PHiLiP_DNS_result_as_reference):
        labels_store.append("DNS ($256^3$ DOFs, P$7$)")
        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/brillon"
        filename=path_to_reference_result+"/"+"turbulent_quantities_256dofs_p7"+".txt"
        time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        # -- compute dissipation
        dissipation = get_dissipation_discrete(time,kinetic_energy)
        dissipation_store.append(dissipation)
        enstrophy_store.append(enstrophy)
        vorticity_based_dissipation_store.append(vorticity_based_dissipation)
        pressure_dilatation_based_dissipation_store.append(pressure_dilatation_based_dissipation)
        strain_rate_based_dissipation_store.append(strain_rate_based_dissipation)
        deviatoric_strain_rate_based_dissipation_store.append(deviatoric_strain_rate_based_dissipation)
        eps_K_minus_eps_S_plus_eps_p_store.append(dissipation - strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_K_minus_eps_Sd_plus_eps_p_store.append(dissipation - deviatoric_strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_K_minus_eps_Sd_store.append(dissipation - deviatoric_strain_rate_based_dissipation)
        eps_K_minus_eps_S_store.append(dissipation - strain_rate_based_dissipation)
        eps_S_minus_eps_Sd_store.append(strain_rate_based_dissipation - deviatoric_strain_rate_based_dissipation)
        eps_S_plus_eps_p_store.append(strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_p_store.append(pressure_dilatation_based_dissipation)
        eps_S_store.append(strain_rate_based_dissipation)
        numerical_viscosity_store.append(dissipation/vorticity_based_dissipation)
        which_lines_black_input.append(i_curve)
        # which_lines_dashed_input.append(i_curve) # uncomment for dashed DNS result
        i_curve += 1
        if(clr_input!=[]):
            clr_input_store.append('k')
        if(mrkr_input!=[]):
            mrkr_input_store.append('None')
        if(lnstl_input!=[]):
            lnstl_input_store.append('solid') # supported values are '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
    if(plot_filtered_dns):
        labels_store.append("Projected DNS ($96^3$ DOFs, P$2$)")
        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/brillon/filtered_dns"
        filename=path_to_reference_result+"/"+"turbulent_quantities_96dofs_p2"+".txt"
        time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        # -- compute dissipation
        dissipation = get_dissipation_discrete(time,kinetic_energy)
        dissipation_store.append(dissipation)
        enstrophy_store.append(enstrophy)
        vorticity_based_dissipation_store.append(vorticity_based_dissipation)
        pressure_dilatation_based_dissipation_store.append(pressure_dilatation_based_dissipation)
        strain_rate_based_dissipation_store.append(strain_rate_based_dissipation)
        deviatoric_strain_rate_based_dissipation_store.append(deviatoric_strain_rate_based_dissipation)
        eps_K_minus_eps_S_plus_eps_p_store.append(dissipation - strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_K_minus_eps_Sd_plus_eps_p_store.append(dissipation - deviatoric_strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_K_minus_eps_Sd_store.append(dissipation - deviatoric_strain_rate_based_dissipation)
        eps_K_minus_eps_S_store.append(dissipation - strain_rate_based_dissipation)
        eps_S_minus_eps_Sd_store.append(strain_rate_based_dissipation - deviatoric_strain_rate_based_dissipation)
        eps_S_plus_eps_p_store.append(strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_p_store.append(pressure_dilatation_based_dissipation)
        eps_S_store.append(strain_rate_based_dissipation)
        numerical_viscosity_store.append(dissipation/vorticity_based_dissipation)
        which_lines_black_input.append(i_curve)
        which_lines_dashed_input.append(i_curve) # for dashed filtered DNS result
        i_curve += 1
        if(clr_input!=[]):
            clr_input_store.append('k')
        if(mrkr_input!=[]):
            mrkr_input_store.append('None')
        if(lnstl_input!=[]):
            lnstl_input_store.append('.') # supported values are '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
    # PHiLiP Results
    number_of_result_curves = len(filenames)
    for i in range(0,number_of_result_curves):
        if(black_line_flag[i]):
            which_lines_black_input.append(i_curve)
        if(dashed_line_flag[i]):
            which_lines_dashed_input.append(i_curve)
        i_curve += 1
        labels_store.append(labels[i])
        # load file
        filename = data_directory_base+"/"+subdirectories[i]+"/"+filenames[i]
        time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        # store data
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        # -- compute dissipation
        if(dissipation_rate_smoothing!=[]):
            dissipation = get_dissipation_discrete(time,kinetic_energy,dissipation_rate_smoothing[i])
            dissipation_store.append(dissipation)
        else:
            dissipation = get_dissipation_discrete(time,kinetic_energy)
            dissipation_store.append(dissipation)
        # -- store other quantities
        if(dissipation_rate_smoothing!=[]):
            # smooth this too for the numerical dissipation plots
            enstrophy = splrep(time,enstrophy,k=5,s=0.0000001)
            enstrophy = splev(time,enstrophy)
            vorticity_based_dissipation = splrep(time,vorticity_based_dissipation,k=5,s=0.0000001)
            vorticity_based_dissipation = splev(time,vorticity_based_dissipation)
        enstrophy_store.append(enstrophy)
        vorticity_based_dissipation_store.append(vorticity_based_dissipation)
        pressure_dilatation_based_dissipation_store.append(pressure_dilatation_based_dissipation)
        strain_rate_based_dissipation_store.append(strain_rate_based_dissipation)
        deviatoric_strain_rate_based_dissipation_store.append(deviatoric_strain_rate_based_dissipation)
        eps_K_minus_eps_S_plus_eps_p_store.append(dissipation - strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_K_minus_eps_Sd_plus_eps_p_store.append(dissipation - deviatoric_strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_K_minus_eps_Sd_store.append(dissipation - deviatoric_strain_rate_based_dissipation)
        eps_K_minus_eps_S_store.append(dissipation - strain_rate_based_dissipation)
        eps_S_minus_eps_Sd_store.append(strain_rate_based_dissipation - deviatoric_strain_rate_based_dissipation)
        eps_S_plus_eps_p_store.append(strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_p_store.append(pressure_dilatation_based_dissipation)
        eps_S_store.append(strain_rate_based_dissipation)
        numerical_viscosity_store.append(dissipation/vorticity_based_dissipation)
        # store inputted line color, markers, and linestyles
        if(clr_input!=[]):
            clr_input_store.append(clr_input[i])
        if(mrkr_input!=[]):
            mrkr_input_store.append(mrkr_input[i])
        if(lnstl_input!=[]):
            lnstl_input_store.append(lnstl_input[i])

    # line parameters if doing solid and dashed lines
    if(solid_and_dashed_lines):
        clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green','tab:orange','tab:orange']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        mrkr_input_store = ['None','None','None','None','None','None','None','None']
        lnstl_input_store = ['solid','solid','dashed','solid','dashed','solid','dashed','solid','dashed']
        if(plot_filtered_dns):
            clr_input_store.insert(1,'k')
            mrkr_input_store.insert(1,'None')
            lnstl_input_store.insert(1,'dashed')
    elif(dashed_and_solid_lines):
        clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green','tab:orange','tab:orange']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        mrkr_input_store = ['None','None','None','None','None','None','None','None']
        lnstl_input_store = ['dashed','dashed','solid','dashed','solid','dashed','solid','dashed','solid']
        if(plot_filtered_dns):
            clr_input_store.insert(1,'k')
            mrkr_input_store.insert(1,'None')
            lnstl_input_store.insert(1,'solid')

    # zoom box location
    zoom_box_origin_and_extent_input=[0.65, 0.65, 0.32, 0.32]
    if(tmax<20.0):
        zoom_box_origin_and_extent_input=[0.60, 0.01, 0.32, 0.32]

    #-----------------------------------------------------
    # evolution of kinetic energy:
    #-----------------------------------------------------
    if (plot_kinetic_energy):
        # kolmogorov_slope = (-5.0/3.0)*time/10000.0+kinetic_energy[0]
        qp.plotfxn(xdata=time_store,#[time,time],
                ydata=kinetic_energy_store,#[kinetic_energy,kolmogorov_slope],
                ylabel='Nondimensional Kinetic Energy, $K^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'kinetic_energy_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                ylimits=[0.0,0.14],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=True,#transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input,
                legend_location="lower left")#,
                # legend_anchor=[0.025,0.3])
    
    #-----------------------------------------------------
    # evolution of kinetic energy dissipation rate:
    #-----------------------------------------------------
    if(plot_reference_result and reference_result_author=="Vermeire"):
        # DNS - dissipation
        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/vermiere"
        filename=path_to_reference_result+"/"+"dissipation"+".txt"
        time, dissipation = np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
        time_store[0] = time # replace it -- this is a hack

    if(plot_dissipation_rate):
        qp.plotfxn(xdata=time_store,
                ydata=dissipation_store,
                # ylabel='$\\varepsilon=-\\frac{\\mathrm{d} K^{*}}{\\mathrm{d}t^{*}}$',
                ylabel='Nondimensional Dissipation Rate, $\\varepsilon^{*}$',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'dissipation_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                ylimits=[0.0,0.018],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input,
                legend_location="upper left",
                plot_zoomed_section=plot_zoomed_section_dissipation_rate,
                x_limits_zoom=[8.0, 10.5],y_limits_zoom=[0.010, 0.0135],
                zoom_box_origin_and_extent=zoom_box_origin_and_extent_input)

    # numerical dissipation plot - can do a max of 4 different results (3 curves per result) -- need a custom legend for this -- can hack the indexing
    if(plot_numerical_dissipation):
        KE_molecular_and_numerical_dissipation_y_store = []
        KE_molecular_and_numerical_dissipation_x_store = []
        lnstl_input_dummy=['solid','solid','dashed','dotted','dashdot']
        # mrkr_input_dummy=['None','None','None','None','None']
        clr_input_dummy = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        leg_elements_input=[]
        second_leg_elements_input=[]
        clr_input_store_numerical_dissipation=[]
        mrkr_input_store_numerical_dissipation=[]
        lnstl_input_store_numerical_dissipation=[]
        # legend_components_input = []
        if(number_of_result_curves>10):
            print("ERROR: Can only plot numerical dissipation for 10 result curves. Aborting...")
            exit()
        
        # explain the linestyles
        ls='solid'
        # "$\\varepsilon=-\\frac{\\mathrm{d} K^{*}}{\\mathrm{d}t^{*}}$"
        second_leg_elements_input.append(Line2D([0],[0], label="$\\varepsilon\\left(K^{*}\\right)$", color='grey', marker='None', markersize=6, mfc='None', linestyle=ls))
        ls='dashed'
        second_leg_elements_input.append(Line2D([0],[0], label="$\\varepsilon\\left(\\zeta^{*}\\right)$", color='grey', marker='None', markersize=6, mfc='None', linestyle=ls))
        ls='dotted'
        second_leg_elements_input.append(Line2D([0],[0], label="$\\varepsilon\\left(K^{*}\\right)-\\varepsilon\\left(\\zeta^{*}\\right)$", color='grey', marker='None', markersize=6, mfc='None', linestyle=ls))
        
        index_shift_for_ref_result = 0 # initialize
        if(plot_reference_result and reference_result_author!="Dairay et al."):
            index_shift_for_ref_result = -1 # minus 1 because no reference result

        second_leg_anchor_input=[]
        if(tmax!=20.0):
            second_leg_anchor_input=[0.0,0.5]
        elif(plot_zoomed_section_numerical_dissipation_components):
            second_leg_anchor_input=[0.0,0.7]

        # results
        number_of_refernce_curves_numerical_dissipation=1 # default has one
        if(plot_filtered_dns):
            number_of_refernce_curves_numerical_dissipation+=1

        for i in range(0,number_of_result_curves+number_of_refernce_curves_numerical_dissipation):
            # ls=lnstl_input_dummy[i]
            ls='solid'
            # mk=mrkr_input_dummy[i]
            mk='None'
            if(plot_filtered_dns and i==1):
                clr_input_dummy.insert(i,"tab:gray")
            lc=clr_input_dummy[i]
            leg_elements_input.append(Line2D([0],[0], label=labels_store[i], color=lc, marker=mk, markersize=6, mfc='None', linestyle=ls))

            KE_molecular_and_numerical_dissipation_y_store.append(dissipation_store[i])
            KE_molecular_and_numerical_dissipation_x_store.append(time_store[i])
            clr_input_store_numerical_dissipation.append(lc)
            mrkr_input_store_numerical_dissipation.append(mk)
            lnstl_input_store_numerical_dissipation.append(ls)

            if(((plot_reference_result and reference_result_author!="Dairay et al.") and i>0) or ((plot_PHiLiP_DNS_result_as_reference) or (plot_reference_result and reference_result_author=="Dairay et al."))):
                # molecular dissipation
                KE_molecular_and_numerical_dissipation_y_store.append(vorticity_based_dissipation_store[i+index_shift_for_ref_result])
                KE_molecular_and_numerical_dissipation_x_store.append(time_store[i])
                ls='dashed'
                clr_input_store_numerical_dissipation.append(lc)
                mrkr_input_store_numerical_dissipation.append(mk)
                lnstl_input_store_numerical_dissipation.append(ls)
                # numerical dissipation
                KE_molecular_and_numerical_dissipation_y_store.append(dissipation_store[i] - vorticity_based_dissipation_store[i+index_shift_for_ref_result])
                KE_molecular_and_numerical_dissipation_x_store.append(time_store[i])
                ls='dotted'
                clr_input_store_numerical_dissipation.append(lc)
                mrkr_input_store_numerical_dissipation.append(mk)
                lnstl_input_store_numerical_dissipation.append(ls)

        x_limits_zoom=[8, 10.5]
        y_limits_zoom=[0.010, 0.0135]
        if(plot_zoomed_section_enstrophy):
            if(dofs_for_zoomed_section==96):
                y_limits_zoom=[0.008, 0.0135]
                x_limits_zoom=[7.5, 11]

        qp.plotfxn(xdata=KE_molecular_and_numerical_dissipation_x_store,
                    ydata=KE_molecular_and_numerical_dissipation_y_store,
                    ylabel='Nondimensional Dissipation Components',
                    # ylabel='$\\varepsilon\\left(\\zeta^{*}\\right)$',
                    xlabel='Nondimensional Time, $t^{*}$',
                    figure_filename=figure_subdirectory+'numerical_dissipation_vs_time'+figure_filename_postfix,
                    title_label=figure_title,
                    markers=False,
                    # legend_labels_tex=labels_store,
                    leg_elements_input=leg_elements_input,
                    black_lines=False,
                    xlimits=[0,tmax],
                    ylimits=[0,0.018],
                    log_axes=log_axes_input,
                    which_lines_black=which_lines_black_input,
                    which_lines_dashed=which_lines_dashed_input,
                    legend_on=True,
                    legend_inside=legend_inside_input,
                    nlegendcols=nlegendcols_input,
                    figure_size=(6,6),
                    transparent_legend=transparent_legend_input,
                    legend_border_on=False,
                    grid_lines_on=False,
                    fig_directory=figure_directory_base,
                    clr_input=clr_input_store_numerical_dissipation,
                    mrkr_input=mrkr_input_store_numerical_dissipation,
                    lnstl_input=lnstl_input_store_numerical_dissipation,
                    legend_fontSize=legend_fontSize_input,
                    legend_location="upper left",
                    second_leg_elements_input=second_leg_elements_input,
                    second_leg_anchor=second_leg_anchor_input,
                    plot_zoomed_section=plot_zoomed_section_numerical_dissipation_components,
                    x_limits_zoom=x_limits_zoom,y_limits_zoom=y_limits_zoom,
                    zoom_box_origin_and_extent=zoom_box_origin_and_extent_input)

    if(plot_reference_result and reference_result_author=="Vermeire"):
        # DNS - enstrophy
        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/vermiere"
        filename=path_to_reference_result+"/"+"enstrophy"+".txt"
        time, enstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
        time_store[0] = time # replace it -- this is a hack

    if(plot_enstrophy):
        # entrophy
        x_limits_zoom=[8, 10.5]
        y_limits_zoom=[8.0, 10.5]
        if(plot_zoomed_section_enstrophy):
            if(dofs_for_zoomed_section==96):
                y_limits_zoom=[6.5, 8.5]
                x_limits_zoom=[8, 11]
                
        qp.plotfxn(xdata=time_store,
                ydata=enstrophy_store,
                ylabel='Nondimensional Enstrophy, $\\zeta^{*}$',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'enstrophy_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                ylimits=[0,14],
                xlimits=[0,tmax],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input,
                legend_location="upper left",
                plot_zoomed_section=plot_zoomed_section_enstrophy,
                x_limits_zoom=x_limits_zoom,y_limits_zoom=y_limits_zoom,
                zoom_box_origin_and_extent=zoom_box_origin_and_extent_input)

    # Remove the reference result for the lists
    if(plot_reference_result and reference_result_author!="Dairay et al."):
        # No reference result for the following plots so remove all the DNS data
        # and reset the which_lines_black_input and which_lines_dashed_input
        time_store.pop(0)
        kinetic_energy_store.pop(0)
        dissipation_store.pop(0)
        labels_store.pop(0)
        which_lines_dashed_input = []
        which_lines_black_input = []
        i_curve = 0 # reset the which_lines_black_input and which_lines_dashed_input
        for i in range(0,number_of_result_curves):
            if(black_line_flag[i]):
                which_lines_black_input.append(i_curve)
            if(dashed_line_flag[i]):
                which_lines_dashed_input.append(i_curve)
            i_curve += 1
        if(clr_input!=[]):
            clr_input_store.pop(0)
        if(mrkr_input!=[]):
            mrkr_input_store.pop(0)
        if(lnstl_input!=[]):
            lnstl_input_store.pop(0)

    if(plot_numerical_viscosity):
        qp.plotfxn(xdata=time_store,
                ydata=numerical_viscosity_store,
                # ylabel='$\\varepsilon=-\\frac{\\mathrm{d} K^{*}}{\\mathrm{d}t^{*}}$',
                ylabel='Nondimensional Numerical Viscosity, $\\frac{\\mathrm{Re}_{\\infty}\\varepsilon^{*}}{2\\zeta^{*}}$',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'numerical_viscosity_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[0.0,0.018],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input,
                legend_location="best")

    if(plot_dissipation_components):
        # vorticity component
        qp.plotfxn(xdata=time_store,
                ydata=vorticity_based_dissipation_store,
                ylabel='$\\varepsilon\\left(\\zeta^{*}\\right)$',
                xlabel='$t^{*}$',
                figure_filename=figure_subdirectory+'vorticity_based_dissipation_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[0,0.008],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)

        # Remove the reference result for the lists
        if(plot_reference_result and reference_result_author=="Dairay et al."):
            # No reference result for the following plots so remove all the DNS data
            # and reset the which_lines_black_input and which_lines_dashed_input
            time_store.pop(0)
            kinetic_energy_store.pop(0)
            dissipation_store.pop(0)
            vorticity_based_dissipation_store.pop(0)
            numerical_viscosity_store.pop(0)
            labels_store.pop(0)
            which_lines_dashed_input = []
            which_lines_black_input = []
            i_curve = 0 # reset the which_lines_black_input and which_lines_dashed_input
            for i in range(0,number_of_result_curves):
                if(black_line_flag[i]):
                    which_lines_black_input.append(i_curve)
                if(dashed_line_flag[i]):
                    which_lines_dashed_input.append(i_curve)
                i_curve += 1
            if(clr_input!=[]):
                clr_input_store.pop(0)
            if(mrkr_input!=[]):
                mrkr_input_store.pop(0)
            if(lnstl_input!=[]):
                lnstl_input_store.pop(0)

        # strain rate and pressure dilatation components
        qp.plotfxn(xdata=time_store,
                ydata=eps_S_plus_eps_p_store,
                ylabel='$\\varepsilon\\left(\\mathbf{S}^{d}\\right)$ + $\\epsilon\\left(p\\right)$',
                xlabel='$t^{*}=\\frac{tV_{\\infty}}{L}$',
                figure_filename=figure_subdirectory+'theoretical_dissipation_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[0,0.014],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)
        # strain rate component
        qp.plotfxn(xdata=time_store,
                ydata=eps_S_store,
                ylabel='$\\varepsilon\\left(\\mathbf{S}\\right)$',
                xlabel='$t^{*}$',
                figure_filename=figure_subdirectory+'strain_rate_dissipation_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[0,0.014],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)

        # deviatoric strain rate component
        qp.plotfxn(xdata=time_store,
                ydata=eps_S_store,
                ylabel='$\\varepsilon\\left(\\mathbf{S}^{d}\\right)$',
                xlabel='$t^{*}$',
                figure_filename=figure_subdirectory+'deviatoric_strain_rate_dissipation_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[0,0.014],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)

        # pressure dilatation component
        qp.plotfxn(xdata=time_store,
                ydata=eps_p_store,
                ylabel='$\\varepsilon\\left(p\\right)$',
                xlabel='$t^{*}$',
                figure_filename=figure_subdirectory+'pressure_dilatation_dissipation_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[-1e-1,1e-1],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)

        qp.plotfxn(xdata=time_store,
                ydata=eps_K_minus_eps_Sd_plus_eps_p_store,
                ylabel='$\\varepsilon^{*}-\\left[-\\varepsilon\\left(p\\right)+\\varepsilon\\left(\\mathbf{S}^{d}\\right)\\right]$',
                xlabel='$t^{*}$',
                figure_filename=figure_subdirectory+'error_in_KE_budget_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[-1e-1,1e-1],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)

        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/"
        filename=path_to_reference_result+"/"+"debonis_et_al_pressure_dilatation_kedr"+".txt"
        time, pressure_dilatation_based_dissipation = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
        # time_store.insert(0,time)
        # eps_K_minus_eps_Sd_store.insert(0,pressure_dilatation_based_dissipation)
        # labels_store.insert(0,"DRP $512^2$")
        time_store.append(time)
        eps_K_minus_eps_Sd_store.append(pressure_dilatation_based_dissipation)
        labels_store.append("DRP [DeBonis]\n $512^3$ DOFs")
        which_lines_black_input = [number_of_result_curves]
        which_lines_dashed_input = [number_of_result_curves]

        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/"
        filename=path_to_reference_result+"/"+"chapelier_2012_64p3_pressure_dilatation_kedr"+".txt"
        time, pressure_dilatation_based_dissipation = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
        # time_store.insert(0,time)
        # eps_K_minus_eps_Sd_store.insert(0,pressure_dilatation_based_dissipation)
        # labels_store.insert(0,"DRP $512^2$")
        time_store.append(time)
        eps_K_minus_eps_Sd_store.append(pressure_dilatation_based_dissipation)
        labels_store.append("DG [Chapelier et al.]\n $256^3$ DOFs ($64$p$3$)")
        which_lines_black_input = [number_of_result_curves,number_of_result_curves+1]
        which_lines_dotted_input = [number_of_result_curves+1]

        qp.plotfxn(xdata=time_store,
                ydata=eps_K_minus_eps_Sd_store,
                ylabel='$\\varepsilon^{*}-\\left[\\varepsilon\\left(\\mathbf{S}^{d}\\right)\\right]$',
                xlabel='$t^{*}$',
                figure_filename=figure_subdirectory+'observed_pressure_dissipation_rate_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                # ylimits=[-1e-1,1e-1],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                which_lines_dotted=which_lines_dotted_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)
