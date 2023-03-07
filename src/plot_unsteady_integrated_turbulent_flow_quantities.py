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
#-----------------------------------------------------
# define functions
#-----------------------------------------------------
def get_dissipation_discrete(time,kinetic_energy):
    return -np.gradient(kinetic_energy,time)
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
    path_to_reference_result,
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
    clr_input=[],mrkr_input=[],lnstl_input=[],
    legend_fontSize_input=16,
    tmax=20.0,
    solid_and_dashed_lines=False,
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
    eps_K_minus_eps_S_minus_eps_p_store = []
    eps_S_plus_eps_p_store = []
    eps_p_store = []
    eps_S_store = []
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
        # DNS - Kinetic Energy
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
        dissipation = get_dissipation_discrete(time,kinetic_energy)
        dissipation_store.append(dissipation)
        # -- store other quantities
        enstrophy_store.append(enstrophy)
        vorticity_based_dissipation_store.append(vorticity_based_dissipation)
        pressure_dilatation_based_dissipation_store.append(pressure_dilatation_based_dissipation)
        strain_rate_based_dissipation_store.append(strain_rate_based_dissipation)
        deviatoric_strain_rate_based_dissipation_store.append(deviatoric_strain_rate_based_dissipation)
        eps_K_minus_eps_S_minus_eps_p_store.append(dissipation - strain_rate_based_dissipation - pressure_dilatation_based_dissipation)
        eps_S_plus_eps_p_store.append(strain_rate_based_dissipation + pressure_dilatation_based_dissipation)
        eps_p_store.append(pressure_dilatation_based_dissipation)
        eps_S_store.append(strain_rate_based_dissipation)
        # store inputted line color, markers, and linestyles
        if(clr_input!=[]):
            clr_input_store.append(clr_input[i])
        if(mrkr_input!=[]):
            mrkr_input_store.append(mrkr_input[i])
        if(lnstl_input!=[]):
            lnstl_input_store.append(lnstl_input[i])

    # line parameters if doing solid and dashed lines
    if(solid_and_dashed_lines):
        clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        mrkr_input_store = ['None','None','None','None','None','None','None']
        lnstl_input_store = ['solid','dashed','solid','dashed','solid','dashed','solid']
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
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)
    
    #-----------------------------------------------------
    # evolution of kinetic energy dissipation rate:
    #-----------------------------------------------------
    if(plot_reference_result):
        # DNS - dissipation
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
                ylimits=[0.0,0.016],
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

    if(plot_reference_result):
        # DNS - enstrophy
        filename=path_to_reference_result+"/"+"enstrophy"+".txt"
        time, enstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
        time_store[0] = time # replace it -- this is a hack

    if(plot_enstrophy):
        # entrophy
        qp.plotfxn(xdata=time_store,
                ydata=enstrophy_store,
                ylabel='Nondimensional Enstrophy, $\\zeta^{*}$',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'enstrophy_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                ylimits=[0,12],
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
                legend_fontSize=legend_fontSize_input)

    # Remove the reference result for the lists
    if(plot_reference_result):
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
                figure_size=(8,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)

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
                ylimits=[0,0.014],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(8,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)
        # strain rate component
        qp.plotfxn(xdata=time_store,
                ydata=eps_S_store,
                ylabel='$\\varepsilon\\left(\\mathbf{S}^{d}\\right)$',
                xlabel='$t^{*}$',
                figure_filename=figure_subdirectory+'strain_rate_dissipation_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                xlimits=[0,tmax],
                ylimits=[0,0.014],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(8,6),
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
                ylimits=[-1e-1,1e-1],
                log_axes=log_axes_input,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(8,6),
                transparent_legend=transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                fig_directory=figure_directory_base,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input)
