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
    plot_enstrophy=True,
    plot_palinstrophy=True,
    clr_input=[],mrkr_input=[],lnstl_input=[],
    legend_fontSize_input=16,
    tmax=1.0,
    solid_and_dashed_lines=False,
    dashed_and_solid_lines=False,
    reference_result_author="Keetels et al.",
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
    which_lines_only_markers_input=[]
    # data store
    time_store = []
    kinetic_energy_store = []
    enstrophy_store = []
    palinstrophy_store = []
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
        if(reference_result_author=="Keetels et al."):
            # DNS - Kinetic Energy
            path_to_reference_result=CURRENT_PATH+"../cases/dipole_wall_collision/data/keetels"
            filename=path_to_reference_result+"/"+"kinetic_energy"+".txt"
            time, kinetic_energy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
            time_store.append(time)
            kinetic_energy_store.append(kinetic_energy)
            labels_store.append("Keetels et al.")
            filename=path_to_reference_result+"/"+"enstrophy"+".txt"
            time, enstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
            enstrophy_store.append(enstrophy)
            filename=path_to_reference_result+"/"+"palinstrophy"+".txt"
            time, palinstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
            palinstrophy_store.append(palinstrophy)
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
        which_lines_only_markers_input.append(i_curve)
        # which_lines_dashed_input.append(i_curve) # uncomment for dashed DNS result
        i_curve += 1
        if(clr_input!=[]):
            clr_input_store.append('k')
        if(mrkr_input!=[]):
            mrkr_input_store.append('o')
        if(lnstl_input!=[]):
            lnstl_input_store.append('None') # supported values are '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
    elif(plot_PHiLiP_DNS_result_as_reference):
        labels_store.append("DNS ($256^3$ DOFs, P$7$)")
        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/brillon"
        filename=path_to_reference_result+"/"+"turbulent_quantities_256dofs_p7"+".txt"
        time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        enstrophy_store.append(enstrophy)
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
        labels_store.append("Filtered DNS, $96^3$ P$2$")
        path_to_reference_result=CURRENT_PATH+"../cases/taylor_green_vortex/data/brillon/filtered_dns"
        filename=path_to_reference_result+"/"+"turbulent_quantities_96dofs_p2"+".txt"
        time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        enstrophy_store.append(enstrophy)
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
        time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, incompressible_kinetic_energy, incompressible_enstrophy, incompressible_palinstrophy = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        # store data
        time_store.append(time)
        domain_area = 4.0
        kinetic_energy_store.append(domain_area*incompressible_kinetic_energy)
        # -- store other quantities
        if(dissipation_rate_smoothing!=[]):
            # smooth this too for the numerical dissipation plots
            enstrophy = splrep(time,incompressible_enstrophy,k=5,s=0.0000001)
            enstrophy = splev(time,enstrophy)
            vorticity_based_dissipation = splrep(time,vorticity_based_dissipation,k=5,s=0.0000001)
            vorticity_based_dissipation = splev(time,vorticity_based_dissipation)
            incompressible_enstrophy = 1.0*enstrophy
        enstrophy_store.append(domain_area*incompressible_enstrophy)
        palinstrophy_store.append(domain_area*incompressible_palinstrophy)
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
                ylimits=[0.6,2.0],
                log_axes=None,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                which_lines_only_markers=which_lines_only_markers_input,
                legend_on=legend_on_input,
                legend_inside=legend_inside_input,
                nlegendcols=nlegendcols_input,
                figure_size=(6,6),
                transparent_legend=True,#transparent_legend_input,
                legend_border_on=False,
                grid_lines_on=False,
                clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
                legend_fontSize=legend_fontSize_input,
                legend_location="lower left",
                marker_size=3)#,
                # legend_anchor=[0.025,0.3])

    if(plot_reference_result and reference_result_author=="Keetels et al."):
        # enstrophy
        path_to_reference_result=CURRENT_PATH+"../cases/dipole_wall_collision/data/keetels"
        filename=path_to_reference_result+"/"+"enstrophy"+".txt"
        time, enstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
        time_store[0] = time # replace it -- this is a hack

    if(plot_enstrophy):
        qp.plotfxn(xdata=time_store,
                ydata=enstrophy_store,
                ylabel='Nondimensional Enstrophy, $\\zeta^{*}$',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'enstrophy_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                ylimits=[200,1800],
                xlimits=[0,tmax],
                log_axes=None,
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                which_lines_only_markers=which_lines_only_markers_input,
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
                legend_location="upper right",
                marker_size=3)
    
    if(plot_reference_result and reference_result_author=="Keetels et al."):
        # palinstrophy
        path_to_reference_result=CURRENT_PATH+"../cases/dipole_wall_collision/data/keetels"
        filename=path_to_reference_result+"/"+"palinstrophy"+".txt"
        time, palinstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
        time_store[0] = time # replace it -- this is a hack
    
    if(plot_palinstrophy):
        qp.plotfxn(xdata=time_store,
                ydata=palinstrophy_store,
                ylabel='Nondimensional Palinstrophy',
                xlabel='Nondimensional Time, $t^{*}$',
                figure_filename=figure_subdirectory+'palinstrophy_vs_time'+figure_filename_postfix,
                title_label=figure_title,
                markers=False,
                legend_labels_tex=labels_store,
                black_lines=False,
                ylimits=[1e5,1e8],
                xlimits=[0,tmax],
                log_axes="y",
                which_lines_black=which_lines_black_input,
                which_lines_dashed=which_lines_dashed_input,
                which_lines_only_markers=which_lines_only_markers_input,
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
                legend_location="upper right",
                marker_size=3)
