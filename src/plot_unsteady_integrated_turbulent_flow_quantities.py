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
import sys
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
# define functions
#-----------------------------------------------------
def get_spectra(kinetic_energy_input,dt):
    signal = 1.0*kinetic_energy_input
    fourier = np.fft.rfft(signal)
    n = signal.size
    sample_rate = dt
    print("dt: ")
    print(dt)
    freq = np.fft.rfftfreq(n, d=1./sample_rate)
    return [fourier,freq]
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
    transparent_legend_input=True):
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
    kinetic_energy_spectra_store = []
    wavenumber_spectra_store = []
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
        labels_store.append("DNS ($260^{3}$ DOFs)\n [Vermeire, 2014]")
        # DNS - dissipation
        filename=data_directory_base+"/"+"dns"+"/"+"dissipation"+".txt"
        time, dissipation = np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
        dissipation_store.append(dissipation)
        # black_line_flag.append(True) # inputs
        # dashed_line_flag.append(True) # inputs
        # if(black_line_flag[i_curve]):
        #   which_lines_black_input.append(i_curve)
        # if(dashed_line_flag[i_curve]):
        #   which_lines_dashed_input.append(i_curve)
        which_lines_black_input.append(i_curve)
        which_lines_dashed_input.append(i_curve)
        i_curve += 1
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
    
    #-----------------------------------------------------
    # evolution of kinetic energy:
    #-----------------------------------------------------
    # kolmogorov_slope = (-5.0/3.0)*time/10000.0+kinetic_energy[0]
    qp.plotfxn(xdata=time_store,#[time,time],
            ydata=kinetic_energy_store,#[kinetic_energy,kolmogorov_slope],
            ylabel='$K^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='$t^{*}$',
            figure_filename=figure_subdirectory+'kinetic_energy_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels_store,
            black_lines=False,
            xlimits=[0,20.0],
            ylimits=[0.01,0.14],
            log_axes=log_axes_input,
            which_lines_black=which_lines_black_input,
            which_lines_dashed=which_lines_dashed_input,
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(8,6),
            transparent_legend=transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False)
    #-----------------------------------------------------
    # evolution of dissipation energy:
    #-----------------------------------------------------
    if(plot_reference_result):
        # DNS - dissipation
        filename=path_to_reference_result+"/"+"dissipation"+".txt"
        time, dissipation = np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
        time_store[0] = time # replace it -- this is a hack
        # dissipation_store.append(dissipation)

    qp.plotfxn(xdata=time_store,
            ydata=dissipation_store,
            # ylabel='$\\varepsilon=-\\frac{\\mathrm{d} K^{*}}{\\mathrm{d}t^{*}}$',
            ylabel='$\\varepsilon=-\\mathrm{d} K^{*}/\\mathrm{d}t^{*}$',
            xlabel='$t^{*}$',
            figure_filename=figure_subdirectory+'dissipation_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels_store,
            black_lines=False,
            xlimits=[0,20.0],
            ylimits=[0.0,0.016],
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
            fig_directory=figure_directory_base)

    # PLOT KE SPECTRA HERE BEFORE REMOVING DNS

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

    for i in range(0,number_of_result_curves):
        kinetic_energy_current = kinetic_energy_store[i]
        dissipation_current = dissipation_store[i]
        time_current = time_store[i]
        index_of_max_dissipation = np.argmax(dissipation_current)
        timestep_current = np.mean(np.diff(time_current)) # average time step
        kinetic_energy_spectra_current, wavenumber_current = get_spectra(kinetic_energy_current[index_of_max_dissipation:],timestep_current)
        kinetic_energy_spectra_store.append(kinetic_energy_spectra_current)
        wavenumber_spectra_store.append(wavenumber_current)

    qp.plotfxn(xdata=time_store,
            ydata=enstrophy_store,
            ylabel='Nondimensional Enstrophy, $\\zeta^{*}$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'enstrophy_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels_store,
            black_lines=False,
            ylimits=[0,11],
            xlimits=[0,20],
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
            fig_directory=figure_directory_base)

# WARNING: If below is uncommented, add the argument: fig_directory=figure_directory_base
    # qp.plotfxn(xdata=time_store,
    #         ydata=vorticity_based_dissipation_store,
    #         ylabel='$\\varepsilon\\left(\\zeta^{*}\\right)$',
    #         xlabel='$t^{*}$',
    #         figure_filename=figure_subdirectory+'vorticity_based_dissipation_vs_time'+figure_filename_postfix,
    #         title_label=figure_title,
    #         markers=False,
    #         legend_labels_tex=labels_store,
    #         black_lines=False,
    #         xlimits=[0,20.0],
    #         ylimits=[0,0.008],
    #         log_axes=log_axes_input,
    #         which_lines_black=which_lines_black_input,
    #         which_lines_dashed=which_lines_dashed_input,
    #         legend_on=legend_on_input,
    #         legend_inside=legend_inside_input,
    #         nlegendcols=nlegendcols_input,
    #         figure_size=(8,6),
    #         transparent_legend=transparent_legend_input,
    #         legend_border_on=False,
    #         grid_lines_on=False,
            # fig_directory=figure_directory_base)

    # qp.plotfxn(xdata=wavenumber_spectra_store,
    #         ydata=kinetic_energy_spectra_store,
    #         ylabel='$E(\\kappa)$',
    #         xlabel='$\\kappa$',
    #         figure_filename=figure_subdirectory+'spectra'+figure_filename_postfix,
    #         title_label=figure_title,
    #         legend_labels_tex=labels_store,
    #         #xlimits=[1e0,1e3],
    #         ylimits=[1e-2,1e2],
    #         log_axes="both",
    #         which_lines_black=which_lines_black_input,
    #         which_lines_dashed=which_lines_dashed_input,
    #         legend_on=legend_on_input,legend_inside=legend_inside_input,nlegendcols=nlegendcols_input)

    # qp.plotfxn(xdata=time_store,
    #         ydata=eps_S_plus_eps_p_store,
    #         ylabel='$\\varepsilon\\left(\\mathbf{S}^{d}\\right)$ + $\\epsilon\\left(p\\right)$',
    #         xlabel='Nondimensional Time, $t^{*}=\\frac{tV_{\\infty}}{L}$',
    #         figure_filename=figure_subdirectory+'theoretical_dissipation_vs_time'+figure_filename_postfix,
    #         title_label=figure_title,
    #         markers=False,
    #         legend_labels_tex=labels_store,
    #         black_lines=False,
    #         ylimits=[0,0.014],
    #         log_axes=log_axes_input,
    #         which_lines_black=which_lines_black_input,
    #         which_lines_dashed=which_lines_dashed_input,
    #         legend_on=legend_on_input,legend_inside=legend_inside_input,nlegendcols=nlegendcols_input)

    # qp.plotfxn(xdata=time_store,
    #         ydata=eps_S_store,
    #         ylabel='$\\varepsilon\\left(\\mathbf{S}^{d}\\right)$',
    #         xlabel='Nondimensional Time, $t^{*}$',
    #         figure_filename=figure_subdirectory+'strain_rate_dissipation_vs_time'+figure_filename_postfix,
    #         title_label=figure_title,
    #         markers=False,
    #         legend_labels_tex=labels_store,
    #         black_lines=False,ylimits=[0,0.014],
    #         log_axes=log_axes_input,
    #         which_lines_black=which_lines_black_input,
    #         which_lines_dashed=which_lines_dashed_input,
    #         legend_on=legend_on_input,legend_inside=legend_inside_input,nlegendcols=nlegendcols_input)

    # qp.plotfxn(xdata=time_store,
    #         ydata=eps_p_store,
    #         ylabel='$\\varepsilon\\left(p\\right)$',
    #         xlabel='Nondimensional Time, $t^{*}$',
    #         figure_filename=figure_subdirectory+'pressure_dilatation_dissipation_vs_time'+figure_filename_postfix,
    #         title_label=figure_title,
    #         markers=False,
    #         legend_labels_tex=labels_store,
    #         black_lines=False,ylimits=[-1e-1,1e-1],
    #         log_axes=log_axes_input,
    #         which_lines_black=which_lines_black_input,
    #         which_lines_dashed=which_lines_dashed_input,
    #         legend_on=legend_on_input,legend_inside=legend_inside_input,nlegendcols=nlegendcols_input)

# #-----------------------------------------------------
# # SPECTRA CALCULATION TEST -- To be reviewed 
# #-----------------------------------------------------
# filename = "data/taylor_green_vortex/sip/p3_n16.txt"
# time, kinetic_energy = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
# figure_filename_postfix = "spectra_test"
# figure_subdirectory = 'taylor_green_vortex/TIMS/'

# def get_spectra(kinetic_energy_input,dt):
#     signal = 1.0*kinetic_energy_input
#     fourier = np.fft.rfft(signal)
#     n = signal.size
#     sample_rate = dt
#     print("dt: ")
#     print(dt)
#     freq = np.fft.rfftfreq(n, d=1./sample_rate)
#     return [fourier,freq]

# filename=data_directory_base+"/"+"dns"+"/"+"dissipation"+".txt"
# time_dns, dissipation_dns=np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
# # index_of_max_dissipation_dns = np.argmin(np.abs(time_per_block - desired_tau_value))
# index_of_max_dissipation_dns = np.argmax(dissipation_dns)
# time_of_max_dissipation_dns = time_dns[index_of_max_dissipation_dns]

# index_of_max_dissipation = np.argmin(np.abs(time - time_of_max_dissipation_dns))

# dt = time[1]-time[0]
# kinetic_energy_spectra, wavenumber = get_spectra(kinetic_energy[index_of_max_dissipation:],dt)

# wavenumber_store = []
# kinetic_energy_spectra_store = []
# wavenumber_store.append(wavenumber)
# wavenumber_store.append(wavenumber)
# kinetic_energy_spectra_store.append(kinetic_energy_spectra)
# kinetic_energy_spectra_store.append(kinetic_energy_spectra)

# plotfxn(xdata=wavenumber_store,
#         ydata=kinetic_energy_spectra_store,
#         ylabel='Kinetic Energy Spectra, $E(\\kappa)$',
#         xlabel='$\\kappa$',
#         figure_filename=figure_subdirectory+'spectra'+figure_filename_postfix,
#         title_label="Viscous TGV, $64^{3}$, P3, CFL=0.025",
#         legend_labels_tex=["Weak DG","Weak DG"],
#         #xlimits=[1e0,1e3],
#         #ylimits=[1e-15,1e1],
#         log_axes="both")
