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
    plotting_subsonic_result=False):
    
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
    #-----------------------------------------------------
    if(plotting_subsonic_result):
        time, kinetic_energy = np.loadtxt("./data/chapelier2024/subsonic/kinetic_energy.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    else:
        time, kinetic_energy = np.loadtxt("./data/chapelier2024/kinetic_energy.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    time_store.append(time)
    kinetic_energy_store.append(kinetic_energy)
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
    labels.append("Chapelier et al. ($2048^3$ DOFs)")
    black_line_flag.append(True)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    clr_input_store = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    mrkr_input_store = ['o','None','None','None','None','None','None','None']
    lnstl_input_store = ['None','solid','solid','solid','solid','dashed','solid','dashed','solid']
    if(plotting_subsonic_result):
        lnstl_input_store = ['None','solid','dashed','solid','solid','dashed','solid','dashed','solid']
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
        filename = data_directory_base+"/"+subdirectories[i]+"/"+filenames[i]
        if(i==2):
            time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation, corrected_pressure_dilatation_based_dissipation, corrected_dilatational_dissipation, uncorrected_pressure_dilatation_based_dissipation, uncorrected_dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
            pressure_dissipation_store.append(corrected_pressure_dilatation_based_dissipation)
            dilatational_dissipation_store.append(uncorrected_dilatational_dissipation)
            print(dilatational_dissipation[-1])
            print(uncorrected_dilatational_dissipation[-1])
            print(corrected_dilatational_dissipation[-1])
        else:
            time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
            pressure_dissipation_store.append(pressure_dilatation_based_dissipation)
            dilatational_dissipation_store.append(dilatational_dissipation)
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        solenoidal_dissipation_store.append(solenoidal_dissipation)
        # dilatational_dissipation_calc = 2.0*(strain_rate_based_dissipation - deviatoric_strain_rate_based_dissipation)
        # print(np.linalg.norm(dilatational_dissipation_calc-dilatational_dissipation))
        
        

    final_time_for_plot = 20.0
    if(plotting_subsonic_result):
        final_time_for_plot = 10.0

    ylimits_for_plot = [0.0,0.14]
    if(plotting_subsonic_result):
        ylimits_for_plot = [0.04,0.13]
    '''
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
    ylimits_for_plot = [0.0,0.016]
    if(plotting_subsonic_result):
        ylimits_for_plot = [0.0,0.012]
    time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/solenoidal_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    if(plotting_subsonic_result):
        time, solenoidal_dissipation = np.loadtxt("./data/chapelier2024/subsonic/solenoidal_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    time_store[0] = time # replace it -- this is a hack
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
            legend_location="best")
    '''
    #-----------------------------------------------------
    ylimits_for_plot = [0.0,0.002]
    if(plotting_subsonic_result):
        # ylimits_for_plot = []
        ylimits_for_plot = [0.0,5.0e-6]
    time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/dilatational_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    if(plotting_subsonic_result):
        time, dilatational_dissipation = np.loadtxt("./data/chapelier2024/subsonic/dilatational_dissipation.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")    
    time_store[0] = time # replace it -- this is a hack
    qp.plotfxn(xdata=time_store,
            ydata=dilatational_dissipation_store,
            ylabel='Nondimensional Dilational Dissipation, $\\varepsilon_{d}^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'dilational_dissipation_vs_time'+figure_filename_postfix,
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
            legend_border_on=True,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")

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
    figure_filename_postfix = "_128"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512",\
    # "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p15_procs128",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $16$p$7$\n ($128^3$ DOF) CFL=0.1",\
    "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$3$\n ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$7$\n ($256^3$ DOF) CFL=0.1",\
    "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $8$p$15$\n ($128^3$ DOF) CFL=0.01",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot)

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
    "subsonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512_no_limiter",\
    "subsonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512_corrected_quantities",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$7$\n ($256^3$ DOF) CFL=0.1",\
    "$c_{DG}$ NSFR.CH$_{RA}$+Roe $32$p$7$\n ($256^3$ DOF) CFL=0.1",\
    "with correction",\
    # "$c_{+}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$3$\n ($128^3$ DOF) CFL=0.1",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $8$p$15$\n ($128^3$ DOF) CFL=0.01",\
    # "$c_{DG}$ NSFR.CH$_{RA}$+Roe+PPL $32$p$7$\n ($256^3$ DOF) CFL=0.1",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,True,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,plotting_subsonic_result=True)

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
