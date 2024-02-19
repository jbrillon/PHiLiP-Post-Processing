#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
# load tools
sys.path.append(CURRENT_PATH+"../../src/tools");
from generate_spectra_files import get_tke_spectra, get_fluctuating_velocity_field
# load submodules
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
# figure_filetype_input='png'
figure_filetype_input='pdf'
nCPUs=16
p_min=1
p_max=30 #15
NSFR_poly_degree=np.arange(p_min,p_max+1)
NSFR_poly_degree=np.arange(p_min,p_max+1)
std_sDG_poly_degree=np.arange(p_min,p_max+1)
NSFR_num_quad_int_strength = NSFR_poly_degree+1
std_sDG_num_quad_int_strength = 2*(std_sDG_poly_degree+1)
# std_sDG_poly_degree=np.array([1,2,3,4,5,6,7,9,10,11]) # failed: 8,12,13,14,15
# note: as soon as nquad becomes >=16 it fails
# std_sDG_poly_degree=np.arange(p_min,p_max+1) 
number_of_elements_per_direction=4
# base_directories=["/Users/Julien/NarvalFiles/2023_JCP/cpu_time_advantage/"]
# base_directories=["/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_2/"]# UPDATE THIS
strong_DG_vs_NSFR=True
if(strong_DG_vs_NSFR):
    base_directories=[\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_01/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_02/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_03/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_04/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_05/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_06/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_07/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_08/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_09/",\
    "NarvalFiles/2023_JCP/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_10/"\
    ]
else:
    base_directories=[\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_01/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_02/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_03/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_04/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_05/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_06/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_07/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_08/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_09/",\
    "NarvalFiles/2023_JCP/cpu_time_coll_vs_uncoll_runs_for_averaging/cpu_time_coll_vs_uncoll_run_10/"\
    ]


n_sets_of_runs_for_averaging=len(base_directories)
NSFR_n_poly_degree=len(NSFR_poly_degree)
std_sDG_n_poly_degree=len(std_sDG_poly_degree)
NSFR_store_cpu_time_per_step = np.zeros((n_sets_of_runs_for_averaging,NSFR_n_poly_degree))
std_sDG_store_cpu_time_per_step = np.zeros((n_sets_of_runs_for_averaging,std_sDG_n_poly_degree))

if(strong_DG_vs_NSFR):
    for i,base_dir in enumerate(base_directories):
        for j,p in enumerate(std_sDG_poly_degree):
            dofs = number_of_elements_per_direction*(p+1)
            oi = p+1 # overintegration
            std_sDG_job_name = "viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-%i_dofs0%i_p%i_CFL-0.1_procs%i" % (oi,dofs,p,nCPUs)
            std_sDG_cpu_time_per_step = np.loadtxt(filesystem+base_dir+std_sDG_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
            std_sDG_store_cpu_time_per_step[i,j] = std_sDG_cpu_time_per_step
else:
    for i,base_dir in enumerate(base_directories):
        for j,p in enumerate(std_sDG_poly_degree):
            dofs = number_of_elements_per_direction*(p+1)
            std_sDG_job_name = "viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-%i_dofs0%i_p%i_CFL-0.1_procs%i" % (0,dofs,p,nCPUs)
            std_sDG_cpu_time_per_step = np.loadtxt(filesystem+base_dir+std_sDG_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
            std_sDG_store_cpu_time_per_step[i,j] = std_sDG_cpu_time_per_step

for i,base_dir in enumerate(base_directories):
    for j,p in enumerate(NSFR_poly_degree):
        dofs = number_of_elements_per_direction*(p+1)
        NSFR_job_name = "viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-%i_dofs0%i_p%i_CFL-0.1_procs%i" % (0,dofs,p,nCPUs)
        NSFR_cpu_time_per_step = np.loadtxt(filesystem+base_dir+NSFR_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
        NSFR_store_cpu_time_per_step[i,j] = NSFR_cpu_time_per_step

# if(strong_DG_vs_NSFR):
#     for i,base_dir in enumerate(base_directories):
#         for j,p in enumerate(std_sDG_poly_degree):
#             if(p==15):
#                 print("Run %i, CPU Time %1.5e" % (i+1,std_sDG_store_cpu_time_per_step[i,j])) 

# average the values for all sets of runs
avg_NSFR_store_cpu_time_per_step = np.zeros(NSFR_n_poly_degree)
avg_std_sDG_store_cpu_time_per_step = np.zeros(std_sDG_n_poly_degree)

for j,p in enumerate(std_sDG_poly_degree):
    avg_std_sDG_store_cpu_time_per_step[j] = np.average(std_sDG_store_cpu_time_per_step[:,j])

for j,p in enumerate(NSFR_poly_degree):
    avg_NSFR_store_cpu_time_per_step[j] = np.average(NSFR_store_cpu_time_per_step[:,j])

#-----------------------------------------------------
# Plot 1
#-----------------------------------------------------
x_store = []
y_store = []
labels_store = []
# - Strong DG
x_store.append(std_sDG_poly_degree)
y_store.append(avg_std_sDG_store_cpu_time_per_step)
if(strong_DG_vs_NSFR):
    figure_filename_input = "cpu_timing_vs_poly_averaged_log_0"
    labels_store.append("Strong DG-Roe-GL-OI")
else:
    figure_filename_input = "cpu_timing_coll_vs_uncoll_averaged_log_0"
    labels_store.append("$c_{DG}$ NSFR.IR-GLL")
# # - NSFR
# x_store.append(NSFR_poly_degree)
# y_store.append(avg_NSFR_store_cpu_time_per_step)
# labels_store.append("$c_{DG}$ NSFR.IR-GL")
title_label="TGV at Re$_{\\infty}=1600$ with $%i^3$ Elements on %i CPUs" % (4,nCPUs)
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Polynomial Degree, P",ylabel="CPU Time for One Time Step [s]",
            # title_label=title_label,
            fig_directory="figures/2023_JCP",
            figure_filename=figure_filename_input,
            # figure_filename="cpu_timing_vs_poly",
            log_axes="both",figure_filetype=figure_filetype_input,
            figure_size=(6,6),
            nlegendcols=1,
            xlimits=[1e0,3e1],
            ylimits=[1e0,2e3],
            markers=True,legend_on=True,legend_labels_tex=labels_store,
            # which_lines_black=[0],
            # which_lines_markers=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
            legend_fontSize=14,
            # legend_location="upper left",
            # legend_anchor=[0.025,0.3]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )

#-----------------------------------------------------
# Plot 1
#-----------------------------------------------------
x_store = []
y_store = []
labels_store = []
# - Strong DG
x_store.append(std_sDG_poly_degree)
y_store.append(avg_std_sDG_store_cpu_time_per_step)
if(strong_DG_vs_NSFR):
    figure_filename_input = "cpu_timing_vs_poly_averaged_log_1"
    labels_store.append("Strong DG-Roe-GL-OI")
else:
    figure_filename_input = "cpu_timing_coll_vs_uncoll_averaged_log_1"
    labels_store.append("$c_{DG}$ NSFR.IR-GLL")

# - NSFR
x_store.append(NSFR_poly_degree)
y_store.append(avg_NSFR_store_cpu_time_per_step)
labels_store.append("$c_{DG}$ NSFR.IR-GL")
title_label="TGV at Re$_{\\infty}=1600$ with $%i^3$ Elements on %i CPUs" % (4,nCPUs)

# reference curve
shift = 6.55
reference_curve_x = np.linspace(0.1,50.0,100)
order_for_ref_curve = 3.0+1.0
reference_curve = (reference_curve_x**(order_for_ref_curve))/np.exp(shift)
x_store.append(reference_curve_x)
y_store.append(reference_curve) # p^(d+1)
labels_store.append("P$^{d+1}$")
# to do: (1) See task (2) first; write these x and y data to files so that I can delete the entire cpu_time_advantage_runs_for_averaging directory and plot faster
#        (2) split this code into the assembly that generates the file for (1); and the other that takes that assembled file and plots it

# loglog
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Polynomial Degree, P",ylabel="CPU Time for One Time Step [s]",
            # title_label=title_label,
            fig_directory="figures/2023_JCP",
            figure_filename=figure_filename_input,
            # figure_filename="cpu_timing_vs_poly",
            log_axes="both",figure_filetype=figure_filetype_input,
            nlegendcols=1,
            figure_size=(6,6),
            xlimits=[1e0,3e1],
            ylimits=[1e0,2e3],
            # markers=True,
            legend_on=True,legend_labels_tex=labels_store,
            # which_lines_black=[0],
            which_lines_markers=[0,1],
            which_lines_dashed=[2],
            which_lines_black=[2],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
            legend_fontSize=14,
            # legend_location="upper left",
            # legend_anchor=[0.025,0.3]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )
# semilog
if(strong_DG_vs_NSFR):
    figure_filename_input = "cpu_timing_vs_poly_averaged_semilog_1"
else:
    figure_filename_input = "cpu_timing_coll_vs_uncoll_averaged_semilog_1"
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Polynomial Degree, P",ylabel="CPU Time for One Time Step [s]",
            # title_label=title_label,
            fig_directory="figures/2023_JCP",
            figure_filename=figure_filename_input,
            # figure_filename="cpu_timing_vs_poly",
            log_axes=None,figure_filetype=figure_filetype_input,
            nlegendcols=1,
            figure_size=(6,6),
            xlimits=[1e0,3e1],
            ylimits=[1e0,2e3],
            markers=True,legend_on=True,legend_labels_tex=labels_store,
            # which_lines_black=[0],
            # which_lines_markers=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
            legend_fontSize=14,
            # legend_location="upper left",
            # legend_anchor=[0.025,0.3]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )

if(strong_DG_vs_NSFR):
    #-----------------------------------------------------
    # Plot 2: ratio of cpu time
    #-----------------------------------------------------
    x_store = []
    y_store = []
    labels_store = []
    # x_store.append(NSFR_poly_degree)
    x_store.append(std_sDG_poly_degree)
    y_store.append(avg_std_sDG_store_cpu_time_per_step/avg_NSFR_store_cpu_time_per_step)
    figure_filename_input = "cpu_timing_vs_poly_averaged_ratio"
    labels_store.append("Strong DG-Roe-GL-OI / $c_{DG}$ NSFR.IR-GLL")

    title_label="TGV at Re$_{\\infty}=1600$ with $%i^3$ Elements on %i CPUs" % (4,nCPUs)

    # loglog
    qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Polynomial Degree, P",ylabel="Ratio of CPU Time for One Time Step",
                # title_label=title_label,
                fig_directory="figures/2023_JCP",
                figure_filename=figure_filename_input,
                # figure_filename="cpu_timing_vs_poly",
                log_axes=None,figure_filetype=figure_filetype_input,
                nlegendcols=1,
                figure_size=(6,6),
                xlimits=[1e0,3e1],
                # ylimits=[1e0,2e3],
                # markers=True,
                legend_on=True,legend_labels_tex=labels_store,
                # which_lines_black=[0],
                which_lines_markers=[0],
                # which_lines_dashed=[2],
                transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
                legend_fontSize=14,
                # legend_location="upper left",
                # legend_anchor=[0.025,0.3]
                # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
                )

# #-----------------------------------------------------
# # Plot 2
# #-----------------------------------------------------
# x_store = []
# y_store = []
# labels_store = []
# # - Strong DG
# x_store.append(std_sDG_num_quad_int_strength)
# y_store.append(avg_std_sDG_store_cpu_time_per_step)
# labels_store.append("Strong DG-Roe-GL-OI")
# # - NSFR
# x_store.append(NSFR_num_quad_int_strength)
# y_store.append(avg_NSFR_store_cpu_time_per_step)
# labels_store.append("$c_{DG}$ NSFR.IR-GL")
# title_label="TGV at Re$_{\\infty}=1600$ with $%i^3$ Elements on %i CPUs" % (4,nCPUs)
# qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Numerical Quad. Int. Strength",ylabel="CPU Time for One Time Step [s]",
#             title_label=title_label,
#             fig_directory="figures/2023_JCP",
#             figure_filename="cpu_timing_vs_nquad_averaged_2",
#             # figure_filename="cpu_timing_vs_nquad",
#             log_axes="y",figure_filetype="pdf",
#             nlegendcols=1,
#             # xlimits=[2e0,10e2],
#             # ylimits=[1e2,1e5],
#             markers=True,legend_on=True,legend_labels_tex=labels_store,
#             # which_lines_black=[0],
#             # which_lines_markers=[0],
#             transparent_legend=True,legend_border_on=False,grid_lines_on=True,#lnstl_input=['solid','dashed','dotted'],
#             legend_fontSize=14,
#             # legend_location="upper left",
#             # legend_anchor=[0.025,0.3]
#             # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#             )