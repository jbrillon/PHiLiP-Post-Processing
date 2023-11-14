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
nCPUs_store=np.array([16,32,64,128,256,512])
number_of_elements_per_direction=16
poly_degree=5
dofs = number_of_elements_per_direction*(poly_degree+1)
base_directories=[\
"NarvalFiles/2023_JCP/strong_scaling_p5_96_runs_for_averaging/strong_scaling_TGV_p5_96_run_01/",\
"NarvalFiles/2023_JCP/strong_scaling_p5_96_runs_for_averaging/strong_scaling_TGV_p5_96_run_02/",\
"NarvalFiles/2023_JCP/strong_scaling_p5_96_runs_for_averaging/strong_scaling_TGV_p5_96_run_03/",\
"NarvalFiles/2023_JCP/strong_scaling_p5_96_runs_for_averaging/strong_scaling_TGV_p5_96_run_04/",\
"NarvalFiles/2023_JCP/strong_scaling_p5_96_runs_for_averaging/strong_scaling_TGV_p5_96_run_05/",\
]

n_sets_of_runs_for_averaging=len(base_directories)
n_data_points_per_run=len(nCPUs_store)
NSFR_store_cpu_time_per_step = np.zeros((n_sets_of_runs_for_averaging,n_data_points_per_run))

for i,base_dir in enumerate(base_directories):
    for j,nCPUs in enumerate(nCPUs_store):
        NSFR_job_name = "viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-%i_dofs0%i_p%i_CFL-0.1_procs%i" % (0,dofs,poly_degree,nCPUs)
        NSFR_cpu_time_per_step = np.loadtxt(filesystem+base_dir+NSFR_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
        NSFR_store_cpu_time_per_step[i,j] = NSFR_cpu_time_per_step

print(NSFR_store_cpu_time_per_step)

# average the values for all sets of runs
avg_NSFR_store_cpu_time_per_step = np.zeros(n_data_points_per_run)

for j in range(0,n_data_points_per_run):
    avg_NSFR_store_cpu_time_per_step[j] = np.average(NSFR_store_cpu_time_per_step[:,j])

#-----------------------------------------------------
# Plot 1
#-----------------------------------------------------
x_store = []
y_store = []
labels_store = []
# - Strong DG
x_store.append(nCPUs_store)
y_store.append(avg_NSFR_store_cpu_time_per_step)
labels_store.append("$c_{DG}$ NSFR.IR-GL")
title_label="TGV at Re$_{\\infty}=1600$ with $%i^{3}$ DOFs (P$%i$, $N_{el}=%i^3$)" % (dofs,poly_degree,number_of_elements_per_direction)
# figure_filename_input = "strong_scaling_p5_96_averaged"
# qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Number of CPUs",ylabel="CPU Time for One Time Step [s]",
#             title_label=title_label,
#             fig_directory="figures/2023_JCP",
#             figure_filename=figure_filename_input,
#             log_axes="both",figure_filetype=figure_filetype_input,
#             figure_size=(6,6),
#             nlegendcols=1,
#             # xlimits=[1e0,3e1],
#             # ylimits=[1e0,2e3],
#             markers=True,legend_on=True,legend_labels_tex=labels_store,
#             # which_lines_black=[0],
#             # which_lines_markers=[0],
#             transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
#             legend_fontSize=14,
#             # legend_location="upper left",
#             # legend_anchor=[0.025,0.3]
#             # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#             )
# figure_filename_input = "strong_scaling_p5_96_averaged_semilog"
# qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Number of CPUs",ylabel="CPU Time for One Time Step [s]",
#             title_label=title_label,
#             fig_directory="figures/2023_JCP",
#             figure_filename=figure_filename_input,
#             log_axes="x",figure_filetype=figure_filetype_input,
#             figure_size=(6,6),
#             nlegendcols=1,
#             # xlimits=[1e0,3e1],
#             # ylimits=[1e0,2e3],
#             markers=True,legend_on=True,legend_labels_tex=labels_store,
#             # which_lines_black=[0],
#             # which_lines_markers=[0],
#             transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
#             legend_fontSize=14,
#             # legend_location="upper left",
#             # legend_anchor=[0.025,0.3]
#             # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#             )
figure_filename_input = "strong_scaling_p5_96_averaged"
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Number of CPUs",ylabel="CPU Time for One Time Step [s]",
            title_label=title_label,
            fig_directory="figures/2023_JCP",
            figure_filename=figure_filename_input,
            log_axes=None,figure_filetype=figure_filetype_input,
            figure_size=(6,6),
            nlegendcols=1,
            # xlimits=[1e0,3e1],
            # ylimits=[1e0,2e3],
            markers=True,legend_on=True,legend_labels_tex=labels_store,
            # which_lines_black=[0],
            # which_lines_markers=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
            legend_fontSize=14,
            # legend_location="upper left",
            # legend_anchor=[0.025,0.3]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )
