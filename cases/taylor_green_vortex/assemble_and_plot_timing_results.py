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
#=====================================================
# Helper functions
#=====================================================
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
base_directories=[\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_01/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_02/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_03/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_04/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_05/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_06/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_07/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_08/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_09/",\
"/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage_runs_for_averaging/cpu_time_advantage_run_10/"\
]


n_sets_of_runs_for_averaging=len(base_directories)
NSFR_n_poly_degree=len(NSFR_poly_degree)
std_sDG_n_poly_degree=len(std_sDG_poly_degree)
NSFR_store_cpu_time_per_step = np.zeros((n_sets_of_runs_for_averaging,NSFR_n_poly_degree))
std_sDG_store_cpu_time_per_step = np.zeros((n_sets_of_runs_for_averaging,std_sDG_n_poly_degree))

for i,base_dir in enumerate(base_directories):
    for j,p in enumerate(std_sDG_poly_degree):
        dofs = number_of_elements_per_direction*(p+1)
        oi = p+1 # overintegration
        std_sDG_job_name = "viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-%i_dofs0%i_p%i_CFL-0.1_procs%i" % (oi,dofs,p,nCPUs)
        std_sDG_cpu_time_per_step = np.loadtxt(base_dir+std_sDG_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
        std_sDG_store_cpu_time_per_step[i,j] = std_sDG_cpu_time_per_step

for i,base_dir in enumerate(base_directories):
    for j,p in enumerate(NSFR_poly_degree):
        dofs = number_of_elements_per_direction*(p+1)
        NSFR_job_name = "viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-%i_dofs0%i_p%i_CFL-0.1_procs%i" % (0,dofs,p,nCPUs)
        NSFR_cpu_time_per_step = np.loadtxt(base_dir+NSFR_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
        NSFR_store_cpu_time_per_step[i,j] = NSFR_cpu_time_per_step

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
labels_store.append("Strong DG-Roe-GL-OI")
# - NSFR
x_store.append(NSFR_poly_degree)
y_store.append(avg_NSFR_store_cpu_time_per_step)
labels_store.append("$c_{DG}$ NSFR.IR-GL")
title_label="TGV at Re$_{\\infty}=1600$ with $%i^3$ Elements on %i CPUs" % (4,nCPUs)
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Polynomial Degree",ylabel="CPU Time for One Time Step [s]",
            title_label=title_label,
            fig_directory="figures",
            figure_filename="cpu_timing_vs_poly_averaged_log",
            # figure_filename="cpu_timing_vs_poly",
            log_axes="both",figure_filetype="pdf",
            nlegendcols=1,
            # xlimits=[2e0,10e2],
            # ylimits=[1e2,1e5],
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
# Plot 2
#-----------------------------------------------------
x_store = []
y_store = []
labels_store = []
# - Strong DG
x_store.append(std_sDG_num_quad_int_strength)
y_store.append(avg_std_sDG_store_cpu_time_per_step)
labels_store.append("Strong DG-Roe-GL-OI")
# - NSFR
x_store.append(NSFR_num_quad_int_strength)
y_store.append(avg_NSFR_store_cpu_time_per_step)
labels_store.append("$c_{DG}$ NSFR.IR-GL")
title_label="TGV at Re$_{\\infty}=1600$ with $%i^3$ Elements on %i CPUs" % (4,nCPUs)
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Numerical Quad. Int. Strength",ylabel="CPU Time for One Time Step [s]",
            title_label=title_label,
            fig_directory="figures",
            figure_filename="cpu_timing_vs_nquad_averaged_2",
            # figure_filename="cpu_timing_vs_nquad",
            log_axes="y",figure_filetype="pdf",
            nlegendcols=1,
            # xlimits=[2e0,10e2],
            # ylimits=[1e2,1e5],
            markers=True,legend_on=True,legend_labels_tex=labels_store,
            # which_lines_black=[0],
            # which_lines_markers=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=True,#lnstl_input=['solid','dashed','dotted'],
            legend_fontSize=14,
            # legend_location="upper left",
            # legend_anchor=[0.025,0.3]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )