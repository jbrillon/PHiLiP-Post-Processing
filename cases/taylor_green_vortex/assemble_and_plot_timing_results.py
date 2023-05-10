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
# p_min=1
# p_max=15
# poly_degree=np.arange(p_min,p_max+1)
# poly_degree=np.array([1,2,3,4,5,6,7,9,10,11])
poly_degree=np.array([1,2,3,5,6,7,9,10,11]) # for nicer plot
number_of_elements_per_direction=4
base_directories=["/Users/Julien/julien_phd/cluster-scripts/outputs/jcp/cpu_time_advantage/"]# UPDATE THIS


n_sets_of_runs_for_averaging=len(base_directories)
n_poly_degree=len(poly_degree)
NSFR_store_cpu_time_per_step = np.zeros((n_sets_of_runs_for_averaging,n_poly_degree))
std_sDG_store_cpu_time_per_step = np.zeros((n_sets_of_runs_for_averaging,n_poly_degree))

# for i,base_dir in base_directories:
i=0
base_dir=base_directories[0]
for j,p in enumerate(poly_degree):
    dofs = number_of_elements_per_direction*(p+1)
    oi = p+1 # overintegration
    std_sDG_job_name = "viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-%i_dofs0%i_p%i_CFL-0.1_procs512" % (oi,dofs,p)
    NSFR_job_name = "viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-%i_dofs0%i_p%i_CFL-0.1_procs512" % (0,dofs,p)
    std_sDG_cpu_time_per_step = np.loadtxt(base_dir+std_sDG_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
    NSFR_cpu_time_per_step = np.loadtxt(base_dir+NSFR_job_name+"/timer_values.txt",usecols=(2),skiprows=1)
    print(NSFR_cpu_time_per_step)
    std_sDG_store_cpu_time_per_step[i,j] = std_sDG_cpu_time_per_step
    NSFR_store_cpu_time_per_step[i,j] = NSFR_cpu_time_per_step

# average the values for all sets of runs
avg_NSFR_store_cpu_time_per_step = np.zeros(n_poly_degree)
avg_std_sDG_store_cpu_time_per_step = np.zeros(n_poly_degree)

for j,p in enumerate(poly_degree):
    avg_std_sDG_store_cpu_time_per_step[j] = np.average(std_sDG_store_cpu_time_per_step[:,j])
    avg_NSFR_store_cpu_time_per_step[j] = np.average(NSFR_store_cpu_time_per_step[:,j])

x_store = []
y_store = []
labels_store = []

x_store.append(poly_degree)
y_store.append(avg_std_sDG_store_cpu_time_per_step)
labels_store.append("Strong DG-Roe-GL-OI")

x_store.append(poly_degree)
y_store.append(avg_NSFR_store_cpu_time_per_step)
labels_store.append("$c_{DG}$ NSFR.IR-GL")

title_label="TGV at Re$_{\\infty}=1600$ with $%i^3$ Elements" % 4
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Polynomial Degree",ylabel="CPU Time for One Time Step [s]",
            title_label=title_label,
            fig_directory="figures",figure_filename="cpu_timing",log_axes="y",figure_filetype="pdf",
            nlegendcols=1,
            # xlimits=[2e0,10e2],ylimits=[1e-7,3e-2],
            markers=True,legend_on=True,legend_labels_tex=labels_store,
            # which_lines_black=[0],
            # which_lines_markers=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,#lnstl_input=['solid','dashed','dotted'],
            legend_fontSize=14,
            # legend_location="upper left",
            # legend_anchor=[0.025,0.3]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )