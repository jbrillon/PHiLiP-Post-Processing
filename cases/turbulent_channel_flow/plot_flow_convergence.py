#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
# Import personal libraries
# from finite_difference_library import first_derivative, fd_non_uniform_grid
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
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
# Input variables for plotting
#=====================================================

# filename=filesystem+"NarvalFiles/2023_JCP/"
# # filename+="flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512"
# # filename+="time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.34_procs512"
# # filename+="time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.15_procs512"
# filename+="time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.26_procs512"
# filename+="/turbulent_quantities.txt"
# time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
# kinetic_energy = np.loadtxt(filename,skiprows=1,usecols=1,dtype=np.float64,unpack=True)

x_store=[]
y_store=[]
labels_store=[]

# current
filename="/home/julien/Codes/2023-07-24/PHiLiP/build_release/tau_wall_first_run.txt"
time,tau_wall = np.loadtxt(filename,skiprows=1,usecols=(0,1),dtype=np.float64,unpack=True)
x_store.append(time)
y_store.append(tau_wall)
labels_store.append("$\\tau_{w}(t)$")

qp.plotfxn(x_store,y_store,
    figure_filename="flow_convergence",
    figure_size=(6,6),
    legend_labels_tex=labels_store,
    figure_filetype="pdf",
    title_label="Flow Satistical Convergence",
    xlabel="$t$",
    ylabel="Bulk Flow Property")