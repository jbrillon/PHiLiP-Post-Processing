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
# from sys import platform
# if platform == "linux" or platform == "linux2":
#     # linux
#     filesystem="/media/julien/Samsung_T5/"
# elif platform == "darwin":
#     # OS X
#     filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
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

# og
filename="./turbulent_quantities-og.txt"
time = np.loadtxt(filename,skiprows=1,usecols=0,dtype=np.float64,unpack=True)
time_step = time[1:] - time[:-1]
x_store.append(time[:-1])
y_store.append(time_step/time_step[0])
labels_store.append("og")
# 19070503
filename="./turbulent_quantities-19070503.txt"
time = np.loadtxt(filename,skiprows=1,usecols=0,dtype=np.float64,unpack=True)
time_step = time[1:] - time[:-1]
x_store.append(time[:-1])
y_store.append(time_step/time_step[0])
labels_store.append("19070503")
# current
filename="./turbulent_quantities.txt"
time = np.loadtxt(filename,skiprows=1,usecols=0,dtype=np.float64,unpack=True)
time_step = time[1:] - time[:-1]
x_store.append(time[:-1])
y_store.append(time_step/time_step[0])
labels_store.append("current")

qp.plotfxn(x_store,y_store,
    figure_filename="time_step_strong_DG_DNS",
    figure_size=(6,6),
    legend_labels_tex=labels_store,
    figure_filetype="pdf",
    title_label="Time Step vs Time",
    xlabel="$t$",
    ylabel="$\\Delta t_{n}/\\Delta t_{0}$")