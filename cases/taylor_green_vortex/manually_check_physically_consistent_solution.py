#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
filename=filesystem+"NarvalFiles/2023_JCP/"
# filename+="flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512"
# filename+="time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.34_procs512"
# filename+="time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.15_procs512"
filename+="time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.26_procs512"
filename+="/turbulent_quantities.txt"
# time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
kinetic_energy = np.loadtxt(filename,skiprows=1,usecols=1,dtype=np.float64,unpack=True)

print("Testing file: %s\n" % filename)

# Check that KE never increases
ke_previous = np.nan # init as nan
for step,ke in enumerate(kinetic_energy):
    if(ke > ke_previous):
        print("ERROR: Kinetic energy increased (from %1.15e to %1.15e) at step %i. Aborting...\n" % (ke_previous,ke,step))
        quit()
    else:
        ke_previous = ke
print("Kinetic energy did not increase.\n")
