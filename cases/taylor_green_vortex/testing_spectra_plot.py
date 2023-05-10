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


# velocity field
# file="/Users/Julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/flow_field_files/velocity_vorticity-0_reordered.dat"
file="/Users/Julien/NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/flow_field_files/velocity_vorticity-0_reordered.dat"
print("Generating spectra file from flow field file...")
print(" - Loading file: %s" % file)
data = np.loadtxt(file,skiprows=1,usecols=(3,4,5),dtype=np.float64)
print(" - done.")
velocity_field = [data[:,0],data[:,1],data[:,2]] # u,v,w
velocity_field = get_fluctuating_velocity_field(velocity_field)
# TurboGen no smoothing
spectra = get_tke_spectra(velocity_field,True)
file_out="/Users/Julien/NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra_no-smoothing.dat"
np.savetxt(file_out,spectra)
# FN
spectra = get_tke_spectra(velocity_field,False)
file_out="/Users/Julien/NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra_FN.dat"
np.savetxt(file_out,spectra)

x_store = []
y_store = []
labels_store = []

# x,y=np.loadtxt("/Users/Julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/flow_field_files/velocity_vorticity-0_reordered_spectra.dat",unpack=True)
# x_store.append(x);y_store.append(y);labels_store.append("TG-smoothing")
# x,y=np.loadtxt("/Users/Julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/flow_field_files/velocity_vorticity-0_reordered_spectra_no-smoothing.dat",unpack=True)
# x_store.append(x);y_store.append(y);labels_store.append("TG-no-smoothing")
# x,y=np.loadtxt("/Users/Julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/flow_field_files/velocity_vorticity-0_reordered_spectra_FN.dat",unpack=True)
# x_store.append(x);y_store.append(y);labels_store.append("FN")

x,y=np.loadtxt("/Users/Julien/NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra.dat",unpack=True)
x_store.append(x);y_store.append(y);labels_store.append("TG-smoothing")
x,y=np.loadtxt("/Users/Julien/NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra_no-smoothing.dat",unpack=True)
x_store.append(x);y_store.append(y);labels_store.append("TG-no-smoothing")
x,y=np.loadtxt("/Users/Julien/NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra_FN.dat",unpack=True)
x_store.append(x);y_store.append(y);labels_store.append("FN")


title_label="$%i^3$ DOFs, P7, TGV at $t=8.0$" % 64
qp.plotfxn(xdata=x_store,ydata=y_store,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
            title_label=title_label,
            fig_directory="figures",figure_filename="spectra_64p7_test",log_axes="both",figure_filetype="pdf",
            nlegendcols=1,
            xlimits=[2e0,10e2],ylimits=[1e-7,3e-2],
            markers=False,legend_on=True,legend_labels_tex=labels_store,
            which_lines_black=[0],
            # which_lines_markers=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,lnstl_input=['solid','dashed','dotted'],
            legend_fontSize=14,
            legend_location="upper left",
            legend_anchor=[0.025,0.3]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )