#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
# load submodules
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
# Check reference DNS
#-----------------------------------------------------
figure_directory_base="figures"

#-----------------------------------------------------
#                  DISSIPATION RATE
#-----------------------------------------------------
labels_store = []
time_store = []
dissipation_store = []
labels_store.append("Vermeire")
path_to_reference_result="/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/vermiere"
filename=path_to_reference_result+"/"+"dissipation"+".txt"
time, dissipation = np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
time_store.append(time)
dissipation_store.append(dissipation)

labels_store.append("Dairay et al.")
path_to_reference_result="/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data"
filename=path_to_reference_result+"/"+"TGV_Re1600"+".dat"
time, dissipation = np.loadtxt(filename,skiprows=43,dtype=np.float64,unpack=True,usecols=(0,2))
time_store.append(time)
dissipation_store.append(dissipation)

qp.plotfxn(xdata=time_store,
        ydata=dissipation_store,
        ylabel='Dissipation Rate',
        xlabel='Time',
        figure_filename='dissipation_rate_dns_ref',
        figure_filetype="pdf",
        title_label='Reference DNS Check',
        legend_labels_tex=labels_store,
        black_lines=False,
        xlimits=[0,20.0],
        #ylimits=[0.0,0.018],
        which_lines_dashed=[1])

#-----------------------------------------------------
#                   KINETIC ENERGY
#-----------------------------------------------------
labels_store = []
time_store = []
kinetic_energy_store = []
labels_store.append("Vermeire")
path_to_reference_result="/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/vermiere"
filename=path_to_reference_result+"/"+"kinetic_energy"+".txt"
time, kinetic_energy = np.loadtxt(filename,skiprows=0,dtype=np.float64,unpack=True)
time_store.append(time)
kinetic_energy_store.append(kinetic_energy)
# labels_store.append("Spectral DNS (P4, $260^{3}$ DOFs)")
# labels_store.append("DNS ($260^{3}$ DOFs)\n [Vermeire, 2014]")

labels_store.append("Dairay et al.")
path_to_reference_result="/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data"
filename=path_to_reference_result+"/"+"TGV_Re1600"+".dat"
time, kinetic_energy = np.loadtxt(filename,skiprows=43,dtype=np.float64,unpack=True,usecols=(0,1))
time_store.append(time)
kinetic_energy_store.append(kinetic_energy)

qp.plotfxn(xdata=time_store,
        ydata=kinetic_energy_store,
        ylabel='Kinetic Energy',
        xlabel='Time',
        figure_filename='kinetic_energy_dns_ref',
        figure_filetype="pdf",
        title_label='Reference DNS Check',
        fig_directory=figure_directory_base,
        legend_labels_tex=labels_store,
        black_lines=False,
        xlimits=[0,20.0],
        #ylimits=[0.0,0.018],
        which_lines_dashed=[1])


#-----------------------------------------------------
#                     ENSTROPHY
#-----------------------------------------------------
labels_store = []
time_store = []
enstrophy_store = []
labels_store.append("Vermeire")
path_to_reference_result="/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/vermiere"
filename=path_to_reference_result+"/"+"enstrophy"+".txt"
time, enstrophy = np.loadtxt(filename,skiprows=1,delimiter=",",dtype=np.float64,unpack=True)
time_store.append(time)
enstrophy_store.append(enstrophy)

labels_store.append("Dairay et al.")
path_to_reference_result="/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data"
filename=path_to_reference_result+"/"+"TGV_Re1600"+".dat"
time, enstrophy = np.loadtxt(filename,skiprows=43,dtype=np.float64,unpack=True,usecols=(0,4))
time_store.append(time)
enstrophy_store.append(enstrophy)

qp.plotfxn(xdata=time_store,
        ydata=enstrophy_store,
        ylabel='Enstrophy',
        xlabel='Time',
        figure_filename='enstrophy_dns_ref',
        figure_filetype="pdf",
        title_label='Reference DNS Check',
        legend_labels_tex=labels_store,
        black_lines=False,
        xlimits=[0,20.0],
        #ylimits=[0.0,0.018],
        which_lines_dashed=[1])