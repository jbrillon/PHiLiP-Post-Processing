#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
# load tools
sys.path.append(CURRENT_PATH+"../../src/tools");
from assemble_mpi_flow_field_files_and_reorder import assemble_mpi_flow_field_files_and_reorder
from generate_spectra_files import generate_spectra_file_from_flow_field_file
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
#=====================================================
# Helper functions
#=====================================================
#-----------------------------------------------------
def append_to_plot(x_,y_,label_):
    global x,y,labels
    labels.append(label_);x.append(x_);y.append(y_)
#-----------------------------------------------------
def batch_append_to_plot(paths_,labels_,filename):
    global x,y,labels
    for i,path in enumerate(paths_):
        spectra = np.loadtxt(filesystem+path+filename)
        append_to_plot(spectra[:,0],spectra[:,1],labels_[i])
#-----------------------------------------------------
def batch_plot_spectra(nDOF_,figure_filename_post_fix,batch_paths,batch_labels,
    solid_and_dashed_lines=False,
    title_off=False,
    figure_directory="figures",
    legend_fontSize_input=14,
    lnstl_input_store=['solid','solid','solid','solid','solid','solid','solid','solid'],
    plot_PHiLiP_DNS_result_as_reference=False,
    plot_reference_result=True,
    title_postfix_input=""):
    # TO DO: Move this function to its own file
    global x,y,labels
    x=[];y=[];labels=[];

    if(solid_and_dashed_lines):
        clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        mrkr_input_store = ['None','None','None','None','None','None','None']
        lnstl_input_store = ['solid','solid','dashed','solid','dashed','solid','dashed','solid']

    if(nDOF_=="all"):
        title_label = "TKE Spectra at $t^{*}=8.0$"
    else:
        title_label = "TKE Spectra at $t^{*}=8.0$ ($%i^3$ DOFs)" % nDOF_
    if(title_postfix_input!=""):
        title_label += title_postfix_input
    figure_filename = "spectra_t8_%s" % (str(nDOF_))
    if(figure_filename_post_fix!=""):
        figure_filename += "_%s" % (figure_filename_post_fix)
    x=[];y=[];labels=[];
    # reference result
    which_lines_black=[]
    if(plot_PHiLiP_DNS_result_as_reference):
        which_lines_black.append(0)
        filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256-0_reordered_spectra.dat"
        spectra = np.loadtxt(filepath_to_reference_result)
        append_to_plot(spectra[:,0],spectra[:,1],"DNS ($256^3$ DOFs, P$7$)")
    elif(plot_reference_result):
        which_lines_black.append(0)
        spectra = np.loadtxt(CURRENT_PATH+"data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
        append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
    batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra_no_smoothing.dat")
    if(title_off):
        title_label = " "
    if(solid_and_dashed_lines):
        qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
        title_label=title_label,
        fig_directory=figure_directory,figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
        nlegendcols=1,
        xlimits=[2e0,5.5e1],ylimits=[2.45e-5,3e-2],
        markers=False,legend_on=True,legend_labels_tex=labels,
        which_lines_black=which_lines_black,
        # which_lines_markers=[0],
        transparent_legend=True,legend_border_on=False,grid_lines_on=False,
        clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
        legend_fontSize=legend_fontSize_input,
        legend_location="upper left",
        legend_anchor=[0.0,0.45]
        # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
        )
    else:
        qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
            title_label=title_label,
            fig_directory=figure_directory,figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
            nlegendcols=1,
            # xlimits=[2e0,5.5e1],ylimits=[1.0e-5,3e-2],
            xlimits=[2.0e0,2.0e2],ylimits=[1.0e-9,5e-2],
            markers=False,legend_on=True,legend_labels_tex=labels,
            which_lines_black=which_lines_black,
            # which_lines_markers=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,lnstl_input=lnstl_input_store,
            legend_fontSize=legend_fontSize_input,
            # legend_location="upper left",
            # legend_anchor=[0.0,0.45]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )
    # qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    #     title_label=title_label,
    #     fig_directory=figure_directory,figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
    #     nlegendcols=1,
    #     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
    #     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
    #     markers=False,legend_on=True,legend_labels_tex=labels,
    #     which_lines_black=[0],
    #     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
    #     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    #     )
    if(nDOF_=="all"):
        title_label = "TKE Spectra at $t^{*}=9.0$"
    else:
        title_label = "TKE Spectra at $t^{*}=9.0$ ($%i^3$ DOFs)" % nDOF_
    if(title_postfix_input!=""):
        title_label += title_postfix_input
    figure_filename = "spectra_t9_%s" % (str(nDOF_))
    if(figure_filename_post_fix!=""):
        figure_filename += "_%s" % (figure_filename_post_fix)
    x=[];y=[];labels=[];
    # reference result
    if(plot_PHiLiP_DNS_result_as_reference):
        filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256-1_reordered_spectra.dat"
        spectra = np.loadtxt(filepath_to_reference_result)
        append_to_plot(spectra[:,0],spectra[:,1],"DNS ($256^3$ DOFs, P$7$)")
    elif(plot_reference_result):
        spectra = np.loadtxt(CURRENT_PATH+"data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
        append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
    # - results
    # same as above
    batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra_no_smoothing.dat")
    if(title_off):
        title_label = " "
    if(solid_and_dashed_lines):
        qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
        title_label=title_label,
        fig_directory=figure_directory,figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
        nlegendcols=1,
        xlimits=[2e0,7e1],ylimits=[2.45e-5,3e-2],
        markers=False,legend_on=True,legend_labels_tex=labels,
        which_lines_black=which_lines_black,
        # which_lines_markers=[0],
        transparent_legend=True,legend_border_on=False,grid_lines_on=False,
        clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
        legend_fontSize=legend_fontSize_input,
        legend_location="upper left",
        legend_anchor=[0.0,0.45]
        # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
        )
    else:
        qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
            title_label=title_label,
            fig_directory=figure_directory,figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
            nlegendcols=1,
            # xlimits=[2e0,7e1],ylimits=[2.45e-5,3e-2],
            xlimits=[2.0e0,2.0e2],ylimits=[1.0e-9,5e-2],
            markers=False,legend_on=True,legend_labels_tex=labels,
            which_lines_black=which_lines_black,
            transparent_legend=True,legend_border_on=False,grid_lines_on=False,lnstl_input=lnstl_input_store,
            legend_fontSize=legend_fontSize_input,
            # legend_location="upper left",
            # legend_anchor=[0.0,0.45]
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )
        # qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
        #     title_label=title_label,
        #     fig_directory=figure_directory,figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
        #     nlegendcols=1,
        #     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
        #     markers=False,legend_on=True,legend_labels_tex=labels,
        #     # which_lines_black=[0],
        #     transparent_legend=True,legend_border_on=False,grid_lines_on=False
        #     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
        #     )
    return

#=====================================================
# Global variables
#=====================================================
global x,y,labels
x=[];y=[];labels=[];
title_off_input=True
# fig_dir_input="figures"
# fig_dir_input="/Users/Julien/julien_phd/presentations/slides/wip_paper/20221205-AIAA/figures"
fig_dir_input="./figures/2023_JCP"
# =====================================================
# =====================================================
# =====================================================
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.18$", \
    "Smag. $C_{S}=0.10$", \
    "WALE $C_{W}=0.50$", \
    "VRMN $C_{V}=0.081$", \
    "SI.Smag. $C_{S}=0.10$", \
    ]
    batch_plot_spectra(96,"p5_sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.10$", \
    "WALE $C_{W}=0.50$", \
    "VRMN $C_{V}=0.081$", \
    "Smag.LRNC $C_{S}=0.10$", \
    "WALE.LRNC $C_{W}=0.10$", \
    "VRMN.LRNC $C_{V}=0.10$", \
    ]
    batch_plot_spectra(96,"p5_lrnc_sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG.LRNC_CLIPMC-0.01-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.10$", \
    "Smag.LRNC $C_{S}=0.10$", \
    "SI.Smag.LRNC $C_{S}=0.10$", \
    "HPF.Smag.LRNC $C_{S}=0.10$ $P_{L}=3$", \
    "HPF.SI.Smag.LRNC $C_{S}=0.10$ $P_{L}=3$", \
    "Dyn.Smag.LRNC ($P_{TF}=3$, $C_{max}=0.1$)", \
    ]
    batch_plot_spectra(96,"p5_lrnc_advanced_sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    # "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL2_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    # "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL4_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    # "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SS.VMS_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.10$", \
    #HPF.Smag. $C_{S}=0.10$ $P_{L}=2$", \
    "HPF.Smag. $C_{S}=0.10$ $P_{L}=3$", \
    #HPF.Smag. $C_{S}=0.10$ $P_{L}=4$", \
    "SI.Smag. $C_{S}=0.10$", \
    # "$SS.VMS $P_{L}=3$", \
    "HPF.SI.Smag. $C_{S}=0.10$ $P_{L}=3$", \
    "Dyn.Smag. ($P_{TF}=3$, $C_{max}=0.1$)", \
    "Smag.LRNC $C_{S}=0.10$", \
    ]
    batch_plot_spectra(96,"p5_advanced_sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs024_p5_procs16/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs012_p5_procs16/",\
    ]
    batch_labels = [ \
    "$96^{3}$ DOFs, P$5$, Strong DG-Roe-GL-OI", \
    "$64^{3}$ DOFs, P$7$, Strong DG-Roe-GL-OI", \
    "$48^{3}$ DOFs, P$5$, Strong DG-Roe-GL-OI",\
    "$24^{3}$ DOFs, P$5$, Strong DG-Roe-GL-OI",\
    "$12^{3}$ DOFs, P$5$, Strong DG-Roe-GL-OI",\
    ]
    # title_postfix_input=" using $c_{DG}$ NSFR.IR-GL"
    # batch_plot_spectra("all","cDG_NSFR_convergence",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(4,len(batch_paths)):
        figure_filename_postfix_input="strong_DG_convergence_%i"%i
        batch_plot_spectra("all",figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)

# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs012_p5_procs16/",\
    ]
    batch_labels = [ \
    "$96^{3}$ DOFs, P$5$, $c_{DG}$ NSFR.IR-GL", \
    "$64^{3}$ DOFs, P$7$, $c_{DG}$ NSFR.IR-GL", \
    "$48^{3}$ DOFs, P$5$, $c_{DG}$ NSFR.IR-GL",\
    "$24^{3}$ DOFs, P$5$, $c_{DG}$ NSFR.IR-GL",\
    "$12^{3}$ DOFs, P$5$, $c_{DG}$ NSFR.IR-GL",\
    ]
    # title_postfix_input=" using $c_{DG}$ NSFR.IR-GL"
    # batch_plot_spectra("all","cDG_NSFR_convergence",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(4,len(batch_paths)):
        figure_filename_postfix_input="cDG_NSFR_convergence_%i"%i
        batch_plot_spectra("all",figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GLL", \
    ]
    batch_plot_spectra(96,"p5_flux_nodes",batch_paths,batch_labels,solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GL-LxF", \
    "$c_{DG}$ NSFR.IR-GL-Roe", \
    "$c_{DG}$ NSFR.IR-GL-L2R", \
    ]
    # batch_plot_spectra(96,"upwind_gl",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(3,len(batch_paths)):
        figure_filename_postfix_input="p5_upwind_gl_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/two_point_flux/viscous_TGV_ILES_NSFR_cDG_KG_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/two_point_flux/viscous_TGV_ILES_NSFR_cDG_CH_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/two_point_flux/viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.KG-GL", \
    "$c_{DG}$ NSFR.CH-GL", \
    "$c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GL", \
    ]
    # batch_plot_spectra(96,"two_point_flux",batch_paths,batch_labels,solid_and_dashed_lines=True,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(3,len(batch_paths)):
        figure_filename_postfix_input="p5_two_point_flux_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/correction_parameter/viscous_TGV_ILES_NSFR_cSD_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/correction_parameter/viscous_TGV_ILES_NSFR_cHU_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{SD}$ NSFR.IR-GL", \
    "$c_{HU}$ NSFR.IR-GL", \
    "$c_{+}$ NSFR.IR-GL", \
    ]
    for i in range(3,len(batch_paths)):
        figure_filename_postfix_input="p5_correction_parameter_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/",\
    ]
    batch_labels = [ \
    "$256^3$ P$7$ $c_{DG}$ NSFR.IR-GL",\
    "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL",\
    "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL-Roe",\
    "$256^3$ P$3$ Strong DG-Roe-GL-OI",\
    ]
    lnstl_input=['solid','solid','solid','solid','dashed']
    for i in range(3,len(batch_paths)):# change to 0 if OK
        figure_filename_postfix_input="verification_%i"%i
        batch_plot_spectra(256,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,lnstl_input_store=lnstl_input,plot_PHiLiP_DNS_result_as_reference=False)
exit()
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/", \
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/", \
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI", \
    "$c_{DG}$ NSFR.IR-GL", \
    ]
    lnstl_input=['solid','solid']
    for i in range(1,len(batch_paths)):
        figure_filename_postfix_input="p3_OI_vs_SF_%i"%i
        batch_plot_spectra(256,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,lnstl_input_store=lnstl_input,plot_PHiLiP_DNS_result_as_reference=False,plot_reference_result=False)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    # "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL2_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    # "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL4_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.10$", \
    # "$c_{DG}$ NSFR.IR-GL-HPF.Smag. $C_{S}=0.10$ $P_{L}=2$", \
    "$c_{DG}$ NSFR.IR-GL-HPF.Smag. $C_{S}=0.10$ $P_{L}=3$", \
    # "$c_{DG}$ NSFR.IR-GL-HPF.Smag. $C_{S}=0.10$ $P_{L}=4$", \
    "$c_{DG}$ NSFR.IR-GL-SI.Smag. $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GL-HPF.SI.Smag. $C_{S}=0.10$ $P_{L}=3$", \
    "$c_{DG}$ NSFR.IR-GL-Dynamic.Smag. ($P_{TF}=3$, $C_{max}=0.1$)", \
    "$c_{DG}$ NSFR.IR-GL-Smag.LRNC $C_{S}=0.10$", \
    ]
    # batch_plot_spectra(96,"sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=True,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(0,len(batch_paths)):
        figure_filename_postfix_input="advanced_sgs_models_gl_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True,legend_fontSize_input=12)
exit()
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI", \
    "$c_{DG}$ NSFR.IR-GL", \
    ]
    for i in range(0,len(batch_paths)):
        figure_filename_postfix_input="p5_OI_vs_SF_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/over_integration_accuracy_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs096_p5_CFL-0.10_procs512/", \
    "NarvalFiles/2023_JCP/over_integration_accuracy_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-2_dofs096_p5_CFL-0.10_procs512/", \
    "NarvalFiles/2023_JCP/filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI-6", \
    "Strong DG-Roe-GL-OI-4", \
    "Strong DG-Roe-GL-OI-2", \
    "Strong DG-Roe-GL", \
    ]
    batch_plot_spectra(96,"overintegration_accuracy_strong_DG",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI", \
    "Strong DG-Roe-GL", \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GL-Roe", \
    ]
    batch_plot_spectra(96,"overintegration_accuracy",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL-Roe", \
    "Strong DG-Roe-GL-OI", \
    ]
    batch_plot_spectra(96,"cDG_roe_vs_sDG",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/", \
    "NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/", \
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI", \
    "$c_{DG}$ NSFR.IR-GL", \
    ]
    batch_plot_spectra(64,"high_poly_degree",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GLL", \
    "Strong DG-Roe-GL-OI", \
    "Strong DG-Roe-GLL-OI", \
    ]
    batch_plot_spectra(96,"flux_nodes",batch_paths,batch_labels,solid_and_dashed_lines=True,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GLL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GLL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GLL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GLL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GLL", \
    "$c_{DG}$ NSFR.IR-GLL-Smag. $C_{S}=0.18$", \
    "$c_{DG}$ NSFR.IR-GLL-Smag. $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GLL-WALE $C_{W}=0.50$", \
    "$c_{DG}$ NSFR.IR-GLL-VRMN $C_{V}=0.081$", \
    ]
    batch_plot_spectra(96,"sgs_models_gll",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    # "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512_filter_width_flad_and_gassner/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.18$", \
    "$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GL-WALE $C_{W}=0.50$", \
    "$c_{DG}$ NSFR.IR-GL-VRMN $C_{V}=0.081$", \
    # "$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$", \
    ]
    # batch_plot_spectra(96,"sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=True,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(0,len(batch_paths)):
        figure_filename_postfix_input="sgs_models_gl_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True,legend_fontSize_input=12)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GL-WALE $C_{W}=0.50$", \
    "$c_{DG}$ NSFR.IR-GL-VRMN $C_{V}=0.081$", \
    "$c_{DG}$ NSFR.IR-GL-Smag.LRNC $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GL-WALE.LRNC $C_{W}=0.50$", \
    "$c_{DG}$ NSFR.IR-GL-VRMN.LRNC $C_{V}=0.081$", \
    ]
    # batch_plot_spectra(96,"sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=True,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(0,len(batch_paths)):
        figure_filename_postfix_input="lrnc_sgs_models_gl_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True,legend_fontSize_input=12)
# =====================================================
if(True):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG.LRNC_CLIPMC-0.01-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GL-Smag.LRNC $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GL-SI.Smag.LRNC $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GL-HPF.Smag.LRNC $C_{S}=0.10$ $P_{L}=3$", \
    "$c_{DG}$ NSFR.IR-GL-HPF.SI.Smag.LRNC $C_{S}=0.10$ $P_{L}=3$", \
    "$c_{DG}$ NSFR.IR-GL-Dyn.Smag.LRNC ($P_{TF}=3$, $C_{max}=0.1$)", \
    ]
    # batch_plot_spectra(96,"sgs_models_gl",batch_paths,batch_labels,solid_and_dashed_lines=True,title_off=title_off_input,figure_directory=fig_dir_input)
    for i in range(0,len(batch_paths)):
        figure_filename_postfix_input="lrnc_advanced_sgs_models_gl_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True,legend_fontSize_input=12)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/upwind_dissipation_GLL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/upwind_dissipation_GLL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/upwind_dissipation_GLL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GLL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GLL", \
    "$c_{DG}$ NSFR.IR-GLL-LxF", \
    "$c_{DG}$ NSFR.IR-GLL-Roe", \
    "$c_{DG}$ NSFR.IR-GLL-L2R", \
    ]
    batch_plot_spectra(96,"upwind_gll",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512_filter_width_flad_and_gassner/", \
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI", \
    "Strong DG-Roe-GL-Smag. $C_{S}=0.18$", \
    "Strong DG-Roe-GL-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$", \
    ]
    batch_plot_spectra(96,"filter_width_stabilization",batch_paths,batch_labels,solid_and_dashed_lines=True,title_off=title_off_input,figure_directory=fig_dir_input)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.32_procs512/", \
    "NarvalFiles/2023_JCP/time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.24_procs512/", \
    # "NarvalFiles/2023_JCP/time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.14_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL, CFL=$0.10$",\
    "$c_{+}$ NSFR.IR-GL, CFL=$0.10$",\
    "$c_{+}$ NSFR.IR-GL, CFL=$0.32$",\
    "$c_{DG}$ NSFR.IR-GL, CFL=$0.24$",\
    #"Strong DG-Roe-GL-OI, CFL=$0.14$",\
    ]
    for i in range(0,len(batch_paths)):
        figure_filename_postfix_input="correction_parameter_cfl_advantage_%i"%i
        batch_plot_spectra(96,figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],solid_and_dashed_lines=False,title_off=title_off_input,figure_directory=fig_dir_input,plot_PHiLiP_DNS_result_as_reference=True)
# =====================================================
if(False):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64/", \
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs048_p5_procs64/", \
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs024_p5_procs16/", \
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs012_p5_procs16/", \
    ]
    batch_labels = [ \
    "$96^{3}$ DOFs, Strong DG-Roe-GL-OI", \
    "$48^{3}$ DOFs, Strong DG-Roe-GL-OI", \
    "$48^{3}$ DOFs, Strong DG-Roe-GL", \
    "$24^{3}$ DOFs, Strong DG-Roe-GL-OI", \
    "$12^{3}$ DOFs, Strong DG-Roe-GL-OI", \
    ]
    batch_plot_spectra("all","strong_DG_convergence",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)
# # =====================================================
# batch_paths = [ \
# "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
# "NarvalFiles/2023_JCP//", \
# "NarvalFiles/2023_JCP//", \
# "NarvalFiles/2023_JCP//", \
# ]
# batch_labels = [ \
# "$c_{DG}$ NSFR-GL", \
# "", \
# "", \
# "", \
# ]
# batch_plot_spectra(96,"insert_post_fix_here",batch_paths,batch_labels,title_off=title_off_input,figure_directory=fig_dir_input)