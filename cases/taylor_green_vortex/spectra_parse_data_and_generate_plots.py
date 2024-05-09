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
# for testing quickplotlib changes; uncomment the lines below and comment the line above
# from sys import platform
# if platform == "linux" or platform == "linux2":
#     # linux
#     sys.path.append("/home/julien/Codes/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes
# elif platform == "darwin":
#     # OS X
#     sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes

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
def get_grid_cutoff_wavenumber(number_of_elements_per_direction):
    grid_cutoff_wavenumber = 0.5*number_of_elements_per_direction
    return grid_cutoff_wavenumber
#-----------------------------------------------------
def get_cutoff_wavenumber(poly_degree,number_of_elements_per_direction,truncate_spectra_at_effective_DOFs):
    nDOF = (poly_degree+1)*number_of_elements_per_direction
    effective_nDOF = (poly_degree)*number_of_elements_per_direction
    cutoff_wavenumber = 0.5*effective_nDOF
    if(truncate_spectra_at_effective_DOFs==False):
        cutoff_wavenumber = 0.5*nDOF
    return cutoff_wavenumber
#-----------------------------------------------------
def get_truncated_spectra_from_cutoff_wavenumber_and_spectra(spectra, cutoff_wavenumber):
    idx = (np.abs(spectra[:,0] - cutoff_wavenumber)).argmin()
    return spectra[:(idx+1),:]
#-----------------------------------------------------
def get_truncated_spectra_from_DOFs_information(spectra, poly_degree, number_of_elements_per_direction,truncate_spectra_at_effective_DOFs):
    cutoff_wavenumber = get_cutoff_wavenumber(poly_degree,number_of_elements_per_direction,truncate_spectra_at_effective_DOFs)
    idx = (np.abs(spectra[:,0] - cutoff_wavenumber)).argmin()
    return spectra[:(idx+1),:]
#-----------------------------------------------------
def append_to_plot(x_,y_,label_):
    global x,y,labels
    labels.append(label_);x.append(x_);y.append(y_)
#-----------------------------------------------------
def batch_append_to_plot(paths_,labels_,filename,list_of_poly_degree_,list_of_number_of_elements_per_direction_,truncate_spectra_at_effective_DOFs):
    global x,y,labels
    for i,path in enumerate(paths_):
        spectra_ = np.loadtxt(filesystem+path+filename)
        poly_degree = list_of_poly_degree_[i]
        number_of_elements_per_direction = list_of_number_of_elements_per_direction_[i]
        spectra = get_truncated_spectra_from_DOFs_information(spectra_, poly_degree, number_of_elements_per_direction,truncate_spectra_at_effective_DOFs)
        # spectra = 1.0*spectra_ # uncomment for no truncation
        append_to_plot(spectra[:,0],spectra[:,1],labels_[i])
#-----------------------------------------------------
def batch_plot_spectra(nDOF_,figure_filename_post_fix,batch_paths,batch_labels,
    solid_and_dashed_lines=False,
    dashed_and_solid_lines=False,
    title_off=False,
    figure_directory="figures",
    legend_fontSize_input=14,
    plot_PHiLiP_DNS_result_as_reference=False,
    plot_reference_result=True,
    plot_filtered_dns=False,
    title_postfix_input="",
    plot_zoomed_section=False,
    which_lines_dashed=[],
    x_limits_zoom_input=[],
    y_limits_zoom_input=[],
    plot_cutoff_wavenumber_asymptote=False,
    list_of_poly_degree_input=[],
    list_of_number_of_elements_per_direction_input=[],
    truncate_spectra_at_effective_DOFs=True): # Modify this here
    # TO DO: Move this function to its own file
    global x,y,labels
    x=[];y=[];labels=[];
    if(solid_and_dashed_lines):
        clr_input_store = ['tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        mrkr_input_store = ['None','None','None','None','None','None','None']
        lnstl_input_store = ['solid','dashed','solid','dashed','solid','dashed','solid']
    elif(dashed_and_solid_lines):
        clr_input_store = ['tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green','tab:orange','tab:orange']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        mrkr_input_store = ['None','None','None','None','None','None','None','None']
        lnstl_input_store = ['dashed','solid','dashed','solid','dashed','solid','dashed','solid']
    else:
        clr_input_store = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        mrkr_input_store = []
        lnstl_input_store = ['solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']

    cutoff_wavenumber_store = []
    grid_cutoff_wavenumber_store = []
    if(plot_cutoff_wavenumber_asymptote):
        cutoff_wavenumber_store = []
        grid_cutoff_wavenumber_store = []
        if(list_of_poly_degree_input==[] or list_of_number_of_elements_per_direction_input==[]):
            print("ERROR: list_of_poly_degree_input or list_of_number_of_elements_per_direction_input is empty with plot_cutoff_wavenumber_asymptote==True.")
            print("Aborting...")
            exit()
        for i,poly_degree in enumerate(list_of_poly_degree_input):
            number_of_elements_per_direction = list_of_number_of_elements_per_direction_input[i]
            cutoff_wavenumber = get_cutoff_wavenumber(poly_degree,number_of_elements_per_direction,True) # True for effective DOFs
            grid_cutoff_wavenumber = get_grid_cutoff_wavenumber(number_of_elements_per_direction)
            if((cutoff_wavenumber in cutoff_wavenumber_store)==False):
                cutoff_wavenumber_store.append(cutoff_wavenumber)
            if((grid_cutoff_wavenumber in grid_cutoff_wavenumber_store)==False):
                grid_cutoff_wavenumber_store.append(grid_cutoff_wavenumber)

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
    i_curve=0
    if(plot_PHiLiP_DNS_result_as_reference):
        # which_lines_black.append(0)
        clr_input_store.insert(i_curve,"k")
        filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256-0_reordered_spectra.dat"
        spectra = np.loadtxt(filepath_to_reference_result)
        append_to_plot(spectra[:,0],spectra[:,1],"DNS ($256^3$ DOFs, P$7$)")
        i_curve += 1
    elif(plot_reference_result):
        # which_lines_black.append(0)
        clr_input_store.insert(i_curve,"k")
        spectra = np.loadtxt(CURRENT_PATH+"data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
        append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
        i_curve += 1
    # batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra_no_smoothing.dat")

    if(plot_PHiLiP_DNS_result_as_reference or plot_reference_result):
        if(solid_and_dashed_lines or dashed_and_solid_lines):
            mrkr_input_store.insert(i_curve-1,'None')
            lnstl_input_store.insert(i_curve-1,'solid')

    if(plot_filtered_dns):
        # filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256_projected_to_p2_dofs096-0_reordered_spectra_oversampled_nquad12.dat"
        # spectra = np.loadtxt(filepath_to_reference_result)
        # append_to_plot(spectra[:,0],spectra[:,1],"Projected DNS ($96^3$ DOFs, P$2$)")
        clr_input_store.insert(i_curve,"k")
        # which_lines_black.append(i_curve)
        if(solid_and_dashed_lines or dashed_and_solid_lines):
            mrkr_input_store.insert(i_curve,'None')
        lnstl_input_store.insert(i_curve,'dashed') # for dashed filtered DNS result
        # which_lines_dashed.append() 
        i_curve += 1

    # compute reference curve 1
    index_of_reference_curve = len(batch_labels)+i_curve
    x_ref_curve = np.linspace(2.0e0,2.0e2,100)
    order_for_ref_curve = -5.0/3.0
    ref_curve_label = "$\\left(k^{*}\\right)^{-5/3}$"
    shift = 2.0
    y_ref_curve = (x_ref_curve**(order_for_ref_curve))/np.exp(shift)
    clr_input_store.insert(index_of_reference_curve,"k")#"tab:gray"
    lnstl_input_store.insert(index_of_reference_curve,"dotted")

    # add reference curve
    # append_to_plot(x_ref_curve,y_ref_curve,ref_curve_label)
    if(nDOF_==256):
        x_limits_zoom=[30, 60]
        y_limits_zoom=[1.0e-4, 3.0e-4]
    elif(nDOF_==96):
        x_limits_zoom=[20, 40]
        y_limits_zoom=[3.5e-5, 5.0e-4]
    elif(nDOF_=="all"):
        x_limits_zoom=[30, 60]
        y_limits_zoom=[1.0e-4, 3.0e-4]
    # set manually if desired
    if(x_limits_zoom_input!=[]):
       x_limits_zoom=x_limits_zoom_input
    if(y_limits_zoom_input!=[]):
       y_limits_zoom=y_limits_zoom_input

    if(title_off):
        title_label = " "

    if(which_lines_dashed!=[]):
        for i in which_lines_dashed:
            if(plot_PHiLiP_DNS_result_as_reference or plot_reference_result):
                lnstl_input_store[i+1] = "dashed"
            else:
                lnstl_input_store[i] = "dashed"

    # UNCOMMENT FOR t=8
    # qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    #     title_label=title_label,
    #     fig_directory=figure_directory,figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    #     nlegendcols=1,
    #     # xlimits=[2e0,5.5e1],ylimits=[1.0e-5,3e-2],
    #     xlimits=[2.0e0,2.0e2],ylimits=[1.0e-8,5e-2],
    #     markers=False,legend_on=True,legend_labels_tex=labels,
    #     which_lines_black=which_lines_black,
    #     # which_lines_markers=[0],
    #     transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    #     clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
    #     legend_fontSize=legend_fontSize_input,
    #     # legend_location="upper left",
    #     # legend_anchor=[0.0,0.45]
    #     # which_lines_only_markers=[1,2,3],
    #     which_lines_dashed=which_lines_dashed,
    #     plot_zoomed_section=plot_zoomed_section,
    #     x_limits_zoom=x_limits_zoom,y_limits_zoom=y_limits_zoom,
    #     zoom_box_origin_and_extent=[0.65, 0.65, 0.32, 0.32],
    #     vertical_lines=cutoff_wavenumber_store,
    #     secondary_vertical_lines=grid_cutoff_wavenumber_store
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
    i_curve=0
    if(plot_PHiLiP_DNS_result_as_reference):
        filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256-1_reordered_spectra_oversampled_nquad16.dat"
        spectra = np.loadtxt(filepath_to_reference_result)
        append_to_plot(spectra[:,0],spectra[:,1],"DNS ($256^{3}$p$7$)")
        i_curve += 1
    elif(plot_reference_result):
        spectra = np.loadtxt(CURRENT_PATH+"data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
        append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
        i_curve += 1
    if(plot_filtered_dns):
        filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256_projected_to_p2_dofs096-1_reordered_spectra_oversampled_nquad12.dat"
        spectra_ = np.loadtxt(filepath_to_reference_result)
        # spectra = get_truncated_spectra_from_DOFs_information(spectra_, 2, 32)
        spectra = get_truncated_spectra_from_DOFs_information(spectra_, 5, 16, truncate_spectra_at_effective_DOFs) # to match cut-off for 96P5
        append_to_plot(spectra[:,0],spectra[:,1],"Projected DNS ($96^{3}$p$2$)")
        # clr_input_store.insert(i_curve,"k")
        # which_lines_black.append(i_curve)
        # mrkr_input_store.insert(i_curve,'None')
        # lnstl_input_store.insert(i_curve,'dashed') # for dashed filtered DNS result
        # which_lines_dashed.append() 
        i_curve += 1
        # TO DO MAKE THIS DASHED
    # - results
    # same as above
    batch_append_to_plot(batch_paths, batch_labels, 
        "flow_field_files/velocity_vorticity-1_reordered_spectra_no_smoothing.dat",
        list_of_poly_degree_input, list_of_number_of_elements_per_direction_input,truncate_spectra_at_effective_DOFs)
    # adjust the reference curve shift
    shift = 2.2
    y_ref_curve = (x_ref_curve**(order_for_ref_curve))/np.exp(shift)
    # add reference curve
    append_to_plot(x_ref_curve,y_ref_curve,ref_curve_label)
    if(title_off):
        title_label = " "

    qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
        title_label=title_label,
        fig_directory=figure_directory,figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
        nlegendcols=1,
        # xlimits=[2e0,7e1],ylimits=[2.45e-5,3e-2],
        xlimits=[2.0e0,112],ylimits=[4.0e-7,5e-2],
        markers=False,legend_on=True,legend_labels_tex=labels,
        which_lines_black=which_lines_black,
        transparent_legend=False,legend_border_on=False,grid_lines_on=False,
        clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
        legend_fontSize=legend_fontSize_input,
        legend_location="lower left",
        # legend_anchor=[0.0,0.45]
        # which_lines_only_markers=[1,2,3],
        # which_lines_dashed=which_lines_dashed,
        plot_zoomed_section=plot_zoomed_section,
        x_limits_zoom=x_limits_zoom,y_limits_zoom=y_limits_zoom,
        zoom_box_origin_and_extent=[0.65, 0.65, 0.32, 0.32],
        vertical_lines=cutoff_wavenumber_store,
        secondary_vertical_lines=grid_cutoff_wavenumber_store
        )
    
    return

#=====================================================
# Global variables
#=====================================================
global x,y,labels
x=[];y=[];labels=[];
regenerate_all_plots=False
title_off_input=True
# fig_dir_input="figures"
# fig_dir_input="/Users/Julien/julien_phd/presentations/slides/wip_paper/20221205-AIAA/figures"
fig_dir_input="./figures/2023_JCP/oversampled_spectra"
# =====================================================
# =====================================================
# =====================================================
# =====================================================

# =====================================================
# MISSING FIG 3
# =====================================================

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/",\
    ]
    batch_labels = [ \
    # "$96^{3}$ DOFs, P$5$", \
    # "$64^{3}$ DOFs, P$7$", \
    # "$48^{3}$ DOFs, P$5$",\
    # "$24^{3}$ DOFs, P$5$",\
    "$96^{3}$p$5$", \
    "$64^{3}$p$7$", \
    "$48^{3}$p$5$",\
    "$24^{3}$p$5$",\
    ]
    list_of_poly_degree=[5,7,5,5]
    list_of_number_of_elements_per_direction=[16,8,8,4]
    
    batch_plot_spectra("all","cDG_NSFR_convergence",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64/",\
    ]
    batch_labels = [ \
    "$96^{3}$p$5$, $c_{DG}$", \
    "$96^{3}$p$5$, sDG", \
    "$64^{3}$p$7$, $c_{DG}$", \
    "$64^{3}$p$7$, sDG", \
    "$48^{3}$p$5$, $c_{DG}$",\
    "$48^{3}$p$5$, sDG",\
    ]
    list_of_poly_degree=[5,5,7,7,5,5]
    list_of_number_of_elements_per_direction=[16,16,8,8,8,8]
    
    batch_plot_spectra("all","convergence_comparison",batch_paths,batch_labels,
        solid_and_dashed_lines=True,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-3_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/over_integration_accuracy_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-2_dofs096_p5_CFL-0.10_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL",\
    "$c_{DG}$ NSFR.IR-GL-OI-3",\
    "sDG-OI-6", \
    "sDG-OI-2",\
    "sDG-OI-0",\
    ]
    list_of_poly_degree=[5,5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16,16]
    
    batch_plot_spectra(96,"OI_stability_GL",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        which_lines_dashed=[4],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512/",\
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI",\
    "Strong DG-Roe-GL",\
    "$c_{DG}$ NSFR.IR-GL",\
    "$c_{DG}$ NSFR.IR-GL-Roe",\
    ]
    list_of_poly_degree=[5,5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16,16]
    
    batch_plot_spectra(96,"sDG_with_and_without_OI_vs_NSFR",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/correction_parameter/viscous_TGV_ILES_NSFR_cSD_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/correction_parameter/viscous_TGV_ILES_NSFR_cHU_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$", \
    "$c_{SD}$", \
    "$c_{HU}$", \
    "$c_{+}$", \
    ]
    list_of_poly_degree=[5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16]
    
    batch_plot_spectra(96,"p5_correction_parameter",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.26_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.36_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.14_procs512/",\
    ] 
    batch_labels = [ \
    "$c_{DG}$: CFL=$0.10$",\
    "$c_{DG}$: CFL=$0.26$",\
    "$c_{+}$: CFL=$0.10$",\
    "$c_{+}$: CFL=$0.36$",\
    "sDG: CFL=$0.10$",\
    "sDG: CFL=$0.14$",\
    ]
    list_of_poly_degree=[5,5,5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16,16,16]
    
    batch_plot_spectra(96,"p5_correction_parameter_cfl_advantage",batch_paths,batch_labels,
        solid_and_dashed_lines=True,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-3_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-3_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "Strong DG-Roe-GL-OI", \
    "GL", \
    "GLL", \
    "GL-OI-3", \
    "GLL-OI-3", \
    ]
    list_of_poly_degree=[5,5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16,16]
    
    batch_plot_spectra(96,"p5_flux_nodes",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        which_lines_dashed=[5],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/two_point_flux/viscous_TGV_ILES_NSFR_cDG_KG_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/two_point_flux/viscous_TGV_ILES_NSFR_cDG_CH_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/two_point_flux/viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.KG-GL", \
    "$c_{DG}$ NSFR.CH-GL", \
    "$c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GL", \
    ]
    list_of_poly_degree=[5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16]
    
    batch_plot_spectra(96,"p5_two_point_flux",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)
# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "No Upwinding", \
    "LxF", \
    "Roe", \
    "L$^2$Roe", \
    ]
    list_of_poly_degree=[5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16]
    
    batch_plot_spectra(96,"p5_upwind_gl",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG.LRNC_CLIPMC-0.01-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16_corrected/",\
    ]
    batch_labels = [ \
    "No Model", \
    "SM", \
    "SI.SM", \
    "HPF.SM", \
    "HPF.SI.SM", \
    "DSM", \
    ]
    list_of_poly_degree=[5,5,5,5,5,5]
    list_of_number_of_elements_per_direction=[16,16,16,16,16,16]
    
    batch_plot_spectra(96,"p5_selected_sgs_models_gl",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        y_limits_zoom_input=[2.7e-5, 5.0e-4],
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG.LRNC_MC-0.10_strong_DG_Roe_GLL_OI-0_dofs096_p5_CFL-0.1_procs512/", \
    ]
    batch_labels = [ \
    "sDG-GLL-OI", \
    "sDG-GLL-HPF.SI.SM", \
    ]
    list_of_poly_degree=[5,5]
    list_of_number_of_elements_per_direction=[16,16]
    
    batch_plot_spectra(96,"sDG_gll_sgs_model_stabilization",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        y_limits_zoom_input=[4e-5, 5.0e-4],
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)