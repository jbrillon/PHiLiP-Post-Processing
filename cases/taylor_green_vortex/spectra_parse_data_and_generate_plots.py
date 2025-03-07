#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
# load tools
sys.path.append(CURRENT_PATH+"../../src");
from plot_unsteady_integrated_turbulent_flow_quantities import get_dissipation_discrete
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
from scipy import integrate
import matplotlib;from matplotlib.lines import Line2D
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
    return get_truncated_spectra_from_cutoff_wavenumber_and_spectra(spectra,cutoff_wavenumber)
#-----------------------------------------------------
def get_unresolved_spectra_from_DOFs_information(spectra, poly_degree, number_of_elements_per_direction,truncate_spectra_at_effective_DOFs):
    cutoff_wavenumber = get_cutoff_wavenumber(poly_degree,number_of_elements_per_direction,truncate_spectra_at_effective_DOFs)
    idx = (np.abs(spectra[:,0] - cutoff_wavenumber)).argmin()
    # return spectra[(idx+1):,:]
    return spectra[(idx):,:]
#-----------------------------------------------------
def append_to_plot(x_,y_,label_):
    global x,y,labels
    labels.append(label_);x.append(x_);y.append(y_)
#-----------------------------------------------------
def batch_append_to_plot(paths_,labels_,filename,list_of_poly_degree_,list_of_number_of_elements_per_direction_,truncate_spectra_at_effective_DOFs,
    append_unresolved_wavenumber_range=False):
    global x,y,labels
    for i,path in enumerate(paths_):
        spectra_ = np.loadtxt(filesystem+path+filename)
        poly_degree = list_of_poly_degree_[i]
        number_of_elements_per_direction = list_of_number_of_elements_per_direction_[i]
        if(append_unresolved_wavenumber_range):
            spectra = get_unresolved_spectra_from_DOFs_information(spectra_, poly_degree, number_of_elements_per_direction,truncate_spectra_at_effective_DOFs)
            x.append(spectra[:,0]);y.append(spectra[:,1])
            labels.append("unresolved")
        else:
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
    truncate_spectra_at_effective_DOFs=True, # Modify this here
    fix_legend_location_for_presentation=False,
    markers_on=False,
    plot_unresolved_wavenumber_range_as_dashed=False,
    plot_full_wavenumber_range_of_reference_DNS=False,
    extend_y_max_limit=False,
    zoom_box_on_left=False):
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
        mrkr_input_store = ['None','None','None','None','None','None','None','None']
        lnstl_input_store = ['solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']
    if(markers_on):
        lnstl_input_store =['solid','dashed','dashdot','dotted','dotted','dotted']
        mrkr_input_store = ['+','d','s','o','s','^','d','v','>','<']

    cutoff_wavenumber_store = []
    grid_cutoff_wavenumber_store = []
    if(plot_cutoff_wavenumber_asymptote):
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
        append_to_plot(spectra[:,0],spectra[:,1],"DNS ($256^3$ DOF, p$7$)")
        i_curve += 1
    elif(plot_reference_result):
        # which_lines_black.append(0)
        clr_input_store.insert(i_curve,"k")
        spectra = np.loadtxt(CURRENT_PATH+"data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
        append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
        i_curve += 1
    # batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra_no_smoothing.dat")

    if(plot_PHiLiP_DNS_result_as_reference or plot_reference_result):
        if(solid_and_dashed_lines or dashed_and_solid_lines or markers_on):
            mrkr_input_store.insert(i_curve-1,'None')
            lnstl_input_store.insert(i_curve-1,'solid')

    if(plot_filtered_dns):
        # filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256_projected_to_p2_dofs096-0_reordered_spectra_oversampled_nquad12.dat"
        # spectra = np.loadtxt(filepath_to_reference_result)
        # append_to_plot(spectra[:,0],spectra[:,1],"Projected DNS ($96^3$ DOFs, p$2$)")
        clr_input_store.insert(i_curve,"k")
        # which_lines_black.append(i_curve)
        if(solid_and_dashed_lines or dashed_and_solid_lines or markers_on):
            mrkr_input_store.insert(i_curve,'None')
        lnstl_input_store.insert(i_curve,'dashed') # for dashed filtered DNS result
        # which_lines_dashed.append() 
        i_curve += 1

    # compute reference curve 1
    index_of_reference_curve = len(batch_labels)+i_curve
    x_ref_curve = np.linspace(2.0e0,3.0e2,100)
    order_for_ref_curve = -5.0/3.0
    ref_curve_label = "$\\left(\\kappa^{*}\\right)^{-5/3}$"
    shift = 2.0
    y_ref_curve = (x_ref_curve**(order_for_ref_curve))/np.exp(shift)
    clr_input_store.insert(index_of_reference_curve,"k")#"tab:gray"
    lnstl_input_store.insert(index_of_reference_curve,"dotted")
    mrkr_input_store.insert(index_of_reference_curve,'None')

    # add reference curve
    # append_to_plot(x_ref_curve,y_ref_curve,ref_curve_label)
    if(nDOF_==256):
        # x_limits_zoom=[35, 70]
        # y_limits_zoom=[1.0e-5, 1.2e-4]
        x_limits_zoom=[80, 112]
        y_limits_zoom=[4.5e-7, 7.0e-6]
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

    # t=9

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
        append_to_plot(spectra[:,0],spectra[:,1],"DNS ($256^{3}$ DOF, p$7$)")
        i_curve += 1
    elif(plot_reference_result):
        # Note: This is the 512^3 spectral reference result in their paper
        spectra = np.loadtxt(CURRENT_PATH+"data/carton2014_spectral_512_spectra_t9.txt",skiprows=1,delimiter=',')
        append_to_plot(spectra[:,0],spectra[:,1],"$512^{3}$ [Carton de Wiart et al.]")
        i_curve += 1
    if(plot_filtered_dns):
        filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256_projected_to_p2_dofs096-1_reordered_spectra_oversampled_nquad12.dat"
        spectra_ = np.loadtxt(filepath_to_reference_result)
        # spectra = get_truncated_spectra_from_DOFs_information(spectra_, 2, 32)
        spectra = get_truncated_spectra_from_DOFs_information(spectra_, 5, 16, truncate_spectra_at_effective_DOFs) # to match cut-off for 96P5
        append_to_plot(spectra[:,0],spectra[:,1],"Projected DNS\n ($96^{3}$ DOF, p$2$)")
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

    legend_anchor_input=[]
    legend_location_input="lower left"
    if(fix_legend_location_for_presentation):
        legend_anchor_input=[0.0,0.45]
        legend_location_input="upper left"

    # for plotting unresolved spectra component
    leg_elements_input=[]
    clr_input_store_=[]
    mrkr_input_store_=[]
    lnstl_input_store_=[]
    number_of_legend_entries = len(labels)
    for i in range(0,number_of_legend_entries):
        leg_elements_input.append(Line2D([0],[0], label=labels[i], color=clr_input_store[i],
                                 marker=mrkr_input_store[i], markersize=6, mfc='None', linestyle=lnstl_input_store[i]))
        clr_input_store_.append(clr_input_store[i])
        mrkr_input_store_.append(mrkr_input_store[i])
        lnstl_input_store_.append(lnstl_input_store[i])

    if(plot_unresolved_wavenumber_range_as_dashed):
        index_where_result_curves_start = labels.index(batch_labels[0])
        number_of_result_curves = len(batch_labels)
        # append the unresolved wavenumber range to the plotted x,y arrays
        batch_append_to_plot(batch_paths, batch_labels, 
            "flow_field_files/velocity_vorticity-1_reordered_spectra_no_smoothing.dat",
            list_of_poly_degree_input, list_of_number_of_elements_per_direction_input,truncate_spectra_at_effective_DOFs,append_unresolved_wavenumber_range=True)
        # set the line styles
        for i in range(0,number_of_result_curves):
            index_for_line_stuff = i+index_where_result_curves_start
            clr_input_store_.append(clr_input_store[index_for_line_stuff])
            mrkr_input_store_.append(mrkr_input_store[index_for_line_stuff])
            lnstl_input_store_.append("dotted")

    # regular range in paper
    xlimits_=[4.0e0,112]; ylimits_=[4.0e-7,1e-2]
    if(plot_full_wavenumber_range_of_reference_DNS):
        # xlimits_=[2e0,7e1]; ylimits_=[2.45e-5,3e-2]
        xlimits_=[4.0e0,240]; ylimits_=[1.0e-9,1e-2]
    if(extend_y_max_limit):
        ylimits_[1]=2e-2

    transparent_legend_=False
    zoom_box_origin_and_extent_=[0.65, 0.65, 0.32, 0.32]
    if(zoom_box_on_left):
        zoom_box_origin_and_extent_=[0.03, 0.47, 0.32, 0.32]
        leg_elements_input[0], leg_elements_input[-1] = leg_elements_input[-1], leg_elements_input[0]
        transparent_legend_=True

    qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $\\kappa^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(\\kappa^{*},t^{*})$",
        title_label=title_label,
        leg_elements_input=leg_elements_input,
        fig_directory=figure_directory,figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
        nlegendcols=1,
        xlimits=xlimits_,ylimits=ylimits_,
        markers=False,
        legend_on=True,
        # legend_labels_tex=labels,
        which_lines_black=which_lines_black,
        transparent_legend=transparent_legend_,
        legend_border_on=False,
        grid_lines_on=False,
        clr_input=clr_input_store_,
        mrkr_input=mrkr_input_store_,
        lnstl_input=lnstl_input_store_,
        legend_fontSize=legend_fontSize_input,
        legend_location=legend_location_input,
        legend_anchor=legend_anchor_input,
        # which_lines_only_markers=[1,2,3],
        # which_lines_dashed=which_lines_dashed,
        plot_zoomed_section=plot_zoomed_section,
        x_limits_zoom=x_limits_zoom,y_limits_zoom=y_limits_zoom,
        zoom_box_origin_and_extent=zoom_box_origin_and_extent_,
        # vertical_lines=cutoff_wavenumber_store,
        # secondary_vertical_lines=grid_cutoff_wavenumber_store
        )
    
    return
#-----------------------------------------------------
def get_total_turbulent_kinetic_energy_from_spectra(spectra):
    wavenumbers = spectra[:,0]
    turbulent_kinetic_energy = spectra[:,1] 
    total_turbulent_kinetic_energy = integrate.trapezoid(turbulent_kinetic_energy,x=wavenumbers)
    return total_turbulent_kinetic_energy
#-----------------------------------------------------
def get_total_dissipation_rate_from_spectra(spectra,reynolds_number):
    wavenumbers = spectra[:,0]
    turbulent_kinetic_energy = spectra[:,1]
    dissipation_rate_components = wavenumbers*wavenumbers*turbulent_kinetic_energy
    total_dissipation_rate = 2.0*np.sum(dissipation_rate_components)/reynolds_number # non-dimensional form
    return total_dissipation_rate
#-----------------------------------------------------
def batch_compute_resolved_turbulent_kinetic_energy(
    paths_,
    labels_,
    list_of_poly_degree_,
    list_of_number_of_elements_per_direction_):
    
    #----------------------------------------------------------------
    # Reference result
    #----------------------------------------------------------------
    # "DNS ($256^{3}$p$7$)"
    filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256-1_reordered_spectra_oversampled_nquad16.dat"
    reference_spectra = np.loadtxt(filepath_to_reference_result)
    
    # Parameters:
    truncate_spectra_at_effective_DOFs=True
    filename="flow_field_files/velocity_vorticity-1_reordered_spectra_no_smoothing.dat"
    unsteady_filename="turbulent_quantities.txt"
    spectra_time = 9.0
    
    print("Label \t\t\t| P, Nel, Total DOFs, Eff.-DOFs\t | Total KE from Spectra\t| Inst. Total KE\t| Rel. Percentage")
    print("---------------------------------------------------------------------------------")
    # Loop through all the paths
    for i,path in enumerate(paths_):
        spectra_ = np.loadtxt(filesystem+path+filename)
        poly_degree = list_of_poly_degree_[i]
        number_of_elements_per_direction = list_of_number_of_elements_per_direction_[i]
        cutoff_wavenumber = get_cutoff_wavenumber(poly_degree,number_of_elements_per_direction,truncate_spectra_at_effective_DOFs) # True for effective DOFs
        spectra = get_truncated_spectra_from_cutoff_wavenumber_and_spectra(spectra_, cutoff_wavenumber)
        total_tke = get_total_turbulent_kinetic_energy_from_spectra(spectra)
        # load the KE from the unsteady quantity
        if("spectra_fix/" in path):
            path = path.replace("spectra_fix/","")
        elif(path=="NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/"):
            print("here")
            path="NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/"

        time, kinetic_energy = np.loadtxt(filesystem+path+unsteady_filename,usecols=(0,1),skiprows=1,dtype=np.float64,unpack=True) # 8.5737827439557435e-02
        # get instantaneous physical KE
        idx_of_KE_at_spectra_time = (np.abs(time - spectra_time)).argmin()
        inst_total_KE = kinetic_energy[idx_of_KE_at_spectra_time]
        # spectra = 1.0*spectra_ # uncomment for no truncation
        truncated_reference_spectra = get_truncated_spectra_from_cutoff_wavenumber_and_spectra(reference_spectra, cutoff_wavenumber)
        # reference_total_tke = get_total_turbulent_kinetic_energy_from_spectra(truncated_reference_spectra)
        reference_total_tke = 1.0*inst_total_KE
        # percentage_of_tke_captured = 100.0*total_tke/reference_total_tke
        percentage_of_tke_captured = 100.0*np.abs(total_tke-reference_total_tke)/reference_total_tke
        total_nDOFs_per_dim = (poly_degree+1)*number_of_elements_per_direction
        effective_nDOFs_per_dim = (poly_degree)*number_of_elements_per_direction
        print("%s \t| %i, %i, %i^3, %i^3 \t\t| %.6e\t\t\t| %.6e\t\t| %3.6f %%" % (labels_[i], poly_degree, number_of_elements_per_direction, total_nDOFs_per_dim, effective_nDOFs_per_dim, total_tke, reference_total_tke, percentage_of_tke_captured))
        # labels_[i]
    print("---------------------------------------------------------------------------------")
    return
#-----------------------------------------------------
def batch_compute_resolved_dissipation_rate(
    paths_,
    labels_,
    list_of_poly_degree_,
    list_of_number_of_elements_per_direction_,
    reynolds_number_=1600.0):
    
    #----------------------------------------------------------------
    # Reference result
    #----------------------------------------------------------------
    # "DNS ($256^{3}$p$7$)"
    filepath_to_reference_result=CURRENT_PATH+"data/brillon/flow_field_files/velocity_vorticity_p7_dofs256-1_reordered_spectra_oversampled_nquad16.dat"
    reference_spectra = np.loadtxt(filepath_to_reference_result)
    
    # Parameters:
    truncate_spectra_at_effective_DOFs=True
    filename="flow_field_files/velocity_vorticity-1_reordered_spectra_no_smoothing.dat"
    unsteady_filename="turbulent_quantities.txt"
    spectra_time = 9.0
    
    print("Label \t\t\t| P, Nel, Total DOFs, Eff.-DOFs\t | Total DR from Spectra\t| Inst. Total DR\t| Rel. Percentage")
    print("---------------------------------------------------------------------------------")
    # Loop through all the paths
    for i,path in enumerate(paths_):
        spectra_ = np.loadtxt(filesystem+path+filename)
        poly_degree = list_of_poly_degree_[i]
        number_of_elements_per_direction = list_of_number_of_elements_per_direction_[i]
        cutoff_wavenumber = get_cutoff_wavenumber(poly_degree,number_of_elements_per_direction,truncate_spectra_at_effective_DOFs) # True for effective DOFs
        spectra = get_truncated_spectra_from_cutoff_wavenumber_and_spectra(spectra_, cutoff_wavenumber)
        total_DR = get_total_dissipation_rate_from_spectra(spectra,reynolds_number_)
        # load the KE from the unsteady quantity
        if("spectra_fix/" in path):
            path = path.replace("spectra_fix/","")
        elif(path=="NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/"):
            print("here")
            path="NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/"

        time, kinetic_energy = np.loadtxt(filesystem+path+unsteady_filename,usecols=(0,1),skiprows=1,dtype=np.float64,unpack=True) # 8.5737827439557435e-02
        dissipation_rate = get_dissipation_discrete(time, kinetic_energy)
        # get instantaneous physical dissipation_rate (DR)
        idx_of_DR_at_spectra_time = (np.abs(time - spectra_time)).argmin()
        inst_total_DR = dissipation_rate[idx_of_DR_at_spectra_time]
        # spectra = 1.0*spectra_ # uncomment for no truncation
        truncated_reference_spectra = get_truncated_spectra_from_cutoff_wavenumber_and_spectra(reference_spectra, cutoff_wavenumber)
        # reference_total_DR = get_total_dissipation_rate_from_spectra(truncated_reference_spectra,reynolds_number_)
        reference_total_DR = 1.0*inst_total_DR
        # percentage_of_tke_captured = 100.0*total_tke/reference_total_tke
        percentage_of_DR_captured = 100.0*np.abs(total_DR-reference_total_DR)/reference_total_DR
        total_nDOFs_per_dim = (poly_degree+1)*number_of_elements_per_direction
        effective_nDOFs_per_dim = (poly_degree)*number_of_elements_per_direction
        print("%s \t| %i, %i, %i^3, %i^3 \t\t| %.6e\t\t\t| %.6e\t\t| %3.6f %%" % (labels_[i], poly_degree, number_of_elements_per_direction, total_nDOFs_per_dim, effective_nDOFs_per_dim, total_DR, reference_total_DR, percentage_of_DR_captured))
        # labels_[i]
    print("---------------------------------------------------------------------------------")
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
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/",\
    "NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/",\
    ]
    batch_labels = [ \
    "$96^{3}$ DOF, p$5$", \
    "$64^{3}$ DOF, p$7$", \
    "$48^{3}$ DOF, p$5$",\
    "$24^{3}$ DOF, p$5$",\
    # "$96^{3}$", \
    # "$64^{3}$", \
    # "$48^{3}$",\
    # "$24^{3}$",\
    "$32^{3}$p$7$",\
    "$64^{3}$p$3$",\
    ]
    list_of_poly_degree=[5,7,5,5,7,3]
    list_of_number_of_elements_per_direction=[16,8,8,4,32,64]
    
    batch_plot_spectra("all","cDG_NSFR_convergence",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction,
        plot_unresolved_wavenumber_range_as_dashed=True,
        extend_y_max_limit=True)
    batch_compute_resolved_turbulent_kinetic_energy(batch_paths,batch_labels,list_of_poly_degree,list_of_number_of_elements_per_direction)
    batch_compute_resolved_dissipation_rate(batch_paths,batch_labels,list_of_poly_degree,list_of_number_of_elements_per_direction)
# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/",\
    "NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/",\
    ]
    batch_labels = [ \
    "$96^{3}$ DOF, p$5$", \
    "$64^{3}$ DOF, p$7$", \
    "$48^{3}$ DOF, p$5$",\
    "$24^{3}$ DOF, p$5$",\
    # "$96^{3}$", \
    # "$64^{3}$", \
    # "$48^{3}$",\
    # "$24^{3}$",\
    "$32^{3}$p$7$",\
    "$64^{3}$p$3$",\
    ]
    list_of_poly_degree=[5,7,5,5,7,3]
    list_of_number_of_elements_per_direction=[16,8,8,4,32,64]
    
    batch_plot_spectra("all","cDG_NSFR_convergence_without_oversampling",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction,
        plot_unresolved_wavenumber_range_as_dashed=True,
        extend_y_max_limit=True)
    batch_compute_resolved_turbulent_kinetic_energy(batch_paths,batch_labels,list_of_poly_degree,list_of_number_of_elements_per_direction)
    batch_compute_resolved_dissipation_rate(batch_paths,batch_labels,list_of_poly_degree,list_of_number_of_elements_per_direction)
exit()
# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    # "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-3_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-3_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    # "Strong DG-Roe-GL-OI", \
    "GL", \
    "GLL", \
    "GL-OI.3", \
    "GLL-OI.3", \
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
        which_lines_dashed=[4],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

#=====================================================
# DOFs: ALL | NSFR CONVERGENCE VERSION FOR WCCM
#-----------------------------------------------------
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/",\
    ]
    batch_labels = [ \
    "$256^{3}$ ($32^{3}$p$7$)",\
    "$96^{3}$ ($16^{3}$p$5$)", \
    "$64^{3}$ ($8^{3}$p$7$)", \
    "$48^{3}$ ($8^{3}$p$5$)",\
    "$24^{3}$ ($4^{3}$p$7$)",\
    ]
    list_of_poly_degree=[7,5,7,5,5]
    list_of_number_of_elements_per_direction=[32,16,8,8,4]

    for i in range(0,len(batch_paths)):
        figure_filename_postfix_input="cDG_NSFR_overview_for_wccm_%i"%i
        batch_plot_spectra("all",figure_filename_postfix_input,batch_paths[:(i+1)],batch_labels[:(i+1)],
            solid_and_dashed_lines=False,
            title_off=title_off_input,figure_directory=fig_dir_input,
            plot_cutoff_wavenumber_asymptote=False,
            truncate_spectra_at_effective_DOFs=True,
            plot_PHiLiP_DNS_result_as_reference=False,
            plot_filtered_dns=False,
            which_lines_dashed=[],
            plot_zoomed_section=False,
            list_of_poly_degree_input=list_of_poly_degree,
            list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction,
            fix_legend_location_for_presentation=True)

# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/",\
    "NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs1024/",\
    ]
    batch_labels = [ \
    "p$7$ $c_{DG}$ NSFR.IR-GL",\
    "p$3$ $c_{DG}$ NSFR.IR-GL",\
    "p$3$ $c_{DG}$ NSFR.IR-GL-Roe",\
    "p$3$ Strong DG-Roe-GL-OI",\
    "p$7$ $c_{DG}$ NSFR.IR-GLL",\
    ]
    lnstl_input=['solid','solid','solid','solid','dashed','dashed']
    list_of_poly_degree=[7,3,3,3,7]
    list_of_number_of_elements_per_direction=[32,64,64,64,32]
    
    batch_plot_spectra(256,"verification_before_oversampling",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        truncate_spectra_at_effective_DOFs=True,
        plot_PHiLiP_DNS_result_as_reference=False,
        plot_filtered_dns=False,
        which_lines_dashed=[3,4],
        plot_zoomed_section=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction,
        plot_unresolved_wavenumber_range_as_dashed=True,
        plot_full_wavenumber_range_of_reference_DNS=True,
        zoom_box_on_left=True)

# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs1024/",\
    ]
    batch_labels = [ \
    "p$7$ $c_{DG}$ NSFR.IR-GL",\
    "p$3$ $c_{DG}$ NSFR.IR-GL",\
    "p$3$ $c_{DG}$ NSFR.IR-GL-Roe",\
    "p$3$ Strong DG-Roe-GL-OI",\
    "p$7$ $c_{DG}$ NSFR.IR-GLL",\
    ]
    lnstl_input=['solid','solid','solid','solid','dashed','dashed']
    list_of_poly_degree=[7,3,3,3,7]
    list_of_number_of_elements_per_direction=[32,64,64,64,32]
    
    batch_plot_spectra(256,"verification",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        truncate_spectra_at_effective_DOFs=True,
        plot_PHiLiP_DNS_result_as_reference=False,
        plot_filtered_dns=False,
        which_lines_dashed=[3,4],
        plot_zoomed_section=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs1024/",\
    ]
    batch_labels = [ \
    "p$7$ $c_{DG}$ NSFR.IR-GL",\
    "p$3$ $c_{DG}$ NSFR.IR-GL",\
    "p$3$ $c_{DG}$ NSFR.IR-GL-Roe",\
    "p$3$ Strong DG-Roe-GL-OI",\
    "p$7$ $c_{DG}$ NSFR.IR-GLL",\
    ]
    lnstl_input=['solid','solid','solid','solid','dashed','dashed']
    list_of_poly_degree=[7,3,3,3,7]
    list_of_number_of_elements_per_direction=[32,64,64,64,32]
    
    batch_plot_spectra(256,"verification_full_range",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        truncate_spectra_at_effective_DOFs=True,
        plot_PHiLiP_DNS_result_as_reference=False,
        plot_filtered_dns=False,
        which_lines_dashed=[3,4],
        plot_zoomed_section=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction,
        plot_unresolved_wavenumber_range_as_dashed=True,
        plot_full_wavenumber_range_of_reference_DNS=True,
        zoom_box_on_left=True)
# =====================================================
# MISSING FIG 3
# =====================================================
'''
batch_paths = [ \
"NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
"NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
"NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
"NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
"NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/",\
"NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
"NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64/",\
"NarvalFiles/2023_JCP/spectra_fix/over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-3_dofs096_p5_procs512/",\
"NarvalFiles/2023_JCP/spectra_fix/over_integration_accuracy_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-2_dofs096_p5_CFL-0.10_procs512/",\
"NarvalFiles/2023_JCP/spectra_fix/filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512/",\
"NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512/",\
"NarvalFiles/2023_JCP/spectra_fix/correction_parameter/viscous_TGV_ILES_NSFR_cSD_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
"NarvalFiles/2023_JCP/spectra_fix/correction_parameter/viscous_TGV_ILES_NSFR_cHU_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
"NarvalFiles/2023_JCP/spectra_fix/correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
"NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GL_OI-0_dofs096_p5_procs512/", \
"NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GL_OI-0_dofs096_p5_procs512/", \
]
batch_labels = [ \
"$c_{DG}$ NSFR.IR-GL",\
"$c_{DG}$ NSFR.IR-GLL",\
"sDG-GL-Roe-OI", \
"$c_{DG}$ NSFR.IR-GL", \
"sDG-GL-Roe-OI", \
"$c_{DG}$ NSFR.IR-GL",\
"sDG-GL-Roe-OI",\
"$c_{DG}$ NSFR.IR-GL-OI-3",\
"sDG-GL-Roe-OI-2",\
"sDG-GL-Roe-OI-0",\
"$c_{DG}$ NSFR.IR-GL-Roe",\
"$c_{SD}$ NSFR.IR-GL",\
"$c_{HU}$ NSFR.IR-GL",\
"$c_{+}$ NSFR.IR-GL",\
"$c_{DG}$ NSFR.IR-GL-L$^2$Roe",\
"$c_{DG}$ NSFR.IR-GL-LxF",\
]
list_of_poly_degree=[5,5,5,7,7,5,5,5,5,5,5,5,5,5,5,5]
list_of_number_of_elements_per_direction=[16,16,16,8,8,8,8,16,16,16,16,16,16,16,16,16]
batch_compute_resolved_turbulent_kinetic_energy(batch_paths,batch_labels,list_of_poly_degree,list_of_number_of_elements_per_direction)
'''
# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "p$5$ $c_{DG}$ NSFR.IR-GL",\
    "p$5$ $c_{DG}$ NSFR.IR-GLL",\
    ]
    lnstl_input=['solid','solid']
    list_of_poly_degree=[5,5]
    list_of_number_of_elements_per_direction=[16,16]
    
    batch_plot_spectra(96,"p5",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        truncate_spectra_at_effective_DOFs=True,
        plot_PHiLiP_DNS_result_as_reference=False,
        plot_filtered_dns=False,
        which_lines_dashed=[3,4],
        plot_zoomed_section=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/",\
    "NarvalFiles/2023_JCP/spectra_fix/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs1024/",\
    ]
    batch_labels = [ \
    "p$7$ $c_{DG}$ NSFR.IR-GL",\
    "p$7$ $c_{DG}$ NSFR.IR-GLL",\
    ]
    lnstl_input=['solid','solid']
    list_of_poly_degree=[7,7]
    list_of_number_of_elements_per_direction=[32,32]
    
    batch_plot_spectra(256,"p7_verification",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        truncate_spectra_at_effective_DOFs=True,
        plot_PHiLiP_DNS_result_as_reference=False,
        plot_filtered_dns=False,
        which_lines_dashed=[3,4],
        plot_zoomed_section=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/",\
    "NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_2refinements/",\
    # "NarvalFiles/2023_JCP/verification_tke_fix/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024_3refinements/",\
    ]
    batch_labels = [ \
    "$256^{3}$p$7$, $c_{DG}$", \
    "$256^{3}$p$7$, $c_{DG}$, oversampled", \
    ]
    list_of_poly_degree=[7,7]
    list_of_number_of_elements_per_direction=[32,32]
    
    batch_plot_spectra(256,"csme_oversampling_256",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=False,
        plot_filtered_dns=False,
        plot_zoomed_section=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    ]
    batch_labels = [ \
    "$96^{3}$ DOF, $c_{DG}$, without oversampling", \
    "$96^{3}$ DOF, $c_{DG}$, oversampled", \
    ]
    list_of_poly_degree=[5,5,7,7,5,5]
    list_of_number_of_elements_per_direction=[16,16,8,8,8,8]
    
    batch_plot_spectra(96,"csme_oversampling_96",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(False or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/",\
    ]
    batch_labels = [ \
    # "$96^{3}$ DOFs, p$5$", \
    # "$64^{3}$ DOFs, p$7$", \
    # "$48^{3}$ DOFs, p$5$",\
    # "$24^{3}$ DOFs, p$5$",\
    "$96^{3}$", \
    "$64^{3}$", \
    "$48^{3}$",\
    "$24^{3}$",\
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
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/",\
    "NarvalFiles/2023_JCP/spectra_fix/robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64/",\
    ]
    batch_labels = [ \
    "$96^{3}$ DOF, p$5$, $c_{DG}$", \
    "$96^{3}$ DOF, p$5$, sDG-OI.6", \
    "$64^{3}$ DOF, p$7$, $c_{DG}$", \
    "$64^{3}$ DOF, p$7$, sDG-OI.8", \
    "$48^{3}$ DOF, p$5$, $c_{DG}$",\
    "$48^{3}$ DOF, p$5$, sDG-OI.6",\
    ]
    list_of_poly_degree=[5,5,7,7,5,5]
    list_of_number_of_elements_per_direction=[16,16,8,8,8,8]
    
    batch_plot_spectra("all","convergence_comparison",batch_paths,batch_labels,
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
if(True or regenerate_all_plots):
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
if(True or regenerate_all_plots):
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
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=True,
        which_lines_dashed=[2,4,6],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction)

# =====================================================
if(True or regenerate_all_plots):
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
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction,
        markers_on=True)
# =====================================================
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512/", \
    "NarvalFiles/2023_JCP/spectra_fix/upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GL_OI-0_dofs096_p5_procs512/", \
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL\n(No Upwinding)", \
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
if(True or regenerate_all_plots):
    batch_paths = [ \
    "NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16/",\
    "NarvalFiles/2023_JCP/spectra_fix/sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG.LRNC_CLIPMC-0.01-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16_corrected/",\
    ]
    batch_labels = [ \
    "$c_{DG}$ NSFR.IR-GL\n(No SGS Model)", \
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
    "sDG-OI.6", \
    "sDG-OI.0-HPF.SI.SM", \
    ]
    list_of_poly_degree=[5,5]
    list_of_number_of_elements_per_direction=[16,16]
    
    batch_plot_spectra(96,"sDG_gll_sgs_model_stabilization",batch_paths,batch_labels,
        solid_and_dashed_lines=False,
        title_off=title_off_input,figure_directory=fig_dir_input,
        plot_cutoff_wavenumber_asymptote=True,
        plot_PHiLiP_DNS_result_as_reference=True,
        plot_filtered_dns=True,
        plot_zoomed_section=False,
        y_limits_zoom_input=[4e-5, 5.0e-4],
        which_lines_dashed=[],
        list_of_poly_degree_input=list_of_poly_degree,
        list_of_number_of_elements_per_direction_input=list_of_number_of_elements_per_direction,
        extend_y_max_limit=True)