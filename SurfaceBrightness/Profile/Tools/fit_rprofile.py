'''
Script to fit a single King Beta Model using a pre-created radial SB profile
'''
from sherpa.plot import DataPlot
from sherpa.astro.all import *
from sherpa.astro.ui import *
from pychips.all import *
from sherpa.all import *
import os



def profiles(scaling,flux=False,model_type='single'):
    # Create basic profile
    if flux == False:
        load_data(1,'rprofile_r_data.fits', 3, ["R_NEW","SUR_BRI","SUR_BRI_ERR"])
    elif flux == True:
        load_data(1,'rprofile_r_data.fits', 3, ["R_NEW","SUR_FLUX","SUR_FLUX_ERR"])
    if model_type == 'single':
        set_source("beta1d.src1")
    elif model_type == 'double':
        set_source("beta1d.src1+beta1d.src2")
    data = get_data()
    data.x = data.x*scaling #now in kpc
    plot_data()
    set_plot_xlabel("R (kpc)")
    set_plot_title("Surface Brightness Profile")
    set_plot_ylabel("photons s^{-1} cm^{-2} arcsec^{-2}")
    log_scale(XY_AXIS)
    print_window("SurfaceBrightnessProfile.png",['clobber','true'])
    #stats, accept, params = get_draws(niter=1e4)
    # King Beta Plots
    plot_fit_delchi()
    set_current_plot("plot1")
    log_scale(X_AXIS)
    print_window("KingBetaFit.png",['clobber','true'])
    log_scale(XY_AXIS)
    print_window("KingBetaFitLOG.png",['clobber','true'])

    return None
