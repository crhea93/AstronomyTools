'''
Clean script to fit SpARCS1049 using McDonald Background Model
'''
from sherpa.astro.xspec import *
from sherpa.astro.all import *
from sherpa.astro.ui import *
from pychips.all import *
from sherpa.all import *
#---------------------------INPUTS------------------------------#
energy_min = 0.1
energy_max = 2.4
flux_emin = 0.01
flux_emax = 50.0
redshift_cluster = 0.3755
group_count = 5
stat = 'chi2gehrels'
method = 'levmar'
guess_temp = 2
region = 'R_total'
save_name = 'MACS_total'
#---------------------------------------------------------------#
from sherpa_contrib.xspec.xsconvolve import load_xscflux
load_xscflux("cflux")
def set_log_sherpa():
    p = get_data_plot_prefs()
    p["xlog"] = True
    p["ylog"] = True
    return None
#--------------------------READING IN---------------------------#
#Read in pha files
pha1 = load_pha(1, "10481/repro/"+region+".pi",use_errors=True)
pha2 = load_pha(2, "17233/repro/"+region+".pi",use_errors=True)
pha3 = load_pha(3, "18825/repro/"+region+".pi",use_errors=True)
#Ignore certain energies
set_stat(stat)
set_method(method)
ignore(0,energy_min)
ignore(energy_max,)
#-------------------------SOURCE MODEL--------------------------#
set_source(1, xsphabs.abs1 * (xsapec.apec1))
set_source(2 , abs1*(xsapec.apec2))
set_source(3 , abs1*(xsapec.apec3))
apec1.kT = guess_temp ;apec2.kT = apec1.kT; apec3.kT = apec1.kT
#Previously calculated abundance
apec1.Abundanc = 0.3; apec2.Abundanc = apec1.Abundanc;
apec3.Abundanc = apec1.Abundanc
apec1.redshift = redshift_cluster
apec2.redshift = apec1.redshift
apec3.redshift = apec1.redshift
abs1.nH = 2.27e-2
freeze(abs1.nH)
group_counts(1,group_count);group_counts(2,group_count)
group_counts(3,group_count)
notice(energy_min,energy_max)
#-------------------------BACKGROUND MODEL--------------------------#
set_bkg_model(1, xsapec.bkgapec1+abs1*xsbremss.brem1)
set_bkg_model(2, xsapec.bkgapec2+abs1*xsbremss.brem2)
set_bkg_model(3, xsapec.bkgapec3+abs1*xsbremss.brem3)
bkgapec1.kT = 0.18;freeze(bkgapec1.kT);
brem1.kT = 40;freeze(brem1.kT)
bkgapec2.kT = 0.18;freeze(bkgapec2.kT);
brem2.kT = 40;freeze(brem2.kT)
bkgapec3.kT = 0.18;freeze(bkgapec3.kT);
brem3.kT = 40;freeze(brem3.kT)
#-------------------------FIT MODEL----------------------------------#
#freeze(apec1.kT)
fit()
#-------------------------PLOT---------------------------------------#
set_log_sherpa()
plot("fit", 1, "fit", 2)#, "bkgfit", 1, "bkgfit", 2)
print_window("Spectra_1_2.png",['clobber', True])
plot("bkgfit", 1, "bkgfit", 2)
print_window("Spectra_1_2_bkg.png",['clobber', True])
#-------------------------ERROR CALCULATIONS-------------------------#
#covar(apec1.kT)
#-------------------------FLUX CALCULATIONS--------------------------#
print("Now calculating Flux...")
freeze(apec1.kT)
set_source(1, abs1*cflux(apec1))
set_source(2, abs1*cflux(apec2))
set_source(3, abs1*cflux(apec3))
cflux.lg10Flux.val = -11. #From previous runs
set_method('levmar')
cflux.Emin = flux_emin
cflux.Emax = flux_emax
fit()
Flux = cflux.lg10Flux.val
#-------------------------PLOT FLUX------------------------#
set_log_sherpa()
plot("fit", 1, "fit", 2)#, "bkgfit", 1, "bkgfit", 2)
print_window("Flux_Spectra_1_2.png",['clobber', True])
plot("bkgfit", 1, "bkgfit", 2)
print_window("Flux_Spectra_1_2_bkg.png",['clobber', True])
#-------------------------CLEANUP--------------------------#
#save(save_name+'.shp')
reset(get_model())
reset(get_source())
clean()
