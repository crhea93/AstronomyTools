# AstronomyTools
A set of tools written for astronomy applications

If you wish to use any of my programs feel free, but please let me first if it is for published material (if you don't mind that is!). I am curious to know how these are being used. This will be updated regularly with new tools :)

I advise cloning the entire repository since several programs require functions in GeneralUse/ToolBox.py


Astrometry Folder:

	LScalc -- Calculate the linear distance given redshift and angular seperation

	AScalc -- Calculate the angular seperation given redshift and linear distance

GeneralUse Folder:

	Precursor.py -- Run specextract given folder and region of interest	
	
	Merge.py -- Merge several observations

	ToolBox.py -- Collection of subroutines used in other scripts (i.e. surface brightness calculations)

	ScaledImage.py -- Create scaled image from image.fits

SurfaceBrightness Folder:

	SurfBright.py -- Calculate Surface Brightness for 40 and 400 kpc using srcflux

	CSB.py -- Calculate Surface Brightness Coefficient from srcflux

	SBCalc_complete.py -- Calculates Surface Brightness in a rigorous manner using aprates

	SB_BoundsCalc.py -- To be used in conjunction with SBCalc_complete.py to generate error bounds

	read_input.py -- Read input files 

TemperatureMapPipeline Folder:

	Bin_data.py -- Bins data given binning file (See WVT repository)

	Xspec.py -- Calculates Temperature for each bin

	Cones_paraview.py -- Generates Temperature Map Graphic

WeightedVoronoiTessellations Folder:
	
	WVT.py -- Weighted Voronoi Tesselations Algorithm
