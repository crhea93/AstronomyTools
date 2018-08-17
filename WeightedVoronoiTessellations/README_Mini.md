# WVT
Weighted Voronoi Tessellations in Python for X-ray image analysis

Additional modules to download (using pip install module_name OR pip3 install module_name)

	tqdm

	sklearn


IF YOU BELIEVE ANY DOCUMENTATION IS MISSING PLEASE LET ME KNOW :)

Please read full documentation pdf...

The data used in testing came from a fits file provided by Chandra (OBSID 12833)

TO RUN:

	python WVT.py input.i

If you are like me, then you may have to run:
	
	python3 WVT.py input.i

input.i contains all necessary input information (see example file)

-----------------------------------------------------------------
Created fits by using the following prescription:
Be sure to have ciao installed and running (%unix ciao)
http://cxc.harvard.edu/ciao/download/

####

1 - find region of interest! Please save in IMAGE Coordinates when saving the .reg file
	Also create a _bkg.reg file for the background region (background of interest)

####

These steps are outlined in detail with generic commands in the documentation
WVT program WILL NOT RUN if you have ciao bash initiated in the terminal
-----------------------------------------------------------------

-----------------------------------------------------------------
ADDITIONAL NOTES:
Verification: Contains 2 programs

	2 - Verification.ipynb --> Jupyter Notebook to check that
	residual is zero
-----------------------------------------------------------------
