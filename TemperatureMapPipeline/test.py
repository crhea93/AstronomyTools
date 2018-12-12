from astropy.io import fits
hdu_list = fits.open("../Merged/carre.fits", memmap=True)
counts = hdu_list[1].data
print(counts)
