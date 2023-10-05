# plst-fitsregs
Plot Fits Files and overlay (Circular) DS9 Regions from files.
	Given a list of fits file images, this will iterate over them and: 
	1. Attempt to convert to MJy/sr unless units in Jy/beam
	2. Reproject the image so that pixels are aligned and North is up if, needed.
	3. Read in a list of DS9 region files, and from each region file get the first region whose title contains a text specified by the user.
	4. Make a nice WCS aligned plot.
	5. The title of the plot will contain a name specified by the user, and from the fits header the telescope name and the wavelength in microns.
	6. Display or save as pdf, the plot.
So far works with WISE, Herschel, SCUBA1 and 2, MIPS data.
