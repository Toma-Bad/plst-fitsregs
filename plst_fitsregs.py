#!/usr/bin/env python3

import reproject
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.io import fits,ascii
import astropy.units as u
from astropy.io.fits import ImageHDU
#from regions import CRTFParser,read_ds9
from astropy.coordinates import SkyCoord
from regions.shapes.ellipse import EllipseSkyRegion
from astropy.table import Table,hstack
import re
import sys
import warnings
import pandas as pd
import pyregion
import copy
import numpy as np
import pkg_resources
import photutils as pht
from photutils.aperture import SkyCircularAperture, ApertureStats
from astropy.wcs import WCS
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse
import montage_wrapper as montage

wise_dn_to_jy = {1:1.9350E-06,
		2:2.7048E-06,
		3:1.8326e-06,
		4:5.2269E-05}

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)



def get_files_regs_data(fitsfilename,regionfilelist,unc_fitsfilename = None,beam=None,beam_params = None,rname='on'):
	print("   ")
	print("################")
	print(fitsfilename)
	with fits.open(fitsfilename) as fitsfile:
		try:
			hdulist = fitsfile
			hdu = copy.copy(fitsfile[0])
			data = np.squeeze(fitsfile[0].data)
			gooddata = data==data
			header = fitsfile[0].header
			wcs = WCS(header)
		except Exception as e:
			print(e, "reading fits")
		if len(header["INSTR*"]) != 0:
			instr = header["INSTR*"][0]
			print("Instrument:",instr)
		else:
			instr = None

	
		if "TELESCOP" in header:
			telname = header.get("TELESCOP")
			if "Herschel" in telname:
			
				hdu = copy.copy(fitsfile['image'])
				data = fitsfile['image'].data
				gooddata = data==data
				header = fitsfile['PRIMARY'].header + fitsfile['image'].header
				wcs = WCS(fitsfile['image'].header)
				if 'error' in fitsfile:
					unc_hdu = fitsfile['error']
					unc_header = fitsfile['error'].header
					unc_data = np.squeeze(fitsfile['error'].data)
					unc_gooddata = unc_data == unc_data
				if 'stDev' in fitsfile:
					unc_hdu = fitsfile['stDev']
					unc_header = fitsfile['stDev'].header
					unc_data = np.squeeze(fitsfile['stDev'].data)
					unc_gooddata = unc_data == unc_data
			if 'JCMT' in telname:
				hdu = copy.copy(fitsfile['PRIMARY'])
				data = np.squeeze(fitsfile['PRIMARY'].data)
				gooddata = data==data
				header = fitsfile['PRIMARY'].header 
				wcs = WCS(header)
				if 'VARIANCE' in fitsfile:
					header = fitsfile['PRIMARY'].header + fitsfile['VARIANCE'].header
					unc_hdu = fitsfile['VARIANCE']
					unc_header = fitsfile['VARIANCE'].header
					unc_data = np.sqrt(np.squeeze(fitsfile['VARIANCE'].data))
					print(unc_data.shape,"<<<<<<<<<<<<<")
					unc_gooddata = unc_data == unc_data
					
	print(telname)
	if unc_fitsfilename is not None:
		print("unc_file present!")
		with fits.open(unc_fitsfilename) as fitsfile:
			unc_hdu = copy.copy(fitsfile[0])
			unc_data = np.squeeze(fitsfile[0].data)
			unc_gooddata = unc_data==unc_data
			unc_header = fitsfile[0].header
	

	if len(header["WAVE*"]) != 0:
		wlen_val = header["WAVE*"][0]
		s = str(header["WAVE*"].comments[0])
		print("wavelength header data:",wlen_val,s)
		try:
			wlen_unit = u.Unit(s)
		except:
			print("Wavelength unit not present directly, searching")
		if '[' in s:
			unit_text= s[s.find('[')+1:s.find(']')]#re.search(r"\[([A-Za-z]+)\]",header["WAVE*"].comments[0]).group(1)
		elif '(' in s:
			unit_text= s[s.find('(')+1:s.find(')')] #re.search(r"\(([A-Za-z]+)\)",header["WAVE*"].comments[0]).group(1)
		if 'micron' in unit_text:
			unit_text = "micron"

		wlen_unit = u.Unit(unit_text)
		wlen = (wlen_val * wlen_unit).to(u.micron)
		print("retrieved WAVE*:",wlen)

	elif len(header["*WAVE*"]) != 0:
		wlen_val = header["*WAVE*"][0]
		s = str(header["*WAVE*"].comments[0])
		print("wavelength header data:",wlen_val,s)
		try:
			wlen_unit = u.Unit(s)
		except:
			print("Wavelength unit not present directly, searching")
		if '[' in s:
			unit_text= s[s.find('[')+1:s.find(']')]#re.search(r"\[([A-Za-z]+)\]",header["WAVE*"].comments[0]).group(1)
		elif '(' in s:
			unit_text= s[s.find('(')+1:s.find(')')] #re.search(r"\(([A-Za-z]+)\)",header["WAVE*"].comments[0]).group(1)
		if 'micron' in unit_text:
			unit_text = "micron"
		wlen_unit = u.Unit(unit_text)
		wlen = wlen_val * wlen_unit
		print("retrieved WAVE*:",wlen,wlen_unit)

	else:
		print("NO WLEN FOUND!")


	if "BAND" in header:
		band = header.get("BAND")
		print("BAND:",band)
	elif "DETECTOR" in header:
		band = header.get("DETECTOR")
	else:
		band = None


	if len(header["INSTR*"]) != 0:
		instr = header["INSTR*"][0]
		print("Instrument:",instr)
	else:
		instr = None
	if 'CDELT1' in header:
		CD1 = abs(header.get('CDELT1'))
		CD2 = abs(header.get('CDELT2'))
	elif 'CD1_1' in header:
		CD1 = abs(header.get('CD1_1'))
		CD2 = abs(header.get('CD2_2'))
	
	if 'CUNIT1' in header:
		print("CUNIT in header")
		unit_text = header.get("CUNIT1")
		CU = u.Unit(unit_text)
		if 'pix' in unit_text:
			CU = CU * u.pixel

	else:
		CU = u.deg
	print("Space unit:",CD1,CU,(CD1*CU).to(u.arcsec))
	pixarea = abs(CD1*CD2) * CU**2
	print("Pixarea:",pixarea)
	#check for beam:
	
	if "BUNIT" in header:
		try:
			str_unit = header.get("BUNIT")	
			real_unit = u.Unit(str_unit)
			if "pix" in str_unit:
				conv_factor = 1. * u.pixel
			elif "beam" in str_unit:
				if 'BMAJ' in header and 'BMIN' in header:
					bmaj = float(header['BMAJ']) * u.deg
					bmin = float(header['BMIN']) * u.deg
					conv_factor = 1./((2*np.pi*bmin*bmaj / pixarea / (8*np.log(2))).decompose() * 1. / u.beam)
					print("Beam present, bmin, bmaj, conv:",bmin,bmaj)
				elif 'JCMT' in telname and 'SCUBA-2' in instr:
					print("Beam inferred",telname,instr,wlen)
					if  430 *u.micron < wlen < 460*u.micron:
						conv_factor = 1./( np.pi*9.8*9.8 / 4 / np.log(2) *  1. /u.beam)
					if  840 *u.micron < wlen < 870*u.micron:
						conv_factor = 1./( np.pi*14.6*14.6 / 4 / np.log(2) * 1. /u.beam)
				elif 'JCMT' in telname and 'SCUBA' in instr:
					print("Beam inferred,",telname,instr,wlen)
					if  430 *u.micron < wlen < 460*u.micron:
						conv_factor = 1./( np.pi*7.8*7.8 / 4 / np.log(2) *  1. /u.beam)
					if  840 *u.micron < wlen < 870*u.micron:
						conv_factor = 1./( np.pi*13.8*13.8 / 4 / np.log(2) * 1. /u.beam)
				
	
			else:
				conv_factor = pixarea #/ u.pixel
			if "WISE" in telname:
				conv_factor = wise_dn_to_jy[band] * u.Jy/u.DN
				
	
		except Exception as e:
			print(e,"cant read unit")
			str_unit = " "
			real_unit = 1
	elif "ZUNITS" in header:
		try:
			str_unit = header.get("ZUNITS")	
			real_unit = u.Unit(str_unit)
			if "pix" in str_unit:
				conv_factor = 1. * u.pixel
			else:
				conv_factor = pixarea #/ u.pixel
			if "WISE" in telname:
				conv_factor = wise_dn_to_jy[band] * u.Jy/u.DN
		except Exception as e:
			print(e,"cant read unit")
			str_unit = " "
			real_unit = 1

	else:
		str_unit = " "
		real_unit = 1
	print("fl_unit:",str_unit,real_unit)
	
	if beam_params is not None:
		print("Beam override,",beam_params)
		bmaj,bmin = beam_params
		conv_factor = 1./((2*np.pi*bmin*bmaj / pixarea / (8*np.log(2))).decompose() * 1. / u.beam)

	#print("conv_factor:",conv_factor)
	
	resframe = pd.DataFrame(index=["sum","mean","stdv","rms","npix"])
	pht_to_return = []
	for regionfile in regionfilelist:
		regs = pyregion.open(regionfile)
		#print("regs")
	
		for reg in regs:
			regl = pyregion.ShapeList([reg])
			pht_aper = SkyCircularAperture(SkyCoord(reg.coord_list[0],reg.coord_list[1],frame='fk5',unit = (u.deg,u.deg)),reg.coord_list[2]*u.deg)
			aper_title = reg.attr[1]['text']
			if rname in aper_title:
				print(pht_aper)
				pht_to_return.append(pht_aper)
	cc = 1
	if "WISE" in telname:
		cc = conv_factor / u.pix
	
	data_w_units = data * (real_unit * cc)
	print(data_w_units.unit,"<<<<<<<<<<<")
	if data_w_units.unit == "Jy/pix":
		data_w_units = (data_w_units * u.pix / pixarea).to(u.MJy/u.sr)
		print(data_w_units.unit,"<<<<<<<<<<<<")
	return data_w_units,pht_to_return,{'wcs':wcs,'wlen':wlen,'telname':telname,'pixarea':pixarea},header



if __name__ == "__main__":
	description = """Plot Fits Files and overlay (Circular) DS9 Regions from files.
	Given a list of fits file images, this will iterate over them and: 
	1. Attempt to convert to MJy/sr unless units in Jy/beam
	2. Reproject the image so that pixels are aligned and North is up if, needed.
	3. Read in a list of DS9 region files, and from each region file get the first region whose title contains a text specified by the user.
	4. Make a nice WCS aligned plot.
	5. The title of the plot will contain a name specified by the user, and from the fits header the telescope name and the wavelength in microns.
	6. Display or save as pdf, the plot.

	So far works with WISE, Herschel, SCUBA1 and 2, MIPS data.
	"""
	parser = argparse.ArgumentParser(description="PlotNiceRegions")
	parser.add_argument('filelist',type=str,nargs='+',help="The list of fitsfiles to be plotted.")
	parser.add_argument('--rlist',type=str,nargs='+',required = True,help="The list of DS9 region files.")
	parser.add_argument('--rname',type=str,nargs='?',const='on',default='on',required=False,help="A string that is in the title of the region inside the region file. The first region to match is used. Default is 'on'")
	parser.add_argument('--gname',type=str,nargs=1,required=True,help="The name of the (galaxy) object, to be placed in the plot title.")	
	parser.add_argument('--pdf',action='store_true',required=False,help="If set, save plots as pdf files instead of showing them.")
	args = parser.parse_args()
	fits_filelist = args.filelist
	regfilelist = args.rlist
	rname = args.rname
	gname = args.gname[0]

#	fits_filelist = ["1728p590_ac51-w3-int-3_ra172.1345833_dec58.56194444_asec600.000.fits",
#	"1728p590_ac51-w4-int-3_ra172.1345833_dec58.56194444_asec600.000.fits",
#	"hpacs_30HPPJSMAPB_blue_1129_p5834_00_v1.0_1458303662690.fits",
#	"hpacs_30HPPJSMAPG_green_1129_p5834_00_v1.0_1458306378033.fits",
#	"hpacs_30HPPJSMAPR_1129_p5834_00_v1.0_1458306378736.fits",
#	"hspirepsw411_30pxmp_1128_p5833_1456975755847.fits",
#	"hspirepmw411_30pxmp_1128_p5834_1456975755620.fits"]
#	
#	reg_file1 = '25as.reg'
#	reg_file2 = '27as.reg'
#	regfilelist = [reg_file1,reg_file2]
#	
#	results=[get_files_regs_data(ff,regfilelist) for ff in fits_filelist]
	plt.rc('xtick',labelsize=18)
	plt.rc('ytick',labelsize=18)
	
	
	for ff in fits_filelist:
		rr = list(get_files_regs_data(ff,regfilelist,rname=rname))
		print(rr[1:])
		centralcoords = rr[1][0]
		wcs = rr[2]['wcs']
		if "Herschel" in r[2]['telname']:
			hdu = rr[-1]
			hdu = fits.PrimaryHDU(data = rr[0].value,header=rr[-1],do_not_scale_image_data=True)
			hdu = montage.reproject_hdu(hdu, north_aligned=True)
			rr[0] = hdu.data*rr[0].unit
			header = hdu.header
			wcs = WCS(header)
		fig = plt.figure(figsize = (8,8))
		ax = fig.add_subplot(111,projection = wcs)
		implot = ax.imshow(rr[0].value,origin='lower',cmap='jet')
		cenpix =  np.round(wcs.wcs_world2pix([[centralcoords.positions.ra.value,centralcoords.positions.dec.value]],0)[0])
		pixsize = np.sqrt(rr[2]['pixarea'])
		apradius = rr[1][0].r
		apradiuspix = apradius / pixsize
		ax.set_xlim(cenpix[0] - 4*apradiuspix,cenpix[0] + 4*apradiuspix)
		ax.set_ylim(cenpix[1] - 4*apradiuspix,cenpix[1] + 4*apradiuspix)
		ax.set_xlabel("")
		ax.set_ylabel("")
		ax.set_title(gname+" "+rr[2]["telname"].split(" ")[0]+" {:.0f}".format(rr[2]["wlen"].to(u.micron).value)+r"$\,\mathrm{\mu m}$",fontsize=22)
		for reg in rr[1]:
			regp = reg.to_pixel(wcs)
			regp.plot(ax=ax,color="white",lw=1.5)
	
		cbar = fig.colorbar(implot,ax=ax)
		cbar.ax.ticklabel_format(style='sci',scilimits=(-2,2))
		if rr[0].unit == "MJy/sr":
			cbar.set_label(r"Intensity [$\mathrm{MJy\,sr^{-1}}$]",size = 18)
		if rr[0].unit == "Jy/beam":
			cbar.set_label(r"Flux Density [$\mathrm{Jy\,beam^{-1}}$]",size = 18)
		
		ax = plt.gca()
		ra = ax.coords[0]
		ra.set_ticks(number = 3)
		ra.set_axislabel(" ")
		dec = ax.coords[1]
		dec.set_ticks(number=4)
		dec.set_ticklabel(rotation = "vertical",exclude_overlapping=True,va='top')
		dec.set_axislabel(" ")
		ax.tick_params(axis='both', which='both', direction='in',width = 2, color = 'white')
		if args.pdf:
			plt.savefig(gname+rr[2]["telname"].replace(" ","")+"{:.0f}".format(rr[2]["wlen"].to(u.micron).value)+".pdf")
		else:
			plt.show()
		plt.close()
