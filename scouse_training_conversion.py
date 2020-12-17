import astropy.table 
from astropy.io import fits
from astropy.table import Table 
from spectral_cube import SpectralCube
import pickle 

tbl = Table.read('/orange/adamginsburg/sgrb2/madelinehall/refine_grid/cutout/HC3N_TP_7m_12m_feather_cutout/stage_6/best_fit_solutions.dat', format = 'ascii')
cube = SpectralCube.read('/home/madelinehall/sgrb2files/bin/fitsfiles/HC3N_TP_7m_12m_feather_cutout.fits')

#create header 
header = fits.getheader('/home/madelinehall/sgrb2files/bin/fitsfiles/HC3N_TP_7m_12m_feather_cutout.fits', 0) 

#create a guasspy+ training dictionary 
data = {'data_list':[], 'location':[], 'amplitudes':[], 'index':[], 'error':[],
	 'best_fit_rchi2':[], 'pvalue':[],'fwhms':[], 'means':[], 
	 'signal_ranges':[], 'x_values':[], 'header':[]}

header = cube.header

header_dict = {'SIMPLE': 'T', 'BITPIX': -64, 'NAXIS': 3, 'NAXIS1': 100, 'NAXIS2': 100, 'NAXIS3': 201, 
               'BMAJ': 0.001031441357401, 'BMIN': 0.0005507485071818, 'BPA': 70.47904968262, 'BTYPE': 'Intensity',
               'OBJECT': 'SgrB2', 'BUNIT': 'K', 'ALTRVAL': -200000.0000769, 'ALTRPIX': 1.0, 'VELREF': 257,
               'TELESCOP': 'ALMA', 'OBSERVER': 'kelfavich', 'OBSRA': 266.8308333333, 'OBSDEC': -28.39138888889, 
               'DATE': '2015-09-14T06:47:47.344000', 'ORIGIN': 'CASA 4.4.0-REL (r33623)',
               'BEAM': 'Beam: BMAJ=3.7131888866436 arcsec BMIN=1.98269462585448 arcsec &', 'CONTINUE': 'BPA=70.47904968262 deg',
               'SLICE':'[[(None, None, None), (263, 363, None), (439, 539, None)]]', 'MJDREFI': 0.0, 'MJDREFF': 0.0 ,
               'WCSAXES': 3, 'CRPIX1': 52.0, 'CRPIX2': 213.0, 'CRPIX3': 1.0, 'CDELT1': -0.000125, 'CDELT2': 0.000125, 
               'CDELT3': 1.99999999997, 'CUNIT1': 'deg', 'CUNIT2': 'deg', 'CUNIT3': 'km s-1', 'CTYPE1': 'RA---SIN',
               'CTYPE2': 'DEC---SIN','CTYPE3': 'VRAD', 'CRVAL1': 266.830833333, 'CRVAL2': -28.3913888889, 'CRVAL3': -200.00000007701,
               'PV2_1': 0.0, 'PV2_2': 0.0, 'LONPOLE': 180, 'LATPOLE': -28.3913888889, 'RESTFRQ': 90979020000.0, 
               'TIMESYS':'UTC', 'MJDREF': 0.0, 'DATE-OBS': '2014-07-01T06:00:41.183999', 'MJD-OBS': 56839.2504767,
               'OBSGEO-X': 2225142.18027, 'OBSGEO-Y': -5440307.37035, 'OBSGEO-Z': -2481029.85187, 'RADESYS':'FK5', 
               'EQUINOX': 2000.0, 'SPECSYS': 'LSRK'}  

#location, amplitude, header
for row in tbl:
    location = (row['y'], row['x'])
    data['location'].append(location)
    data['header'].append(header_dict)
    if location in data['location']:
        index = data['location'].index(location)
        data['amplitudes'].append(tbl['amplitude'])
    else: 
        data['location'].append(location)
        data['amplitudes'].append(tbl['amplitude'])

#index
for row in tbl:
    location = (row['y'], row['x'])
    index = row['y']*row['x']
    if location in data['location']:
        data['index'].append(index)
    else: 
        data['index'].append(index)

#error, best fit rchi2, pvalue, fwhms, mean, x_values
for row in tbl:
    location = (row['y'], row['x'])
    index = row['y']*row['x']
    pvalue = 0.5
    x_values = cube.spectral_axis
    data['location'].append(location)
    if location in data['location']:
        data['amplitudes'].append(tbl['amplitude'])
        data['data_list'].append(cube[:,int(row['y']),int(row['x']) ])
        data['error'].append(tbl['rms']) 
        data['best_fit_rchi2'].append(tbl['chi2'])
        data['pvalue'].append(pvalue)
        data['fwhms'].append(tbl['width'])
        data['means'].append(tbl['shift'])
        data['x_values'].append(x_values)
        data['signal_ranges'].append(len(x_values))
    else:
        data['amplitudes'].append(tbl['amplitude'])
        data['data_list'].append(cube[:,int(row['y']),int(row['x']) ])
        data['error'].append(tbl['rms'])
        data['best_fit_rchi2'].append(tbl['chi2']) 
        data['pvalue'].append(pvalue)
        data['fwhms'].append(tbl['width'])
        data['means'].append(tbl['shift'])
        data['x_values'].append(x_values)
        data['signal_ranges'].append(len(x_values))
        
pickle.dump( data, open( "cutout_data_dictionary.p", "wb" ) )
