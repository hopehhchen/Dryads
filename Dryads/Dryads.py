import numbers
import warnings

#
import numpy as np

#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import styles

#
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.modeling import models, fitting
from astrodendro import Dendrogram

#
from stat_tools import *



class Dryads(object):
    '''
    A toolkit for exploring the dendrogram.

    Inputs
    ------
    data: numpy.ndarray
          This is used to build the dendrogram and the "contour masks."  In
          theory, the array can be from 1D to 3D, but the package at the moment
          is built with a 2D array in mind.

    header: fits.header.Header
            The header that includes astrometry information for data.

    distance: numbers.Number or astropy.Quantity
              The physical distance to the supposed astronomical object.  If no
              unit, it is assumed in pc.

    cbins: int or 1D numpy.ndarray
           This is used to determine the contour bins, which in turns determine
           the contour levels where physical quantities are calculated.

    cscale: str; the default is 'linear'.
            The contour scale.  This has to be either 'linear' or 'log', and is
            used when cbins is an integer, to calculate the contour bins.

    *args, **kargs
    '''


    def __init__(self, data, header, distance, cbins, cscale = 'linear', *args, **kargs):


        # dendro
        ## data (used to calculate the dendro and the contours)
        self.data = data  ## push the type checking to astrodendro
        ## calculate the dendrogram
        dendro = Dendrogram.compute(data, *args, **kargs)
        self.dendro = dendro
        ## make the dendrogram mask list (ordered by the structure id)
        mask0 = np.zeros(data.shape,
                         dtype = bool)
        maskl_dendro = []
        for i in range(len(dendro)):
            mask = mask0.copy()
            mask[dendro[i].indices(subtree = False)] = True ###
            maskl_dendro.append(mask)
        maskl_dendro = np.array(maskl_dendro)
        self.maskl_dendro = maskl_dendro



        # contours
        ## cbins
        if isinstance(cbins, numbers.Number) and (cbins%1 == 0):
            if isinstance(cscale, str) and (cscale.lower() == 'linear'):
                cbins = int(cbins)
                cbins = np.linspace(max(np.nanmin(data),
                                        dendro.params['min_value']),
                                    np.nanmax(data),
                                    cbins)
            elif isinstance(cscale, str) and (cscale.lower() == 'log'):
                cbins = int(cbins)
                cbins = np.logspace(np.log10(max(np.nanmin(data),
                                                 dendro.params['min_value'])),
                                    np.log10(np.nanmax(data)),
                                    cbins)
            else:
                raise ValueError('Choose cscale from "linear" and "log".')
        elif isinstance(cbins, np.ndarray) and\
        issubclass(cbins.dtype.type, numbers.Number) and\
        (cbins.ndim == 1):
            warnings.warn('Custom cbins; cscale is not used.')
            cbins = cbins
        else:
            raise ValueError('Input an 1D numpy array or a number for cbins.')
        self.cbins = cbins
        ## make the contour mask list (ascending)
        maskl_contours = np.array([(data >= cbins[i]) for i in range(len(cbins))])
        self.maskl_contours = maskl_contours


        # header
        ## header
        if isinstance(header, fits.header.Header):
            self.header = header
        else:
            raise ValueError('The header should be compatible with astropy.io.fits.')
        ## image scale (angular)
        if ('CDELT1' in header) and (abs(header['CDELT1']) == abs(header['CDELT2'])):
            imgscale_ang = abs(header['CDELT1'])
            self.imgscale_ang = imgscale_ang  ## assuming degrees
        elif ('CD1_1' in header) and (abs(header['CD1_1']) == abs(header['CD2_2'])):
            imgscale_ang = abs(header['CD1_1'])
            self.imgscale_ang = imgscale_ang  ## assuming degrees
        else:
            imgscale_ang = np.nan
            self.imgscale_ang = imgscale_ang
            #raise ValueError('The header should include the angular pix size. No irregular proj.')


        # distance
        ## distance
        if issubclass(distance.dtype.type, numbers.Number):

            ## distance unit
            if hasattr(distance, 'unit'):
                print distance.value
                distance, unit_distance = distance.value, distance.unit
            else:
                distance, unit_distance = distance, u.pc

        self.distance = distance
        self.unit_distance = unit_distance
        ## image scale (physical)
        imgscale = np.radians(imgscale_ang)*distance*unit_distance
        self.imgscale = imgscale.value
        self.unit_imgscale = imgscale.unit


        # prepare for climbing
        self.results = None
        self.method = None



    def climber(self, data_sigma, unit_coldens = u.cm**-2.):

        #
        #stat = statBasic2D(mapTpeak[mask], (meshy[mask], meshx[mask]))
        #stat.calculate()



        # data_coldens, data_sigma
        ## data_coldens
        self.data_coldens = self.data ## allowing data_coldens == data for now.
        ## coldens unit
        if isinstance(unit_coldens, u.core.CompositeUnit):
            self.unit_coldens = unit_coldens
        else:
            raise ValueError('Use astropy.units for unit_coldens.')


        ## data_sigma
        if isinstance(data_sigma, np.ndarray) and (data_sigma.shape == self.data.shape):

            ## sigma unit
            if hasattr(data_sigma, 'unit'):
                self.data_sigma, self.unit_sigma = data_sigma.value, data_sigma.unit
            else:
                self.data_sigma = data_sigma
                self.unit_sigma = u.km/u.s

        else:
            raise ValueError('The inputs should be numpy arrays with the same shape.')


        # nobs
        nobs = np.tensordot(self.maskl_dendro.astype(float),
                            self.maskl_contours.astype(float),
                            axes = ([1, 2], [1, 2]))  ## number of pix
        self.nobs = nobs


        # radius
        radius = np.sqrt(nobs*(self.imgscale*self.unit_imgscale)**2./np.pi)
        self.radius = radius.value
        self.unit_radius = radius.unit


        # mass
        coldens = self.data_coldens.copy()  ## replace nan and inf with 0; different from nan_to_num.
        coldens[~np.isfinite(coldens)] = 0.
        mass = np.tensordot(self.maskl_dendro*coldens,
                            self.maskl_contours,
                            axes = ([1, 2], [1, 2]))
        mass_baseline = np.meshgrid(self.cbins, range(len(self.dendro)))[0]
        mass = mass-mass_baseline*nobs  ## baseline adjusted
        try:  ## try whether unit_coldens is in [mass]/[length]^2 or 1/[length]^2
            mass = (mass*self.unit_coldens*(2.37*u.u)*(self.imgscale*self.unit_imgscale)**2.).to(u.Msun)
        except UnitConversionError:
            mass = (mass*self.unit_coldens*(self.imgscale*self.unit_imgscale)**2.).to(u.Msun)
        self.mass = mass.value
        self.unit_mass = mass.unit
        self.mass_baseline = mass_baseline



        # sigma
        sigma = self.data_sigma.copy()  ## replace nan and inf with 0; different from nan_to_num.
        sigma[~np.isfinite(sigma)] = 0.
        sigma = np.tensordot(self.maskl_dendro*sigma,
                             self.maskl_contours,
                             axes = ([1, 2], [1, 2]))
        sigma = sigma/nobs  ## take the average
        self.avgsigma = sigma
        self.unit_avgsigma = self.unit_sigma


        # pressure
        pressure = (self.mass*self.unit_mass)\
                   /(4./3.*np.pi*(self.radius*self.unit_radius)**3.)\
                   *(self.avgsigma*self.unit_sigma)**2.
        pressure = pressure.to(u.Ba)
        self.pressure = pressure.value
        self.unit_pressure = pressure.unit


        # virial
        virial = 5.*(self.avgsigma*self.unit_sigma)**2.*(self.radius*self.unit_radius)\
                 /(1.*c.G*(self.mass*self.unit_mass))  ## assuming constant density profile (see Bertoldi+McKee92)
        virial = virial.decompose().value
        self.virial = virial


        return self


    def climber3D(self):
        '''
        Climber Feature for parsing through 3D cubes.  The properties ave
        calculated using the statistics provided by the astrodendro package.
        The constants used for converting flux to density/mass are stored
        in `constants.py`.

        The feature is in development.  For convenience at the current stage,
        the cube is assumed to have the third dimension in velocity units
        (intead of frequency units).  The cubes is also assumed to be in the
        main beam temperature units.  The main beam efficiency is given in
        `constants.py`.  (For cubes that are already corrected, please set
        the efficiency to 1 in `constants.py`.)
        '''

        if self.data.ndim != 3:
            raise ValueError('`climber3D` only applies to 3D data cubes.')

        # velocity scale from the header if it is a cube.
        if 'CDELT3' in self.header:
            velocity_scale = abs(self.header['CDELT3'])*observation['vel_unit']
        elif 'CD3_3' in self.header:
            velocity_scale = abs(self.header['CD3_3'])*observation['vel_unit']
        else:
            velocity_scale = np.nan
            #raise ValueError('The header needs to contain info on the 3rd axis.')

        # initiate PPVStatistic with pixel values. ##

        # calculate dendro masks with full extent.
        maskl_dendro_full = []
        for i in range(len(self.dendro)):
            #mask = mask0.copy()
            #mask[dendro[i].indices(subtree = False)] = True ###
            maskl_dendro_full.append(self.dendro[i].get_mask())
        maskl_dendro_full = np.array(maskl_dendro_full)
        self.maskl_dendro_full = maskl_dendro_full

        # loop through all combinations.
        results = {'nobs': np.zeros([len(self.maskl_dendro),
                                     len(self.maskl_contours)])*np.nan,
                   'radius': np.zeros([len(self.maskl_dendro),
                                       len(self.maskl_contours)])*np.nan,
                   'unit_radius': None,
                   'mass': np.zeros([len(self.maskl_dendro),
                                     len(self.maskl_contours)])*np.nan,
                   'unit_mass': None,
                   'sigma': np.zeros([len(self.maskl_dendro),
                                      len(self.maskl_contours)])*np.nan,
                   'unit_sigma': None,
                   'pressure': np.zeros([len(self.maskl_dendro),
                                         len(self.maskl_contours)])*np.nan,
                   'unit_pressure': None,
                   'virial': np.zeros([len(self.maskl_dendro),
                                       len(self.maskl_contours)])*np.nan}
        # looping
        for i in range(len(self.maskl_dendro)):
            for j in range(len(self.maskl_contours)):

                # create the intersection mask.
                mask = (self.maskl_dendro_full[i] & self.maskl_contours[j])
                mask_noLeaves = (self.maskl_dendro[i] & self.maskl_contours[j])

                # start calculation only if the two masks intersect.
                if (np.sum(mask_noLeaves) != 0):

                    # initiate with indices and values.
                    meshgrid = np.meshgrid(np.arange(mask.shape[0]),
                                           np.arange(mask.shape[1]),
                                           np.arange(mask.shape[2]),
                                           indexing = 'ij') ## this will be z, y, x (0, 1, 2)
                    ## indices
                    indices = (meshgrid[0][mask],
                               meshgrid[1][mask],
                               meshgrid[2][mask])
                    ## values
                    values = self.data[indices]
                    ## statBasic3D
                    spatial_scale = self.imgscale*self.unit_imgscale
                    #velocity_scale = velocity_scale
                    stat = statBasic3D(values, indices,
                                       metadata = {'spatial_scale': None,
                                                   'velocity_scale': None})
                    stat.calculate()

                    # nobs
                    nobs = np.sum(mask)
                    results['nobs'][i, j] = nobs

                    # radius
                    radius = stat.radius.value
                    results['radius'][i, j] = radius ## pix

                    # mass (total flux)
                    mass = np.sum(values)
                    results['mass'][i, j] = mass

                    # sigma
                    sigma = stat.v_rms.value
                    results['sigma'][i, j] = sigma ## pix in the v-direction

                    # pressure (flux density * sigma**2)
                    #pressure = mass/(4./3.*np.pi*radius**3.) * sigma**2.
                    #results['pressure'][i, j] = pressure

                    # virial
                    virial = 5.*sigma**2.*radius/mass
                    results['virial'][i, j] = virial

                else:

                    continue

        # pressure (flux density * sigma**2)
        pressure_result = results['mass']/(4./3.*np.pi*results['radius']**3.)\
                          *results['sigma']**2.
        results['pressure'] = pressure_result

        # virial
        virial_result = 5.*results['sigma']**2.*results['radius']/results['mass']
        results['virial'] = virial_result

        # return
        self.nobs = results['nobs']
        self.radius = results['radius']
        self.unit_radius = results['unit_radius']
        self.mass = results['mass']
        self.unit_mass = results['unit_mass']
        self.sigma = results['sigma']
        self.unit_sigma = results['unit_sigma']
        self.pressure = results['pressure']
        self.unit_pressure = results['unit_pressure']
        self.virial = results['virial']


        self.method = 'clipping'
        self.results = {'clipping': results}


    def extrapolator3D(self):
        '''
        Based on the results from clipping, the extrapolator extrapolates
        properties to the zero point.
        '''

        if (self.method != 'clipping') or ('extrapolation' in self.results):

            raise ValueError('Need to run climber first/use switch to change scheme.')

        # loop through all combinations.
        results = {'nobs': np.zeros([len(self.maskl_dendro),
                                     len(self.maskl_contours)])*np.nan,
                   'radius': np.zeros([len(self.maskl_dendro),
                                       len(self.maskl_contours)])*np.nan,
                   'unit_radius': None,
                   'mass': np.zeros([len(self.maskl_dendro),
                                     len(self.maskl_contours)])*np.nan,
                   'unit_mass': None,
                   'sigma': np.zeros([len(self.maskl_dendro),
                                      len(self.maskl_contours)])*np.nan,
                   'unit_sigma': None,
                   'pressure': np.zeros([len(self.maskl_dendro),
                                         len(self.maskl_contours)])*np.nan,
                   'unit_pressure': None,
                   'virial': np.zeros([len(self.maskl_dendro),
                                       len(self.maskl_contours)])*np.nan}
        badfits = {'flag_radius': np.zeros([len(self.maskl_dendro),
                                            len(self.maskl_contours)],
                                           dtype = bool),
                   'flag_mass': np.zeros([len(self.maskl_dendro),
                                          len(self.maskl_contours)],
                                         dtype = bool),
                   'flag_sigma': np.zeros([len(self.maskl_dendro),
                                           len(self.maskl_contours)],
                                          dtype = bool)}
        # looping
        for i in range(len(self.maskl_dendro)):

            # grab the dendro feature
            dendro_feature = self.dendro[i]

            for j in range(len(self.maskl_contours)):

                # masking
                mask_threshold = max(dendro_feature.vmin, self.cbins[j])
                mask_section = (self.cbins >= mask_threshold)&\
                               (self.nobs[i] >= 4.)

                if np.sum(mask_section) >= 3.:
                    # fitting
                    ## radius
                    xfit, yfit = self.cbins[mask_section], self.radius[i][mask_section]
                    wfit = self.nobs[i][mask_section]
                    t_init = models.Polynomial1D(1)
                    fit_t = fitting.LinearLSQFitter()
                    t_fit = fit_t(t_init, xfit, yfit, weights = wfit)
                    ##
                    radius = t_fit.c0.value
                    if radius < yfit.max():
                        badfits[i, j] = True
                    results['radius'][i, j] = radius

                    ## sigma
                    xfit, yfit = self.cbins[mask_section], self.sigma[i][mask_section]
                    wfit = self.nobs[i][mask_section]
                    t_init = models.Polynomial1D(1)
                    fit_t = fitting.LinearLSQFitter()
                    t_fit = fit_t(t_init, xfit, yfit, weights = wfit)
                    ##
                    sigma = t_fit.c0.value
                    if sigma < yfit.max():
                        badfits[i, j] = True
                    results['sigma'][i, j] = sigma

                    ## mass
                    xfit, yfit = self.cbins[mask_section], self.mass[i][mask_section]
                    wfit = self.nobs[i][mask_section]
                    t_init = models.Polynomial1D(2)
                    fit_t = fitting.LinearLSQFitter()
                    t_fit = fit_t(t_init, xfit, yfit, weights = wfit)
                    ##
                    mass = t_fit.c0.value
                    if mass < yfit.max():
                        t_init = models.Polynomial1D(1)
                        fit_t = fitting.LinearLSQFitter()
                        t_fit = fit_t(t_init, xfit, yfit, weights = wfit)

                        mass = t_fit.c0.value
                        if sigma < yfit.max():
                            badfits[i, j] = True
                    results['mass'][i, j] = mass

                else:
                    continue

        # pressure
        pressure_result = results['mass']/(4./3.*np.pi*results['radius']**3.)\
                          *results['sigma']**2.
        results['pressure'] = pressure_result

        # virial
        virial_result = 5.*results['sigma']**2.*results['radius']/results['mass']
        results['virial'] = virial_result

        # return
        self.radius = results['radius']
        self.unit_radius = results['unit_radius']
        self.mass = results['mass']
        self.unit_mass = results['unit_mass']
        self.sigma = results['sigma']
        self.unit_sigma = results['unit_sigma']
        self.pressure = results['pressure']
        self.unit_pressure = results['unit_pressure']
        self.virial = results['virial']


        self.method = 'extrapolation'
        self.results['extrapolation'] = results


    def switch(self):
        '''
        The feature that allows user to switch between schemes for plotting.
        '''

        # decision
        if self.method == 'clipping':
            results = self.results['extrapolation']
            self.method = 'extrapolation'
        elif self.method == 'extrapolation':
            results = self.results['clipping']
            self.method = 'clipping'
        else:
            raise ValueError('Need to run climber and extrapolator.')

        # switch
        self.radius = results['radius']
        self.unit_radius = results['unit_radius']
        self.mass = results['mass']
        self.unit_mass = results['unit_mass']
        self.sigma = results['sigma']
        self.unit_sigma = results['unit_sigma']
        self.pressure = results['pressure']
        self.unit_pressure = results['unit_pressure']
        self.virial = results['virial']



    def plotter(self):

        from plot import DryadsPlotter

        return DryadsPlotter(self)
