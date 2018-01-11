import numbers
import warnings

#
import numpy as np

#
import matplotlib.pyplot as plt
import styles

#
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astrodendro import Dendrogram



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
            mask[dendro[i].indices(subtree = False)] = True
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
            raise ValueError('The header should include the angular pix size. No irregular proj.')


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



    def climber(self, data_sigma, unit_coldens = u.cm**-2.):

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


    def plotter(self):

        from plot import DryadsPlotter

        return DryadsPlotter(self)
