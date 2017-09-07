import warnings

#
import numpy as np

#
import astropy.units as u
import astropy.constants as c

def feature_builder(dryads, feature):

    # feature to plot
    if isinstance(feature, str) and (feature.lower() == 'nobs'):
        feature = 'nobs'
        data_feature = dryads.nobs.copy()
        data_feature[dryads.data_feature == 0.] = np.nan
        unit_feature = u.dimensionless_unscaled
    elif isinstance(feature, str) and (feature.lower() == 'radius'):
        feature = 'radius'
        data_feature = dryads.radius.copy()
        data_feature[data_feature == 0.] = np.nan
        unit_feature = dryads.unit_radius
    elif isinstance(feature, str) and (feature.lower() == 'mass'):
        feature = 'mass'
        data_feature = dryads.mass.copy()
        data_feature[data_feature == 0.] = np.nan
        unit_feature = dryads.unit_mass
    elif isinstance(feature, str) and (feature.lower() == 'sigma'):
        feature = 'sigma'
        data_feature = dryads.avgsigma.copy()
        data_feature[data_feature == 0.] = np.nan
        unit_feature = dryads.unit_avgsigma
    elif isinstance(feature, str) and (feature.lower() == 'pressure'):
        feature = 'pressure'
        data_feature = dryads.pressure.copy()
        data_feature[data_feature == 0.] = np.nan
        unit_feature = dryads.unit_pressure
    elif isinstance(feature, str) and (feature.lower() == 'virial'):
        feature = 'virial'
        data_feature = dryads.virial.copy()
        data_feature[data_feature == 0.] = np.nan
        unit_feature = u.dimensionless_unscaled
    elif isinstance(feature, str) and (feature.lower() == 'contours'):  ## cbins
        feature = 'contours'
        data_feature = dryads.mass_baseline.copy()
        data_feature[data_feature == 0.] = np.nan
        unit_feature = dryads.unit_mass
    else:
        raise ValueError('Choose to plot one of "nobs", "radius", "mass", "sigma", "pressure", "virial", and "contours".')


    return feature, data_feature, unit_feature
