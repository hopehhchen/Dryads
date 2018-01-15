#
import numpy as np

#
import astropy.units as u
import astropy.constants as c

#
import styles

# observations
## main beam eff of FCRAO is 0.45 for 12CO and 0.49 for 13CO.
observation = {'mb_eff': 0.49}

# conversion
## B13 and A_W13 are from the equation in Pineda+08: AV = B13*W13+A_W13.
## N_over_AV converts AV magnitudes to cm^-2: N(H2) = N_overAV * A_V.
## magnitudes are not expressed explicitly here.
conversion = {'B13': 0.37*(u.K*u.km/u.s)**-1., 'A_W13': 1.45,
              'N_over_AV': 9.4e20*u.cm**-2.}
