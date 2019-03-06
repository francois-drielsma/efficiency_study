import math
import numpy as np
import scipy as sp
import scipy.stats

class AmplitudeQuantile(object):

    # Initializes the amplitude quantile object
    #  -dim     Dimension of the space
    #  -amps    List of amplitudes (array of scalars)
    #  -frac    Fraction of the beam below the quantile
    #  -norm    Normalisation factor (i.e. transmission, default to 1.)
    #  -bsr     Use bootstrap resampling to evaluate uncertainties
    #  -bsr_n   Number of iterations of the bootstrap resampling
    def __init__(self, dim, norm=1., bsr=False, bsr_n=100):

        # Set the characteristics
        self.n = 0
        self.dim = float(dim)
        self.norm = norm
        self.bsr = bsr
        self.bsr_n = bsr_n

        # Set the error
        self.val = 0.
        self.err = 0.

    # Resets the quantile fraction
    #  -amps    List of amplitudes (array of scalars)
    #  -frac    Fraction of the beam below the quantile
    def set_quantile(self, amps, frac):
        self.n = len(amps)
        self.frac = frac
        if self.frac > self.norm:
            raise ValueError('The requested fraction cannot be larger than the transmission')

        self.set_value(amps)
        self.set_uncertainty(amps)

    # Returns the amplitude quantile
    def value(self):
        return self.val

    # Returns the statistical uncertainty on the quantile
    def error(self):
        return self.err

    # Returns the volume of the ellipse corresponding to the alpha-amplitude
    #  -mass    Particle mass in MeV/c (default to muon mass)
    def volume(self, mass=105.66):
        return pow(math.pi*mass*self.val, self.dim/2)/sp.special.gamma(self.dim/2+1)

    # Returns the theoretical statistical uncertainty on the volume
    #  -mass    Particle mass in MeV/c (default to muon mass)
    def volume_error(self, mass=105.66):
        return self.volume(mass)*(self.dim/2)*self.err/self.val

    # Returns the squared radius below which a fraction frac of amplitudes lie
    def chi2_quantile(self):
        return sp.stats.chi2.ppf(self.frac, self.dim)

    # Returns the value of the chi squared distribution in x
    def chi2_dist(self, x):
        return sp.stats.chi2.pdf(x, self.dim)

    # Sets the value of the quantile
    def set_value(self, amps):
        self.val = np.quantile(amps, self.frac/self.norm)

    # Sets the statistical uncertainty on the quantile
    #  -amps    List of amplitudes (array of scalars)
    def set_uncertainty(self, amps):
        if not self.bsr:
            R2 = self.chi2_quantile()
            self.err = self.val*math.sqrt(self.frac*(1-self.frac)/self.n)/R2/self.chi2_dist(R2)
        else:
            quantiles = np.empty(self.bsr_n)
            for k in range(self.bsr_n):
                sample = np.empty(self.n)
                for i in range(self.n):
                    sample[i] = np.random.choice(amps)

                quantiles[k] = np.quantile(sample, self.frac/self.norm)
        
            self.err = np.std(quantiles)
   
