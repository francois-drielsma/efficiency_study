from amplitude_quantile import AmplitudeQuantile
from knn_density_estimator import kNNDensityEstimator
import numpy as np
import math
import ROOT
import array

# Generates a nDarray of points sampled from an uncorrelated beam distribution
#  -dim		Dimension of the space
#  -mom		Total momentum of particles [MeV/c]
#  -mass	Particle mass [MeV/c^2]
#  -eps		Normalised emittance [mm]
#  -beta	Betatron function [mm]
#  -seed	Random number generator seed
def beam_generator(dim, mom, mass, eps, beta, seed=None):

    # Only 2D and 4D transverse beams supported
    if dim != 2 and dim != 4:
	raise ValueError('Dimension not supported')

    # Set widths
    geoeps = mass*eps/mom
    sigx = math.sqrt(geoeps*beta)
    sigp = mass*eps/sigx

    # Set seed
    np.random.seed(seed)

    # Return a Gaussian beam of n particles
    data = np.ndarray((n, dim))
    for i in range(n):
	# Sample from a Gaussian distribution
	data[i, 0] = np.random.normal(0, sigx)
	data[i, 1] = np.random.normal(0, sigp)
	if dim == 4:
	    data[i, 2] = np.random.normal(0, sigx)
	    data[i, 3] = np.random.normal(0, sigp)
	
    return data

# Returns the array of amplitudes associated with the data sample
#  -data	Input data
#  -mass	Particle mass [MeV/c^2]
#  -eps		Normalised emittance [mm]
def amplitudes(data, mass, eps):

    # Calculate the covariance matrix, normalised emittance
    dim = len(data[0])
    cov = np.cov(data.transpose())
    eps = pow(np.linalg.det(cov), 1./dim)/mass

    # Calculate the amplitudes
    n = len(data)
    amps = np.empty(n)
    invcov = np.linalg.inv(cov)
    for i in range(n):
	v = np.matrix(data[i])
	amps[i] = eps*(v*invcov*v.transpose())[0, 0]

    return amps

# Prints the amplitude quantile
#  -quantile	Amplitude quantile object
#  -bsr		Flag to use the boostrap resampling to evaluate uncertainties
def print_quantile(quantile):
    print "%d%%-amplitude: %0.3f +/- %0.3f mm"\
	%(100*frac, quantile.value(), quantile.error())
    print "%d%%-emittance: %0.3f +/- %0.3f x 10^6 mm^2(MeV/c)^2"\
	%(100*frac, quantile.volume()/1e6, quantile.volume_error()/1e6)

# Draws the amplitude distribution as a histogram
# -amps		Array of amplitudes [mm]
def print_amplitudes(amps):

    # Initialize the histogram
    amp_hist = ROOT.TH1F("amp", ";A_{#perp}  [mm];PDF [/5mm]", 12, 0, 60)
    amp_hist.FillN(n, amps, np.full(n, 1.))
    amp_hist.Scale(1./amp_hist.GetEntries())

    # Print the canvas
    c = ROOT.TCanvas("c", "c", 1200, 800)
    amp_hist.Draw()
    c.SaveAs("amps.pdf")
    del c

# Prints the poincare section of the kNN density estimator (n*(n+1)/2 sections)
#  -knn		kNN density estimator object
def print_poincare(knn):

    # Create a (d-1)x(d-1) canvas
    ROOT.gStyle.SetOptStat(0)
    canv = ROOT.TCanvas('c', 'c', 1200, 1200)
    canv.Divide(knn.dim-1, knn.dim-1, 0., 0.)
    sections = []
    for i in range (knn.dim-1):
	for j in range(i+1, knn.dim):
	    cid = (j-1)*3+i+1
	    canv.cd(cid)
	    sections.append(knn.section((i, j)))
	    sections[-1].Draw("CONTZ")

    canv.SaveAs("poincare.pdf")
    del canv

# Prints the kNN density profile (i.e. the quantiles of the inverse PDF)
#  -knn		kNN density estimator object
def print_profile(knn, bsr_flag):

    # Print the density profile
    profile = knn.profile(bsr=bsr_flag)
    c = ROOT.TCanvas("c", "c", 1200, 800)
    profile.SetFillStyle(1000)
    profile.SetFillColorAlpha(ROOT.kBlack, .2)
    profile.Draw("ALE3")
    c.SaveAs("density_profile.pdf")
    del c

# Main part of the code, calls the subroutines
if __name__ == "__main__":

    # General parameters
    ROOT.gROOT.SetBatch()
    dim = 4			# Transverse space dimension
    mom = 140			# [MeV/c], Total momentum
    mass = 105.66		# [MeV/c^2], Particle mass
    eps = 6			# [mm], Normalised emittance
    beta = 500			# [mm], Betatron function
    n = 100000			# Number of particles to generate
    frac = 0.09			# Fraction of the beam to look at

    # Build the input beam
    data = beam_generator(dim, mom, mass, eps, beta)

    # Calculate the amplitudes
    amps = amplitudes(data, mass, eps)

    # Draw the amplitude PDF
    print_amplitudes(amps)

    # Get the quantile and the corresponding fractional emittance
    #  - Set the flag to True to use the bootsrapped uncertainties
    quantile = AmplitudeQuantile(dim, amps, frac, bsr=False)
    print_quantile(quantile)

    # Initialize the kNN estimator on the beam
    #  - Set the flag to True to use the distances in the metric of the covariance matrix.
    #  - Specify the amount of available cores, the profile estimate 
    #    is easily parallelisable (KD tree is thread safe). The time scales
    #    down linearly with the amount of available cores.
    nthreads = 2
    knn = kNNDensityEstimator(data, True, nthreads)

    # Represent the Poincare sections
    print_poincare(knn)

    # Print the density profile
    #  - Set the flag to True to use the bootsrapped uncertainties
    print_profile(knn, bsr_flag=False)

