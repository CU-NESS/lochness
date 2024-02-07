'''
Description: Some utility functions for rotating maps using healpy.
			 These have been taken as functions from CU-NESS/perses so that 
			 installation of the latter is not necessary to run
			 LOCHNESS.
'''
import healpy as hp

def spherical_rotator(theta, phi, psi, deg=True):
	"""
	Generates a healpy-based Rotator object which rotates points such that the
	zenith direction before rotation becomes the (lat,lon)=(90-theta,phi)
	direction after rotation.

	theta: the colatitude (units determined by deg argument) of the direction
		   to which zenith (before rotation) should be rotated
	phi: the longitude (units determined by deg argument) of the direction to
		 which zenith (before rotation) should be rotated
	psi: angle (units determined by deg argument) through which sphere should
		 be rotated after zenith is rotated to (theta, phi)
	deg: if True (default), angles should be given in degrees

	returns: healpy.rotator.Rotator object capable of performing rotation which
		     puts the zenith to (theta, phi) and applies a rotation of psi
		     about that direction
	"""
	rot_zprime = hp.rotator.Rotator(rot=(-phi, 0, 0), deg=deg, eulertype='y')
	rot_yprime = hp.rotator.Rotator(rot=(0, theta, 0), deg=deg, eulertype='y')
	rot_z = hp.rotator.Rotator(rot=(psi, 0, 0), deg=deg, eulertype='y')
	return rot_zprime * rot_yprime * rot_z
    
def rotator_for_spinning(angle, degrees=True):
	"""
	Generates a healpy Rotator object that would spin around the current pole
	by given angle.

	angle: the angle through which to spin, equivalent to angle in spin_maps
		   function
	degrees: boolean describing whether angle is in degrees or not

	returns: healpy Rotator object which implements the spinning
	"""
	return spherical_rotator(0, 0, -angle, deg=degrees)
