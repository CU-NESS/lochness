"""
File: lochness/LOCHNESS.py
Authors: Joshua J. Hibbard/Valerie Wong with code from Neil Bassett and Keith Tauscher.
Date: Jan 2024

Description: The Lunar Observatory Code in Healpix by the NESS team. Useful for 
			 determining the altitude and azimuth of celestial bodies at various times
			 from landing sites on the lunar surface. Can also be used to rotate healpy
			 maps of the galaxy into the beam frame of the landing site on the lunar surface.
			 By default LOCHNESS currently calculates altitude and azimuth for the following
			 celestial bodies: Sun, Jupiter, and Earth.
"""

import numpy as np
import healpy as hp
import spiceypy as spice
from datetime import datetime
from lochness.map_utilities import spherical_rotator,\
	rotator_for_spinning

MOON_RAD = 1736.0     # km (polar radius)  # avg rad 1737.4 km
del_angle = 1e-6 # arbitrary small angle
# Define naif numbers
sun_naif = 10
venus_naif = 2
earth_naif = 399
moon_naif = 301
jupiter_naif = 5 #599
saturn_naif = 6  #699

class LOCHNESS(object):

	def __init__(self, path_to_spice_kernels, time_list, landing_site_coordinates,\
		MOON_ALT=0., filename_to_save=None, antenna_psi_angle=0., galaxy_map=None,\
		left_handed_mollview_maps=True, run_test_galaxy_map=False):
		'''
		path_to_spice_kernels: file path to the kernels downloaded for use in 
							   spicepy. These must include 'naif0012.tls',
							   'de440.bsp', 'pck00011.tpc' for celestial bodies
							   and 'moon_080317.tf', 'moon_assoc_pa.tf',
							   and 'moon_pa_de421_1900-2050.bpc' among them.
		time_list: a list of lists with format [[year,month,day,hour,minute,second]],[t2],...[tn]]
				   for each observation in UTC.
		landing_site_coordinates: tuple of form (lon, lat) of lunar landing site. This is the 
								  Lunar frame coordinates of the observatory, in degrees.
		MOON_ALT: the altitude (in km) of the landing site. This is the elevation of the landing site
				  with respect to the average radius of the Moon from the poles, given by MOON_RAD.
				  For most cases this can safely be left at 0, as elevation changes are very small compared
				  to the radius of the Moon.
		filename_to_save: filename to save all data generated from this object into. If None
						  no data is saved. Default: None.
		antenna_psi_angle: angle that the antenna is rotated relative to North in the lunar
						   beam frame, given in degrees. Default: 0 degrees.
		galaxy_map: healpy map of the galaxy to be convolved with the antenna. LOCHNESS will rotate
					galaxy_map so that it is in the beam frame, a left-handed coordinate system,
					with North pointing at (0,0) and East at (-90,0) in Mollview projection.
		left_handed_mollview_maps: Boolean determining whether to use a left-handed coordinate system when
								   generating the altitude and azimuth of the celestial bodies, 
								   which is the same as the system that the galaxy_map is rotated into. 
								   Default: True.
								   WARNING: If you set this Boolean to False, then the galaxy_map
								   and celestial bodies will be in different coordinate systems
								   and therefore cannot be plotted together!
		run_test_galaxy_map: this function will use the N and E vectors generated from spicepy
							 at each ephem time to define the beam frame, along with the zenith
							 vector. This function first rotates into the beam frame by rotating
							 so that the zenith is overhead (at (0,90) in Mollview projection).
							 Next it spins the map from the pole so that the N vector is at (0,0),
							 defining North, and E is at (-90,0), defining East. This is a left-
							 handed coordinate system. Note that this function takes less
							 rotations than the usual one, but is difficult to include longer
							 integration times (and thus map smearing) in this framework.
		'''
		self.path_to_spice_kernels = path_to_spice_kernels
		self.time_list = time_list
		self.landing_site_coordinates = landing_site_coordinates
		self.MOON_ALT = MOON_ALT
		self.filename_to_save = filename_to_save
		self.antenna_psi_angle = antenna_psi_angle
		self.galaxy_map = galaxy_map
		self.left_handed_mollview_maps = left_handed_mollview_maps
		self.run_test_galaxy_map = run_test_galaxy_map
			
		print('Furnsh-ing spice kernels...')
		spice.furnsh(self.path_to_spice_kernels + 'naif0012.tls') #leap seconds
		spice.furnsh(self.path_to_spice_kernels + 'de440.bsp') #positions of planets
		spice.furnsh(self.path_to_spice_kernels + 'pck00011.tpc') #planetary constants

		# loading moon_pa kernels
		spice.furnsh(self.path_to_spice_kernels + 'moon_080317.tf')
		##### We want to use MOON_PA for better accuracies #####
		spice.furnsh(self.path_to_spice_kernels + 'moon_assoc_pa.tf')
		spice.furnsh(self.path_to_spice_kernels + 'moon_pa_de421_1900-2050.bpc')
		
		zenith_gal_coordinates_all_times = []
		lunar_north_pole_gal_coordinates_all_times = []
		ephemTime_list = []
		norm_north_vector_list = []
		norm_east_vector_list = []
		norm_north_vector_gal_coord_list = []
		norm_east_vector_gal_coord_list = []
		sun_alt_azi_all_times = []
		jupiter_alt_azi_all_times = []
		earth_alt_azi_all_times = []
		galaxy_map_all_times = []
		jupiter_earth_ang_distance_list = []
		sun_earth_ang_distance_list = []
		jupiter_sun_ang_distance_list = []
		test_galaxy_maps_from_zenith_vectors_list = []
		for time_snap in self.time_list:
			print('Computing lunar simulation for time ', time_snap)
		
			time = datetime(year=int(time_snap[0]), month=int(time_snap[1]), \
				day=int(time_snap[2]), hour=time_snap[3], minute=time_snap[4], \
				second=time_snap[5])
			ephemTime = spice.datetime2et(time)
			ephemTime_list.append(ephemTime)

			landing_site_zenith_lon, landing_site_zenith_lat, zenith_vector_gal_coord = \
				self.getLunarZenithGalLonLat(ephemTime, self.landing_site_coordinates[0],\
					self.landing_site_coordinates[1])
			zenith_gal_coordinates_all_times.append([landing_site_zenith_lon,\
				landing_site_zenith_lat])
				
			lunar_north_pole_lon, lunar_north_pole_lat, lunar_north_pole_vector_gal_coord = \
				self.getLunarZenithGalLonLat(ephemTime, 0.,90.)
			lunar_north_pole_gal_coordinates_all_times.append([lunar_north_pole_lon,\
				lunar_north_pole_lat])
				
			norm_north, norm_east = self.getNorthandEastVectorsFromZenith(ephemTime, zenith_vector_gal_coord)
			
			north_vec_dist, north_vec_galLon, north_vec_galLat = spice.reclat(norm_north)
			norm_north_vector_gal_coord_list.append([north_vec_galLon,north_vec_galLat])
			
			east_vec_dist, east_vec_galLon, east_vec_galLat = spice.reclat(norm_east)
			norm_east_vector_gal_coord_list.append([east_vec_galLon, east_vec_galLat])

			norm_north_vector_list.append(norm_north)
			norm_east_vector_list.append(norm_east)
				
			sun_alt, sun_azi, sun_vec_from_moon = self.getAltAzFromZenith(str(sun_naif), ephemTime,\
				zenith_vector_gal_coord,\
				norm_north, norm_east)
			jupiter_alt, jupiter_azi, jupiter_vec_from_moon = self.getAltAzFromZenith(str(jupiter_naif),\
				ephemTime, zenith_vector_gal_coord,\
				norm_north, norm_east)
			earth_alt, earth_azi, earth_vec_from_moon = self.getAltAzFromZenith(str(earth_naif), ephemTime,\
				zenith_vector_gal_coord,\
				norm_north, norm_east)
				
			jupiter_sun_ang_distance = self.getAngDistanceBetweenBodiesinRadians(jupiter_vec_from_moon,\
				sun_vec_from_moon)
			jupiter_sun_ang_distance_list.append(jupiter_sun_ang_distance)
			jupiter_earth_ang_distance = self.getAngDistanceBetweenBodiesinRadians(jupiter_vec_from_moon,\
				earth_vec_from_moon)
			jupiter_earth_ang_distance_list.append(jupiter_earth_ang_distance)
			sun_earth_ang_distance = self.getAngDistanceBetweenBodiesinRadians(sun_vec_from_moon,\
				earth_vec_from_moon)
			sun_earth_ang_distance_list.append(sun_earth_ang_distance)
				
			sun_alt_azi_all_times.append([sun_alt,sun_azi])
			jupiter_alt_azi_all_times.append([jupiter_alt,jupiter_azi])
			earth_alt_azi_all_times.append([earth_alt,earth_azi])
			
			if self.galaxy_map is not None:
				ndims = self.galaxy_map.ndim
				if ndims == 1:
					lunar_frame_galaxy_map = self.getBeamFrameGalaxyMap(self.galaxy_map, \
						lunar_north_pole_lon, lunar_north_pole_lat,\
						landing_site_zenith_lon, landing_site_zenith_lat)
					
					galaxy_map_all_times.append(lunar_frame_galaxy_map)
				else:
					nmaps = self.galaxy_map.shape[0]
					galaxy_by_map = []
					for imap in range(nmaps):
						lunar_frame_galaxy_map = self.getBeamFrameGalaxyMap(self.galaxy_map[imap], \
							lunar_north_pole_lon, lunar_north_pole_lat,\
							landing_site_zenith_lon, landing_site_zenith_lat)
						galaxy_by_map.append(lunar_frame_galaxy_map)
						
					galaxy_map_all_times.append(galaxy_by_map)
				if self.run_test_galaxy_map:
					lunar_frame_galaxy_map_from_vec = \
					self.getBeamFrameGalaxyMapFromNEZVectors(self.galaxy_map, \
					landing_site_zenith_lon, landing_site_zenith_lat, \
					north_vec_galLon, north_vec_galLat)
					
					test_galaxy_maps_from_zenith_vectors_list.append(lunar_frame_galaxy_map_from_vec)
				
		self.all_ephemTime = np.array(ephemTime_list)
		self.zenith_galactic_coordinates = np.array(zenith_gal_coordinates_all_times)
		self.lunar_north_pole_galactic_coordinates = np.array(lunar_north_pole_gal_coordinates_all_times)
		self.normalized_north_vector_from_landing_site = np.array(norm_north_vector_list)
		self.normalized_east_vector_from_landing_site = np.array(norm_east_vector_list)
		self.sun_altitude_azimuth_from_landing_site = np.array(sun_alt_azi_all_times)
		self.jupiter_altitude_azimuth_from_landing_site = np.array(jupiter_alt_azi_all_times)
		self.earth_altitude_azimuth_from_landing_site = np.array(earth_alt_azi_all_times)
		self.normalized_north_vector_galactic_coordinates = np.array(norm_north_vector_gal_coord_list)
		self.normalized_east_vector_galactic_coordinates = np.array(norm_east_vector_gal_coord_list)
		self.lunar_frame_galaxy_maps = np.array(galaxy_map_all_times)
		self.jupiter_sun_angular_distance_radians = np.array(jupiter_sun_ang_distance_list)
		self.jupiter_earth_angular_distance_radians = np.array(jupiter_earth_ang_distance_list)
		self.sun_earth_angular_distance_radians = np.array(sun_earth_ang_distance_list)
		
		if self.run_test_galaxy_map:
			self.test_galaxy_maps_from_zenith_vectors = np.array(test_galaxy_maps_from_zenith_vectors_list)

		if self.filename_to_save is not None:
			try:
				import h5py
			except ImportError:
				print('You must have h5py installed to use the save simulations method.')
			print('Saving lunar simulation...')
			with h5py.File(self.filename_to_save, 'w') as save_file:
				save_file.create_dataset(name='all_ephemTime',data=self.all_ephemTime)
				save_file.create_dataset(name='zenith_galactic_coordinates',\
					data=self.zenith_galactic_coordinates)
				save_file.create_dataset(name='lunar_north_pole_galactic_coordinates',\
					data=self.lunar_north_pole_galactic_coordinates)
				save_file.create_dataset(name='normalized_north_vector_from_landing_site',\
					data=self.normalized_north_vector_from_landing_site)
				save_file.create_dataset(name='normalized_east_vector_from_landing_site',\
					data=self.normalized_east_vector_from_landing_site)
				save_file.create_dataset(name='sun_altitude_azimuth_from_landing_site',\
					data=self.sun_altitude_azimuth_from_landing_site)
				save_file.create_dataset(name='jupiter_altitude_azimuth_from_landing_site',\
					data=self.jupiter_altitude_azimuth_from_landing_site)
				save_file.create_dataset(name='earth_altitude_azimuth_from_landing_site',\
					data=self.earth_altitude_azimuth_from_landing_site)
				save_file.create_dataset(name='normalized_north_vector_galactic_coordinates',\
					data=self.normalized_north_vector_galactic_coordinates)
				save_file.create_dataset(name='normalized_east_vector_galactic_coordinates',\
					data=self.normalized_east_vector_galactic_coordinates)
				save_file.create_dataset(name='jupiter_sun_angular_distance_radians',\
					data=self.jupiter_sun_angular_distance_radians)
				save_file.create_dataset(name='jupiter_earth_angular_distance_radians',\
					data=self.jupiter_earth_angular_distance_radians)
				save_file.create_dataset(name='sun_earth_angular_distance_radians',\
					data=self.sun_earth_angular_distance_radians)
				if self.galaxy_map is not None:
					save_file.create_dataset(name='lunar_frame_galaxy_maps',\
						data=self.lunar_frame_galaxy_maps)
				if self.run_test_galaxy_map:
					save_file.create_dataset(name='lunar_frame_galaxy_maps_from_NEZ_vectors',\
						data=self.test_galaxy_maps_from_zenith_vectors)
					
	def getAngDistanceBetweenBodiesinRadians(self, body1_vec, body2_vec):
		normed_dot_product = (np.dot(body1_vec, body2_vec))/ \
			(np.linalg.norm(body1_vec) * np.linalg.norm(body2_vec))
			
		return np.arccos(normed_dot_product)
				
	def getLunarZenithGalLonLat(self, ephemTime, lunar_suface_lon, lunar_surface_lat):
		"""
		Get the Galactic Coordinates of the lunar zenith pointing vector
		given a certain time and location on the lunar surface.
		
		Parameters
		----------
		ephemTime: the ephemTime corresponding to the observation of the desired body.
		lunar_surface_lon: longitude of the landing site on the lunar surface in degrees.
		lunar_surface_lat: latitude of the landing site on the lunar surface in degrees.
		
		Returns
		-------
		Zenith_Lon : float
			lunar zenith pointing galactic longitude in radians
		Zenith_Lat : float
			lunar zenith pointing galactic latitude in radians
		moonZenithGalRec: the vector in rectangular coordinates of the zenith from the
						  landing site on the Moon at the given ephemTime.
		"""
		# Define rotation matrix to rotate frame to the specified time (we want 'GALACTIC' frame)
		rotationMatrix = spice.pxform("MOON_PA", "GALACTIC", ephemTime)
		# Transform moon lon, lat, alt+rad (in rad) --> rectangular coordinates (vector)
		moonZenithRec = spice.latrec(MOON_RAD+self.MOON_ALT, \
			lunar_suface_lon * np.pi/180., lunar_surface_lat * np.pi/180.)
		
		# Get the zenith vector (or the vector pointing from the center of the moon to the landing site)
		moonZenithGalRec = spice.mxv(rotationMatrix, moonZenithRec)
		
		# GalLonLat of the zenith vector
		Zenith_Dist, Zenith_Lon, Zenith_Lat = spice.reclat(moonZenithGalRec)

		return Zenith_Lon, Zenith_Lat, moonZenithGalRec
		
	def getNorthandEastVectorsFromZenith(self, ephemTime, moonZenithGalacticRec):
		'''
		ephemTime: the ephemTime corresponding to the observation of the desired body.
		moonZenithGalacticRec: the vector in rectangular coordinates of the zenith from the
							   landing site on the Moon at the given ephemTime.

		returns: normalized vectors pointing Northward and Eastward according
				 to the given zenith vector moonZenithGalacticRec and ephemTime,
				 on the lunar surface.
		'''
		rotationMatrix = spice.pxform("MOON_PA", "GALACTIC", ephemTime)
		MOON_LON = self.landing_site_coordinates[0]
		MOON_LAT = self.landing_site_coordinates[1]
		
		# Normalized vector -> zenith vector
		norm_moonZenithGalacticRec = moonZenithGalacticRec / np.linalg.norm(moonZenithGalacticRec)
		
		# arbitrary vector 1 = moonpoint w lat : lat+del North
		delta_vec_1 = spice.latrec(MOON_RAD + self.MOON_ALT,\
			MOON_LON * np.pi / 180., (MOON_LAT+del_angle) * np.pi / 180.)
		delta_vec_1_gal = spice.mxv(rotationMatrix, delta_vec_1)
		north = delta_vec_1_gal - moonZenithGalacticRec
		norm_north = north / np.linalg.norm(north)
		#print(np.linalg.norm(norm_north))

		# arbitrary vector 2 = moonpoint w lon+del East
		delta_vec_2 = spice.latrec(MOON_RAD + self.MOON_ALT,\
			(MOON_LON+del_angle) * np.pi / 180., MOON_LAT * np.pi / 180.)
		delta_vec_2_gal = spice.mxv(rotationMatrix, delta_vec_2)
		east = delta_vec_2_gal - moonZenithGalacticRec
		norm_east = east / np.linalg.norm(east)
		
		return norm_north, norm_east
		
	def getAltAzFromZenith(self, body, ephemTime, moonZenithGalacticRec,\
		norm_north, norm_east):
		'''
		body: a string of the name of the body for which to calculate altitude and azimuth 
			  from the landing site. Supported are the 'NAIF' types in Spice.
		ephemTime: the ephemTime corresponding to the observation of the desired body.
		moonZenithGalacticRec: the vector in rectangular coordinates of the zenith from the
							   landing site on the Moon at the given ephemTime.
		norm_north: the unit vector in rectangular coordinates of the northward pointing
					direction from the landing site on the Moon at the given ephemTime.
		norm_east: the unite vector in rectangular coordinates of the eastward pointing
				   direction from the landing site on the Moon at the given ephemTime.
		
		returns: altitude and azimuth (in degrees) of the given body at
				 the given ephemTime from the lunar surface determined by
				 self.landing_site_coordinates,the zenith pointing 
				 moonZenithGalacticRec, and the normalized northward and 
				 eastward pointing vectors on the Lunar surface, determining
				 the North and East directions.
		'''
		# Normalized vector -> zenith vector
		norm_moonZenithGalacticRec = moonZenithGalacticRec / np.linalg.norm(moonZenithGalacticRec)

		# Get the body vector from the center of the moon
		bodyFromMoonCenterGalacticRec = spice.spkpos(body, ephemTime, \
			'GALACTIC', 'NONE', 'MOON')[0]
		# Get the body vector from the landing site
		bodyFromMoonPointGalacticRec = bodyFromMoonCenterGalacticRec - \
			(MOON_RAD+self.MOON_ALT)*norm_moonZenithGalacticRec

		# vec_A dot vec_B / |A||B|
		normedDotProduct = np.dot(moonZenithGalacticRec, bodyFromMoonPointGalacticRec) /\
			(np.linalg.norm(moonZenithGalacticRec) * np.linalg.norm(bodyFromMoonPointGalacticRec))

		# Altitude angle of the object
		bodyAlt = np.pi/2 - np.arccos(normedDotProduct)   # in rad

		# vector projected on N-E plane
		body_proj_z = np.dot(bodyFromMoonPointGalacticRec, norm_moonZenithGalacticRec) /\
			np.dot(norm_moonZenithGalacticRec, norm_moonZenithGalacticRec) * norm_moonZenithGalacticRec
		body_proj_NE = bodyFromMoonPointGalacticRec - body_proj_z

		body_N_norm_dot = np.dot(body_proj_NE, norm_north) /\
			(np.linalg.norm(body_proj_NE) * np.linalg.norm(norm_north))
		body_N = np.arccos(body_N_norm_dot)   # in rad

		body_E_norm_dot = np.dot(body_proj_NE, norm_east) /\
			(np.linalg.norm(body_proj_NE) * np.linalg.norm(norm_east))
		body_E = np.arccos(body_E_norm_dot)   # in rad

		if body_E <= np.pi/2:
			bodyAz = body_N
		else:
			bodyAz = 2*np.pi-body_N  # 0 to 360 degree from N (->E->S->W)
		
		if self.left_handed_mollview_maps:
			bodyAz = -1*bodyAz

		return bodyAlt*180./np.pi, bodyAz*180./np.pi, bodyFromMoonPointGalacticRec
		
	def getBeamFrameGalaxyMapFromNEZVectors(self, galaxy_map, zenith_at_ephemTime_lon, \
		zenith_at_ephemTime_lat, north_at_ephemTime_lon, north_at_ephemTime_lat):
		'''
		galaxy_map: healpy map of the galaxy to be convolved with the antenna.
		zenith_at_ephemTime_lon: the Galactic longitude of the zenith vector from the landing site
								 at ephemTime in radians.
		zenith_at_ephemTime_lat: the Galactic latitude of the zenith vector from the landing site
								 at ephemTime in radians.
		north_at_ephemTime_lon: the Galactic longitude of the north-pointing vector from the landing
								site at ephemTime in radians.
		north_at_ephemTime_lat: the Galactic latitude of the north-pointing vector from the landing
								site at ephemTime in radians.

		returns: a 1D healpy map containing the input galaxy_map at the ephemTime from the 
				 beam frame defined by the zenith, lunar North Pole, and landing site
				 coordinates.
		'''
		#Change all coordinates to degrees.
		zenith_at_ephemTime_lon = zenith_at_ephemTime_lon*(180/np.pi)
		zenith_at_ephemTime_lat = zenith_at_ephemTime_lat*(180/np.pi)
		
		north_at_ephemTime_lon = north_at_ephemTime_lon*(180/np.pi)
		north_at_ephemTime_lat = north_at_ephemTime_lat*(180/np.pi)
		
		east_lon = east_lon*(180/np.pi)
		east_lat = east_lat*(180/np.pi)
		
		first_rotation_theta = 90 - zenith_at_ephemTime_lat
		first_rotation_phi = zenith_at_ephemTime_lon
		first_rotation_psi = 0.
		#Create a rotator which takes the galaxy_map and other vectors into the 
		#beam frame, where the zenith is overhead at (0,90).
		first_lunar_zenith_rotator = spherical_rotator(first_rotation_theta,\
			first_rotation_phi, first_rotation_psi).get_inverse()
			
		rotated_north_vec_lon, rotated_north_vec_lat = \
			first_lunar_zenith_rotator(north_at_ephemTime_lon, \
			north_at_ephemTime_lat, lonlat=True)
		#Create a rotator which finds the longitude of the northward pointing
		#vector and rotates so that it is at (0,0) corresponding to North
		#in the beam frame.
		north_rotator = rotator_for_spinning(-rotated_north_vec_lon)
		
		full_vec_rotator = north_rotator * first_lunar_zenith_rotator
		#Rotate the galaxy map.
		rotated_galaxy_map_from_vec = full_vec_rotator.rotate_map_alms(galaxy_map)
		
		return rotated_galaxy_map_from_vec
		
	def getBeamFrameGalaxyMap(self, galaxy_map, north_pole_lon, north_pole_lat,\
		zenith_at_ephemTime_lon, zenith_at_ephemTime_lat):
		'''
		galaxy_map: healpy map of the galaxy to be convolved with the antenna.
		north_pole_lon: the Galactic longitude of the Lunar North Pole at ephemTime in radians.
		north_pole_lat: the Galactic latitude of the Lunar North Pole at ephemTime in radians.
		zenith_at_ephemTime_lon: the Galactic longitude of the zenith vector from the landing site
								 at ephemTime in radians.
		zenith_at_ephemTime_lat: the Galactic latitude of the zenith vector from the landing site
								 at ephemTime in radians.
								 
		returns: a 1D healpy map containing the input galaxy_map at the ephemTime from the 
				 beam frame defined by the zenith, lunar North Pole, and landing site
				 coordinates.
		'''
		north_pole_lon = north_pole_lon*(180/np.pi)
		north_pole_lat = north_pole_lat*(180/np.pi)
		zenith_at_ephemTime_lon = zenith_at_ephemTime_lon*(180/np.pi)
		zenith_at_ephemTime_lat = zenith_at_ephemTime_lat*(180/np.pi)
		
		first_rotation_theta = 90 - north_pole_lat
		first_rotation_phi = north_pole_lon
		first_rotation_psi = 0.

		north_pole_lunar_rotator = spherical_rotator(first_rotation_theta,\
			first_rotation_phi, first_rotation_psi).get_inverse()

		(rotated_lunar_NP_lon, rotated_lunar_NP_lat) =\
			north_pole_lunar_rotator(north_pole_lon, north_pole_lat,\
			lonlat=True)
		(rotated_lunar_zen_lon, rotated_lunar_zen_lat) =\
			north_pole_lunar_rotator(zenith_at_ephemTime_lon,zenith_at_ephemTime_lat,\
			lonlat=True)

		lunar_amount_to_spin = self.landing_site_coordinates[0] -\
			(rotated_lunar_zen_lon)

		spinning_lunar_rotator = rotator_for_spinning(lunar_amount_to_spin)

		first_full_lunar_rotator = spinning_lunar_rotator * north_pole_lunar_rotator

		rot_lunar_sky_map = \
			first_full_lunar_rotator.rotate_map_alms(galaxy_map)

		#Smearing is almost certainly unnecessary for the Moon, given that the
		#Moon spins so slowly and most integration times will be less than
		#a minute long or so. So let's omit smearing for now. However, it can be
		#easily incorporated using small modifications of existing functions within
		#perses.

		xpart = np.cos(np.radians(90 - self.landing_site_coordinates[1])) * \
			np.cos(np.radians(self.landing_site_coordinates[0]))
		ypart = np.cos(np.radians(90 - self.landing_site_coordinates[1])) *\
			np.sin(np.radians(self.landing_site_coordinates[0]))
		zpart = -np.sin(np.radians(90 - self.landing_site_coordinates[1]))
		lunar_northhat = (-1. * np.array([xpart, ypart, zpart]))

		lunar_beam_frame_rotator = spherical_rotator(np.radians(90 - self.landing_site_coordinates[1]),\
			np.radians(self.landing_site_coordinates[0]), 0, deg=False).get_inverse()
		rot_lunar_northhat = lunar_beam_frame_rotator(lunar_northhat)

		#The displacement angle 
		lunar_displacement =\
			np.degrees(np.arctan2(rot_lunar_northhat[1], rot_lunar_northhat[0]))

		lunar_antenna_psi_rotator =\
			rotator_for_spinning(self.antenna_psi_angle - lunar_displacement)
		
		lunar_final_full_rotator = lunar_antenna_psi_rotator *\
			lunar_beam_frame_rotator
		
		final_rot_lunar_sky_maps = \
			lunar_final_full_rotator.rotate_map_alms(rot_lunar_sky_map)

		return final_rot_lunar_sky_maps
		
	@property
	def all_ephemTime(self):
		if not hasattr(self, '_all_ephemTime'):
			raise AttributeError('all_ephemTime referenced before set!')
		return self._all_ephemTime
		
	@all_ephemTime.setter
	def all_ephemTime(self, value):
		self._all_ephemTime = value
		
	@property
	def zenith_galactic_coordinates(self):
		if not hasattr(self, '_zenith_galactic_coordinates'):
			raise AttributeError('zenith_galactic_coordinates referenced before set!')
		return self._zenith_galactic_coordinates
		
	@zenith_galactic_coordinates.setter
	def zenith_galactic_coordinates(self, value):
		self._zenith_galactic_coordinates = value
		
	@property
	def lunar_north_pole_galactic_coordinates(self):
		if not hasattr(self, '_lunar_north_pole_galactic_coordinates'):
			raise AttributeError('lunar_north_pole_galactic_coordinates referenced before set!')
		return self._lunar_north_pole_galactic_coordinates
		
	@lunar_north_pole_galactic_coordinates.setter
	def lunar_north_pole_galactic_coordinates(self, value):
		self._lunar_north_pole_galactic_coordinates = value
		
	@property
	def normalized_north_vector_from_landing_site(self):
		if not hasattr(self, '_normalized_north_vector_from_landing_site'):
			raise AttributeError('normalized_north_vector_from_landing_site referenced before set!')
		return self._normalized_north_vector_from_landing_site
		
	@normalized_north_vector_from_landing_site.setter
	def normalized_north_vector_from_landing_site(self, value):
		self._normalized_north_vector_from_landing_site = value
		
	@property
	def normalized_east_vector_from_landing_site(self):
		if not hasattr(self, '_normalized_east_vector_from_landing_site'):
			raise AttributeError('normalized_east_vector_from_landing_site referenced before set!')
		return self._normalized_east_vector_from_landing_site
		
	@normalized_east_vector_from_landing_site.setter
	def normalized_east_vector_from_landing_site(self, value):
		self._normalized_east_vector_from_landing_site = value
		
	@property
	def sun_altitude_azimuth_from_landing_site(self):
		if not hasattr(self, '_sun_altitude_azimuth_from_landing_site'):
			raise AttributeError('sun_altitude_azimuth_from_landing_site referenced before set!')
		return self._sun_altitude_azimuth_from_landing_site
		
	@sun_altitude_azimuth_from_landing_site.setter
	def sun_altitude_azimuth_from_landing_site(self, value):
		self._sun_altitude_azimuth_from_landing_site = value
		
	@property
	def jupiter_altitude_azimuth_from_landing_site(self):
		if not hasattr(self, '_jupiter_altitude_azimuth_from_landing_site'):
			raise AttributeError('jupiter_altitude_azimuth_from_landing_site referenced before set!')
		return self._jupiter_altitude_azimuth_from_landing_site
		
	@jupiter_altitude_azimuth_from_landing_site.setter
	def jupiter_altitude_azimuth_from_landing_site(self, value):
		self._jupiter_altitude_azimuth_from_landing_site = value
		
	@property
	def earth_altitude_azimuth_from_landing_site(self):
		if not hasattr(self, '_earth_altitude_azimuth_from_landing_site'):
			raise AttributeError('earth_altitude_azimuth_from_landing_site referenced before set!')
		return self._earth_altitude_azimuth_from_landing_site
		
	@earth_altitude_azimuth_from_landing_site.setter
	def earth_altitude_azimuth_from_landing_site(self, value):
		self._earth_altitude_azimuth_from_landing_site = value
		
	@property
	def normalized_north_vector_galactic_coordinates(self):
		if not hasattr(self, '_normalized_north_vector_galactic_coordinates'):
			raise AttributeError('normalized_north_vector_galactic_coordinates referenced before set!')
		return self._normalized_north_vector_galactic_coordinates
		
	@normalized_north_vector_galactic_coordinates.setter
	def normalized_north_vector_galactic_coordinates(self, value):
		self._normalized_north_vector_galactic_coordinates = value
		
	@property
	def normalized_east_vector_galactic_coordinates(self):
		if not hasattr(self, '_normalized_east_vector_galactic_coordinates'):
			raise AttributeError('normalized_east_vector_galactic_coordinates referenced before set!')
		return self._normalized_east_vector_galactic_coordinates
		
	@normalized_east_vector_galactic_coordinates.setter
	def normalized_east_vector_galactic_coordinates(self, value):
		self._normalized_east_vector_galactic_coordinates = value
		
	@property
	def jupiter_sun_angular_distance_radians(self):
		if not hasattr(self, '_jupiter_sun_angular_distance_radians'):
			raise AttributeError('jupiter_sun_angular_distance_radians referenced before set!')
		return self._jupiter_sun_angular_distance_radians
		
	@jupiter_sun_angular_distance_radians.setter
	def jupiter_sun_angular_distance_radians(self, value):
		self._jupiter_sun_angular_distance_radians = value
		
	@property
	def jupiter_earth_angular_distance_radians(self):
		if not hasattr(self, '_jupiter_earth_angular_distance_radians'):
			raise AttributeError('jupiter_earth_angular_distance_radians referenced before set!')
		return self._jupiter_earth_angular_distance_radians
		
	@jupiter_earth_angular_distance_radians.setter
	def jupiter_earth_angular_distance_radians(self, value):
		self._jupiter_earth_angular_distance_radians = value
		
	@property
	def sun_earth_angular_distance_radians(self):
		if not hasattr(self, '_sun_earth_angular_distance_radians'):
			raise AttributeError('sun_earth_angular_distance_radians referenced before set!')
		return self._sun_earth_angular_distance_radians
		
	@sun_earth_angular_distance_radians.setter
	def sun_earth_angular_distance_radians(self, value):
		self._sun_earth_angular_distance_radians = value
		
	@property
	def lunar_frame_galaxy_maps(self):
		if not hasattr(self, '_lunar_frame_galaxy_maps'):
			raise AttributeError('lunar_frame_galaxy_maps referenced before set!')
		return self._lunar_frame_galaxy_maps
		
	@lunar_frame_galaxy_maps.setter
	def lunar_frame_galaxy_maps(self, value):
		self._lunar_frame_galaxy_maps = value
		
	@property
	def test_galaxy_maps_from_zenith_vectors(self):
		if not hasattr(self, '_test_galaxy_maps_from_zenith_vectors'):
			raise AttributeError('test_galaxy_maps_from_zenith_vectors referenced before set!')
		return self._test_galaxy_maps_from_zenith_vectors
		
	@test_galaxy_maps_from_zenith_vectors.setter
	def test_galaxy_maps_from_zenith_vectors(self, value):
		self._test_galaxy_maps_from_zenith_vectors = value

