
# v2: we don't want to use the coastline - we want to use the 25m isobath as the inshore wall

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
import scipy.interpolate as spint


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

isobath_inshore_depth = 25

isobath_offshore_depth = 300
#isobath_offshore_depth = 250


grid_path_in = '/home/blaughli/tracking_project/grid_data/wc15n_grd.nc'

base_dir = '/home/blaughli/jordan_project/'
output_dir= base_dir + 'y_setup_output/'

coastline_file_in = output_dir + 'coastline_coords_psi_wc15n.npz'

dset = netCDF4.Dataset(grid_path_in, 'r')

lon_field = np.array(dset['lon_rho'])
lat_field = np.array(dset['lat_rho'])
h = np.array(dset['h'])

dset.close


output_file = output_dir + 'polygon_seed_isobath_{0:04}_lonlat_wc15n'.format(isobath_offshore_depth)


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------



RGI = spint.RegularGridInterpolator
# create interpolator to get lat/lon at isoline points
x = np.arange(np.shape(lon_field)[0])
y = np.arange(np.shape(lon_field)[1])
rgi_lon = RGI([x,y],lon_field)
rgi_lat = RGI([x,y],lat_field)


contours_offshore = measure.find_contours(h, isobath_offshore_depth)

# Assume that the longest contour returned is the one we want....
max_length = 0 
contour_index = 0 
dex = 0 
for contour in contours_offshore:
    if np.shape(contour)[0] > max_length:
        max_length = np.shape(contour)[0]
        contour_index = dex 
    dex += 1

isoline_ij = contours_offshore[contour_index]

isoline_offshore_lonlat = np.zeros(np.shape(isoline_ij))
for ii in range(np.shape(isoline_ij)[0]):
    isoline_offshore_lonlat[ii,0] = rgi_lon((isoline_ij[ii,0], isoline_ij[ii,1]))
    isoline_offshore_lonlat[ii,1] = rgi_lat((isoline_ij[ii,0], isoline_ij[ii,1]))


contours_inshore = measure.find_contours(h, isobath_inshore_depth)

# Assume that the longest contour returned is the one we want....
max_length = 0 
contour_index = 0 
dex = 0 
for contour in contours_inshore:
    if np.shape(contour)[0] > max_length:
        max_length = np.shape(contour)[0]
        contour_index = dex 
    dex += 1

isoline_ij = contours_inshore[contour_index]

isoline_inshore_lonlat = np.zeros(np.shape(isoline_ij))
for ii in range(np.shape(isoline_ij)[0]):
    isoline_inshore_lonlat[ii,0] = rgi_lon((isoline_ij[ii,0], isoline_ij[ii,1]))
    isoline_inshore_lonlat[ii,1] = rgi_lat((isoline_ij[ii,0], isoline_ij[ii,1]))




# make the polygon

poly_lon = list(isoline_inshore_lonlat[:,0])
poly_lat = list(isoline_inshore_lonlat[:,1])

poly_lon = poly_lon + list(reversed(isoline_offshore_lonlat[:,0]))
poly_lat = poly_lat + list(reversed(isoline_offshore_lonlat[:,1]))

poly_lon.append(isoline_inshore_lonlat[0,0])
poly_lat.append(isoline_inshore_lonlat[0,1])

plt.plot(poly_lon,poly_lat)
plt.show()

#poly_lon.append(d

#d['isoline_offshore_lonlat'] = isoline_offshore_lonlat

#np.savez(output_file, **d)








