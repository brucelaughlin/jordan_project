# v2: stop using pickle

# walking the psi_mask_bl grid
# point, assuming the prevoius point was blow, and checking in a
# clockwise manner

# Use "psi_bl", as defined by Chris!

import pickle
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import ast

#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

base_path = '/home/blaughli/jordan_project/y_setup_output/'

grid_path_in = '/home/blaughli/tracking_project/grid_data/wc15n_grd.nc'
psi_bl_path_in = base_path + 'mask_psi_bl_continent.p'

coastline_file_out = base_path + 'coastline_coords_psi_wc15n'

dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_psi = 'psi'
lon_psi = np.array(dset['lon_{}'.format(points_type_psi)])
lat_psi = np.array(dset['lat_{}'.format(points_type_psi)])

mask_u = np.array(dset['mask_u'])
mask_v = np.array(dset['mask_v'])

dset.close

#---------------------------------------------------------------------
#---------------------------------------------------------------------


# Need to modify this to work for islands, ie no longer have the condition
# that the coast ends at the edge of the grid

# load psi_bl mask
file = open(psi_bl_path_in,'rb')
psi = pickle.load(file)
file.close

# size of domains i and j
n_i = np.shape(psi)[0]
n_j = np.shape(psi)[1]

# begin algorithm...
# walk around clockwise, coast on the right, looking left first then clockwise

coast_lon = []
coast_lat = []

# If I got i/j of the grid dims backwards, switch the order of the starting coordinates
# Start in lower right of grid (for wc12, this is land near scb)
# This assumes we start on land.  need to adjust if that's not the case

ii = 0 
jj = n_j - 1

current_psi = psi[ii,jj]

while current_psi == 0:
    jj -= 1;
    current_psi = psi[ii,jj]
 
jj = jj + 1    
 
# current point = "cp"
cp = [ii,jj]



# Make fake last_point, to the south of starting point 
# (may be better choice, depending on where starting point is - want the
# fake point to enforce the desired starting direction of walk/traversal)

# Last point = "lp"
lp = [cp[0]-1,cp[1]] # for islands, starting in a NW corner
#lp = [cp[0],cp[1]+1] # for west coast continent, starting in most SW point of land mask

coast_lon.append(lon_psi[ii,jj])
coast_lat.append(lat_psi[ii,jj])

while True:
    
    # cp left of lp
    if cp[1] < lp[1]:
        try:
            # look down
            if mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            # look left
            elif mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            # look up
            elif mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            # look right
            elif mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            else:
                break
        except:
            break


    # cp above lp
    elif cp[0] > lp[0]:
        try:
            # look left
            if mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            # look up
            elif mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            # look right
            elif mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            # look down
            elif mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            else:
                break
        except:
            break

    # cp right of lp
    elif cp[1] > lp[1]:
        try:
            # look up
            if mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            # look right
            elif mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            # look down
            elif mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            # look left
            elif mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            else:
                break
        except:
            break

    # cp below lp
    elif cp[0] < lp[0]:
        try:
            # look right
            if mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            # look down
            elif mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            # look left
            elif mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            # look up
            elif mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            else:
                break
        except:
            break

    if ii >= n_i or jj >= n_j:
        break
    else:
        coast_lon.append(lon_psi[ii,jj])
        coast_lat.append(lat_psi[ii,jj])
        lp = cp
        cp = [ii,jj]
        #print(cp)
        #print([n_i,n_j])
        #print('\n')


final_coordinates = np.zeros((len(coast_lon),2))
final_coordinates[:,0] = coast_lon
final_coordinates[:,1] = coast_lat

d = {}
d['final_coordinates'] = final_coordinates

np.savez(coastline_file_out, **d)







