#!/usr/bin/env python

from __future__ import print_function

from math import floor
import numpy as np
import tracebacK

def F(r, phi0, phi1):
  # find which bins it belongs to
  r_bin = floor( (r - rs[0]) / delta_r )
  if (r_bin < 0):
    print("r index >>>", r_bin, "<<< out of bounds")
    traceback.print_stack()
    quit()
  phi0_bin = floor( (phi0 - phi0s[0]) / delta_phi )
  phi1_bin = floor( (phi1 - phi1s[0]) / delta_phi )

  # now interpolate
  
  # On a periodic and cubic lattice, let r_d, phi0_d, and phi1_d be the differences between each of r, phi0, phi1 and the smaller coordinate related, that is:

  r_d = (r - rs[r_bin])/delta_r
  phi0_d = (phi0 - phi0s[phi0_bin])/delta_phi
  phi1_d = (phi1 - phi1s[phi1_bin])/delta_phi

  #First we interpolate along r (imagine we are pushing the front face of the cube to the back), giving:

  # c_00 = V[x_0,y_0, z_0] * (1 - x_d) + V[x_1, y_0, z_0] * x_d 
  c_00 = pmft[r_bin][phi0_bin][phi1_bin] * (1 - r_d) + pmft[r_bin + 1][phi0_bin][phi1_bin] * r_d 

  # c_10 = V[x_0,y_1, z_0] * (1 - x_d) + V[x_1, y_1, z_0] * x_d 
  c_10 = pmft[r_bin][phi0_bin + 1][phi1_bin] * (1 - r_d) + pmft[r_bin + 1][phi0_bin + 1][phi1_bin] * r_d 
  # c_01 = V[x_0,y_0, z_1] * (1 - x_d) + V[x_1, y_0, z_1] * x_d 
  c_01 = pmft[r_bin][phi0_bin][phi1_bin + 1] * (1 - r_d) + pmft[r_bin + 1][phi0_bin][phi1_bin + 1] * r_d 
  #c_11 = V[x_0,y_1, z_1] * (1 - x_d) + V[x_1, y_1, z_1] * x_d 
  c_11 = pmft[r_bin][phi0_bin + 1][phi1_bin + 1] * (1 - r_d) + pmft[r_bin + 1][phi0_bin + 1][phi1_bin + 1] * r_d 

  # Then along phi0

  #c_0 = c_{00}(1 - y_d) + c_{10}y_d
  c_0 = c_00 * (1 - phi0_d) + c_10 * phi0_d
  #c_1 = c_{01}(1 - y_d) + c_{11}y_d
  c_1 = c_01 * (1 - phi0_d) + c_11 * phi0_d

  # And finally we interpolate along phi1

  # c = c_0(1 - z_d) + c_1z_d 
  c = c_0 * (1 - phi1_d) + c_1 * phi1_d 

  return c

def Fr(r, phi0, phi1):
  delta = 0.001
  return ((F(r+delta, phi0, phi1) - F(r-delta, phi0, phi1))/(2*delta))

def Fphi0(r, phi0, phi1):
  delta = 0.001
  return ((F(r, phi0 + delta, phi1) - F(r, phi0 - delta, phi1))/(2*delta))

def Fphi1(r, phi0, phi1):
  delta = 0.001
  return ((F(r, phi0, phi1 + delta) - F(r, phi0, phi1 + delta))/(2*delta))

##################################################################################################################
##################################################################################################################
##################################################################################################################

dat = np.load('avg_0.0_0.6.npz')

min_inds = dat['min_inds']
adj_list = dat['adj_list']
sdl_inds = dat['sdl_inds']
pmft     = dat['pmft']
rs       = dat['rs']
phi0s    = dat['phi0s']
phi1s    = dat['phi1s']

delta_r = rs[1] - rs[0]
delta_phi = phi0s[1] - phi0s[0]


#for r in range (1200, 2400):
#  print(r/1000.0, "\t",F(r/1000.0, 2, 3), "\t", Fr(r/1000.0, 2, 3))
  
#print(rs[0], "\t", rs[len(rs) -1])


##################################################################################################################
# badass variational vector calc stuff starts here..
## grab the loop stuff from other code

'''
for i in range(sdl_inds.__len__()):
  # draw the basin
  print(" sphere{<", phi1s[min_inds[i][0]], ",", phi0s[min_inds[i][1]], ",", rs[min_inds[i][2]], ">, 0.100000 texture{ pigment {color rgb <0, 0, 1> transmit 0.000000} finish {phong 1.000000}}}")

  for j in range(adj_list[i].__len__()):
#  print("saddles for 0 are at ", phi1s[sdl_inds[0][j][0]], phi0s[sdl_inds[0][j][1]], rs[sdl_inds[0][j][2]]);
#    print("cylinder { <", phi1s[min_inds[0][0]], ",", phi0s[min_inds[0][1]], ",", rs[min_inds[0][2]], ">, <", phi1s[sdl_inds[0][j][0]], ", ", phi0s[sdl_inds[0][j][1]], ", ", rs[sdl_inds[0][j][2]], " >, 0.5 open ", texture, " }")
    print("cylinder { <", phi1s[min_inds[i][0]], ",", phi0s[min_inds[i][1]], ",", rs[min_inds[i][2]], ">, <", phi1s[sdl_inds[i][j][0]], ", ", phi0s[sdl_inds[i][j][1]], ", ", rs[sdl_inds[i][j][2]], " >, 0.05 open ", texture, " }")
    #draw the saddle
    print(" sphere{<", phi1s[sdl_inds[i][j][0]], ", ", phi0s[sdl_inds[i][j][1]], ", ", rs[sdl_inds[i][j][2]], ">, 0.100000 texture{ pigment {color rgb <1, 0, 0> transmit 0.000000} finish {phong 1.000000}}}")

'''

phi_dimension = phi0s[0] - phi0s[len(phi0s) - 1] 
phi_dimension_sq = phi_dimension**2

# for each basin
#for basin in range(sdl_inds.__len__()):
for basin in range(1):

  # to each saddle point
  #for saddle in range(sdl_inds[basin].__len__()):
  for saddle in range(1):

    # create a string of points

    # end points
    basin_r = rs[min_inds[basin][2]]
    saddle_r = rs[sdl_inds[basin][saddle][2]]
    basin_phi0 = phi0s[min_inds[basin][1]]
    saddle_phi0 = phi0s[sdl_inds[basin][saddle][1]]
    basin_phi1 = phi1s[min_inds[basin][0]]
    saddle_phi1 = phi1s[sdl_inds[basin][saddle][0]]

    # apply minimum image convention
    if ((saddle_phi0 - basin_phi0) > +phi_dimension): saddle_phi0 -= phi_dimension
    if ((saddle_phi0 - basin_phi0) < -phi_dimension): saddle_phi0 += phi_dimension
    if ((saddle_phi1 - basin_phi1) > +phi_dimension): saddle_phi1 -= phi_dimension
    if ((saddle_phi1 - basin_phi1) < -phi_dimension): saddle_phi1 += phi_dimension

    # initialize variational points
    points_r=[]
    points_phi0=[]
    points_phi1=[]
    for point in range(256):
      points_r.append(basin_r + (saddle_r - basin_r) * (point/256.0))
      points_phi0.append(basin_phi0 + (saddle_phi0 - basin_phi0) * (point/256.0))
      points_phi1.append(basin_phi1 + (saddle_phi1 - basin_phi1) * (point/256.0))


    # need to wrap a loop around the following to iterate
    for repeat in range(1000):

      # and iterate by nudging these points perpendicular to the vector spanning the adjacent points
      for point in range(1,255):
	# Get the displacement
	
	# Gradient
	g_r = Fr(points_r[point], points_phi0[point], points_phi1[point])
	g_phi0 = Fphi0(points_r[point], points_phi0[point], points_phi1[point])
	g_phi1 = Fphi1(points_r[point], points_phi0[point], points_phi1[point])

	# line between adjacent points
	delta_r = points_r[point + 1] - points_r[point - 1]
	delta_phi0 = points_phi0[point + 1] - points_phi0[point - 1]
	delta_phi1 = points_phi1[point + 1] - points_phi1[point - 1]
	delta_sq = delta_r * delta_r + delta_phi0 * delta_phi0 + delta_phi1 * delta_phi1

	# vector rejection
	g_dot_delta = g_r * delta_r + g_phi0 * delta_phi0 + g_phi1 * delta_phi1 
	g_parallel_r = delta_r * (g_dot_delta / delta_sq) 
	g_parallel_phi0 = delta_phi0 * (g_dot_delta / delta_sq) 
	g_parallel_phi1 = delta_phi1 * (g_dot_delta / delta_sq) 
	g_perp_r = g_r - g_parallel_r
	g_perp_phi0 = g_phi0 - g_parallel_phi0
	g_perp_phi1 = g_phi1 - g_parallel_phi1
	g_perp_sq = g_perp_r * g_perp_r + g_perp_phi0 * g_perp_phi0 + g_perp_phi1 * g_perp_phi1
	g_perp_modulus = sqrt(g_perp_sq)

	# now apply the nudge. First, normalize it:
	g_perp_r = g_perp_r / g_perp_modulus
	g_perp_phi0 = g_perp_phi0 / g_perp_modulus
	g_perp_phi1 = g_perp_phi1 / g_perp_modulus
	
	# Then add a bit of it:
	points_r[point] += g_perp_r * 0.001
	points_phi0[point] += g_perp_phi0 * 0.001  
	points_phi1[point] += g_perp_phi1 * 0.001
	
      # end loop over points
      
    # end repeats

    # write em out
    for point in range(256):
      print(" sphere{<", phi1s[sdl_inds[basin][saddle][0]], ", ", phi0s[sdl_inds[basin][saddle][1]], ", ", rs[sdl_inds[basin][saddle][2]], ">, 0.100000 texture{ pigment {color rgb <1, 0, 0> transmit 0.000000} finish {phong 1.000000}}}")

    ## until there is no lesser energy they can be nudged to for each point 

    # then print out the list of points
  # and move on to the next path
# and the next basin

