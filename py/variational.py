#!/usr/bin/env python

from __future__ import print_function
from random import random
from math import exp
from math import floor
from math import sqrt
import pdb
import numpy as np
import traceback

def F(x,y):
  return exp(-0.025*(x**2 + y**2)) + exp(-0.025*((x-1)**2 + y**2)) + exp(-0.025*(x**2 + (y-1)**2)) - 3*exp(-0.25*((x-1)**2 + (y-1)**2))

def color(F):

  color_min = -4;
  color_max = 2;
  color_mid = (color_min + color_max) * 0.5

  if (F < color_mid):
    green = (F - color_min) / (color_mid - color_min)
    blue  = 1.0 - green
    red   = 0.0
  else:
    blue  = 0.0
    red   = (F - color_mid) / (color_max - color_mid)
    green = 1.0 - red

  color = "<" + str(red) + ", " + str(green) + ", " + str(blue) + ">"
  return color

# Main program... draw a line from the saddle at (0.5, 0.0) to basin at (1.0, 1.0) ###########################################################################################

n_variational_points = 64
n_repeats = 1024
x=[]
y=[]
x_0 = 0.5
y_0 = 0.0
x_1 = 1.0
y_1 = 1.0
step_size = 0.01

for point in range(n_variational_points):
  x.append(x_0 + point * (x_1 - x_0) / n_variational_points)
  y.append(y_0 + point * (y_1 - y_0) / n_variational_points)

for repeat in range(n_repeats):
  for point in range(1, n_variational_points - 1):
    #print(x[point], y[point], F(x[point], y[point]))
    dx = (random() - 0.5) * step_size
    dy = (random() - 0.5) * step_size
    if (F(x[point]+dx, y[point]+dy) < F(x[point],y[point])): 
      x[point] = x[point]+dx
      y[point] = y[point]+dy

for point in range(n_variational_points):
  F_color=F(x[point], y[point])
  texture ="texture{ pigment {color rgb " + color(F_color) + " transmit 0.000000} finish {phong 1.000000}}"
  print(" sphere{<", x[point], ", ", y[point], ", ", 0.0, ">, 0.012500 ", texture, "}")
  if(1+point < n_variational_points): 
    print("cylinder { <", x[point], ",", y[point], ",", 0.0, ">, <", x[point+1], ", ", y[point+1], ", ", 0.0, " >, 0.006250 open ", texture, " }") 

quit()

  
'''
#  print(" sphere{<", phi1s[min_inds[i][0]], ",", phi0s[min_inds[i][1]], ",", rs[min_inds[i][2]], ">, 0.100000 texture{ pigment {color rgb <0, 0, 1> transmit 0.000000} finish {phong 1.000000}}}")
    print("cylinder { <", phi1s[min_inds[i][0]], ",", phi0s[min_inds[i][1]], ",", rs[min_inds[i][2]], ">, <", phi1s[sdl_inds[i][j][0]], ", ", phi0s[sdl_inds[i][j][1]], ", ", rs[sdl_inds[i][j][2]], " >, 0.05 open ", texture, " }")
    #draw the saddle
    print(" sphere{<", phi1s[sdl_inds[i][j][0]], ", ", phi0s[sdl_inds[i][j][1]], ", ", rs[sdl_inds[i][j][2]], ">, 0.100000 texture{ pigment {color rgb <1, 0, 0> transmit 0.000000} finish {phong 1.000000}}}")


  F_color=F(basin_r, basin_phi0, basin_phi1)
  texture ="texture{ pigment {color rgb " + color(F_color) + " transmit 0.000000} finish {phong 1.000000}}"
  print(" sphere{<", basin_phi1, ", ", basin_phi0, ", ", basin_r, ">, 0.025000 ", texture, "}")


    saddle_phi0 = phi0s[sdl_inds[basin][saddle][1]]
    saddle_r = rs[sdl_inds[basin][saddle][2]]
    saddle_phi1 = phi1s[sdl_inds[basin][saddle][0]]
    F_color=F(saddle_r, saddle_phi0, saddle_phi1)
    texture ="texture{ pigment {color rgb " + color(F_color) + " transmit 0.000000} finish {phong 1.000000}}"
    print(" sphere{<", saddle_phi1, ", ", saddle_phi0, ", ", saddle_r, ">, 0.025000 ", texture, "}")

    # apply minimum image convention
    if ((saddle_phi0 - basin_phi0) > +phi_dimension): 
      continue
    if ((saddle_phi0 - basin_phi0) < -phi_dimension):
      continue
    if ((saddle_phi1 - basin_phi1) > +phi_dimension):
      continue
    if ((saddle_phi1 - basin_phi1) < -phi_dimension): 
      continue

    for point in range(n_variational_points):
      points_r.append(basin_r + (saddle_r - basin_r) * (point/(n_variational_points+0.0)))
      points_phi0.append(basin_phi0 + (saddle_phi0 - basin_phi0) * (point/(n_variational_points + 0.0)))
      points_phi1.append(basin_phi1 + (saddle_phi1 - basin_phi1) * (point/(n_variational_points + 0.0)))
      new_points_r.append(basin_r + (saddle_r - basin_r) * (point/(n_variational_points+0.0)))
      new_points_phi0.append(basin_phi0 + (saddle_phi0 - basin_phi0) * (point/(n_variational_points + 0.0)))
      new_points_phi1.append(basin_phi1 + (saddle_phi1 - basin_phi1) * (point/(n_variational_points + 0.0)))


	# vector rejection
	g_dot_delta = g_r * delta_r + g_phi0 * delta_phi0 + g_phi1 * delta_phi1 
	g_parallel_r = delta_r * (g_dot_delta / delta_sq) 
	g_parallel_phi0 = delta_phi0 * (g_dot_delta / delta_sq) 
	g_parallel_phi1 = delta_phi1 * (g_dot_delta / delta_sq) 
	g_perp_r = g_r - g_parallel_r
	g_perp_phi0 = g_phi0 - g_parallel_phi0
	g_perp_phi1 = g_phi1 - g_parallel_phi1
	g_perp_sq = g_perp_r * g_perp_r + g_perp_phi0 * g_perp_phi0 + g_perp_phi1 * g_perp_phi1
#	g_perp_modulus = sqrt(g_perp_sq)

	# now apply the nudge. First, normalize it:
#	g_perp_r = g_perp_r / g_perp_modulus
#	g_perp_phi0 = g_perp_phi0 / g_perp_modulus
#	g_perp_phi1 = g_perp_phi1 / g_perp_modulus
	
	# Then add a bit of it:
	# new_points_r[point]    = points_r[point] + g_perp_r * 0.001
	# new_points_phi0[point] = points_phi0[point] + g_perp_phi0 * 0.001  
	# new_points_phi1[point] = points_phi1[point] + g_perp_phi1 * 0.001

	# Then add a bit of it:
#        print(g_perp_r, g_perp_phi0, g_perp_phi1)
	if (g_perp_r**2 < .0001):    new_points_r[point]    = points_r[point] + g_perp_r
        else: g_perp_r += 0.01
	if (g_perp_phi0**2 < .0001): new_points_phi0[point] = points_phi0[point] + g_perp_phi0
        else: g_perp_phi0 += 0.01  
	if (g_perp_phi1**2 < .0001): new_points_phi1[point] = points_phi1[point] + g_perp_phi1
        else: g_perp_phi1 += 0.01  
	
      # end loop over points

      for point in range(1,n_variational_points - 1):
        points_r[point]    = new_points_r[point]
        points_phi0[point] = new_points_phi0[point]
        points_phi1[point] = new_points_phi1[point]
      
    # end repeats

    # write em out
    for point in range(n_variational_points-1):
      #print(" sphere{<", phi1s[sdl_inds[basin][saddle][0]], ", ", phi0s[sdl_inds[basin][saddle][1]], ", ", rs[sdl_inds[basin][saddle][2]], ">, 0.100000 texture{ pigment {color rgb <1, 0, 0> transmit 0.000000} finish {phong 1.000000}}}")

#      print(" sphere{<", points_phi1[point], ", ", points_phi0[point], ", ", points_r[point], ">, 0.100000 texture{ pigment {color rgb <1, 0, 0> transmit 0.000000} finish {phong 1.000000}}}")
#      texture ="texture{ pigment {color rgb <1, 2, 3> transmit 0.000000} finish {phong 1.000000}}"
      F_color=F(points_r[point], points_phi0[point], points_phi1[point])
      texture ="texture{ pigment {color rgb " + color(F_color) + " transmit 0.000000} finish {phong 1.000000}}"
      print("cylinder { <", points_phi1[point], ",", points_phi0[point], ",", points_r[point], ">, <", points_phi1[point + 1], ", ", points_phi0[point + 1], ", ", points_r[point + 1], " >, 0.025 open ", texture, " }")

    ## until there is no lesser energy they can be nudged to for each point 

    # then print out the list of points
  # and move on to the next path
# and the next basin
'''
