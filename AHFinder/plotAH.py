# GRChombo
# Copyright 2012 The GRChombo collaboration.
# Please refer to LICENSE in GRChombo's root directory.
#

# script to convert the output of AHFinder into .png images for the AH and for the area and spin

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import glob, os
from os import path
import sys #argv

# expected files "coords_AH1_%04d.out" and "stats_AH1.out" and same for "AH1" -> "AH2" and "merged"

scale_mult = 0.0 # add ratio margin to plots
scale_add  = 1.0 # add margin to plots
jump  = 1 	# increase to jump through files and skip some

# Code

def readCoords():
	coords = glob.glob("coords_*")
	coords = [int(name.replace('.','_').split('_')[2]) for name in coords]
	coords = [value for value in coords if value!=0] #remove the 0

	last = np.max(coords)
	step = np.min(coords)

	num_AH = len(glob.glob("stats_*.out"))

	print("Last file found: %d\nFile step: %d\nNumber of stats files %d" % (last, step, num_AH))

	if num_AH<1:
		exit()

	return [num_AH, last, step]

[num_AH, last, step] = readCoords()

names	= [[]]*num_AH
times 	= names[:]
areas 	= names[:]
spins 	= names[:]
c_x		= names[:]
c_y		= names[:]
c_z		= names[:]

def readStats():

	for h in range(0, num_AH):
		names[h] = ("AH%d" % (h+1))

		stats = np.loadtxt("stats_" + names[h] + ".out")
		ncols = np.size(stats,1)

		times[h] = stats[:,0]

		if ncols>4:
			areas[h] = stats[:,1]
			if ncols>5:
				spins[h] = stats[:,2]

		c_x[h] = stats[:,-3]
		c_y[h] = stats[:,-2]
		c_z[h] = stats[:,-1]

def plotAH():

	mini = c_x[0][0]
	maxi = mini

	for i in range(0, last/step, jump):

		fig = plt.figure(figsize=(10,10))
		ax = fig.add_subplot(111, projection='3d', aspect='equal')
		fig.suptitle('Time = %09.4f' % (times[0][1] / step * i))
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z')

		for h in range(0, num_AH):
			# skip if first time is later than i
			if (times[h][0])/(times[h][1]-times[h][0])>i:
				continue
			# skip if no more timesteps for this AH
			if len(c_x[h])<=i:
				continue 

			print("Plotting Horizon %d, time step %d" % (h+1, i*step))

			if not path.exists('coords_%s_%04d.out' % (names[h], i*step)):
				continue

			out = np.loadtxt('coords_%s_%04d.out' % (names[h], i*step), unpack=True)

			if out.size == 0:
				continue

			[u, v, F] = out

			x_c = c_x[h][i]
			y_c = c_y[h][i]
			z_c = c_z[h][i]

			x = F*np.sin(u)*np.cos(v) 	+ x_c;
			y = F*np.sin(u)*np.sin(v) 	+ y_c;
			z = F*np.cos(u)				+ z_c;
		
			if i==0:
				mini = min(min(x),min(y),min(z),mini)
				maxi = max(max(x),max(y),max(z),maxi)
				
			# Plot the surface
			ax.scatter(x,y,z, s=(0.1 if (num_AH > 1) else 5), color = 'black')
			ax.scatter([x_c],[y_c],[z_c], s=(5 if (num_AH > 1) else 30), color = 'red')

		if i==0:
			#add scale to each side
			dx = (maxi-mini)*scale_mult + scale_add
			mini-=dx
			maxi+=dx
			print("Using (min,max) scale = (%f, %f)" % (mini, maxi))

		ax.set_xlim([mini,maxi])
		ax.set_ylim([mini,maxi])
		ax.set_zlim([mini,maxi])

		# plt.savefig("AHs_%04d.png" % ((i+1)*step), bbox_inches = 'tight', dpi=150)
		plt.savefig("AHs_%04d.png" % i, bbox_inches = 'tight', dpi=150)
		plt.close()

	make_movie = raw_input('Make video (y/n)? ')
	if make_movie=='y' or make_movie=='Y':
		os.system('ffmpeg -r 5 -f image2 -s 1920x1080 -i AHs_%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p AHs.mp4')

def plotAreas():

	if(len(areas[0])>0):
		print("Plotting areas")

		plt.figure(figsize=(12,8))

		colors = cm.rainbow(np.linspace(0, 1, num_AH))
		for h in range(0, num_AH):
			plt.scatter(times[h],areas[h], label = "Area AH%d" % (h+1), facecolor=colors[h])
		plt.ylabel("Area ")
		plt.xlabel('t $[m]$')
		plt.legend(loc= "best")
		# #plt.ylim([34,30])
		plt.savefig("areasAHs.png", bbox_inches = 'tight')
		plt.close()

def plotSpins():

	if(len(spins[0])>0):
		print("Plotting spins")

		plt.figure(figsize=(12,8))

		colors = cm.rainbow(np.linspace(0, 1, num_AH))
		for h in range(0, num_AH):
			plt.scatter(times[h],spins[h], label = "Spin AH%d" % (h+1), facecolor=colors[h])
		plt.ylabel("Spin ")
		plt.xlabel('t $[m]$')
		plt.legend(loc= "best")
		# #plt.ylim([34,30])
		plt.savefig("spinsAHs.png", bbox_inches = 'tight')
		plt.close()

readStats()
plotAreas()
plotSpins()
plotAH()
