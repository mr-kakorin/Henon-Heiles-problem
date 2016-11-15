from scipy.integrate import odeint
import numpy as np
import math
import matplotlib.pyplot as plt
import itertools

def parce_arr(path_todata):
	#path_todata - absolute path to your file with entrance-data for sode
	#data have to be in form:
	#a b w
	#q1_0 q2_0 p2_0 E_0
	#step
	#q1 poincare_map_depth
	parc_arr = []
	try:
		with open("H:\\anaconda phython\MyScripts\Henon-Heiles\entrance.txt") as f:
			parc_arr = [list(map(float, line.split())) for line in f]
	except IOError as e:
		print ("I/O error({0}): {1}".format(e.errno, e.strerror))
	except ValueError:
   		print ("Could not convert data to an integer.")
	except:
		print ("Unexpected error:"+str(sys.exc_info()[0]))
	return parc_arr

def read_parced_arr(parced_arr):
	#parced_arr - inputed entrance-data from file
	#return parameters, initial cond, step, poincare map(x, depth)
	a,b,w = parced_arr[0]
	x0,y0,py0,E0 = parced_arr[1]
	px0 = (2.*E0 - py0**2 - (w**2)*y0**2 - (w**2)*x0**2 + (2./3)*b*y0**3)
	if px0 < 0:
		raise ValueError("Incorrect entrance parameters and condtions, p1 cant be calculated")
	px0 = math.sqrt(px0)
	h = parced_arr[2][0]
	x,s = parced_arr[3]
	e = math.sqrt(E0)	
	return [[a,b,w,e],[x0,y0,px0,py0],h,[x,s]];

def ode_system(y,t,a,b,w,e):
	#generalized ode system of Henon-Heiles problem
	#y = [x,y,px,py], t = tval, a,b,w,e = parameters of system
	x,y,px,py = y
	dydt = [px,py,-(w**2)*x-2.*e*a*x*y,-(w**2)*y-e*a*x**2+e*b*y**2]
	return dydt

def poincare_section(c_phase_space, x, pdepth):
	#create poincare map for q2,p2 plane with q1 = x and depth = pdeth
	#c_phase_space - matrix of phase space calculated by solving sode
	psec = [];
	for i in range(len(c_phase_space)-1):
		if (c_phase_space[i,0]<x-pdepth and c_phase_space[i+1,0]>x+pdepth) or (c_phase_space[i,0]>x-pdepth and c_phase_space[i+1,0]<x+pdepth):
			psec.append([.5*(c_phase_space[i,1]+c_phase_space[i+1,1]),.5*(c_phase_space[i,3]+c_phase_space[i+1,3])])
	return psec

def plot_poincare_section(q,p,x,pdepth):
	#plot poincare section for q2,p2 plane with q,p calculated with poincare_section call
	#q1 fixed variable = x, depth of poincare section = pdeth
	plt.scatter(q,p)
	plt.ylabel('p2')
	plt.xlabel('q2')
	plt.title("Poincare section q1 = " + str(x) +", with depth = " + str(pdepth))
	plt.grid()
	return

def solve(y0,t,param,tol=1e-6):
	#solve Henon-Heiles problem with initial conditions = y0
	#t - tval evolution in time
	#param - vector of parametrs for sode
	#tol - tolerance using in calculating sode solution
	a,b,w,e = param
	e=1 #normilize trajectory
	y = odeint(ode_system,y0,t,args=(a,b,w,e),atol=tol)
	return y

