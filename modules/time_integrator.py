



import sys
sys.path.insert(0,'./athena-public-version/vis/python')
sys.path.insert(0,'./.local/lib/python3.8/site-packages')
sys.path.insert(0,'./working')

import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
from scipy import integrate
#IMPORT APPROPRIATE ATHINPUT FILE
#import athinput.hgb as athin


import athena_read
data = athena_read.hst('./working/HGB.hst')

#Load other resolution data sets
data_std_res = athena_read.hst('./working/HGB_std_res.hst')
data_double_res = athena_read.hst('./working/HGB_double_res.hst')

#set x values as time
x_vals=((data_std_res['time']))

#define saturation point
sat=500

#to match HGB, divide by 2 pi, because of difference in code units: 1/omega vs period
x_vals_adj = x_vals/(2*np.pi)
#print(x_vals)

#Calculate volume from athinput variables
#x_len = athin.x1max-athin.x1min
#y_len = athin.x2max-athin.x2min
#z_len = athin.x3max-athin.x3max
#volume = x_len*y_len*z_len
#print(volume)
#-------------------------------------------------------------------------
#Magnetic Stress

#Y_vals are BxBy dsata from HGB, still pending correction factor
y_vals=(data['-BxBy'])
y_vals_std_res=(data_std_res['-BxBy'])
y_vals_double_res=(data_double_res['-BxBy'])

#set length of x_vals for the 4x dimension to be equal to whatever hst length we currently have

x_vals_adj_4x = x_vals_adj[0:len(y_vals)]

#For time integration, assume saturation occurs halfway through the set

x_vals_adj = x_vals_adj[sat:]


#Scale to HGB (change code units for B, divide by volume for volume avg)
volume = np.pi*2

y_vals= y_vals /volume
y_vals_std_res= y_vals_std_res[sat:] /volume
y_vals_double_res= y_vals_double_res[sat:] /volume


print('y set lengths')
print(len(y_vals))
print(len(y_vals_std_res))
print(len(y_vals_double_res))

print('x set lengths')
print(len(x_vals_adj))
print(len(x_vals_adj_4x))


int_mag_stress_std_res = scipy.integrate.simpson(y_vals_std_res,x_vals_adj)
avg_mag_stress_std_res = int_mag_stress_std_res/5
int_mag_stress_double_res = scipy.integrate.simpson(y_vals_double_res,x_vals_adj)
avg_mag_stress_double_res = int_mag_stress_double_res/5

print('time avg mag stress std res')
print(avg_mag_stress_std_res)

print('time avg mag stress double res')
print(avg_mag_stress_double_res)
#-------------------------------------------------------------
#Copy for B^2

#Y_vals are BxBy dsata from HGB, still pending correction factor

#ADJUST FOR B^2 (magnitude of columns 10-12)

mex = data['1-ME']
mey = data['2-ME']
mez = data['3-ME']

#print('Length of b_x')
#print(len(mex))

me_total = (mex+mey+mez)

#print('Length of b_2')
#print(len(me_total))
y_vals=(me_total)
#-------------------------
#For std res
mex = data_std_res['1-ME']
mey = data_std_res['2-ME']
mez = data_std_res['3-ME']

#print('Length of b_x')
#print(len(mex))

me_total = (mex+mey+mez)

#print('Length of b_2')
#print(len(me_total))
y_vals_std_res=(me_total)


#------------------------
#For double res
mex = data_double_res['1-ME']
mey = data_double_res['2-ME']
mez = data_double_res['3-ME']

#print('Length of b_x')
#print(len(mex))

me_total = (mex+mey+mez)

#print('Length of b_2')
#print(len(me_total))
y_vals_double_res=(me_total)


#Scale to HGB (divide by volume for volume avg)
y_vals= y_vals / volume
y_vals_std_res= y_vals_std_res[sat:] / volume
y_vals_double_res= y_vals_double_res[sat:] / volume


int_mag_e_std_res = scipy.integrate.simpson(y_vals_std_res,x_vals_adj)
avg_mag_e_std_res = int_mag_e_std_res/5
int_mag_e_double_res = scipy.integrate.simpson(y_vals_double_res,x_vals_adj)
avg_mag_e_double_res = int_mag_e_double_res/5

print('time avg mag energy std res')
print(avg_mag_e_std_res)

print('time avg mag energy double res')
print(avg_mag_e_double_res)


#-------------------------------------------------
#Copy for VxVy
#Y_vals are dVxVy  dsata from HGB, still pending correction factor
y_vals=(data['dVxVy'])
y_vals_std_res=(data_std_res['dVxVy'])
y_vals_double_res=(data_double_res['dVxVy'])
#Scale to HGB (divide by volume for volume avg)
y_vals= y_vals / volume
y_vals_std_res= y_vals_std_res[sat:] / volume
y_vals_double_res= y_vals_double_res[sat:] / volume

int_rey_stress_std_res = scipy.integrate.simpson(y_vals_std_res,x_vals_adj)
avg_rey_stress_std_res = int_rey_stress_std_res/5
int_rey_stress_double_res = scipy.integrate.simpson(y_vals_double_res,x_vals_adj)
avg_rey_stress_double_res = int_rey_stress_double_res/5

print('time avg rey stress std res')
print(avg_rey_stress_std_res)

print('time avg rey stress double res')
print(avg_rey_stress_double_res)


#-------------------------------------------------
#Copy for v^2
#Y_vals are KE column data from HGB, still pending correction factor

mass = data['mass'][4]

ke_x = data['1-KE']
ke_y = data['2-KE']
ke_z = data['3-KE']
ke_sum = ke_x+ke_y+ke_z
v_2_div_2 = ke_sum

y_vals = v_2_div_2
#copy for other resolutions

ke_x = data_std_res['1-KE']
ke_y = data_std_res['2-KE']
ke_z = data_std_res['3-KE']
ke_sum = ke_x+ke_y+ke_z
v_2_div_2_std_res = ke_sum


y_vals_std_res=v_2_div_2_std_res
#copy for other resolutions
ke_x = data_double_res['1-KE']
ke_y = data_double_res['2-KE']
ke_z = data_double_res['3-KE']
ke_sum = ke_x+ke_y+ke_z
v_2_div_2_double_res = ke_sum


y_vals_double_res=v_2_div_2_double_res

#Scale to HGB (divide by volume for volume avg)
y_vals= y_vals / volume
y_vals_std_res= y_vals_std_res[sat:] / volume
y_vals_double_res= y_vals_double_res[sat:] / volume


int_ke_std_res = scipy.integrate.simpson(y_vals_std_res,x_vals_adj)
avg_ke_std_res = int_ke_std_res/5
int_ke_double_res = scipy.integrate.simpson(y_vals_double_res,x_vals_adj)
avg_ke_double_res = int_ke_double_res/5

print('time avg kinetic energy std res')
print(avg_ke_std_res)

print('time avg kinetic energy double res')
print(avg_ke_double_res)


#----------------------------------

#Separate dt graph, log scaled



# Data for plotting
#dt = data['dt']
#y_vals=np.log10(dt)

#print(dt[2])
#print(y_vals[2])
#dt_std_res = data_std_res['dt']
#y_vals_std_res=np.log10(dt_std_res)

#dt_double_res = data_double_res['dt']
#y_vals_double_res=np.log10(dt_double_res)

        


