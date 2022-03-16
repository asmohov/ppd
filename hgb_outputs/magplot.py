import sys
sys.path.insert(0,'/home5/amohov/athena/vis/python')
import athena_read

import numpy as np
import matplotlib.pyplot as plt


def magplot(filename,col,lbl):
	data = athena_read.hst(filename)

	x_vals = data['time']
	y_vals = data['1-ME']+data['2-ME']+data['3-ME']
	y_vals =(1/(2*np.pi))* y_vals

	x_vals = x_vals/6.283

	plt.semilogy(x_vals,y_vals,col,label=lbl)
#plt.show()

#magplot('./std_res_long.hst','b:','std_res')
magplot('./std_res_long_adv.hst','b','std_res_adv')
#magplot('./2_res_long.hst','r:','2_res')
magplot('./2_res_long_adv.hst','r','2_res_adv')
magplot('./4_res_long_adv.hst','g','4_res_adv')
plt.legend()
plt.xlabel('Orbits')
plt.ylabel('Magnetic Energy Magnitude (B^2)')
plt.title('Resolution Convergence in Magnetic Energy for Orbitally Advected HGB')

plt.savefig('adv_res_conv_mag.png')

