#import sys
#sys.path.append('~/athena-public-version/vis/python/')
#sys.path.append('~/.local/lib/python3.8/site-packages/')
#sys.path.append('~/working')

import numpy as np
#import matplotlib.pyplot as plt


#IMPORT APPROPRIATE ATHINPUT FILE
#import athinput.hgb as athin


import athena_read


def init_hgb():
    #necessary imports
    import matplotlib.pyplot as plt

    fig, ax=plt.subplots(2,2, sharex=False, sharey = False, constrained_layout=True)

    fig.set_figheight(10)
    fig.set_figwidth(10)
    print('Figure initialised!')
    return fig,ax


def hgb_plot(file_name,labl,fig,ax,vol=2*np.pi):

    #necessary imports
    import sys
    sys.path.append('~/athena-public-version/vis/python/')
    #sys.path.append('~/.local/lib/python3.8/site-packages/')
    sys.path.append('~/working')

    import numpy as np
    import matplotlib.pyplot as plt


    data = athena_read.hst(file_name)
    #set x values as time
    x_vals=((data['time']))

    #to match HGB, divide by 2 pi, because of difference in code units: 1/omega vs period
    x_vals_adj = x_vals/(2*np.pi)
    #Magnetic Stress

    #Y_vals are BxBy dsata from HGB, still pending correction factor
    y_vals=(data['-BxBy'])
    
    #set length of x_vals for the 4x dimension to be equal to whatever hst length we currently have

    x_vals_adj_4x = x_vals_adj[0:len(y_vals)]

    #Scale to HGB (change code units for B, divide by volume for volume avg)
    volume = vol

    y_vals= y_vals /volume
    

    ax[1,0].semilogy(x_vals_adj_4x,y_vals,label=labl)
    ax[1,0].set(xlabel='time(periods)',ylabel='-BxBy',title = 'Average Magnetic Stress')
    ax[1,0].grid()
    ax[1,0].legend()
    
    

    
    #-------------------------------------------------------------
    #Copy for B^2
    mex = data['1-ME']
    mey = data['2-ME']
    mez = data['3-ME']
    me_total = (mex+mey+mez)
    y_vals=(me_total)
    y_vals= y_vals / volume
    #figure
    ax[0,0].semilogy(x_vals_adj_4x,y_vals)
    ax[0,0].set(xlabel='time(orbits)',ylabel='B^2',title = 'Magnetic Magnitude')
    ax[0,0].grid()
    
    #-------------------------------------------------
    #Copy for VxVy
    #Y_vals are dVxVy  dsata from HGB, still pending correction factor
    y_vals=(data['dVxVy'])

    #Scale to HGB (divide by volume for volume avg)
    y_vals= y_vals / volume
    #figure
    ax[1,1].plot(x_vals_adj_4x,y_vals)
    ax[1,1].set(xlabel='time(orbits)',ylabel='dVxVy',title = 'Reynolds Stress')
    ax[1,1].grid()
    
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
    #Scale to HGB (divide by volume for volume avg)
    y_vals= y_vals / volume
    
    #figure
    
    ax[0,1].semilogy(x_vals_adj_4x,y_vals)
    ax[0,1].set(xlabel='time(orbits)',ylabel='v^2',title = 'Average Kinetic Energy')
    ax[0,1].grid()