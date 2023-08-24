#take data name as string and a file name to make a surface profile

import sys
sys.path.append('~/athena-public-version/vis/python/')
#sys.path.append('~/.local/lib/python3.8/site-packages/')
sys.path.append('~/working')


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#IMPORT APPROPRIATE ATHINPUT FILE
#import athinput.hgb as athin
import athena_read
from celluloid import Camera
fig2 = plt.figure()
cam = Camera(fig2)
def surface_profile(file_name,data_name):
    #print('current file is :'+file_name)
    title=file_name[10:-20]
    title = data_name + ' surface profile for ' +title
    data = athena_read.athdf(file_name)
    #print(data)
    #for 8x8x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 1/len(data['x3v'])
    #print(side_length,' side length')
    volume = side_length**3
    #print(volume,' volume')
    
    #constants
    omega0 = 1.0
    qshear = 1.5
    Nx = len(data['x1v'])
    Ny = len(data['x2v'])
    Nz = len(data['x3v'])
    
    #assuming 64x256x256, but should work for any 
    overall_length = Nx*Ny*Nz
    data_arr = []
    #array of vshear possible values
    vsh = -qshear*omega0*data['x1v']
    vsh_3d = np.broadcast_to(vsh,(Nz,Ny,Nx))
    vel_y = data['vel2']-vsh_3d
    if data_name == '-BxBy':
        data_arr = -1*data['Bcc1']*data['Bcc2']
    elif data_name == 'rhovxvy':
        data_arr = data['vel1']*data['rho']*vel_y
    elif data_name == 'bmag':
        data_arr = data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']
 
    else:
        data_arr = data[data_name]
    if data_name == 'vel2':
        data_arr = data_arr-vsh_3d
    data_arr = np.sum(data_arr,axis=(0))
    radial_dim = data['x1f'][-1]-data['x1f'][0]
    div_vol = radial_dim*side_length
    #divide by number of cells in vertical slice
    data_arr = data_arr/(Nz)
    #print(np.shape(data_arr))

    #plotting section
    #plt.ioff()
    t=plt.pcolormesh(data_arr,norm=mpl.colors.Normalize(vmin=min(data_arr.flatten()),vmax=max(data_arr.flatten())),shading = 'gouraud',cmap = 'RdBu_r')
    #plt.colorbar(t,label="$\\rho/\\rho_0$")
    plt.xlabel('x cell')
    plt.ylabel('y cell')
    plt.title(title)

    return(t