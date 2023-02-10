import sys
import time
sys.path.append('~/athena-public-version/vis/python/')
#sys.path.append('~/.local/lib/python3.8/site-packages/')
sys.path.append('~/working')


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#IMPORT APPROPRIATE ATHINPUT FILE
#import athinput.hgb as athin
import athena_read

#alpha profiles-------------------------------------------------------------------------------
def oned_alpha_profile(file_name):
    #print('current file is :'+file_name)
    data = []
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
    #adjust vel y minus the shear velocity
    vely = data['vel2']-vsh_3d
    #reynolds stress
    reyn = data['rho']*vely*data['vel1']
    #average across y and z
    reyn = np.sum(reyn,axis=(0,1))/(Nz*Ny)
    #divide by density profile
    rho_prof = np.sum(data['rho'],axis=(0,1))/(Nz*Ny)
    reyn = reyn/rho_prof
    #maxewell stress -BxBy
    maxw=-1*data['Bcc2']*data['Bcc1']
    maxw=np.sum(maxw,axis=(0,1))/(Nz*Ny)
    maxw= maxw/rho_prof
    return(reyn,maxw)

#print(oned_stress_profile('./ad_prof/amp_1/sig_1/HGB.out2.00100.athdf'))
#iterate over last 20 outputs and return reynolds, maxwell, total alpha profiles
def avg_alpha_prof(file_path):
    prof_list_reyn = []
    prof_list_maxw = []
    #x_arr = np.linspace(-4,4,512)
    for i in range(20):
        #print('passing step ',i)
        if i == 0:
            fname = file_path+'/HGB.out2.00100.athdf'
        else:
            fname = file_path+'/HGB.out2.000'+str(100-i)+'.athdf'
        lbl = str(100-i)
        reyn_list,maxw_list=oned_alpha_profile(fname)
        prof_list_reyn.append(reyn_list)
        prof_list_maxw.append(maxw_list)
    #convert to numpy
    prof_list_reyn = np.array(prof_list_reyn)
    prof_list_maxw = np.array(prof_list_maxw)
    #plt.plot(x_arr,prof_list[i],color = 'c')

    prof_avg_reyn = np.sum(prof_list_reyn,axis=0)/len(prof_list_reyn)
    prof_avg_maxw = np.sum(prof_list_maxw,axis=0)/len(prof_list_maxw)
    prof_avg_tot = prof_avg_reyn + prof_avg_maxw
    return prof_avg_reyn,prof_avg_maxw,prof_avg_tot
#stress profiles-------------------------------------------------------------------------------
def oned_stress_profile(file_name):
    #print('current file is :'+file_name)
    data = []
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
    #adjust vel y minus the shear velocity
    vely = data['vel2']-vsh_3d
    #reynolds stress
    reyn = data['rho']*vely*data['vel1']
    #average across y and z
    reyn = np.sum(reyn,axis=(0,1))/(Nz*Ny)
    #maxewell stress -BxBy
    maxw=-1*data['Bcc2']*data['Bcc1']
    maxw=np.sum(maxw,axis=(0,1))/(Nz*Ny)
    return(reyn,maxw)
def avg_stress_prof(file_path):
    prof_list_reyn = []
    prof_list_maxw = []
    #x_arr = np.linspace(-4,4,512)
    for i in range(20):
        #print('passing step ',i)
        if i == 0:
            fname = file_path+'/HGB.out2.00100.athdf'
        else:
            fname = file_path+'/HGB.out2.000'+str(100-i)+'.athdf'
        lbl = str(100-i)
        reyn_list,maxw_list=oned_stress_profile(fname)
        prof_list_reyn.append(reyn_list)
        prof_list_maxw.append(maxw_list)
    #convert to numpy
    prof_list_reyn = np.array(prof_list_reyn)
    prof_list_maxw = np.array(prof_list_maxw)
    #plt.plot(x_arr,prof_list[i],color = 'c')

    prof_avg_reyn = np.sum(prof_list_reyn,axis=0)/len(prof_list_reyn)
    prof_avg_maxw = np.sum(prof_list_maxw,axis=0)/len(prof_list_maxw)
    prof_avg_tot = prof_avg_reyn + prof_avg_maxw
    return prof_avg_reyn,prof_avg_maxw,prof_avg_tot
#density-------------------------------------------------------------------------------
def oned_rho_profile(file_name):
    #print('current file is :'+file_name)
    data = []
    data = athena_read.athdf(file_name)
    #print(data)
    #for 8x8x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 1/len(data['x3v'])
    #print(side_length,' side length')
    volume = side_length**3
    
    Nx = len(data['x1v'])
    Ny = len(data['x2v'])
    Nz = len(data['x3v'])
    
    #assuming 64x256x256, but should work for any 
    overall_length = Nx*Ny*Nz
    data_arr = []

    rho_prof = np.sum(data['rho'],axis=(0,1))/(Nz*Ny)
    return(rho_prof)
def avg_rho_prof(file_path):
    prof_rho = []
    prof_rho_sig = []
    #x_arr = np.linspace(-4,4,512)
    for i in range(20):
        #print('passing step ',i)
        if i == 0:
            fname = file_path+'/HGB.out2.00100.athdf'
        else:
            fname = file_path+'/HGB.out2.000'+str(100-i)+'.athdf'
        lbl = str(100-i)
        rho_list=oned_rho_profile(fname)
        prof_rho.append(rho_list)

    #convert to numpy
    prof_rho = np.array(prof_rho)
    #plt.plot(x_arr,prof_list[i],color = 'c')
    #stddeviation
    prof_rho_sig = np.std(prof_rho,axis=0)
    #average
    prof_rho = np.sum(prof_rho,axis=0)/len(prof_rho)
    return prof_rho,prof_rho_sig
#magnetic field-------------------------------------------------------
def oned_mag_profile(file_name):
    #print('current file is :'+file_name)
    data = []
    data = athena_read.athdf(file_name)
    #print(data)
    #for 8x8x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 1/len(data['x3v'])
    #print(side_length,' side length')
    volume = side_length**3
    
    Nx = len(data['x1v'])
    Ny = len(data['x2v'])
    Nz = len(data['x3v'])
    
    #assuming 64x256x256, but should work for any 
    overall_length = Nx*Ny*Nz
    #print(data['Bcc1'])
    bx = np.sum(data['Bcc1'],axis=(0,1))/(Nz*Ny)
    by = np.sum(data['Bcc2'],axis=(0,1))/(Nz*Ny)
    bz = np.sum(data['Bcc3'],axis=(0,1))/(Nz*Ny)
    bmag = np.sqrt(bx*bx+by*by+bz*bz)
    
    return bx,by,bz,bmag
def avg_mag_prof(file_path):
    prof_bx = []
    prof_by = []
    prof_bz = []
    prof_bmag = []

    #x_arr = np.linspace(-4,4,512)
    for i in range(20):
        #print('passing step ',i)
        if i == 0:
            fname = file_path+'/HGB.out2.00100.athdf'
        else:
            fname = file_path+'/HGB.out2.000'+str(100-i)+'.athdf'
        bx,by,bz,bmag = oned_mag_profile(fname)
        prof_bx.append(bx)
        prof_by.append(by)
        prof_bz.append(bz)
        prof_bmag.append(bmag)


    #convert to numpy
    prof_bx = np.array(prof_bx)
    prof_by = np.array(prof_by)
    prof_bz = np.array(prof_bz)
    prof_bmag = np.array(prof_bmag)

    #average
    prof_bx = np.sum(prof_bx,axis=0)/len(prof_bx)
    prof_by = np.sum(prof_by,axis=0)/len(prof_by)
    prof_bz = np.sum(prof_bz,axis=0)/len(prof_bz)
    prof_bmag = np.sum(prof_bmag,axis=0)/len(prof_bmag)

    return prof_bx,prof_by,prof_bz,prof_bmag