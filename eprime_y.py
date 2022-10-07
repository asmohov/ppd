import sys
sys.path.append('~/athena-public-version/vis/python/')
#sys.path.append('~/.local/lib/python3.8/site-packages/')
sys.path.append('~/working')


import numpy as np
import matplotlib.pyplot as plt

#IMPORT APPROPRIATE ATHINPUT FILE
#import athinput.hgb as athin
import athena_read



def eprime_y(file_name):

    data = athena_read.athdf(file_name)
    
    #for 4x4x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 4/len(data['x1v'])
    #print(side_length,' side length')
    volume = side_length**3
    #print(volume,' volume')
    
    #constants
    omega0 = 1.0
    qshear = 1.5

    
    #arrays for plotting to visualize v_shear
    vy_arr = []
    x_arr = []
    v_shear = []
    counter = 0
    #assuming 64x256x256, but should work for any 
    overall_length = len(data['x1v'])*len(data['x2v'])*len(data['x3v'])
    eprime_arr = np.zeros(overall_length)
    for k in range(len(data['x3v'])):
        #iterate over z, store current z
        z = data['x3v'][k]
        

        for j in range(len(data['x2v'])):
            #iterate over y, store current y
            y = data['x2v'][j]

            for i in range(len(data['x1v'])):
                #iterate over x, retrieve x value from x1v array
                x = data['x1v'][i]
                #retrieve indexed density, y velocity
                rho = data['rho'][k][j][i]
                vy = data['vel2'][k][j][i]
                #calculate cell centered eprime for chosen cell
                vsh = -qshear*omega0*x
                eprime_arr[counter]=(((vy-vsh)**2)*rho)/2
                counter+=1
                #append for plotting purposes for fiducial y,z
                #if y == 0.0078125 and z == 0.0078125:
                    #x_arr.append(x)
                    #vy_arr.append(vy)
                    #v_shear.append(-qshear*omega0*x)
    #Plotting section --------------------------------------------------------
    #Comment out if running code iteratively, i.e. for more than one time step
    #print('x_arr length is ',len(x_arr))
    #print('vy_arr length is ',len(vy_arr))
    #print('eprime_arr length is ',len(eprime_arr))
    #plt.plot(x_arr,vy_arr,label = 'numerical vy')
    #plt.plot(x_arr,v_shear,label = 'shearing velocity')
    #plt.legend()
    #plt.ylabel('vy')
    #plt.xlabel('radial component x')
    #plt.show()
    #print('eprime is ',np.average(eprime_arr))    
    #calculate and plot zy average
    #eprime_arr = np.array(eprime_arr)
    #eprime_arr = np.reshape(eprime_arr,(16384,256))
    #eprime_zyavg = []
    #for i in range(256):
        #eprime_zyavg.append(np.average(eprime_arr[:,i]))
    #plt.plot(x_arr,eprime_zyavg)
    #plt.ylabel('zy avg of eprime')
    #plt.xlabel('radial component x')
    #plt.show()
    #sec_eprime = np.average(eprime_zyavg)
    #print('secondary eprime avg is ',sec_eprime
    #print('eprime min is ',min(eprime_arr))
    #print('eprime length is ',len(eprime_arr))
    #print('zeros array length is ',overall_length)
    #print('entering eprime_y loop')
    return(np.average(eprime_arr))