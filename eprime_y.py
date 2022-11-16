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




def eprime_y(file_name):
    data = athena_read.athdf(file_name)
    
    #for 4x4x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 8/len(data['x1v'])
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
    eprime_arr = []
    for k in range(len(data['x3v'])):
        #iterate over z, store current z
        z = data['x3v'][k]
        

        for j in range(len(data['x2v'])):
            #iterate over y, store current y
            y = data['x2v'][j]

            for i in range(len(data['x1v'])):
                #reset variables
                vsh = 0
                vy = 0
                rho = 0
                #iterate over x, retrieve x value from x1v array
                x = data['x1v'][i]
                #retrieve indexed density, y velocity
                rho = data['rho'][k][j][i]
                vy = data['vel2'][k][j][i]
                #calculate cell centered eprime for chosen cell
                vsh = -qshear*omega0*x
                eprime_arr.append((((vy-vsh)**2)*rho)/2)
                counter+=1
                #append for plotting purposes for fiducial y,z
                #if y == 0.0078125 and z == 0.0078125:
                    #print('appending at fiducial yz')
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



def eprime_x(file_name):
    data = []
    data = athena_read.athdf(file_name)
    
    #for 4x4x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 8/len(data['x1v'])
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
                #reset variables
                vsh = 0
                vx = 0
                rho = 0
                #iterate over x, retrieve x value from x1v array
                x = data['x1v'][i]
                #retrieve indexed density, y velocity
                rho = data['rho'][k][j][i]
                vx = data['vel1'][k][j][i]
                #calculate cell centered eprime for chosen cell
                vsh = -qshear*omega0*x
                eprime_arr[counter]=(((vx)**2)*rho)/2
                counter+=1
    return(np.average(eprime_arr))



def eprime_z(file_name):
    data = []
    data = athena_read.athdf(file_name)
    
    #for 4x4x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 8/len(data['x1v'])
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
                #reset variables
                vsh = 0
                vz = 0
                rho = 0
                #iterate over x, retrieve x value from x1v array
                x = data['x1v'][i]
                #retrieve indexed density, y velocity
                rho = data['rho'][k][j][i]
                vz = data['vel3'][k][j][i]
                #calculate cell centered eprime for chosen cell
                vsh = -qshear*omega0*x
                eprime_arr[counter]=(((vz)**2)*rho)/2
                counter+=1
    return(np.average(eprime_arr))

def v_over_time_plot(file_name,fig,ax):
    data = athena_read.athdf(file_name)
    
    #for 4x4x1 scale height box, with cubic cells, needs to be adjusted for other sizes
    side_length = 8/len(data['x1v'])
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
    eprime_arr = []
    for k in range(len(data['x3v'])):
        #iterate over z, store current z
        z = data['x3v'][k]
        

        for j in range(len(data['x2v'])):
            #iterate over y, store current y
            y = data['x2v'][j]

            for i in range(len(data['x1v'])):
                #reset variables
                vsh = 0
                vy = 0
                rho = 0
                #iterate over x, retrieve x value from x1v array
                x = data['x1v'][i]
                
                #retrieve indexed density, y velocity
                rho = data['rho'][k][j][i]
                vy = data['vel2'][k][j][i]
                #calculate cell centered eprime for chosen cell
                vsh = -qshear*omega0*x
                eprime_arr.append((rho))
                counter+=1
    labl = file_name[file_name.find('00'):]
    #calculate and plot zy average
    eprime_arr = np.array(eprime_arr)
    eprime_arr = np.reshape(eprime_arr,(16384,256))
    eprime_zyavg = []
    for i in range(256):
        eprime_zyavg.append(np.average(eprime_arr[:,i]))
    ax.plot(data['x1v'],eprime_zyavg,label = labl)
    #plt.show()
    return(np.average(eprime_arr))

#AZIMUTHAL AVERAGE OVER RADIAL COORDINATE (AVERAGED IN YZ PLANE FOR INDICIAL X)

def az_avg(file_name,data_name):
    #print('current file is :'+file_name)
    data = athena_read.athdf(file_name)
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
    data_arr = data[data_name]
    if data_name == 'vel2':
        data_arr = data_arr-vsh_3d
    data_arr = np.sum(data_arr,axis=(0,1))
    radial_dim = data['x1f'][-1]-data['x1f'][0]
    div_vol = radial_dim*side_length
    #divide by number of cells in azimuthal slice
    data_arr = data_arr/(Ny*Nz)
    return(data_arr)

#creates .npz file with calculated az_avg arrays for density and vy
def az_avg_rho_vy(file_path,output_name,min_orbits = 20,max_orbits = 100):
    arr=range(min_orbits,max_orbits+1)
    radial_coords = athena_read.athdf(file_path+'/HGB.out2.00002.athdf')['x1v']
    contour_list = []
    radial_dim = len(athena_read.athdf(file_path+'/HGB.out2.00002.athdf')['x1v'])
    print('orbits selected: ',arr)
    for i in arr:
        #print(file_path+'/HGB.out2.000'+str(i)+'.athdf')
        if i<10:
            contour_list.append(az_avg(file_path+'/HGB.out2.0000'+str(i)+'.athdf','rho'))
        elif 9<i<100:
            contour_list.append(az_avg(file_path+'/HGB.out2.000'+str(i)+'.athdf','rho'))
        else:
            contour_list.append(az_avg(file_path+'/HGB.out2.00'+str(i)+'.athdf','rho'))
        if i%25 == 0:
            print('passing step ',i)


    radial_dim = len(athena_read.athdf(file_path+'/HGB.out2.00002.athdf')['x1v'])
    contour_list= np.vstack(contour_list)
    print(np.shape(contour_list))
    az_avg_rho = np.transpose(contour_list)

    #vy-vshear, 8x8 AM = 1
    print('rho calculation done')

    contour_list = []

    for i in arr:
        if i<10:
            contour_list.append(az_avg(file_path+'/HGB.out2.0000'+str(i)+'.athdf','vel2'))
        elif 9<i<100:
            contour_list.append(az_avg(file_path+'/HGB.out2.000'+str(i)+'.athdf','vel2'))
        else:
            #print('current file at 100 orbit is :')
            #print(file_path+'/HGB.out2.00'+str(i)+'.athdf')
            #val = az_avg(file_path+'/HGB.out2.00'+str(i)+'.athdf','vel2')
            #print('current az_avg list at 100 orbit is :')
            #print(val)                         
            contour_list.append(az_avg(file_path+'/HGB.out2.00'+str(i)+'.athdf','vel2'))
        if i%25 == 0:
            print('passing step ',i)

    radial_dim = len(athena_read.athdf(file_path+'/HGB.out2.00002.athdf')['x1v'])
    contour_list= np.vstack(contour_list)
    print(np.shape(contour_list))
    az_avg_vy = np.transpose(contour_list)
    #print('az avg vy at calculation and saving step is:')
    #print(az_avg_vy)
    
    np.savez(output_name,az_avg_vy=az_avg_vy,az_avg_rho = az_avg_rho)
    print('finished calculating az_avg for ',output_name)
    
def az_avg_plotter(npz_name):
    #for vy-vshear
    fig, ax=plt.subplots(sharex=False, sharey = False, constrained_layout=True)

    fig.set_figheight(5)
    fig.set_figwidth(10)
    ax.set_xlabel('orbit number')
    ax.set_ylabel('radial coordinate x')
    title = 'Azimuthal Average Vy - Vshear for '+npz_name 
    fig.suptitle(title)
    #load dataset
    npzfile = np.load(npz_name)
    data = npzfile['az_avg_vy']
    
    
    #axes values
    num_rows,num_col = np.shape(data)

    #for normal runs
    #arr=range(20,20+num_col)
    #for devel runs
    arr = range(0,num_col)
    
    
    
    #hard code for 4x4 box
    radial_coords = athena_read.athdf('./ad_prof/const_am/hundred/HGB.out2.00002.athdf')['x1v']
    
    #hard coded for 8x8 box
    #radial_coords = athena_read.athdf('./ad_prof/const_am_big/1point0/HGB.out2.00044.athdf')['x1v']

    plt.pcolormesh(arr,radial_coords,data,norm=mpl.colors.CenteredNorm(vcenter =0),shading = 'gouraud',cmap = 'RdBu_r')
    plt.colorbar()
    plt.show()

    #same but for rho
    fig, ax=plt.subplots(sharex=False, sharey = False, constrained_layout=True)

    fig.set_figheight(5)
    fig.set_figwidth(10)
    ax.set_xlabel('orbit number')
    ax.set_ylabel('radial coordinate x')
    title = 'Azimuthal Average Rho for '+npz_name 
    fig.suptitle(title)
    #load dataset

    #axes values
    print(data.shape)
    num_rows,num_col = np.shape(data)
    data = npzfile['az_avg_rho']
    npzfile.close()
    #scaling statement, if desired
    #norm=mpl.colors.CenteredNorm(vcenter =1),
    plt.pcolormesh(arr,radial_coords,data,norm=mpl.colors.CenteredNorm(vcenter =1),shading = 'gouraud',cmap = 'RdBu_r')
    plt.colorbar()
    plt.show()