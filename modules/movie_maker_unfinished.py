#MOVIE MAKING SECTION
#turn into a method that takes a file path,data type (see surface profiles) and a range tuple
#using celluloid
import sys
sys.path.append('~/ppd/modules')
from surface_profiles import surface_profile
import matplotlib as mpl
import time
from matplotlib import pyplot as plt
from celluloid import Camera as cam
import ffmpeg
fig = plt.figure()
start = time.time()
outrange = range(0,101)
z=[]
Cam = cam(fig)

file_path = './ad_prof/amp_50/'
for i in outrange:
        if z != []:
            z.remove()
        #plt.clf()
        if i<10:
            fname = file_path+'/HGB.out2.0000'+str(i)+'.athdf'
        elif 9<i<100:
            fname = file_path+'/HGB.out2.000'+str(i)+'.athdf'
        else:
            fname = file_path+'/HGB.out2.00'+str(i)+'.athdf'
        if i%25 == 0:
            print('passing step ',i)
        x=surface_profile(fname,'Bcc3')
        z =plt.colorbar(label="$\\rho/\\rho_0$",location='right')
        Cam.snap()

print('Beginning animation')
anim = Cam.animate(blit=False,interval=300)
anim.save('am_50_flat_video_bz.mp4',dpi=700)
print('Run time is ',(time.time()-start),' seconds' )