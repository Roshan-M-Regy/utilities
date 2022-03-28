# Roshan Mammen Regy
# roshanm.regy@tamu.edu
import numpy as np 
import matplotlib.pyplot as plt 
plt.rcParams.update({'text.usetex':True,"font.family":'serif'})
import gif 


@gif.frame
def plot(y):
    # LAF-1 
    markersize=4
    lw=0.5
    mlw=1.0
    mec='black'
    mew=0.2
    lim=0
    fig,ax = plt.subplots(figsize=[5,5],dpi=300)
    ax=[]
    ax.append(plt.subplot2grid((2,2),(0,0),colspan=1,rowspan=1))
    ax.append(plt.subplot2grid((2,2),(0,1),colspan=1,rowspan=1))
    ax.append(plt.subplot2grid((2,2),(1,0),colspan=1,rowspan=1))
    ax.append(plt.subplot2grid((2,2),(1,1),colspan=1,rowspan=1))
    ax[0].set_ylim(0,600)
    ax[1].set_ylim(0,90)
    ax[2].set_ylim(0,600)
    ax[3].set_ylim(0,100)
    for i in ax:
        i.tick_params(labelsize=12)
    # data = np.loadtxt containing density profile of a frame
    ax[0].plot(data[:,0],data[:,y],'--o',color='lightblue',markersize=markersize,lw=lw,mec=mec,mew=mew)
    ax[0].plot(data[:,0],data[:,-1],'--',color='black',markersize=markersize,lw=mlw)
    # data = np.loadtxt containing density profile of a frame
    ax[1].plot(data[:,0],data[:,y],'--^',color='lightblue',markerfacecolor='red',markersize=markersize,lw=lw,mec=mec,mew=mew)
    ax[1].plot(data[:,0],data[:,-1],'--',color='black',markersize=markersize,lw=mlw)
    # data = np.loadtxt containing density profile of a frame
    ax[2].plot(data[:,0],data[:,y],'--o',color='lightgreen',markersize=markersize,lw=lw,mec=mec,mew=mew)
    ax[2].plot(data[:,0],data[:,-1],'--',color='black',markersize=markersize,lw=mlw)
    # data = np.loadtxt containing density profile of a frame
    ax[3].plot(data[:,0],data[:,y],'--^',color='lightgreen',markerfacecolor='red',markersize=markersize,lw=lw,mec=mec,mew=mew)
    ax[3].plot(data[:,0],data[:,-1],'--',color='black',markersize=markersize,lw=mlw)
    ax[0].set_ylabel('Ylabel',fontsize=14)
    ax[1].set_ylabel('Ylabel',fontsize=14)
    ax[2].set_ylabel('Ylabel',fontsize=14)
    ax[3].set_ylabel('Ylabel',fontsize=14)
    ax[2].set_xlabel('Xlabel',fontsize=14)
    ax[3].set_xlabel('Xlabel',fontsize=14)
    ax[0].set_ylim(-10,600)
    ax[1].set_ylim(-10,100)
    ax[2].set_ylim(-10,600)
    ax[3].set_ylim(-10,100)
    ax[1].annotate(r'%.2f $\mu$s'%(y*0.05),(13,80)) # update with the time step of the frame being plotted
    xx = -0.3
    yy=1.1
    ax[0].text(xx,yy,'A',fontsize=12,transform = ax[0].transAxes,fontweight='bold')
    ax[1].text(xx,yy,'B',fontsize=12,transform = ax[1].transAxes,fontweight='bold')
    ax[2].text(xx,yy+.1,'C',fontsize=12,transform = ax[2].transAxes,fontweight='bold')
    ax[3].text(xx,yy+.1,'D',fontsize=12,transform = ax[3].transAxes,fontweight='bold')
    for i in range(4):
        ax[i].set_xlim(-30,30)
    plt.tight_layout()
frames=[]    

for k,number in enumerate(range(1,101)):
    frame=plot(number)
    frames.append(frame)
print (len(frames))
gif.save(frames, 'density_profiles_movie_mean_front_lw1_nomarkers_newlabels_shiftCD.gif', duration=40, unit="s", between="startend")
