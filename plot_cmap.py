import numpy as np 
import matplotlib.pyplot as plt 
import sys 

axd = plt.figure(figsize=[3,3],dpi=300).subplot_mosaic(
        """
        bAd
        .c.
        """,
        gridspec_kw = {
            "height_ratios": [3,1],
            "width_ratios" : [1,3,0.1],
            },
        )

data = np.load(sys.argv[1])
title = sys.argv[2]
filename = sys.argv[3]
tmean = data.mean(axis=0)
tickskip = 15
mainticklocs = range(0,tmean.shape[0],tickskip)
sideticklocs = np.arange(0.5,tmean.shape[0]+0.5,tickskip)
sideticks = range(1,tmean.shape[0]+1,tickskip)
im = axd['A'].imshow(tmean, aspect='auto',origin='lower',cmap='gist_heat_r',vmin=0,vmax=np.max(tmean))
axd['A'].set_xticks(mainticklocs)
axd['A'].set_yticks(mainticklocs)
axd['A'].set_xticklabels([])
axd['A'].set_yticklabels([])
axd['b'].plot(tmean.mean(axis=0),range(tmean.shape[1]),ls='--',lw=0.5,color='black')
axd['c'].plot(range(tmean.shape[0]),tmean.mean(axis=1),ls='--',lw=0.5,color='black')
axd['b'].set_yticks(sideticklocs)
axd['b'].set_yticklabels(sideticks)
axd['c'].set_xticks(sideticklocs)
axd['c'].set_xticklabels(sideticks)
axd['c'].set_xlim(0,tmean.shape[0])
axd['b'].set_ylim(0,tmean.shape[0])
axd['b'].set_ylabel('Residue no.',fontsize=7)
axd['c'].set_xlabel('Residue no.',fontsize=7)
axd['b'].xaxis.tick_top()
cb = plt.colorbar(im,cax=axd['d'])
cb.ax.tick_params(length=1,pad=0.5,labelsize=5)
cb.ax.set_ylabel('Contacts per chain per frame',fontsize=6)
for i in axd:
    axd[i].tick_params(length=2,pad=0.5,labelsize=4)
axd['A'].tick_params(length=3,pad=0.0,labelsize=0)
plt.subplots_adjust(wspace=0.1, hspace=0.1)
axd['A'].set_title(title,fontsize=7)
plt.savefig('ctmap_%s.png'%filename,dpi=300)
plt.show()
