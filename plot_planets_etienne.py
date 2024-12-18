import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse_single_table
import matplotlib.image as image
import matplotlib.colors
from mw_plot import MWSkyMap
from mw_plot import MWPlot
from astropy import units as u
#plt.style.use('dark_background')
import astropy.coordinates as apycoords

plt.rc('font', family='serif')
plt.rc('xtick', labelsize=25)
plt.rc('ytick', labelsize=25)

plot_instance = MWPlot(radius=12 * u.kpc, unit=u.kpc, coord='galactic')
plot_instance.title = None  # plot title, or it can be None to show no title
plot_instance.fontsize = 35  # fontsize for matplotlib plotting
plot_instance.figsize = (20, 20)  # figsize for matplotlib plotting
plot_instance.dpi = 100  # dpi for matplotlib plotting
plot_instance.cmap = 'plasma'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
plot_instance.clim = (0,10000) # colorbar range
plot_instance.imalpha = 0.5  # alpha value for the milkyway image
plot_instance.s = 100.0  # make the scatter points bigger
plot_instance.tight_layout = True # whether plt.tight_layout() will be run
plot_instance.initialize_mwplot()






table = parse_single_table("PS_2023.09.22_15.59.22.votable")
mm = np.unique(table.array['pl_name'],return_index=True)
table.array = table.array[mm[1]]

breakpoint()

#im = image.imread('earth.jpeg')
#plt.imshow(im, aspect='auto', extent=(0.5/2.7, 1.5/2.7, 0.5, 1.5), zorder=-1)

#im = image.imread('jupit.jpeg')
#plt.imshow(im, aspect='auto', extent=(1/2.7, 10.0/2.7, 0.2*317.6, 1.8*317.6), zorder=-1)
import pdb; pdb.set_trace()
methods = table.array['discoverymethod'].data



unique_method2 = np.unique( table.array['discoverymethod'].data,return_counts=True) 
unique_method = [b'Transit',b'Radial Velocity',b'Microlensing',b'Imaging',]
#unique_method = ['Transit','Radial Velocity',]
colors = np.zeros((len(methods))).astype(str)

cco = ['seagreen','gold','peru','black','grey','lime','grey','k']
shape = ['o','s','p','^']
fig,ax = plt.subplots(1,2)

ax[0].axvline(1,lw=30,c='dodgerblue',alpha=0.5)
ax[0].text(0.8,0.5, 'SNOW LINE', color='dodgerblue',fontsize=25,rotation=90)

ax[1].imshow(plot_instance._img, zorder=0, extent=plot_instance._ext,  rasterized=True,alpha=0.75)
ratio = []
for index,met in  enumerate(unique_method):

	mask = methods == met.decode("utf-8")
	colors[mask] = cco[index]

	ax[0].scatter(table.array['pl_orbsmax'].data[mask]/(2.7*table.array['st_mass'].data[mask]),table.array['pl_bmassj'].data[mask]*317.6,facecolors=cco[index],edgecolors='none',s=100,marker = shape[index],alpha=0.75,label = met.decode("utf-8")+' '+str(len(colors[mask])))

	mask = (methods == met.decode("utf-8")) & (~np.isnan(table.array['sy_dist'].data))	
	c_dr1 = apycoords.SkyCoord(table.array['ra'].data[mask]*u.deg,table.array['dec'].data[mask]*u.deg, distance=table.array['sy_dist'].data[mask]*u.pc, frame='icrs')
	#import pdb; pdb.set_trace()
	colors[mask] = cco[index]
	#plot_instance.ax.scatter(-c_dr1.galactic.cartesian.x[mask], c_dr1.galactic.cartesian.y[mask],c=cco[index],marker=shape[index],s=1000)
	ax[1].scatter(-c_dr1.galactic.cartesian.x/1000, c_dr1.galactic.cartesian.y/1000,facecolors=cco[index],marker=shape[index],s=100,alpha=1,)
	if b'Microlensing' in met:
		ratio.append(table.array['pl_bmassj'].data[mask]*317.6/table.array['st_mass'].data[mask])
	#plt.scatter(table.array['pl_orbsmax'].data[mask],table.array['pl_bmassj'].data[mask]*317.6,facecolors=cco[index],edgecolors='none',s=50,marker = shape[index],alpha=1.0)
	#plt.subplot(122)
	#plt.scatter(table.array['sy_dist'].data[mask],table.array['pl_bmassj'].data[mask]*317.6,facecolors=cco[index],edgecolors='none',s=50,marker = shape[index],alpha=1.0,label = met.decode("utf-8")+' '+str(len(colors[mask])))

	#plt.scatter(table.array['pl_orbper'].data[mask][0],table.array['pl_bmassj'].data[mask][0],facecolors=cco[index],edgecolors='none',label = met,s=150)
	#plt.scatter(table.array['ra'].data[mask],table.array['dec'].data[mask],c=5*np.log10(table.array['sy_dist'].data[mask])-5,edgecolors='none',s=50,marker = shape[index],alpha=1.0)

mask = ~np.isnan(table.array['sy_dist'].data)
#sm=plt.scatter(table.array['ra'].data[mask],table.array['dec'].data[mask],c=5*np.log10(table.array['sy_dist'].data[mask])-5,edgecolors='none',s=50,marker = shape[index],alpha=1.0)#
#plt.colorbar(sm)

ax[0].set_xlabel(r'$a [a_{snow}]$',fontsize=35)
ax[0].set_ylabel(r'$M_{planet} [M_\oplus]$',fontsize=35)
ax[0].tick_params(axis='both', pad=10)

ax[0].set_xlim([8*10**-3,10**3])
ax[0].set_ylim([10**-1,10**4])
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].legend(bbox_to_anchor=(-0.2, 0.9, 3, .102), loc=3, ncol=2, mode="expand", borderaxespad=2,fontsize=35,markerscale=3,frameon=False)
#ax[0].legend()
#plt.legend(scatterpoints=1,loc=4,fontsize=15)
#plt.yscale('log')
#plt.xscale('log')
ax[1].set_xlabel(r'$X [kpc]$',fontsize=35)
ax[1].set_ylabel(r'$Y [kpc]$',fontsize=35)
ax[1].tick_params(axis='both', pad=10)
#plt.xticks(fontsize=15)
#plt.yticks(fontsize=15)



plt.show()
import pdb; pdb.set_trace()

