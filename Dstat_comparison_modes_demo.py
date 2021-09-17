#--------------------------
#    Python script to produce Figure 11 in 
#    "Statistical characteristics of extreme daily precipitation during 1501BCE – 1849CE in the Community Earth System Model", Kim et al. (2021)
#    It compares D statistics (Coles et al, 2001) among the modes of variability, TS-R, and TS-G GPD models. 
#    Written by Woon Mi Kim in September 2021. 
#    
#--------------------------

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl


#-----------------------
#   Load variables 
#-----------------------
#   Transient (Full forcing) simulation
#-----------------------

nllh_trans=xr.open_dataset('CESM122.transient.log-likelihood-GPDmodel-ModesVar.nc')
nllh1=nllh_trans.nllh.values

lat1=nllh_trans.lat.values
lon1=nllh_trans.lon.values

nlat=len(lat1)
nlon=len(lon1)

stat_matrix_tr=nllh1[0,:,:]
modes_matrix_tr=nllh1[1:,:,:]

#-----------------------
#   Orbital-only simulation
#-----------------------

nllh_orb=xr.open_dataset('CESM122.orbital.log-likelihood-GPDmodel-ModesVar.nc')
nllh2=nllh_orb.nllh.values

stat_matrix_ob=nllh2[0,:,:]
modes_matrix_ob=nllh2[1:,:,:]

#---- Precipitation mask file
maskfile=xr.open_dataset('mask_prcp_desert.nc')
mask=np.squeeze(maskfile.PRECT.values)

#-----------------------
#   Functions
#-----------------------

#----- masking
def masking(var_matrix,mask_matrix):
    nlat=mask_matrix.shape[0]
    nlon=mask_matrix.shape[1]
    var_new=var_matrix.copy()
    for i in range(nlat):
        for j in range(nlon):
            if mask_matrix[i,j]==0.:
                var_new[i,j]=np.float('nan')
    return(var_new)

#----- Function that calculates D statistics 
def Dstatistic(l1,l0):
    dval=2.*(-l1-(-l0))
    return(dval)

#---- Function that compares the D statistics of all models. 

modes_names=['EA-WR','NAO','NAM','PDO','PNA','ENSO','SAM','PSA','TS-R','TS-G']

def comparison_D(matrix,threshold=99,ttype='modes'):  
    l=matrix.shape[0]
    #---- with TS-R 
    if ttype=='ts-r':  
        nhpat=[0,1,2,3,4,5,8]
        shpat=[5,6,7,8]
    #---- with TS-G
    elif ttype=='ts-g':  
        nhpat=[0,1,2,3,4,5,9]
        shpat=[5,6,7,9]
    #---- without TS. Only modes
    elif ttype=='modes':
        nhpat=[0,1,2,3,4,5]
        shpat=[5,6,7]
    #---- Choose NLLH for Southern Hem and Northern Hem. 
    msh=matrix[shpat,:,:]
    mnh=matrix[nhpat,:,:]
    nlat=matrix.shape[1]
    nlon=matrix.shape[2]
    if threshold==99:
        Dval=6.634
    elif threshold==95.:
        Dval=3.841
    #---- Comparison starts. Choose the maximum at each grid point. 
    lat_sh=range(0,48)
    lat_nh=range(48,96)
    com=np.zeros((nlat,nlon))
    com[:]=np.float('nan')
    val=np.zeros((l,nlat,nlon))
    val[:]=np.float('nan')
    #---- sh : 
    for m in range(len(shpat)):
        for i in range(0,48):
            for j in range(nlon):
                if msh[m,i,j]>=Dval:
                    val[m,i,j]=msh[m,i,j]
    #---- nh : 
    for m in range(len(nhpat)):
        for i in range(48,96):
            for j in range(nlon):
                if mnh[m,i,j]>=Dval:
                    val[m,i,j]=mnh[m,i,j]
    #------ Replace maximum values. 
    for j in range(nlon):
        for i in range(0,48):
            if np.nansum(val[:,i,j])!=0:              
                I=int(np.nanargmax(val[:,i,j]))
                com[i,j]=shpat[I]+0.5
        #-----------------
        for i in range(48,96):
            if np.nansum(val[:,i,j])!=0:              
                I=int(np.nanargmax(val[:,i,j]))
                com[i,j]=nhpat[I]+0.5
    return(com)
    
    
#---- Function that chooses the regions where the full-forcing and orbital-only are coherent.

def coherence_region(full_value,orb_value):   
    coh_matrix=np.zeros((nlat,nlon))
    coh_matrix[:]=np.float('nan')
    for i in range(nlat):
        for j in range(nlon):
            if (np.isnan(full_value[i,j])==False and np.isnan(orb_value[i,j])==False):
                if full_value[i,j]==orb_value[i,j]:
                    coh_matrix[i,j]=full_value[i,j]
                else:
                    coh_matrix[i,j]=np.float('nan')
    return(coh_matrix)


#---- Plotting function

def plotting_mesh(matrix, minval=1.,maxval=2.,colorlevel=1):
    air_cyclic, lons_cyclic = addcyclic(matrix, lon1)
    air_cyclic, lons_cyclic = shiftgrid(180., air_cyclic, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat1)
    map = Basemap(projection='robin',lon_0=0.,resolution='c')
    map.drawcoastlines(linewidth=0.5)
    map.drawmapboundary()
    x, y = map(lon2d, lat2d)
    if colorlevel==1:
        colors1=[(145./255.,39./255.,180./255.),(50/255.,64./255.,255/255.), (43./255.,55./255.,193./255.),\
                 (150/255.,205/255.,245/255.),(116/255.,245/255.,169/255.),\
                 (213/255.,251/255.,92/255.), (242/255.,189/255.,65/255.),(230/255.,49/255.,34/255.)]
    elif colorlevel==2:
        colors1=[(50/255.,64./255.,255/255.), (43./255.,55./255.,193./255.),\
                (150/255.,205/255.,245/255.),(116/255.,245/255.,169/255.),\
                (213/255.,251/255.,92/255.), (242/255.,189/255.,65/255.),(230/255.,49/255.,34/255.)]
    cmap_name='my_list'
    bounds=np.arange(minval,maxval)
    n=len(bounds)
    cmap=LinearSegmentedColormap.from_list(cmap_name,colors1[::-1],N=(n))
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cs = map.pcolormesh(x, y, air_cyclic, rasterized=True, vmin=minval, vmax=maxval, cmap=cmap, norm=norm)
    return(cs)


#------------------------

Dtr=Dstatistic(modes_matrix_tr,stat_matrix_tr)
Dob=Dstatistic(modes_matrix_ob,stat_matrix_ob)

#---- only modes
max_Dtr=comparison_D(Dtr)
max_Dob=comparison_D(Dob)

#--- with TS-G
max_Dtr_tg=comparison_D(Dtr,ttype='ts-g')
max_Dob_tg=comparison_D(Dob,ttype='ts-g')

#--- with TS-R
max_Dtr_tr=comparison_D(Dtr,ttype='ts-r')
max_Dob_tr=comparison_D(Dob,ttype='ts-r')

#--- coherence regions
coh_all=coherence_region(max_Dtr, max_Dob)
coh_all_tg=coherence_region(max_Dtr_tg, max_Dob_tg)
coh_all_tr=coherence_region(max_Dtr_tr, max_Dob_tr)


#-------------------------
#    Plotting all
#-------------------------

tfont = {'fontsize':6}
my_dpi=100.

#-------------------------
#    1) Each int/ Ts-G / Ts-r
#-------------------------

fig = plt.figure(figsize=(750/my_dpi,150/my_dpi), dpi=my_dpi*3)

#--- Modes
ax0 =plt.subplot2grid((1,3), (0,0))
ax0.set_title('(a)   Modes of variability',**tfont)

cs0=plotting_mesh(masking(coh_all,mask),0,len(modes_names)-1,2)

#--- TS-G
ax1 =plt.subplot2grid((1,3), (0,1))
ax1.set_title('(b)   with TS-G',**tfont)
cs1=plotting_mesh(masking(coh_all_tg,mask),0,len(modes_names),1)

#--- TS-R
ax2 =plt.subplot2grid((1,3), (0,2))
ax2.set_title('(c)   with TS-R',**tfont)
cs2=plotting_mesh(masking(coh_all_tr,mask),0,len(modes_names),1)


fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.05, hspace=0.05)


#--- colorbar 

cb_ax1= fig.add_axes([0.27, 0.13, 0.36, 0.03])
cb1 = fig.colorbar(cs1,cax=cb_ax1,orientation='horizontal')
loc=np.arange(0.5,len(modes_names))
cb1.set_ticks(loc)
cb1.set_ticklabels(modes_names[0:-2]+['TS-G or\n TS-R'])
cb1.ax.tick_params(labelsize=5,length=1)

plt.show()
