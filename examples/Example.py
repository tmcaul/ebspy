#%%
#Example script for exhibiting functionality of EBSP class
#TPM 12/2019 (with BGcorrection derived from AstroEBSD (https://github.com/benjaminbritton/AstroEBSD))

import h5py as h5#
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import ebspy as eb

# Set the directory which the data is stored in
patfolder=r'C:\Users\tpm416\Documents\Data\LR22ht3c+2a'
patfile='yprime1.h5'

patpath=Path(patfolder)
ebspfilepath=patpath/patfile

#%% Import data from h5
#region=[0,60,0,60]
region=[]
Metadata, Data, visual, region = eb.bImportEBSPs(ebspfilepath, True, region)
plt.imshow(visual)
plt.show()

#print radon map
Radon=Metadata['RadonQuality']
Radon=Radon.reshape((Metadata['NROWS'],Metadata['NCOLS']))

fig=plt.figure()
plt.imshow(Radon)#,vmin=0.51,vmax=0.56)
plt.colorbar(fraction=0.075)
plt.axis('off')
plt.show()
#fig.savefig('Radon.tiff',dpi=400)

#%% Convert patterns to EBSP class
fulldataset=Data['Patterns'][:,:,:] #n, patrow, patcol
EBSPs=eb.EBSP(fulldataset,Metadata)
EBSPs.disp(0)

#%% Perform background correction
EBSPsCor=EBSPs.bgcorr(resize=True,normalise=True,screensize=[300,300])
EBSPsCor.disp(0)

#%% Functionality for playing with arrays of EBSPs
EBSPs1=EBSP(fulldataset[1,:,:],Metadata)
EBSPs2=EBSP(fulldataset[2:4,:,:],Metadata)
EBSPs3=EBSP(fulldataset[22:25,:,:],Metadata)

len(EBSPs1) #gives the number in the stack

#can extract by index
EBSPs4=EBSPs3[0:2]

#can replace by index
EBSPs3[1:3]=EBSPs2[0:2]
EBSPs3[0,1]=EBSPs2[0,1]
