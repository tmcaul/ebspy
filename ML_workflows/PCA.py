#%%
import ebspy as eb
import os as o
import h5py as h5#
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

patfolder=r'E:\Tom\GammaPrime_Data\V208C'
patfile=r'yprime3.h5'

#%% Make a saving directory and move to it
patpath=Path(patfolder)
ebspfilepath=patpath/patfile

savingdir=patfile[:-3]+'_analysis'
o.chdir(patfolder)
savingpath=Path(patfolder)/Path(savingdir)
if o.path.isdir(savingdir):
    o.chdir(savingpath)
    
else:
    o.mkdir(savingdir)
    o.chdir(savingpath)

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
fig.savefig('Radon.tiff',dpi=400)


   
#%% Convert to EBSP type
fulldataset=Data['Patterns'][:,:,:] #n, patrow, patcol
EBSPs=eb.EBSP(fulldataset,Metadata)
EBSPs.disp(0)

#%% Perform background correction
EBSPsCor=EBSPs.bgcorr(resize=True,screensize=[300,300])
EBSPsCor.disp(0)

EBSPs=[]
fulldataset=[]
Data=[]
    
#%% Generate spatially weighted testarray
for n in range(3,4):
    
    o.chdir(savingpath)
   
    if o.path.isdir(r'Local'+str(n))==False:
        o.mkdir(r'Local'+str(n))
    
    o.chdir(r'Local'+str(n))

    tA_h=eb.spatial_tA(EBSPsCor.vector,region,n)
    
    #%% Perform PCA
    #testArray=EBSPsCor.vector
    
    from sklearn.decomposition import PCA
    
    X=tA_h
    pca = PCA(n_components=5,whiten=True) #set up PCA object
    pcafit = pca.fit(X.T) #apply pca object to testArray
    scores=pcafit.transform(X.T) #calculate scores
    totalvariance=np.var(X)

    #%
    ncomps=pcafit.n_components_
    coeffs=np.array(pcafit.components_)
    
    #% Transfer variance to scores (for whiten = True)
    # Same result as whiten = False but a bit more transparent as to what is actually going on...
    for i in range(0,scores.shape[1]):
        scores[:,i]=scores[:,i]*pcafit.singular_values_[i]

    #nb: singular values are the lengths (2-norms) of the coefficient vectors. 
    # eigenvalues of the covariance matrix are the explained variances. lambdas = sigmas^2 / (n_objects-1)
    # so we 'transfer variance' by multiplying scores by the lengths of the corresponding coefficient vectors

    #np.linalg.norm(coeffs[0,:]) #=1 for all vectors before variance transfer
    fig=plt.figure()
    plt.plot(range(0,len(pcafit.explained_variance_)),pcafit.explained_variance_,linewidth=1.5,marker='o',color='r')
    plt.xticks([0,1,2,3,4],[1,2,3,4,5])
    plt.xlabel('Component')
    plt.ylabel('Explained variance')
    plt.show()
    fig.savefig('Scree.tiff',dpi=400)
    
    tA_h=None
    #EBSPsCor=None
    X=None
    
    #%%reshape score maps and get label maps
    #Also invert if average coeff < 0
    scoresreshaped=np.zeros((ncomps,region[1]-region[0],region[3]-region[2]))
    for i in range(0,ncomps):
        scoresreshaped[i,:,:]=scores[:,i].reshape((region[1]-region[0],region[3]-region[2]))
        if i==0 or i==1 or i==4 or (i==3 and n!=1): #manual flipping
            scoresreshaped[i,:,:]=-scoresreshaped[i,:,:]
            coeffs[i,:]=-coeffs[i,:]

    coeffEBSPs=eb.EBSP(coeffs,EBSPsCor.Metadata,vector=False)
    coeffEBSPs.disp(0)

    
    labels=np.zeros((region[1]-region[0],region[3]-region[2]))
    for i in range(0,int((region[1]-region[0]))): #loop over rows
        for j in range(0,int((region[3]-region[2]))): #loop over cols
            labels[i,j]=np.argmax(scoresreshaped[:,i,j])
            
    plt.imshow(labels)
    plt.show()
    
    #%% Print scores and coefficients
    for n1 in range(0,ncomps):
        fig=plt.figure()
        scoremap=scoresreshaped[n1,:,:]#+scoresreshaped[n2,:,:]
        pattern=coeffEBSPs[n1].array#+coeffEBSPs[2].array
        ev=pcafit.singular_values_[n1]

        plt.imshow(scoremap,cmap='magma',vmin=-3*ev,vmax=3*ev)
        plt.colorbar(fraction=0.075)
        plt.axis('off')
        plt.show()
        fig.savefig('Scoremap'+str(n1)+'.tiff',dpi=400)
    
        fig=plt.figure()
        plt.imshow(pattern[0,:,:],cmap='viridis',vmin=-0.015,vmax=0.015)
        plt.colorbar(fraction=0.075)
        plt.axis('off')
        plt.show()
        fig.savefig('Pattern'+str(n1)+'.tiff',dpi=400)
    
    #gp point = 20,32
    #g point = 30,28

    #%% Save scores, coefficients, eigenvalues and singular values
    hf=h5.File('scorescoeffs.h5','w')
    hf.create_dataset('scores',data=scoresreshaped)
    hf.create_dataset('coeffs',data=coeffEBSPs.array)
    hf.create_dataset('eigenvalues',data=pcafit.explained_variance_)
    hf.create_dataset('total_var',data=totalvariance)
    #hf.create_dataset('singularvalues',data=pcafit.singular_values_) can calculate from eigenvalues
    hf.close()   
    
    #%% Get total scores
    totalscore=np.sum(scoresreshaped[0:2,:,:],axis=0)
    normscores=np.linalg.norm(scoresreshaped,axis=0)
    
    fig=plt.figure()
    plt.imshow(totalscore,cmap='magma')
    plt.colorbar(fraction=0.075)
    plt.axis('off')
    plt.show()
    fig.savefig('Totalscore.tiff',dpi=400)
    
    fig=plt.figure()
    plt.imshow(normscores,cmap='magma')
    plt.colorbar(fraction=0.075)
    plt.axis('off')
    plt.show()
    fig.savefig('Normscores.tiff',dpi=400)
        
    #o.chdir(savingpath)

    #%% Reconstruct the test array
    s=scoresreshaped.reshape(ncomps,EBSPsCor.number)
    reconstruction = coeffs.T @ s

    s=[]
    #scoresreshaped=[]
    coeffs=[]
    #%%
    plt.imshow(reconstruction[:,2].reshape(EBSPsCor.patheight,EBSPsCor.patheight),cmap='gray')

    #%% Renormalise variance of reconstructed pattern for Xcorrelation
    for i in range(0,EBSPsCor.number):
        std=reconstruction[:,i].std()
        reconstruction[:,i]=reconstruction[:,i]/std
        std=[]
        
    reconEBSPs=eb.EBSP(reconstruction.T,EBSPsCor.Metadata)
    reconstruction=[]

    #%% Export HDF5 file
    # hf=h5.File('reconstruction.h5','w')
    # hf.create_dataset('data',data=reconEBSPs.array)
    # hf.close()


# %% Generate labels and save
ts=np.sum(scoresreshaped[0:3,:,:],axis=0)
labels2=np.zeros_like(labels)
labels2[ts<np.quantile(ts,0.38)]=1

#hist,edges=np.histogram(s.flatten(),10)

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4))
ax1.imshow(labels)
ax2.imshow(ts)

LabelPats=np.zeros((2,300,300))
a=EBSPsCor.array.reshape(72,95,300,300)


#%%
# gamma THEN gamma prime
LabelPats[1,:,:]=np.average(EBSPsCor.array[labels2.reshape(6840)==1,:,:],axis=0)
LabelPats[0,:,:]=np.average(EBSPsCor.array[labels2.reshape(6840)==0,:,:],axis=0)

hf=h5.File('LabelPatterns.h5','w')
hf.create_dataset('PCA',data=LabelPats)
hf.create_dataset('PCA_scorematrix',data=ts)
hf.create_dataset('PCA_region',data=labels2)
hf.close()

# # %%
# labelstest=np.zeros((1,300,300))
# a=EBSPsCor.array.reshape(72,95,300,300)
# for i in range(0,72):
#     for j in range(0,95):
#         if labels2[i,j]==1:
#             labelstest=np.concatenate((labelstest,a[i,j,:,:].reshape(1,300,300)),axis=0)

# %%
