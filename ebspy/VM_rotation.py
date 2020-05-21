# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:02:54 2019

@author: tpm416
"""

#%% VARIMAX rotation
#VM rotation not super applicable to the y/y' case due to near-identity of PC patterns
topPCs=np.array(coeffs)
R=ortho_rotation(topPCs.T)
#R.shape #should be (n,n)
#%
rotarrays=topPCs.T@R
rotPCs=eb.EBSP((topPCs.T@R).T,EBSPsCor.Metadata,vector=True)
rotScores=np.array(scores)@R

#%% Plot VM rotated scores
n=0
rotscoremap=rotScores[:,n].reshape((region[1]-region[0],region[3]-region[2]))
plt.imshow(rotscoremap,cmap='inferno_r')
plt.show()

#%%
n=0
meanEBSP=EBSPsCor.array.mean(0)
coeffs_uc=np.array(coeffs[n,:]).reshape(300,300)+meanEBSP

if coeffs_uc.mean()<0:
    coeffs_uc=coeffs_uc*-1

plt.imshow(coeffs_uc,cmap='inferno')
plt.show()