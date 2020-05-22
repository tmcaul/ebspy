#Create spatially weihgted testArray after Guo et al 
#TPM 14/10/19

def spatial_tA(testArray,region,r=2,normalise=True,normalise_NMF=False):
    
    import numpy as np
    
    #tA_h=np.array(testArray).reshape(int(testArray.shape[0]),int(testArray.shape[1]**0.5),int(testArray.shape[1]**0.5))
    tA_h=np.array(testArray).reshape(int(testArray.shape[0]),int(region[1]-region[0]),int(region[3]-region[2]))
    
    tA_h_original=tA_h
    
    for i in range(0,(region[1]-region[0])):
        for j in range(0,(region[3]-region[2])):
            
            
            d=np.zeros((region[1]-region[0],region[3]-region[2]))
            for itest in range(0,(region[1]-region[0])):
                for jtest in range(region[3]-region[2]):
                    
                    d[itest,jtest]=((itest-i)**2+(jtest-j)**2)**0.5 #get euclidean distances
                    
            kernel=np.nonzero(d<=r)
            kernel_i=kernel[0]
            kernel_j=kernel[1]
            
            w=np.zeros(len(kernel_i))
            for l in range(0,len(kernel_i)):
                w[l]=(1-(d[kernel_i[l],kernel_j[l]]/r)**2)**2
            w=w/sum(w) #normalise
            
            sumpat=np.zeros(tA_h.shape[0]) #reset the sum pattern
            for l in range(0,len(kernel_i)):
                sumpat=sumpat+w[l]*tA_h_original[:,kernel_i[l],kernel_j[l]]
            tA_h[:,i,j]=sumpat #insert sum pattern

            if normalise==True:
                if normalise_NMF==False:
                    #renormalise pattern
                    tA_h[:,i,j]=(tA_h[:,i,j]-np.mean(tA_h[:,i,j]))/np.std(tA_h[:,i,j])
            
            if normalise_NMF==True:
                tA_h[:,i,j]=(tA_h[:,i,j]-np.amin(tA_h[:,i,j]))/np.std(tA_h[:,i,j])
            
    tA_h=tA_h.reshape(testArray.shape[0],testArray.shape[1])
            
    return tA_h