#%% Classes
class EBSP: #feed it n x patheight x patwidth OR n x (patheight*patwidth)
    
    def __init__(self,inputarray1,Metadata,rect=False,vector=False):
        
        import numpy as np
        
        inputarray=np.array(inputarray1) #create a COPY rather than change the thing itself
        
        #attempt to identify if the input is in column or array form
        if (Metadata['PatternHeight']*Metadata['PatternWidth'] in inputarray.shape) or vector==True:
            if inputarray.ndim==1: #single column vector
                array=inputarray.reshape(int(Metadata['PatternHeight']),int(Metadata['PatternWidth']))[:]
            if inputarray.ndim>1: #if it is then find non-pattern index and reshape into array
                ind=inputarray.shape[0]
                array2=inputarray.reshape(ind,int(Metadata['PatternHeight']),int(Metadata['PatternWidth']))[:]
                array=array2[:]
        
        else:
            array=inputarray[:]
            
        #otherwise assume been inputted in shape
        self.array=array[:]
        
        #define some generally useful properties
        self.Metadata=Metadata
        if self.array.ndim==3:
            self.number=array.shape[0]
            self.patheight=array.shape[1]
            self.Metadata['PatternHeight']=self.patheight
            self.patwidth=array.shape[2]
            self.Metadata['PatternWidth']=self.patwidth
            
        if self.array.ndim==2:
            self.number=1
            self.patheight=array.shape[0]
            self.patwidth=array.shape[1]
        
        self.length=self.patheight*self.patwidth
        
        self.vector=self.array.reshape((self.number,self.length)).T
        self.mean=self.vector.mean(axis=0)
        self.std=self.vector.std(axis=0)
        self.norm=np.linalg.norm(self.vector,axis=0)
        
        self.shape=self.array.shape

#%%      
    def __len__(self):
        return(self.number)

#%%    
    # how to deal with ADDING EBSPs (ie. stack them)
    def __add__(self,obj2):
        
        import numpy as np
        
        #get numpy arrays of the objects
        arr1=self.array
        arr2=obj2.array
        
        #reshape to eg. by (1,ph,pw) rather than (ph,pw) if only one
        n_1=self.number
        if n_1==1:
            arr1=arr1.reshape(1,self.patheight,self.patwidth)
        
        #reshape to eg. by (1,ph,pw) rather than (ph,pw) if only one
        n_2=obj2.number
        if n_2==1:
            arr2=arr2.reshape(1,self.patheight,self.patwidth)
        
        #preallocate then overwrite an array of ones
        arr3=np.ones((n_1+n_2,self.patheight,self.patwidth))
        for n in range(0,n_1+n_2):
            if n<n_1:
                arr3[n,:,:]=arr1[n,:,:]
            if n>n_1:
                arr3[n,:,:]=arr2[n-n_1,:,:]
        
        newEBSP=EBSP(arr3,self.Metadata)        
        return newEBSP
    
    #%% grab with EBSP[a:b], for example
    def __getitem__(self,nslice):
        import numpy as np
        
        if type(nslice)==int:
            n=[nslice]
        elif type(nslice)==tuple:
            n=list(nslice)
        else:#it's a slice
            listnums=[]
            for i in range(0,self.number):
                listnums.append(i)
            
            indices=listnums[nslice]     
            n=indices
        
        outputarray=np.ones((len(n),self.patheight,self.patwidth))
        
        for i in range(0,len(n)):
        
            number=n[i]
            
            if number+1>self.number: #ie try to index outside the range of this group of EBSPs
                raise Exception ('EBSP not big enough for this index')
            
            if self.number>1:
                newarray=np.array(self.array[number,:,:])
                newarray=newarray.reshape(1,self.patheight,self.patwidth)
            else:
                newarray=np.array(self.array)
                newarray=newarray.reshape(1,self.patheight,self.patwidth) #so will raise an error if 
                
            outputarray[i,:,:]=newarray
        
        return EBSP(outputarray,self.Metadata)
    
    
    
    #%% change with EBSP[a]=EBSP2[b], for example
    def __setitem__(self,nslice,newEBSP):
        import numpy as np
        
        if type(nslice)==int:
            n=[nslice]
        elif type(nslice)==tuple:
            n=list(nslice)
        else:#it's a slice
            listnums=[]
            for i in range(0,self.number):
                listnums.append(i)
            
            indices=listnums[nslice]     
            n=indices
                
        if newEBSP.number < len(n):
            raise Exception ('EBSP set not big enough to insert here')
        
        if len(n) != newEBSP.number:
            raise Exception ('Wrong exchange size')
        
        outputarray=np.ones((self.number,self.patheight,self.patwidth))
        
        selfarray=np.array(self.array.reshape(self.number,self.patheight,self.patwidth))
        newEBSParray=np.array(newEBSP.array.reshape(newEBSP.number,newEBSP.patheight,newEBSP.patwidth))
        
        for i in range(0,self.number): 
            if i in n:
                outputarray[i,:,:]=newEBSParray[(i-min(n)),:,:]
            else:
                outputarray[i,:,:]=selfarray[i,:,:]
        
        self.array=outputarray
        self.vector=self.array.reshape((self.number,self.length)).T
        self.mean=self.vector.mean(axis=0)
        self=EBSP(outputarray,self.Metadata)

#%%    
    # a method to quickly print the EBSP
    def disp(self,n):
        if n==None:
            n=self.number
            
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots()
        plt.imshow(self.array[n,:,:],cmap='gray')
        ax.set_aspect(1.0)

#%%    
    #Background correction, based off of AstroEBSD bgcorr
    def bgcorr(self,squarecrop=True,radcrop=False,split_BG=False,gaussfilt=True,gaussiansigma=4,line_fix=True,normalise=True,normalise_NMF=False,resize=False,screensize=[200,200]):
        
        
        #will go through and update a new EBSP, newEBSP
        
        from scipy.ndimage import gaussian_filter
        import numpy as np
        
        
        #gaussian background filtering - two options: do by half chip by half chip or full thing
        if gaussfilt==True:
            #preallocate a useful array of ones
            corarray=np.ones((self.number,self.patheight,self.patwidth))
            if split_BG==True:
                for i in range(0,self.number):
                    left,right=np.hsplit(self.array[i,:,:],2)
                    
                    gs=gaussiansigma*self.patheight/100
                    
                    left_g=gaussian_filter(left,sigma=int(gs))
                    left=np.divide(left,left_g)
                    
                    right_g=gaussian_filter(right,sigma=int(gs))
                    right=np.divide(right,right_g)
                    
                    recombined=np.hstack((left,right))
                    corarray[i,:,:]=recombined
                newEBSP=EBSP(corarray,self.Metadata,rect=True)
            
            else:
                for i in range(0,self.number):
                    gs=gaussiansigma*self.patheight/100
                    pat=self.array[i,:,:]
                    pat_g=gaussian_filter(pat,sigma=gs)
                    pat=np.divide(pat,pat_g)
                    corarray[i,:,:]=pat
                newEBSP=EBSP(corarray,self.Metadata,rect=True)
        
        else:
            newEBSP=self
        
        #crop the EBSP to a square           
        if squarecrop==True:
            #preallocate a new array
            corarray=newEBSP.array
            
            #crop
            minsize=min([newEBSP.patheight,newEBSP.patwidth])
            centrepoint=[int(newEBSP.patheight/2),int(newEBSP.patwidth/2)]
            s=np.linspace(0,minsize,minsize,endpoint=False,dtype='int')
            s=s-int(s.mean())
            sy=s+centrepoint[0]-1#+0.5
            sx=s+centrepoint[1]-1#+0.5
            EBSP1=corarray[:,sy[:,None],sx[None,:]]
            
            #update the newEBSP
            newEBSP=EBSP(EBSP1,newEBSP.Metadata,rect=True) #re-run through class to update patheight etc
            
        if line_fix==True:
            
            #find row numbers to convert to noise
            loc1=int(0.015*newEBSP.patheight)
            loc2=int(0.0333*newEBSP.patheight)
            
            for n in range(0,newEBSP.number):
                for i in range(loc1,loc2):
                    line=np.random.rand(newEBSP.patwidth)
                    newEBSP.array[n,i,:]=line*(2*newEBSP.std[n])-1.5*newEBSP.std[n]+newEBSP.mean[n]

        if resize==True:
            
            from PIL import Image
            
            corarray=np.zeros((self.number,screensize[0],screensize[1]))
            for i in range(0,self.number):
                pat=np.array(newEBSP.array[i,:,:])
                
                if squarecrop==True:
                    corarray[i,:,:]=np.array(Image.fromarray(pat).resize((screensize[0],screensize[1]),resample=Image.BICUBIC))
                else:
                    corarray[i,:,:]=np.array(Image.fromarray(pat).resize((screensize[0],screensize[1]),resample=0))
            
                #imresize(array1,(screensize[0],screensize[1]),interp='bicubic')
            newEBSP=EBSP(corarray,newEBSP.Metadata,rect=False)
            
        if normalise==True:
            
            #preallocate a useful array of ones
            corarray=np.ones((newEBSP.number,newEBSP.patheight,newEBSP.patwidth))

            if normalise_NMF==True:
                for n in range(0,newEBSP.number):
                    corarray[n,:,:]=(newEBSP.array[n,:,:]-np.amin(newEBSP.array[n,:,:]))/newEBSP.array[n,:,:].std()
                newEBSP=EBSP(corarray,newEBSP.Metadata,rect=True)

            else:
                for n in range(0,newEBSP.number):
                    corarray[n,:,:]=(newEBSP.array[n,:,:]-newEBSP.array[n,:,:].mean())/newEBSP.array[n,:,:].std()
                newEBSP=EBSP(corarray,newEBSP.Metadata,rect=True)


        return newEBSP
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

