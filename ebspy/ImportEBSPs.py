def bImportEBSPs (ebspfilepath,importdata,region=[]):
    
    import h5py as h5
    from pathlib import Path
    import numpy as np
    import matplotlib.pyplot as plt
    
    #TPM 2019 import all HDF5 values as a dictionary
    #requires Path, h5 and numpy
    
    imp=h5.File(ebspfilepath,'r')
    
    #%
    listofnames=[]
    imp.visit(listofnames.append)
    Metadata={}
    Data={}
    
    
    #First extract metadata
    for i in range(0,len(listofnames)):
        
        dataset=imp[listofnames[i]]  
        
        try: 
            dataset.keys()
        except AttributeError:
            name=listofnames[i]
            
            if r'Phases' in name:
                continue
            
            name2=name[name.rfind('/')+1:]
            val=dataset[()]
              
            if (name2 != r'Counts Corrected') or (name2 != r'Counts Raw') or (name2 != r'Raw Patterns'):
                Metadata[name2]=val            
            
            
    #%now extract actual data
    size1=Metadata['NROWS']
    size2=Metadata['NCOLS']
    patsize1=Metadata['PatternHeight']
    patsize2=Metadata['PatternWidth']
    
    #region: h1 h2 w1 w2
    
    if region==[]:
        region=[0,size1,0,size2]
    
    subset_h1=region[0]
    subset_h2=region[1]
    nrows=subset_h2-subset_h1
    
    subset_w1=region[2]
    subset_w2=region[3]
    ncols=subset_w2-subset_w1
    
    toextract=[]

    for i in range(0,nrows):
        for j in range(0,ncols):
            ind=(subset_h1*size2)+(i*size2)+(subset_w1+j)
            #fully skipped rows + already completed rows + distance along this row
            toextract.append(ind)
            
    testaxis=np.ones(size1*size2)
    toextract=np.array(toextract,dtype='int')
    testaxis2=testaxis
    testaxis2[toextract]=2
    
    #visualisation of output
    visual=np.reshape(testaxis2,(size1,size2))
    
    #import data if requested
    if importdata==True:
        for i in range(0,len(listofnames)):
        
            dataset=imp[listofnames[i]]
            
            try: 
                dataset.keys()
            except AttributeError:
                name=listofnames[i]
                
                name2=name[name.rfind('/')+1:]
            
                
                if name2 == r'Counts Corrected':
                    val=np.array(dataset[()])
                    Data['EDSCor']=val[toextract,:] #spatial, channums
                if name2 == r'Counts Raw':
                    val=np.array(dataset[()])
                    Data['EDSRaw']=val[toextract,:] #spatial, channums
                if name2 == r'RawPatterns':
                    Data['Patterns']=np.zeros((len(toextract),patsize1,patsize2))
                    Data['Patterns'][:,:,:]=(dataset[toextract,:,:])
                    
#                    Data['Patterns']=np.zeros([2,patsize1,patsize2])
#                    for n in range(0,len(toextract)):
#                        val=np.array(dataset[toextract[n],:,:]) #spatial, patheight, pat cols
#                        val=val.reshape(1,patsize1,patsize2)
#                        Data['Patterns']=np.concatenate((Data['Patterns'],val))
#                    Data['Patterns'][0:2,:,:]=None
                    
    return(Metadata,Data,visual,region)