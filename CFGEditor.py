import numpy as np ##### To keep consistency with dynamore


def findAtomsGrain( fileIn):
    
    f = open(fileIn,'r')
#    fout = open(fileOut,'w')
        
    H0 = np.zeros((3,3))
    EdgeLength=np.zeros((3))
    
    header = f.readline()
    header = header.strip()
    ch = header.split()
    
    n_p = int(ch[4])
    
    x = np.zeros((n_p))
    y = np.zeros((n_p))
    z = np.zeros((n_p))
    gId = np.zeros((n_p),int)
    ix = np.zeros((n_p),int)
    iy = np.zeros((n_p),int)
    iz= np.zeros((n_p),int)
    
    
    
    counter = 0
    
    for line in f:
        
        text = line.strip()
        columns = text.split()
        
        if text == '' or text =='.NO_VELOCITY.':
            continue
        
        if columns[0] == 'H0(1,1)':
            H0[0,0]=columns[2]
            EdgeLength[0]=columns[2]
            
            continue
            
        if columns[0] == 'H0(1,2)':
            H0[0,1]=columns[2]
            continue
            
        if columns[0] == 'H0(1,3)':
            H0[0,2]=columns[2]
            continue
            
        if columns[0] == 'H0(2,1)':
            H0[1,0]=columns[2]
            continue
            
        if columns[0] == 'H0(2,2)':
            H0[1,1]=columns[2]
            EdgeLength[1]=columns[2]
            continue
            
        if columns[0] == 'H0(2,3)':
            H0[1,2]=columns[2]
            continue
            
        if columns[0] == 'H0(3,1)':
            H0[2,0]=columns[2]
            continue
            
        if columns[0] == 'H0(3,2)':
            H0[2,1]=columns[2]
            continue
            
        if columns[0] == 'H0(3,3)':
            H0[2,2]=columns[2]
            EdgeLength[1]=columns[2]
            continue
            
        if columns[0] == 'entry_count':
            entryCount = int(columns[2])-3+2
            for i in range(entryCount):
                line = f.__next__()
                line.strip()
                
            continue 
        if counter<=n_p-1:
            x[counter] = columns[0]
            y[counter] = columns[1]
            z[counter] = columns[2]
            gId[counter] = columns[6]
            ix[counter]=columns[7]
            iy[counter]=columns[8]
            iz[counter]=columns[8]
            counter = counter+1
       
    np.resize(x,(counter,1))
    np.resize(y,(counter,1))
    np.resize(z,(counter,1))
    np.resize(gId,(counter,1))
    np.resize(ix,(counter,1))
    np.resize(iy,(counter,1))
    np.resize(iz,(counter,1))
    
    X=np.array([x,y,z]) *
    
       
        #print "%d\n" % counter
        
    return (x,y,z,gId,ix,iy,iz,H0,n_p,EdgeLength)