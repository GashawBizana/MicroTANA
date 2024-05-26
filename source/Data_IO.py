
import numpy as np
import os.path
import time
#Data In

######## Function ###########

def HeaderFile(file):
    # This function reads the header for a CFG file genrated from Grade-A
    with open(file) as myfile:
        head = [next(myfile) for x in range(24)]
    return(head)

def Header(fileIn):
    # This function reads the header for a CFG file 
    f = open(fileIn,'r')
    
    header = f.readline()
    header = header.strip()
    

def ReadGrainData(file):
    # This function reads the positions and other per-atom properties for a CFG file genrated from Grade-A
    GrainData=np.loadtxt(file, skiprows = 24)
    return(GrainData)

def GrainEdgeLength(FileToAnalize):
    # This function extracts the simulation box size from a CFG file
    with open(FileToAnalize) as file:
        for line in file:
            if ('Number of particles ' in line):
                NumberOfAtoms_Str=line.split("= ",1)[1]
                NumberOfAtoms=int(NumberOfAtoms_Str[0:-1])
            elif ('H0(1,1)' in line):
                xl_str=line.split("= ",1)[1]
                xl=float(xl_str[0:-1])
            elif ('H0(2,2)' in line):
                yl_str=line.split("= ",1)[1]
                yl=float(yl_str[0:-1])
            elif ('H0(3,3)' in line):
                zl_str=line.split("= ",1)[1]
                zl=float(zl_str[0:-1])
    return(NumberOfAtoms,xl,yl,zl)


def findAtomsGrain( fileIn):
    '''
    This function extracts the atomic position, periodic image, simulation box siz, header info,
    grain ID and other per atom properties from a CFG file
    '''
    doefileexist=os.path.isfile(fileIn)
    if doefileexist==False:
        time.sleep(10)
        
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
            EdgeLength[2]=columns[2]
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
            if entryCount>9:
                ix[counter]=columns[6]
                iy[counter]=columns[6]
                iz[counter]=columns[6]
            counter = counter+1
       
    np.resize(x,(counter,1))
    np.resize(y,(counter,1))
    np.resize(z,(counter,1))
    np.resize(gId,(counter,1))
    np.resize(ix,(counter,1))
    np.resize(iy,(counter,1))
    np.resize(iz,(counter,1))
    
    X=(np.array([x,y,z]).T)*EdgeLength
    Image=np.array([ix,iy,iz]).T
    
       
        #print "%d\n" % counter
        
    return (X,gId,Image,H0,n_p,EdgeLength)

#Data Out

#######  Function ###################

def OutPutCFGWritter(ModifiedGrainIdListSorted,GrainData,GrainId,file,fname):
    '''
    This functions writes a CFG file conatining atomic positions, grain Id and connectivity information
    for each atoms
    '''
    M=[]
    C=[]
    TCL=[]
    OutputFile='./TestDataIO/MT_'+fname
    for i in range(0,len(ModifiedGrainIdListSorted)):
       M.append(ModifiedGrainIdListSorted[i][1])
       C.append(ModifiedGrainIdListSorted[i][3])
       TL=ModifiedGrainIdListSorted[i][2]
       if (ModifiedGrainIdListSorted[i][1]==-1):
           TL=np.array([TL[0],TL[1],-1,-1])
       elif (ModifiedGrainIdListSorted[i][1]==-2):
           TL=np.array([TL[0],TL[1],TL[2],-1])
       elif (ModifiedGrainIdListSorted[i][1]==-3):
           TL=np.array([TL[0],TL[1],TL[2],TL[3]])
       elif (ModifiedGrainIdListSorted[i][1]>=0):
           TL=np.array([TL[0],TL[0],TL[0],TL[0]])
       else:
           TL=np.array([-1,-1,-1,-1])
       TCL.append(TL)    
                
    MM=np.transpose(np.array(M,ndmin=2))
    CC=np.transpose(np.array(C,ndmin=2))
    
    all_data = np.hstack((GrainData[:,[0,1,2,3,4,5]],GrainId,GrainData[:,[7,8]], MM,CC,TCL))
    Header="".join(map(str,HeaderFile(file)))
    f = open(OutputFile, 'w')
    f.write(Header)
    f.close()
    with open(OutputFile, 'a+') as outfile:
        np.savetxt(outfile,all_data, fmt="%5f")
    
