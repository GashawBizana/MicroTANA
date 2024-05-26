import os
import pandas as pd
import numpy as np
import glob

# Define the directory path and core name
dir_path = "E:/Anneal500+Compression50Al161/GradeAcalcA500C50Al161"
core_name = "dump.anneal500Al161_P50"
nf = 98

# Function to read matrix data from a CSV file and skip the first row and column
def read_matrix(file_path):
    return pd.read_csv(file_path,sep=';', header=None).iloc[1:, 1:].to_numpy()

# Read the matrices
GrainId = read_matrix(f"{dir_path}/TimeEvo/{core_name}__GrainData__TimeEvo_Grain ID.csv")
X = read_matrix(f"{dir_path}/TimeEvo/{core_name}__GrainData__TimeEvo_PosX.csv")
Y = read_matrix(f"{dir_path}/TimeEvo/{core_name}__GrainData__TimeEvo_PosY.csv")
Z = read_matrix(f"{dir_path}/TimeEvo/{core_name}__GrainData__TimeEvo_PosZ.csv")
NumberOfAtoms = read_matrix(f"{dir_path}/TimeEvo/{core_name}__GrainData__TimeEvo_NumAtoms.csv")


# Initialize TrackedId with NaNs
shape = GrainId.shape
TrackedId = np.full(shape, np.nan)

# Create lists to store filenames
C = []
C1 = []

# Generate file lists
for i in range (1,nf):
    
    a=dir_path+"/"+core_name+"__AtomData_" + str(i)+ "_*.cfg"
    b=dir_path+"/"+core_name+"__GrainData_" + str(i)+ "_*.csv"
    C.append(glob.glob(a)[0])
    C1.append(glob.glob(b)[0])
    



# Map grain IDs
for i in range(shape[0]-1):
    GrainData = pd.read_csv( C1[i], skiprows=[0,1],sep=';')
    grainId_GrainData = GrainData.iloc[:, 0].to_numpy()
    X_GrainData = GrainData.iloc[:, 4].to_numpy()
    Y_GrainData = GrainData.iloc[:, 5].to_numpy()
    Z_GrainData = GrainData.iloc[:, 6].to_numpy()
    n_GrainData = GrainData.iloc[:, 1].to_numpy()
    sizeGrainData = GrainData.shape

    for j in range(shape[1]):
        x = X[i, j]
        y = Y[i, j]
        z = Z[i, j]
        n = NumberOfAtoms[i, j]
       
        for k in range(sizeGrainData[0]):
                
            if n == n_GrainData[k] and x == X_GrainData[k] and y == Y_GrainData[k] and z == Z_GrainData[k]:
                
                TrackedId[i, j] = k

# Save the TrackedId array
np.save(f"{dir_path}/{core_name}_grain_id_mapping.npy", TrackedId)
