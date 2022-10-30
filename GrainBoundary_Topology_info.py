# class Microstructure:
      
    






# class Grain:
    
#     def __init__(self, grain_id=0):
#         grain_id=self.grain_id
    
#     def grain_boundary(grain_id,InnerAtom=innerAtom):
        
#         grain_boundary_associated_with_grain=[j[0] for j in  InnerAtom if  (Grain_id[0] in j[0]) ]
        
        
        
    




def Microstructure_topology(InnerAtom,GrainBoundary,TripleLine):
    GrainBoundary_number_of_triple_line=[]
    for i in range(0,len(GrainBoundary)):
    
        Grain_Boundary_id=GrainBoundary[i][0]
        
        triple_lines_associated_with_grain_boundary=[j[0] for j in  TripleLine if  (Grain_Boundary_id[0] in j[0] and Grain_Boundary_id[1] in j[0]) ]
        number_of_triple_lines_associated_with_grain_boundary=len(triple_lines_associated_with_grain_boundary)
        
        GrainBoundary_number_of_triple_line.append([Grain_Boundary_id,number_of_triple_lines_associated_with_grain_boundary, triple_lines_associated_with_grain_boundary])
        
    
        
    Grain_number_of_grain_boundaries=[]
    
    for i in range(0,len(InnerAtom)):
    
        Grain_id=InnerAtom[i][0]
        
        grain_boundary_associated_with_grain=[j[0] for j in  GrainBoundary if  (Grain_id[0] in j[0]) ]
        number_of_grain_boundaries_associated_with_grain=len(grain_boundary_associated_with_grain)
        
        Grain_number_of_grain_boundaries.append([Grain_id,number_of_grain_boundaries_associated_with_grain, grain_boundary_associated_with_grain])
        
    Grain_number_of_triple_lines=[]
    
    for i in range(0,len(InnerAtom)):
    
        Grain_id=InnerAtom[i][0]
        
      
        
        triple_line_associated_with_grain=[j[0] for j in  TripleLine if  (Grain_id[0] in j[0]) ]
        number_of_triple_lines_associated_with_grain=len(triple_line_associated_with_grain)
        
        Grain_number_of_triple_lines.append([Grain_id,number_of_triple_lines_associated_with_grain, triple_line_associated_with_grain])
      