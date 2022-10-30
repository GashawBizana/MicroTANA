def GrainBoundary_distance_computation(EdgeLength_for_all_timesteps,grain_boundary_property_for_all_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps,GrainBoundary_number_of_triple_line_for_all_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps,grain_property_for_all_timesteps,grain_boundary_property_extended_for_all_timesteps,grain_mapping_data,step,Grain_VonNeuman_Data_for_timesteps):
    
    Grain_Boundary_Velocity_Data_for_all_timesteps=[]
    
    for k in range(0,len(grain_boundary_property_for_all_timesteps)-step):
        print(k)
        
        if not Grain_VonNeuman_Data_for_timesteps: 
        Grain_VonNeuman_Data_step_k=Grain_VonNeuman_Data_for_timesteps[k]
        
        
        EdgeLength_step_k= EdgeLength_for_all_timesteps[k]
        EdgeLength_step_k_plus_1=EdgeLength_for_all_timesteps[k+step]
        
        grain_boundary_property_step_k=grain_boundary_property_extended_for_all_timesteps[k]
        grain_boundary_property_step_k_plus_1=grain_boundary_property_extended_for_all_timesteps[k+step]
        
        grain_boundary_total_curvature_step_k=[j[1] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        grain_boundary_total_surface_area_step_k=[j[3] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        grain_boundary_mean_curvature_step_k =[j[2] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        
        
        # grain_boundary_forming_step_k=[j[0] for j in grain_boundary_property_step_k]
        
        grain_boundary_forming_step_k=[j[0] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        
        grain_property_step_k=grain_property_for_all_timesteps[k]
        grain_property_step_k_plus_1=grain_property_for_all_timesteps[k+step]
        
        GrainBoundary_length_of_triple_line_step_k=GrainBoundary_length_of_triple_line_for_all_timesteps[k]
        GrainBoundary_length_of_triple_line_step_k_plus_1=GrainBoundary_length_of_triple_line_for_all_timesteps[k+step]
        GrainBoundary_misorientation_step_k=GrainBoundary_misorientation_for_all_timesteps[k]
        GrainBoundary_misorientation_step_k_plus_1=GrainBoundary_misorientation_for_all_timesteps[k+step]
        GrainBoundary_number_of_triple_line_step_k=GrainBoundary_number_of_triple_line_for_all_timesteps[k]
        GrainBoundary_number_of_triple_line_step_k_plus_1=GrainBoundary_number_of_triple_line_for_all_timesteps[k+step]
        
        Grain_Boundary_Velocity_Data_k=[]
        for i in range(0,len(grain_boundary_total_curvature_step_k)):
            
            # print(i)
            grain_boundary_total_curvature_step_k_i=grain_boundary_total_curvature_step_k[i]
            grain_boundary_total_surface_area_step_k_i=grain_boundary_total_surface_area_step_k[i]
            grain_boundary_mean_curvature_step_k_i =grain_boundary_mean_curvature_step_k[i]
            
            # if grain_boundary_total_surface_area_step_k_i==0:
            #     grain_boundary_mean_curvature_step_k_i=0
            # elif grain_boundary_total_surface_area_step_k_i !=0:
            #     grain_boundary_mean_curvature_step_k_i=grain_boundary_total_curvature_step_k_i/grain_boundary_total_surface_area_step_k_i
            
            
            grain_id_k_i=grain_boundary_forming_step_k[i]
            grain_id_k_i_initial=Initial_Grain_Id(grain_mapping_data, grain_id_k_i ,k)
            
            grain_id_k_plus_1_i=GrainIdMapper(grain_mapping_data, grain_id_k_i ,k,step)
            
            grain_boundary_k_i_triangleId=[j[2] for j in grain_boundary_property_step_k if (grain_id_k_i[0] == j[0][0] and grain_id_k_i[1] == j[0][1])]
            
            grain_boundary_k_i=[j[1] for j in grain_boundary_property_step_k if (grain_id_k_i[0] == j[0][0] and grain_id_k_i[1] == j[0][1])]
            
            grain_boundary_k_plus_1_i=[j[1] for j in grain_boundary_property_step_k_plus_1 if (grain_id_k_plus_1_i[0] == j[0][0] and grain_id_k_plus_1_i[1] == j[0][1])]
            
            
            grain_id_k_i0=grain_id_k_i[0]
            grain_id_k_i1=grain_id_k_i[1]
            grain_id_k_plus_1_i0=grain_id_k_plus_1_i[0]
            grain_id_k_plus_1_i1=grain_id_k_plus_1_i[1]
            
            
            GrainBoundary_length_of_triple_line_step_k_i=[j[1] for j in GrainBoundary_length_of_triple_line_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            GrainBoundary_misorientation_step_k_i= [j[1] for j in GrainBoundary_misorientation_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            GrainBoundary_number_of_triple_line_step_k_i= [j[1] for j in GrainBoundary_number_of_triple_line_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            
            if (np.isnan(grain_id_k_plus_1_i).any()) ==False and ((not grain_boundary_k_plus_1_i)== False) and (len(grain_boundary_k_i[0].faces) >=1) and len(grain_boundary_k_plus_1_i[0].faces) >=1:

                
                
                    

                # velocity,curv_py,vol1,=distanace_between_current_and_next(grain_boundary_k_i[0],grain_boundary_k_plus_1_i[0],EdgeLength_step_k_plus_1)
                velocity,curv_py,vol1,vol2=distance_from_von_neuman(grain_id_k_i_initial,grain_boundary_k_i[0],grain_boundary_k_i_triangleId[0],Grain_VonNeuman_Data_step_k)
                        
                Grain_Boundary_Velocity_Data_k.append([grain_id_k_i,grain_id_k_i_initial,grain_id_k_plus_1_i,velocity,grain_boundary_total_curvature_step_k_i,grain_boundary_mean_curvature_step_k_i,grain_boundary_total_surface_area_step_k_i,np.nan,GrainBoundary_length_of_triple_line_step_k_i,GrainBoundary_number_of_triple_line_step_k_i,GrainBoundary_misorientation_step_k_i,curv_py,vol1,vol2])
        
        '''
        0- grain_id_k_i, 1- grain_id_k_i_initial, 2- grain_id_k_plus_1_i, 3- Velocity, 4- grain_boundary_total_curvature_step_k_i, 5-grain_boundary_mean_curvature_step_k_i, 6- grain_boundary_total_surface_area_step_k_i, 7- grain_boundary_total_volume_step_k_i, 8- GrainBoundary_length_of_triple_line_step_k_i, 9- GrainBoundary_number_of_triple_line_step_k_i, 10- GrainBoundary_misorientation_step_k_i
        '''
        Grain_Boundary_Velocity_Data_for_all_timesteps.append(Grain_Boundary_Velocity_Data_k)
    
    return(Grain_Boundary_Velocity_Data_for_all_timesteps)