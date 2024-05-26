dir_path="/home/uda69/lcz56/GradeA_Al250_500"
% define the file name pattern
core_name = "dump.annelAl250";
nf=299;

% if nargin ~= 3
%     error('Three inputs are required: string1, string2, and number');
% end
% 
% % assign the input values to individual variables
% dir_path = string(varargin{1});
% core_name = string(varargin{2});
% nf= varargin{3};

addpath(dir_path)

GrainId=readmatrix(dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_Grain ID");

GrainId= GrainId(2:end,2:end);
X=readmatrix(dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_PosX");
X=X(2:end,2:end);
Y=readmatrix(dir_path+'/TimeEvo/'+core_name+"__GrainData__TimeEvo_PosY");
Y=Y(2:end,2:end);
Z=readmatrix(dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_PosZ");
Z=Z(2:end,2:end);
NumberOfAtoms=readmatrix(dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_NumAtoms");
NumberOfAtoms=NumberOfAtoms(2:end,2:end);

%%
shape=size(GrainId);
TrackedId=NaN(shape(1),shape(2));

%%
shape=size(GrainId);
TrackedId=NaN(shape(1),shape(2));

C=string.empty;
C1=string.empty;
for i=1:nf
    ii=i;

filename=core_name+"__AtomData_"+ii+"_*.cfg";
filename1=core_name+"__GrainData_"+ii+"_*.csv";

fileList = dir(fullfile(dir_path, filename));
fileList1 = dir(fullfile(dir_path, filename1));


C(i)=fileList.name;
C1(i)=fileList1.name;
end
   

for i=1:shape(1)
    
    GrainData= readmatrix(C1(i));
    GrainData=sortrows(GrainData,[2 3],'descend');
    grainId_GrainData=GrainData(:,1);
    X_GrainData=GrainData(:,5);
    Y_GrainData=GrainData(:,6);
    Z_GrainData=GrainData(:,7);
    n_GrainData=GrainData(:,2);
    sizeGrainData=size(GrainData);
    
    for j=1:shape(2)
   x=X(i,j);
   y=Y(i,j);
   z=Z(i,j);
   n=NumberOfAtoms(i,j);
   for k=1:sizeGrainData(1,1)
       if n== n_GrainData(k) && x== X_GrainData(k) && y== Y_GrainData(k) && z== Z_GrainData(k)
        TrackedId(i,j)= k-1; 
       end
   end
    end
end

save(core_name+"_grain_id_mapping.mat",'TrackedId')
