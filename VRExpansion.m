function [dimen, diam, Exp]=VRExpansion(numpts,edges, edgediam ,tophom)
%numpts is number of 0-simplices
%edges is 2xe matrix of e 1-simplices, lower vertices are on first row.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Zomorodian's Algorithm for Vaetoris-Rips expansion is an algorithm that
%converts the data of a neighborhood graph obtained from point cloud data
%into a VR-Complex. A Neighborhood Graph is an edge-weighted graph where
%weight corresponds to the Euclidean distance between the two points.
%
%The input data is a 3xe matrix, where each column contains the data for
%one edge in the weighted graph. The first two entrees of each column are
%the enumerations of the endpoints, and the third is the weight of the edge
% which corresponds to the distance between the points in the point cloud.
% Typically, the graph will not be a complete graph, as we are interested
% in only studying the filtration of complexes up to a fixed distance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preprocessing: We first make sure that the input Edges Matrix is
%appropriately formatted. We require an 3xe matrix, where e is the number
%of edges in the Neighborhood graph. Edges must be input so that the
%smaller endpoint is above the larger endpoint.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Discretize edgediam into integer values. Assume that no two pairs of
%points have equal distance apart.
[~,perm]=sort(edgediam);
edges=edges(:,perm);
edgediam=1:length(edgediam);

%make sure edges is properly sorted.
edges(1:2,:)=sort(edges(1:2,:)); %Ensure that the lower vertex is always in the first row.
[~,perm]=sort(edges(2,:)); 
edges=edges(:,perm); %Sort columns by the second row
edgediam=edgediam(perm); %keep track of the weights

[~,perm]=sort(edges(1,:));
edges=edges(:,perm); %Sort columns by the first row so that we have dictionary ordering.
edgediam=edgediam(perm); %again, keep track of the weights

%Make sure that the matrix is the right size
[x,~]= size(edges);
if x~=2
    error('The columns of edge matrix are not the correct size.')
end

%Probably need some errors to catch putting in too many edges.

%Make sure that the endpoints of edges are ordered properly.
if prod(edges(1,:)<edges(2,:))==0
   error('The endpoints of edge matrix are not ordered properly.') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization of the Algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Big=100000;
Small =100; %Try to make this more precise.
Exp=zeros(tophom+1,Big,'uint64');
vExp=zeros(tophom+1,Big,'uint64');
Cofaces=zeros(tophom+1,Small,'uint64');% need to calculate a bound of the cofaces

numCofaces=zeros(1,Big,'uint64');
dimen= zeros(1,Big,'uint64');

[lNbrs,Mdim]=lowerNbrs(numpts,edges);
lNbrs=uint64(lNbrs);
Mdim=uint64(Mdim);
inds=zeros(1,numpts,'uint64');
rinds=zeros(1,Big,'uint64');
initial=1;

for i=1:numpts %Cycle through vertices
    
        %Update inds to reference index of vertex i in Exp instead of the
        %enumeration.
        
        M= inds(lNbrs(1:Mdim(i),i));
        Simplex= [0]; %We input each vertex like this so that dimen will be calculated correctly.        
        [inds,rinds, Exp, vExp, Cofaces, dimen, numCofaces, initial ] = AddCofaces(i, lNbrs, Mdim, tophom, inds, rinds, Simplex, initial, Exp, vExp, Cofaces, dimen, numCofaces, M);
end
%Output
Exp=Exp(:,1:initial-1);
dimen=dimen(1:initial-1);
diam=zeros(1,initial-1);

%Computation of weight function
for i=1:initial-1
   if dimen(i)==0
      diam(i)=0; 
   end
   if dimen(i)==1
       diam(i)= edgediam(find((edges(1,:)==rinds(vExp(1,i))) .* (edges(2,:)==rinds(vExp(2,i))))); %This might be slow       
   end
   if dimen(i)>1
      diam(i)=max(diam(Exp(1:dimen(i)+1,i)'));
   end
   
end

end
