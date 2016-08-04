function [lNbrs, Mdim] = lowerNbrs( N, edges )
%The lower-Nbrs function as described in Zomorodian's "Fast Construction of
%the Vaetoris-Rips Complex.
%   The input is a matrix whose colums correspond to endpoints of edges and
%   a weight for each edges in a Neighborhood Graph. For more information
%   regarding the format of edges, see VRExpansion.

%N=number of data points


%Using a matrix
lNbrs=zeros(N,N); 
Mdim=zeros(1,N);

% ordering is a vector that lists the vector out in order of occurence in
% edges. 

 ordering=unique(edges(2,:));
 inputsize=size(ordering);

 for i= ordering %This might be faster as a for loop going through all edges.
     B= edges(2,:)== i;
     Mdim(i)=sum(B);
     lNbrs(1:Mdim(i),i)=edges(1,B);    
 end
end

