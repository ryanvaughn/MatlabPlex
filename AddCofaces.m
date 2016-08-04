function [ inds, rinds, Exp, vExp, Cofaces, dimen, numCofaces, final ] = AddCofaces(currentVertex, lNbrs, Mdim, tophom, inds, rinds, Simplex, initial, Exp, vExp, Cofaces, dimen, numCofaces, M)
%AddCofaces is an implementation of the algorithm of the same name from the
%paper "Fast Construction of the Vaetoris-Rips Complex". 

%VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edges- the input edge matrix data.

% currentVertex- The enumeration of the vertex in the main for loop. This
% is only used to update inds, so is not that important.

% tophom- the highest dimension of interest. Again, not used or updated
% often.

% inds- a vector whose ith entry is the index of the ith vertex in Exp

% Simplex- the simplex being added, in column vector format, where nonzero
% entrees are the indices of the dim-1 cofaces of the simplex

% M- the indices of each vertex which is both adjacent to all vertices in
% the simplex, and whose enumeration is lower than all the vertices of the
% simplex.

%Exp- A matrix whose columns correspond to simplices in the complex. Each
%entry is a listing of the indices of the dimension-1 cofaces of the
%simplex.

%vExp- A matrix whose columns correspond to simplices in the complex. Each
%entry is a listing of the indices of the vertices of the simplex.

%dimen- a vector whose entrees correspond to the dimension of the simplex
%in Exp.

%numCofaces- a vector whose entrees correspond to the number of cofaces
%of the simplex with the given index.

% initial- the index of the next open column in the Exp and/or vExp matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ALGORITHM

%Add in the input simplex to Exp and update inds, vExp, Cofaces, and dimen
%to reflect the new addition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Update Exp
    Exp(1:length(Simplex),initial)=Simplex;

%update dimen
    dimen(initial)=length(Simplex)-1;
    
%update vExp
    if sum(Simplex)==0 
        vExp(1,initial)=initial;
    else
        SimplexVerts=union(vExp(:,Exp(1,initial)),vExp(:,Exp(2,initial)));
        if SimplexVerts(1)==0
            vExp(1:length(Simplex),initial)=SimplexVerts(2:length(SimplexVerts));
        else
            vExp(1:length(Simplex),initial)=SimplexVerts(1:length(SimplexVerts));
        end
    end

        
%update Cofaces
    if sum(Simplex)~=0
        for i=Simplex' %Need to index a for loop by a row vector.
            try
            numCofaces(i)=numCofaces(i)+1;
            Cofaces(numCofaces(i),i)=initial;
            catch
                keyboard;
            end
            end
    end
    
%Update inds if the input simplex is a vertex.
if sum(Simplex)==0
    inds(currentVertex)=initial;
    rinds(initial)=currentVertex;
end

%Final step: move to the next column.
final=initial+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(dimen(initial)>=tophom)%Do nothing if the dimension of the simplex is too high.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for x=M %Cycle through vertices in the input
    if isempty(M)==0
    counter=1;
    
    newSimplex=zeros(length(Simplex)+1,1);%Construct a new simplex
    newSimplex(counter)=initial;%The first entry of the simplex will be the previous simplex.
    counter=counter+1;
    if sum(Simplex)==0 %if the previous simplex was a vertex, construct an edge.
        newSimplex=[initial;x];
        
    else
        for i=Simplex'
            %keyboard;
            for j=Cofaces(1:numCofaces(i),i)'
                if sum(vExp(1:dimen(j)+1,j)==x)==1;
                        newSimplex(counter)=j;
                        counter=counter+1;
                end
            end
        end
    end
    %keyboard;
    vertexnumx=rinds(x);

    Set=intersect(inds(lNbrs(1:Mdim(vertexnumx),vertexnumx)),M);
    [ inds, rinds, Exp, vExp, Cofaces, dimen, numCofaces, final ] = AddCofaces(currentVertex, lNbrs, Mdim, tophom, inds,rinds, newSimplex, final, Exp, vExp, Cofaces, dimen, numCofaces, Set);  
    end
end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end % this closes off the if statement telling us not to go above the tophom.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

