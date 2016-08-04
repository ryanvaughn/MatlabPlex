%Test Script
clear;
n=5;
tophom=2;
sizecompletegraph=6;

%numedges=15;
%A=[1,1,1,1,1,2,2,2,2,3,3,3,4,4,5;2,3,4,5,6,3,4,5,6,4,5,6,5,6,6;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
%translate=[1,2,3,4,5,6];
%B=zeros(3,n*numedges);
%B(1:3,1:numedges)=A;
 
%numedges=6;
 %A=[1,1,1,2,2,3;2,3,4,3,4,4];
 %translate=[1,2,3,4];
 %B=zeros(2,n*numedges);
 %B(1:2,1:numedges)=A;
 
 
 %numpts=n*sizecompletegraph;

%  for i=1:n
%      %keyboard;
%      B(1:2,numedges*(i-1)+1:numedges*i)=translate(A(1:2,1:numedges));
%      translate=translate+sizecompletegraph;
%  end
% rng(0)
% edges=B(1:2,:);
% edgediam=rand(1,n*numedges);


numpts=4;edges=[1 2 3 4 3;2 3 4 1 1];edgediam=[1,2,3,4,5];%square with diagonal
%numpts=8;edges=[1 2 3 4 5 6 7 8 1 2 3 4; 2 3 4 1 6 7 8 5 5 6 7 8]
%numpts=5;edges=[1 2 3 4 1 2 3 4 1; 2 3 4 1 5 5 5 5 3];%pyramid with diagonal on base
%numpts=8;edges=[1 2 3 4 5 6 7 8 1 2 3 4 1 5 1 2 3 4;2 3 4 1 6 7 8 5 5 6 7 8 3 7 6 7 8 5];%cube with triangles on all faces

tic;
edges=uint64(edges);
[dimen, diam, bound]=VRExpansion(numpts,edges, edgediam,tophom+1);
toc;
dimen=double(dimen);
bound=double(bound);
[q,p]=size(bound);
time=diam;
ComputeIntervals
maxtime=max(time);BettiNumbers(L,maxtime+.5,tophom)
PlotBarCodes(L,tophom)

keyboard;
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 tic;
 [dimen,bound]=Edges2VR(numpts,edges,tophom);
 toc;
 
 time=zeros(1,p);
 
  for i=1:p 
     if dimen(i)==1
        test= find((edges(1,:)==bound(1,i)) .* (edges(2,:)==bound(2,i)));
        time(i)=edgediam(test);
     end 
        if dimen(i)==2
            a=bound(:,bound(1,i));
            b= bound(:,bound(2,i));
            c=bound(:,bound(3,i));
            test1= find((edges(1,:)==a(1)) .* (edges(2,:)==a(2)));
            test2= find((edges(1,:)==b(1)) .* (edges(2,:)==b(2)));
            test3= find((edges(1,:)==c(1)) .* (edges(2,:)==c(2)));
            time(i)=max([edgediam(test1),edgediam(test2),edgediam(test3)]);
        end
            
  end
 
 
 [q,p]=size(bound);
 ComputeIntervals
 maxtime=max(time);BettiNumbers(L,maxtime+.5,tophom)
 PlotBarCodes(L,tophom)