function [Cd,D]=Edges2VR(numpts,edges,tophom)
%numpts is number of 0-simplices
%edges is 2xe matrix of e 1-simplices
S=numpts;Big=40000;%to leave space for the maximum number of simplices.
tophom2=tophom+2;
C=zeros(tophom2,Big);D=C; % column i of C are indices of points composing simplex i
C(1,1:numpts)=1:numpts;  % column i of D are indices of faces of simplex i
E=zeros(30,Big); % column i of E are indices of simplices for which i is a face
%E looks up all simplices that contain the simplex
Cd=zeros(1,Big); % Cd(i) is dim of simplex i
Ed=zeros(1,Big); % Ed(i) is length of nonzeros in column i of E 
[~,numedges]=size(edges);% ~ surpresses the output
S=numpts+1;
for i=1:numedges  % while more edges to add
  %if 10*round(i/10)==i %check i in increments of 10
   % disp([num2str(i) ' of ' num2str(numedges) ' edges'])
  %end
  C(1:2,S)=sort(edges(:,i));
  D(1:2,S)=sort(edges(:,i));
  Cd(S)=1;% S is an edge
  Ed(C(1,S))=Ed(C(1,S))+1;E(Ed(C(1,S)),C(1,S))=S; %for each point in the edge, the point is a face of the edge. (first point, next line is for the other point)
  Ed(C(2,S))=Ed(C(2,S))+1;E(Ed(C(2,S)),C(2,S))=S;
  T=S;  %S = number of processed simplices; T = number of simplices found--- an edge might make more than one new simplex, so "processing" means that we have added the single new face out of the possibly many added faces.
 
  %Processing all new simplices
  while S <= T   % while more unprocessed simplices exist
   if (Cd(S) <= tophom) %only add stuff that is low enough.
    p=Cd(S)+1;  % number of points in simplex
    newsimplex=C(1:p,S); % grab next simplex
    %search for p-simplexes that newsimplex is the final face of
    match=zeros(p,1);cand=[];
    for j=1:p %Remove 1 point from the newsimplex and check if there is another point that adds in a simplex
      fj=D(j,S); %fj=index of jth face
      face=C(1:p-1,fj);% need a common pt that forms p-1 simplex with all p faces
      numk(j)=Ed(fj);cand=union(cand,E(1:numk(j),fj));
      g=find(E(1:numk(j),fj)<=S);numg(j)=length(g); %ignore containing simplices newer than S
      for k=1:numg(j)                               %will get them later
        match(j,k)=setxor(face,C(1:p,E(g(k),fj)));
      end
    end
    %The jth row of match is the set of new points such that when you
    %subtract jth point, the new points with all the other is a known
    %simplex that has already been processed.
    
    inter=match(1,1:numg(1));
    for j=2:p
      inter=intersect(inter,match(j,1:numg(j)));
    end % Find any numbers occuring in all rows of match
     for kk=1:length(inter) % if length(inter)>0 a new simplex has been found
      newsimplex1=sort([inter(kk);newsimplex]); 
      C(1:p+1,T+1)=newsimplex1;
      T=T+1;
      for j=1:p+1 %update D and E with then new simplex information
        face=[newsimplex1(1:j-1); newsimplex1(j+1:p+1)];
        [~,~,f]=intersect(face',C(1:p,cand)','rows');
        D(j,T)=cand(f);
        Ed(cand(f))=Ed(cand(f))+1;
        E(Ed(cand(f)),cand(f))=T;
      end
      D(1:p+1,T)=sort(D(1:p+1,T));
      Cd(T)=p;
     end
   end
   S=S+1;
  end
end
D=D(:,1:T);
Cd=Cd(1:T);