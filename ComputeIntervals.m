% Output L is 4xq where q is the number of bars
% Each column has start time, end time, index of generator, dimension of
% generator
T=zeros(p,p);
T1=zeros(1,p);
mark=zeros(1,p);
L=[];
maxtime=max(time);
for j=1:p
  d=RemovePivotRows(j,bound,dimen,T);
  if isempty(d)
    mark(j)=1;
  else
    i=max(d);
    k=dimen(j);
    T1(i)=j;
    T(1:length(d),i)=d;
    L=[L [time(i);time(j);i;dimen(i)]];
  end
end
fm=find(mark);
f=find(T(1,fm)==0);
%Get infinite intervals
for i=1:length(f)
  L=[L [time(fm(f(i)));maxtime+1;fm(f(i));dimen(fm(f(i)))]];
end
%Sort by dimension
[~,r]=sort(L(4,:));
L=L(:,r);


