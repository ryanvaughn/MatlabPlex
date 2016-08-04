function d=RemovePivotRows(sig,bound,dimen,T)
k=dimen(sig);
d=bound(1:dimen(sig)+1,sig);
if length(d)==1
  d=[];
end
while ~isempty(d)
  i=max(d);
  if T(1,i)==0
    break
  end
  d=setxor(d,T(:,i));
  d=d(d>0);
end
d=d(d>0);