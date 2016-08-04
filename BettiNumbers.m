function Betti=BettiNumbers(L,t,tophom)
%Compute Betti numbers at time t
maxdim=max(L(4,:));maxdim=min(maxdim,tophom);
Betti=zeros(maxdim+1,1);
for dim=0:maxdim
  f=find(L(4,:)==dim);num=length(f);L1=L(:,f);
  counter=0;
  for i=1:num
    if L1(1,i)<=t && t<=L1(2,i)
      counter=counter+1;
    end
  end
  Betti(dim+1)=counter;
end