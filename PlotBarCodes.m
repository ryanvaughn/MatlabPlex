function PlotBarCodes(L,tophom)
close all;
maxdim=max(L(4,:));maxdim=min(maxdim,tophom);
maxtime=max(1,max(L(2,:)));
for dim=0:maxdim
  f=find(L(4,:)==dim);num=length(f);L1=L(:,f);
  [~,r]=sort(L1(1,:),'descend');L1=L1(:,r);
  y=linspace(0,1,num+2);y=y(2:end-1);
  figure(dim+1);hold on
  for i=1:num
    plot([L1(1,i) L1(2,i)+.05],[y(i) y(i)],'LineWidth',5)
  end
  hold off
  axis([0 maxtime 0 1]);title(['Homology in dim ' num2str(dim)],'FontSize',20)
end


