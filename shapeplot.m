

load('shapesave2.mat')


figure(99)
hold on
%for np=316:10:346
nplotted=0;
for np=[4	9	14	18	23	28	32	37	42	47	51]
    plot(-shapesavey(np,2:end)+nplotted*1.3,-shapesavex(np,2:end),'-k')
    nplotted=nplotted+1
end
axis equal
axis off




figure(100)
%for np=316:10:346
nplotted=0;
hold on
allnp=[39 46 96 146 196 246 296]%[34 39 44 121 134 146 296];
for np=allnp;
    plot(-shapesavey(np,2:end)+nplotted*12,-shapesavex(np,2:end),'-b')
    nplotted=nplotted+1
end
%for np=[276 286 296]
%    plot(-shapesavey(np,2:end)+(nplotted-7)*12,-shapesavex(np,2:end)-12,'-b')
%    nplotted=nplotted+1
%end
axis equal
axis off
%ylim([-22 0])

%allnp=[21 46 71 96 121 146 171 196 221 246 271 296 321 346];
%allnp=[34 39 44 71 84 96 296];
alltimes=allnp/5+0.8;
nnp=length(allnp);

resusave=csvread('resusave4.csv');
figure(101)
plot(resusave(:,1),-resusave(:,6),'-b')
hold on
box off
for n1=1:nnp
    t=alltimes(n1);
    C=-interp1(resusave(:,1),resusave(:,6),t);
    plot(t,C,'ok','MarkerEdgeColor','b','MarkerFaceColor','w')
end
%xlim([0 70])
ylim([0 9])
xlabel('xlabel')
ylabel('ylabel')