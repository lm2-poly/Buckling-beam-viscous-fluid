
function cool
set(0,'DefaultFigureWindowStyle','docked')
format long
format compact

clear
clc
eta=2;

options=bvpset;


t=1.5

n1=0;
while t<30
    n1=n1+1;
    t0=t;
    t=t0*1.02;
    dt=t-t0;
    
    if exist('solprev')
        clear Sxyp
        Sxyp(:,1)=solprev.x';
        Sxyp(:,2)=solprev.y(5,:)';
        Sxyp(:,3)=solprev.y(6,:)';
        sol = bvp4c(@fsibeam,@beambcs,solprev,options,Sxyp,eta,dt,t);
        solprev=sol;
        
        xcenterold=xcenter;
        xcenter=trapz(sol.x(1,:),sol.y(5,:),2);
        vxcenter=(xcenter-xcenterold)/dt;
    else
        solinit = bvpinit(linspace(0,1,20),@inigshape);
        Sxyp(:,1)=solinit.x';
        Sxyp(:,2)=solinit.y(5,:)';
        Sxyp(:,3)=solinit.y(6,:)';
        sol = bvp4c(@fsibeam,@beambcs,solinit,options,Sxyp,eta,dt,t);
        solprev=sol;
        figure(4)
        plot(Sxyp(:,2)*t,Sxyp(:,3)*t,'-k')
        xlabel('x')
        ylabel('y')
        
        
        xcenterold=trapz(solinit.x(1,:),solinit.y(5,:),2);
        xcenter=trapz(sol.x(1,:),sol.y(5,:),2);
        vxcenter=(xcenter-xcenterold)/dt;
    end
    
    %xint=0:0.1:1
    %Sxint=deval(sol,xint)
    clear Sxysol
    Sxysol(:,1)=sol.x';
    Sxysol(:,2)=sol.y(5,:)';
    Sxysol(:,3)=sol.y(6,:)';
    
    if rem(n1,5)==0
        figure(4)
        hold on
        plot(Sxysol(:,2)*t,Sxysol(:,3)*t,'-k')
        axis equal
    end
    size(Sxysol)
    
    bedingenergy=0.5*trapz(sol.x(1,:),sol.y(1,:).*sol.y(1,:),2);
    reaccomp=sol.y(4,1);
    figure(5)
    subplot(2,1,1)
    hold on
    plot(t,reaccomp,'ok')
    ylabel('reaction')
    subplot(2,1,2)
    hold on
    plot(t,bedingenergy,'ok')
    ylabel('energy')
    xlabel('time')
    
    
    figure(6)
    subplot(2,1,1)
    loglog(t,reaccomp,'ok')
    hold on
    ylabel('reaction')
    subplot(2,1,2)
    loglog(t,bedingenergy,'ok')
    hold on
    ylabel('energy')
    xlabel('time')
    
    figure(7)
    subplot(4,1,1)
    plot(Sxysol(:,1),sol.y(1,:),'-k')
    ylabel('\theta')
    subplot(4,1,2)
    plot(Sxysol(:,1),sol.y(2,:),'-k')
    ylabel('M')
    subplot(4,1,3)
    plot(Sxysol(:,1),sol.y(3,:),'-k')
    ylabel('V')
    subplot(4,1,4)
    plot(Sxysol(:,1),sol.y(4,:),'-k')
    ylabel('T')
    xlabel('S')
    
    figure(8)
    subplot(2,1,1)
    plot(t,xcenter,'ok')
    hold on
    ylabel('x_center')
    subplot(2,1,2)
    plot(t,vxcenter,'ok')
    hold on
    ylabel('vx_center')
    xlabel('time')
    
end

end
% ------------------------------------------------------------
function dds = fsibeam(S,Y,Sxyp,eta,dt,t)

%Y(1)=theta
%Y(2)=M
%Y(3)=V
%Y(4)=T
%Y(5)=x
%Y(6)=y

theta=Y(1);
L=t;

xp=interp1(Sxyp(:,1),Sxyp(:,2),S,'linear');
yp=interp1(Sxyp(:,1),Sxyp(:,3),S,'linear');
dxdt=(Y(5)-xp)/dt;
dydt=(Y(6)-yp)/dt;
vs=L*dxdt*cos(theta)+L*dydt*sin(theta)+1;
vn=-L*dxdt*sin(theta)+L*dydt*cos(theta);

dds=[L*Y(2)
    L*Y(3)
    -eta*L*vn+Y(4)*L*Y(2)
    L*vs-Y(3)*L*Y(2)
    cos(theta)
    sin(theta)];
end
% ------------------------------------------------------------
function [ res ] = beambcs(Yl,Yr,thetap,eta,dt,t)

%Boundary conditions set to zero
res = [  Yl(1)
    Yr(2)
    Yr(3)
    Yr(4)
    Yl(5)
    Yl(6)];
end
% ------------------------------------------------------------
function shap0=inigshape(S)

thetatip=0.001;
t=1.5;
L=t;

shap0=[thetatip*S
    thetatip/L
    0
    L*(S-1)
    1/thetatip*sin(thetatip*S)
    1/thetatip*(1-cos(thetatip*S))];
end