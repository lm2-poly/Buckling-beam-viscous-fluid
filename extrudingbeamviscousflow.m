%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXTRUDINGBEAMVISCOUSFLOW.m SOLVES THE EQUATIONS OF A BEAM EXTRUDING IN A
%HIGHLY VISCOUS FLUID. THIS CODE IS USED TO GENERATE THE RESULTS PUBLISHED
%IN THE REFERENCE BELOW. 
%Gosselin, F. P., Neetzow, P., & Paak, M. (2014). Buckling of a beam
%extruded into highly viscous fluid. Physical Review E, 90(5), 052718.
%IN A NUTSHELL, THE CODE USES A BACKWARD EULER ALGORITHM TO INTEGRATE IN
%TIME THE LARGE-DEFLECTION DYNAMICS OF THE ROD EXTRUDED IN VISCOUS FLUID.
%AT EACH TIME STEP, THE BVP4C MATLAB FUNCTION IS CALLED TO SOLVE IMPLICITLY
%THE EQUATIONS OF THE BEAM ON THE LAGRAGIAN DOMAIN (THE LENGTH OF THE BEAM
%SCALED TO ALWAYS BE OF LENGTH 1).
%THIS CODE WAS DEVELOPPED BY FREDERICK P. GOSSELIN AT ECOLE POLYTECHNIQUE
%DE MONTREAL IN 2014.

function []=extrudingbeamviscousflow
%set(0,'DefaultFigureWindowStyle','docked')
format long
format compact

clear
clc

%VARIABLE INITIALISATION
eta=2;          %RATIO OF NORMAL TO PARALLEL FRICTION FORCE COEFFICIENTS
options=bvpset;
tzero=1;        %TIME (AND BEAM LENGTH) AT THE BEGINNING OF THE SIMULATION
t=tzero;
dt=0.2;         %TIME STEP SIZE
%FINAL TIME OF THE SIMULATION (BETTER TO USE A MULTIPLE OF 10 SINCE RESULTS
%ARE SAVED TO HARD DRIVE EVERY rem(t,10)==0 
tmax=70;        
resupeaks=zeros(1,4);
npts=100;       %NUMBER OF DISCRETIZATION POINTS ALONG THE LENGTH OF THE BEAM
shapesavex=zeros((tmax-tzero)/dt,npts+1);
shapesavey=zeros((tmax-tzero)/dt,npts+1);
tic             %tic FUNCTION USED TO EVALUATE THE DURATION OF THE COMPUTATION
nframes=0;      %NUMBER OF FRAMES IN THE ANIMATION MOVIE

%SIMULATION DURATION WHILE LOOP
n1=0;
while t<tmax
    n1=n1+1;
    t0=t;
    t=tzero+(n1-1)*dt
    %t=t0*1.02;
    %dt=t-t0;
    
    %CONDITIONNAL STATEMENT USED TO APPLY INITIAL CONDITIONS
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
        %FUNCTION BVPINIT ALLOWS FORMATTING THE INITIAL CONDITIONS
        solinit = bvpinit(linspace(0,1,20),@inigshape);
        Sxyp(:,1)=solinit.x';
        Sxyp(:,2)=solinit.y(5,:)';
        Sxyp(:,3)=solinit.y(6,:)';
        sol = bvp4c(@fsibeam,@beambcs,solinit,options,Sxyp,eta,dt,t);
        solprev=sol;
        
        %{
        figure(4)
        axis equal
        set(4,'Position',[50 50 800 480])
        xlim([-10 10])
        ylim([-12 0])
        hold on
        box on
        plot(Sxyp(:,3)*t,-Sxyp(:,2)*t,'color',[0.5 0.5 0.5],'linewidth',3)%'-k')
        xlabel('y')
        ylabel('x')
        
        %set(4,'Position',[50 50 800 480])
        set(findall(4,'type','text'),'fontSize',14,'fontWeight','bold')
        hold on
        %}
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
    
    
    
    x = linspace(0,1,npts);
    y = deval(sol,x);
    shapesavex(n1,1)=t;
    shapesavey(n1,1)=t;
    shapesavex(n1,2:npts+1)=y(5,:)'*t;
    shapesavey(n1,2:npts+1)=y(6,:)'*t;
    
    %figure(5)
    %hold on
    %plot(shapesavex(n1,2:npts+1),shapesavey(n1,2:npts+1),'-k')
    
    %COMMANDS TO CREATE AN ANIMATION MOVIE
    
    if rem(n1,4)==0
         figure(4)
        %axis equal
        %hold on
        %xlim([0 15])
        %ylim([-6 6])
        plot(Sxysol(:,3)*t,-Sxysol(:,2)*t,'color',[0 0 0],'linewidth',3)%'-k')
        %set(4,'Position',[50 50 1000 800])
        %set(findall(4,'type','text'),'fontSize',14,'fontWeight','bold')
        
        nframes=nframes+1;
        abc=text(-9,-10,strcat('time=',num2str(t)),'fontSize',14,'fontWeight','bold');
        M(nframes) = getframe(4);
        delete(abc);
        plot(Sxysol(:,3)*t,-Sxysol(:,2)*t,'color',[0.5 0.5 0.5],'linewidth',3)%'-k')
    end
    
    
    %{
    if rem(t,10)==0
        save('shapesave2.mat','shapesavex','shapesavey')
    end
    %}
    end

%{
for n11=1:15
    nframes=nframes+1
    M(nframes) = getframe(4);
end
movie2avi(M, 'beammodel-theta0=0.001-new4.avi', 'compression', 'none','fps',15);
%}

    
    
    toc

end
% ------------------------------------------------------------
function dds = fsibeam(S,Y,Sxyp,eta,dt,t)



%THE PROBLEM IS FORMULATED AS A SYSTEM OF 6 FIRST ORDER PDES
%DEFINITION OF THE VARIABLES USED
%Y(1)=theta     %LOCAL ANGLE OF THE BEAM
%Y(2)=M         %BENDING MOMENT (OR CURVATURE)
%Y(3)=V         %SHEAR FORCE
%Y(4)=T         %AXIAL TENSION FORCE
%Y(5)=x         %EULERIAN x COORDINATE
%Y(6)=y         %EULERIAN y COORDINATE

theta=Y(1);
L=t;            %INSTANTANEOUS LENGTH OF THE BEAM

xp=interp1(Sxyp(:,1),Sxyp(:,2),S,'pchip');
yp=interp1(Sxyp(:,1),Sxyp(:,3),S,'pchip');
dxdt=(Y(5)-xp)/dt;  %BACKWARD EULEUR SCHEME IS USED TO EVALUATE THE x VELOCITY
dydt=(Y(6)-yp)/dt;  %BACKWARD EULEUR SCHEME IS USED TO EVALUATE THE Y VELOCITY
vs=(L*dxdt)*cos(theta)+(L*dydt)*sin(theta)+1;
vn=-(L*dxdt)*sin(theta)+(L*dydt)*cos(theta);

%COMPUTATION OF THE DERIVATIVES
dds=[L*Y(2)
    L*Y(3)
    -eta*L*vn+Y(4)*L*Y(2)
    L*vs-Y(3)*L*Y(2)
    cos(theta)
    sin(theta)];
end
% ------------------------------------------------------------
function [ res ] = beambcs(Yl,Yr,thetap,eta,dt,t)

%BOUNDARY CONDITIONS SET TO ZERO (DIRICHLET)
%LHS = CLAMP BOUNDARY
%RHS = FREE END
res = [  Yl(1)  %THETA SET TO ZERO AT LHS OF DOMAIN
    Yr(2)       %M SET TO ZERO AT RHS OF DOMAIN
    Yr(3)       %V SET TO ZERO AT RHS OF DOMAIN 
    Yr(4)       %T SET TO ZERO AT RHS OF DOMAIN
    Yl(5)       %x SET TO ZERO AT LHS OF DOMAIN
    Yl(6)];     %y SET TO ZERO AT LHS OF DOMAIN
end
% ------------------------------------------------------------
%FUNCTION FOR APPLYING THE INITIAL CONDITIONS
function shap0=inigshape(S)
%THE INITIAL CONDITION USED IN THE PAPER IS THAT OF A UNIFORM CURVATURE
%ALONG THE LENGTH OF THE BEAM. A SINGLE VALUE OF THETA IS PRESCRIBED AT THE
%FREE END AND THE REST OF THE SHAPE FOLLOWS.
thetatip=.001;
t=1;
L=t;

shap0=[thetatip*S
    thetatip/L
    0
    L*(S-1)
    1/thetatip*sin(thetatip*S)
    1/thetatip*(1-cos(thetatip*S))];
end