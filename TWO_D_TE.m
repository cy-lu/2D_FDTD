%============2D FDTD plane wave  TEz incident===========
%============Total field and scatter field domain=======
%============ PML Boundary outside the scatter field====
%=======================================================
%Total field Ex(1:Nx-1,1:Ny-1),Ey(1:Nx,1:Ny),Hz(1:Nx-1,1:Ny-1)
%========================================================
clear all;
close all;
tic;                       
%=========================Parameters==================== 
nmax=1500;                   % max time steps
epso=8.85*10^(-12);          % permittivity of free space
muo=4*pi*10^(-7);            % permeability of free space
c=1/sqrt(epso*muo);          % speed of light
zo=sqrt(muo/epso);           
      
%=======================parameters of incident wave======
sigma = 0.4e-10;
m_offset = 4*sigma;
tau=2/6*10^(-9);             %reference to taflove textbook
to=tau;

%==============================================================                         
lamda=c/(2/tau);            % Wavelength
ds=lamda/30;                % grid size
dt=ds/(2*c);                % time step
nfreq=1;                    % nfreq
f=(1:nfreq)*1*10^9;         % frequenct range for calculating the optical response in the frequency domain                         % 跑的頻率相差間隔大小

Nx=200;                     %mesh number in the x direction
Ny=200;                     %mesh number in the y direction
%=========================================================
[x,y]=meshgrid(1:Ny,1:Nx);
%=========================================================
%=======matrix of the electric/magnetic field=============
Ex=zeros(Nx,Ny);
Ey=zeros(Nx,Ny);
Hz=zeros(Nx,Ny);
epsx=ones(Nx,Ny)*epso;
epsy=ones(Nx,Ny)*epso;
muz=ones(Nx,Ny)*muo;
cex=ones(Nx,Ny)*dt./(epsx*ds);
cey=ones(Nx,Ny)*dt./(epsy*ds);
dhz=ones(Nx,Ny)*dt./(muz*ds);
%========================material properties================
%=======================a core-shell nanoparticle===========
cent_x=Nx/2;
cent_y=Ny/2;
radius=0.3*c/f(1);   %outer radius
radius2=0.25*c/f(1); %innter radius
  for p = 1:Nx
    for q = 1:Ny
        if ( sqrt( ((p-1/2-(cent_x-1))*ds)^2 +  ((q-1-(cent_y-1))*ds)^2 ) <=  radius )
            cex(p,q)=dt/(epso*ds*4);
     end
        if ( sqrt( ((p-1-(cent_x-1))*ds)^2 +  ((q-1/2-(cent_y-1))*ds)^2 ) <=  radius ) 
            cey(p,q)=dt/(epso*ds*4);
        end
    end
  end
  
    for p = 1:Nx
    for q = 1:Ny
     
       if ( sqrt( ((p-1/2-(cent_x-1))*ds)^2 +  ((q-1-(cent_y-1))*ds)^2 ) <=  radius2 )
           cex(p,q)=dt/(epso*ds);
     end
        if ( sqrt( ((p-1-(cent_x-1))*ds)^2 +  ((q-1/2-(cent_y-1))*ds)^2 ) <=  radius2 ) 
          cey(p,q)=dt/(epso*ds);
        end
        end

  end

%========================PML parameters=====================
pm=12;   % how many mesh as the PML thickness 
m=4;     % PML m paramtere ~3-4 in general
sigmax=0.8*(m+1)/(zo*ds);  % sigmax
ss=Nx+2*(pm-1);            % up & down mesh number in the x direction

%======================SFTF parameter============================
sf=10;             % SF/TF back boundary
sff=Nx-10;         % SF/TF front boundary
sfu=Ny-10;         % SF/TF up boundary
sfd=10;            % SF/TF down boundary

Nxinc=Nx;
Nyinc=Ny;
ssinc=Nxinc+2*(pm-1);  

Exinc=zeros(Nxinc,Nyinc);
Eyinc=zeros(Nxinc,Nyinc);
Hzinc=zeros(Nxinc,Nyinc);
%===========================================================

%====================right PML parameter====================
Eyr=zeros(pm,Ny);
Exr=zeros(pm,Ny);
Hzxr=zeros(pm,Ny);
Hzyr=zeros(pm,Ny);

Eyrinc=zeros(pm,Nyinc);
Exrinc=zeros(pm,Nyinc);
Hzxrinc=zeros(pm,Nyinc);
Hzyrinc=zeros(pm,Nyinc);

caexr=zeros(pm,Ny);
cbexr=zeros(pm,Ny);
caeyr=zeros(pm,Ny);
cbeyr=zeros(pm,Ny);
dahzxr=zeros(pm,Ny);
dbhzxr=zeros(pm,Ny);
dahzyr=zeros(pm,Ny);
dbhzyr=zeros(pm,Ny);
for i= 1:pm
sigmxr=sigmax*(((i-1)/pm)^(m));
sigmxrs=muo/epso*sigmax*(((i-0.5)/pm)^(m));    
caexr(i,1:Ny)=1;
cbexr(i,1:Ny)=dt/(epso*ds);
caeyr(i,1:Ny)=(1-sigmxr*dt/(2*epso))/(1+sigmxr*dt/(2*epso));
cbeyr(i,1:Ny)=dt/(epso*ds)/(1+sigmxr*dt/(2*epso));
dahzxr(i,1:Ny)=(1-sigmxrs*dt/2/muo)/(1+sigmxrs*dt/2/muo);
dbhzxr(i,1:Ny)=dt/muo/ds/(1+sigmxrs*dt/2/muo);
dahzyr(i,1:Ny)=1;
dbhzyr(i,1:Ny)=dt/muo/ds;
end
%=====================left PML parameter====================
Eyl=zeros(pm,Ny);
Exl=zeros(pm,Ny);
Hzxl=zeros(pm,Ny);
Hzyl=zeros(pm,Ny);

Eylinc=zeros(pm,Nyinc);
Exlinc=zeros(pm,Nyinc);
Hzxlinc=zeros(pm,Nyinc);
Hzylinc=zeros(pm,Nyinc);

caexl=zeros(pm,Ny);
cbexl=zeros(pm,Ny);
caeyl=zeros(pm,Ny);
cbeyl=zeros(pm,Ny);
dahzxl=zeros(pm,Ny);
dbhzxl=zeros(pm,Ny);
dahzyl=zeros(pm,Ny);
dbhzyl=zeros(pm,Ny);
for i= 1:pm
sigmxl=sigmax*(((pm-i)/pm)^(m));
sigmxls=muo/epso*sigmax*(((pm-i-0.5)/pm)^(m));    
caexl(i,1:Ny)=1;
cbexl(i,1:Ny)=dt/(epso*ds);
caeyl(i,1:Ny)=(2*epso-sigmxl*dt)/(2*epso+sigmxl*dt);
cbeyl(i,1:Ny)=dt/(epso*ds)/(1+sigmxl*dt/(2*epso));
dahzxl(i,1:Ny)=(1-sigmxls*dt/(2*muo))/(1+sigmxls*dt/(2*muo));
dbhzxl(i,1:Ny)=dt/(muo*ds*(1+sigmxls*dt/2/muo));
dahzyl(i,1:Ny)=1;
dbhzyl(i,1:Ny)=dt/(muo*ds);
end
%=======================up PML parameter====================
Eyu=zeros(ss,pm);
Exu=zeros(ss,pm);
Hzxu=zeros(ss,pm);
Hzyu=zeros(ss,pm);

Eyuinc=zeros(ssinc,pm);
Exuinc=zeros(ssinc,pm);
Hzxuinc=zeros(ssinc,pm);
Hzyuinc=zeros(ssinc,pm);

caexu=zeros(ss,pm);
cbexu=zeros(ss,pm);
caeyu=zeros(ss,pm);
cbeyu=zeros(ss,pm);
dahzxu=zeros(ss,pm);
dbhzxu=zeros(ss,pm);
dahzyu=zeros(ss,pm);
dbhzyu=zeros(ss,pm);
for j= 1:pm
sigmyu=sigmax*(((j-1)/pm)^(m));
sigmyus=muo/epso*sigmax*(((j-0.5)/pm)^(m));    
caexu(1:ss,j)=(1-sigmyu*dt/(2*epso))/(1+sigmyu*dt/(2*epso));
cbexu(1:ss,j)=dt/(ds*epso)/(1+sigmyu*dt/(2*epso));
caeyu(1:ss,j)=1;
cbeyu(1:ss,j)=dt/(epso*ds);
dahzxu(1:ss,j)=1;
dbhzxu(1:ss,j)=dt/(muo*ds);
dahzyu(1:ss,j)=(1-sigmyus*dt/(2*muo))/(1+sigmyus*dt/(2*muo));
dbhzyu(1:ss,j)=dt/muo/ds/(1+sigmyus*dt/(2*muo));
end
for i=pm+Nx-1:ss-1
sigmxr=sigmax*(((i-(pm+Nx-1))/pm)^(m));
sigmxrs=muo/epso*sigmax*(((i+0.5-(pm+Nx-1))/pm)^(m));  
caeyu(i,1:pm)=(1-sigmxr*dt/(2*epso))/(1+sigmxr*dt/(2*epso));
cbeyu(i,1:pm)=dt/(epso*ds)/(1+sigmxr*dt/(2*epso));
dahzxu(i,1:pm)=(1-sigmxrs*dt/2/muo)/(1+sigmxrs*dt/2/muo);
dbhzxu(i,1:pm)=dt/muo/ds/(1+sigmxrs*dt/2/muo);
end
for i=1:pm-1
sigmxl=sigmax*(((pm-i)/pm)^(m));
sigmxls=muo/epso*sigmax*(((pm-i-0.5)/pm)^(m));  
caeyu(i,1:pm)=(1-sigmxl*dt/(2*epso))/(1+sigmxl*dt/(2*epso));
cbeyu(i,1:pm)=dt/(epso*ds)/(1+sigmxl*dt/(2*epso));
dahzxu(i,1:pm)=(1-sigmxls*dt/2/muo)/(1+sigmxls*dt/2/muo);
dbhzxu(i,1:pm)=dt/muo/ds/(1+sigmxls*dt/2/muo);
end
%==
%== We also need to set up a PML boundary for incident light
%==
for j= 1:pm
sigmyu=sigmax*(((j-1)/pm)^(m));
sigmyus=muo/epso*sigmax*(((j-0.5)/pm)^(m));    
caexuinc(1:ssinc,j)=(1-sigmyu*dt/(2*epso))/(1+sigmyu*dt/(2*epso));
cbexuinc(1:ssinc,j)=dt/(ds*epso)/(1+sigmyu*dt/(2*epso));
caeyuinc(1:ssinc,j)=1;
cbeyuinc(1:ssinc,j)=dt/(epso*ds);
dahzxuinc(1:ssinc,j)=1;
dbhzxuinc(1:ssinc,j)=dt/(muo*ds);
dahzyuinc(1:ssinc,j)=(1-sigmyus*dt/(2*muo))/(1+sigmyus*dt/(2*muo));
dbhzyuinc(1:ssinc,j)=dt/muo/ds/(1+sigmyus*dt/(2*muo));
end
for i=pm+Nxinc-1:ssinc-1
sigmxr=sigmax*(((i-(pm+Nxinc-1))/pm)^(m));
sigmxrs=muo/epso*sigmax*(((i+0.5-(pm+Nxinc-1))/pm)^(m));  
caeyuinc(i,1:pm)=(1-sigmxr*dt/(2*epso))/(1+sigmxr*dt/(2*epso));
cbeyuinc(i,1:pm)=dt/(epso*ds)/(1+sigmxr*dt/(2*epso));
dahzxuinc(i,1:pm)=(1-sigmxrs*dt/2/muo)/(1+sigmxrs*dt/2/muo);
dbhzxuinc(i,1:pm)=dt/muo/ds/(1+sigmxrs*dt/2/muo);
end
for i=1:pm-1
sigmxl=sigmax*(((pm-i)/pm)^(m));
sigmxls=muo/epso*sigmax*(((pm-i-0.5)/pm)^(m));  
caeyuinc(i,1:pm)=(1-sigmxl*dt/(2*epso))/(1+sigmxl*dt/(2*epso));
cbeyuinc(i,1:pm)=dt/(epso*ds)/(1+sigmxl*dt/(2*epso));
dahzxuinc(i,1:pm)=(1-sigmxls*dt/2/muo)/(1+sigmxls*dt/2/muo);
dbhzxuinc(i,1:pm)=dt/muo/ds/(1+sigmxls*dt/2/muo);
end

%=======================down PML parameters==========================
Eyd=zeros(ss,pm);
Exd=zeros(ss,pm);
Hzxd=zeros(ss,pm);
Hzyd=zeros(ss,pm);

Eydinc=zeros(ssinc,pm);
Exdinc=zeros(ssinc,pm);
Hzxdinc=zeros(ssinc,pm);
Hzydinc=zeros(ssinc,pm);

caexd=zeros(ss,pm);
cbexd=zeros(ss,pm);
caeyd=zeros(ss,pm);
cbeyd=zeros(ss,pm);
dahzxd=zeros(ss,pm);
dbhzxd=zeros(ss,pm);
dahzyd=zeros(ss,pm);
dbhzyd=zeros(ss,pm);
for j= 1:pm
sigmyd=sigmax*(((j-1)/pm)^(m));
sigmyds=muo/epso*sigmax*(((j-0.5)/pm)^(m));    
caexd(1:ss,j)=(1-sigmyd*dt/(2*epso))/(1+sigmyd*dt/(2*epso));
cbexd(1:ss,j)=dt/(ds*epso)/(1+sigmyd*dt/(2*epso));
caeyd(1:ss,j)=1;
cbeyd(1:ss,j)=dt/(epso*ds);
dahzxd(1:ss,j)=1;
dbhzxd(1:ss,j)=dt/(muo*ds);
dahzyd(1:ss,j)=(1-sigmyds*dt/(2*muo))/(1+sigmyds*dt/(2*muo));
dbhzyd(1:ss,j)=dt/muo/ds/(1+sigmyds*dt/(2*muo));
end
for i=pm+Nx-1:ss-1
sigmxr=sigmax*(((i-(pm+Nx-1))/pm)^(m));
sigmxrs=muo/epso*sigmax*(((i+0.5-(pm+Nx-1))/pm)^(m));  
caeyd(i,1:pm)=(1-sigmxr*dt/(2*epso))/(1+sigmxr*dt/(2*epso));
cbeyd(i,1:pm)=dt/(epso*ds)/(1+sigmxr*dt/(2*epso));
dahzxd(i,1:pm)=(1-sigmxrs*dt/2/muo)/(1+sigmxrs*dt/2/muo);
dbhzxd(i,1:pm)=dt/muo/ds/(1+sigmxrs*dt/2/muo);
end
for i=1:pm-1
sigmxl=sigmax*(((pm-i)/pm)^(m));
sigmxls=muo/epso*sigmax*(((pm-i-0.5)/pm)^(m));  
caeyd(i,1:pm)=(1-sigmxl*dt/(2*epso))/(1+sigmxl*dt/(2*epso));
cbeyd(i,1:pm)=dt/(epso*ds)/(1+sigmxl*dt/(2*epso));
dahzxd(i,1:pm)=(1-sigmxls*dt/2/muo)/(1+sigmxls*dt/2/muo);
dbhzxd(i,1:pm)=dt/muo/ds/(1+sigmxls*dt/2/muo);
end
%== incident light 
for j= 1:pm
sigmyd=sigmax*(((j-1)/pm)^(m));
sigmyds=muo/epso*sigmax*(((j-0.5)/pm)^(m));    
caexdinc(1:ssinc,j)=(1-sigmyd*dt/(2*epso))/(1+sigmyd*dt/(2*epso));
cbexdinc(1:ssinc,j)=dt/(ds*epso)/(1+sigmyd*dt/(2*epso));
caeydinc(1:ssinc,j)=1;
cbeydinc(1:ssinc,j)=dt/(epso*ds);
dahzxdinc(1:ssinc,j)=1;
dbhzxdinc(1:ssinc,j)=dt/(muo*ds);
dahzydinc(1:ssinc,j)=(1-sigmyds*dt/(2*muo))/(1+sigmyds*dt/(2*muo));
dbhzydinc(1:ssinc,j)=dt/muo/ds/(1+sigmyds*dt/(2*muo));
end
for i=pm+Nxinc-1:ssinc-1
sigmxr=sigmax*(((i-(pm+Nxinc-1))/pm)^(m));
sigmxrs=muo/epso*sigmax*(((i+0.5-(pm+Nxinc-1))/pm)^(m));  
caeydinc(i,1:pm)=(1-sigmxr*dt/(2*epso))/(1+sigmxr*dt/(2*epso));
cbeydinc(i,1:pm)=dt/(epso*ds)/(1+sigmxr*dt/(2*epso));
dahzxdinc(i,1:pm)=(1-sigmxrs*dt/2/muo)/(1+sigmxrs*dt/2/muo);
dbhzxdinc(i,1:pm)=dt/muo/ds/(1+sigmxrs*dt/2/muo);
end
for i=1:pm-1
sigmxl=sigmax*(((pm-i)/pm)^(m));
sigmxls=muo/epso*sigmax*(((pm-i-0.5)/pm)^(m));  
caeydinc(i,1:pm)=(1-sigmxl*dt/(2*epso))/(1+sigmxl*dt/(2*epso));
cbeydinc(i,1:pm)=dt/(epso*ds)/(1+sigmxl*dt/(2*epso));
dahzxdinc(i,1:pm)=(1-sigmxls*dt/2/muo)/(1+sigmxls*dt/2/muo);
dbhzxdinc(i,1:pm)=dt/muo/ds/(1+sigmxls*dt/2/muo);
end

%======================PML parameters end============================


%======================near to far field parameters==================

SN=zeros(nfreq,2);                       % cos and sin of discrete fourier transfrom
RD=zeros(nfreq,2);                       % phase and amplitude of light sources
ntff=5;                                  % near to far location

Hzefr=zeros(nfreq,Ny-ntff-1,2);          % record scatter field field Hz
Hzefl=zeros(nfreq,Ny-ntff-1,2);
Hzefu=zeros(nfreq,Nx-ntff,2);
Hzefd=zeros(nfreq,Nx-ntff,2);

Eyefr=zeros(nfreq,Ny-ntff-1,2);          % record scatter field field Ey
Eyefl=zeros(nfreq,Ny-ntff-1,2);
Eyefu=zeros(nfreq,Nx-ntff,2);
Eyefd=zeros(nfreq,Nx-ntff,2);

Exefr=zeros(nfreq,Ny-ntff-1,2);          % record scatter field field Ex
Exefl=zeros(nfreq,Ny-ntff-1,2);
Exefu=zeros(nfreq,Nx-ntff,2);
Exefd=zeros(nfreq,Nx-ntff,2);

RCS=zeros(nfreq,360);                    % save RCS matrix 
tscs=zeros(nfreq,1);                     % save all frequencies at the tscs matrix
%=========================time step===========================
for n=1:nmax
%=======================SF/TF record point====================
Eycc=Ey(sf,sfd+1:sfu-1);     % back
Eyff=Ey(sff,sfd+1:sfu-1);    % front
Exuu=Ex(sf:sff-1,sfu);       % up
Exdd=Ex(sf:sff-1,sfd+1);     % down
%=============================================================
%===========================SF/TF light source================
% plane wave
Hzinc(1,1:Nyinc-1)=exp(-(4*pi*((n*dt-to)/tau)^2));

Hzxuinc(pm,1:pm-1)=exp(-(4*pi*((n*dt-to)/tau)^2));

Hzxdinc(pm,1:pm-1)=exp(-(4*pi*((n*dt-to)/tau)^2));


%Hzinc(1,Ny-5)=exp(-(4*pi*((n*dt-to)/tau)^2));

%=====================Exinc Eyinc update zone=====================
Exinc(1:Nxinc-1,2:Nyinc-1)=Exinc(1:Nxinc-1,2:Nyinc-1)+dt/(epso*ds)*(Hzinc(1:Nxinc-1,2:Nyinc-1)-Hzinc(1:Nxinc-1,1:Nyinc-2));

Eyinc(2:Nxinc-1,1:Nyinc-1)=Eyinc(2:Nxinc-1,1:Nyinc-1)-dt/(epso*ds)*(Hzinc(2:Nxinc-1,1:Nyinc-1)-Hzinc(1:Nxinc-2,1:Nyinc-1));

%========================update zone Eyuinc PML up & Eydinc PML down


Eyuinc(2:ssinc-1,1:pm-1)=caeyuinc(2:ssinc-1,1:pm-1).*Eyuinc(2:ssinc-1,1:pm-1)-cbeyuinc(2:ssinc-1,1:pm-1).*( (Hzxuinc(2:ssinc-1,1:pm-1)+Hzyuinc(2:ssinc-1,1:pm-1) ) - (Hzxuinc(1:ssinc-2,1:pm-1)+Hzyuinc(1:ssinc-2,1:pm-1)) ); 

Eydinc(2:ssinc-1,1:pm-1)=caeydinc(2:ssinc-1,1:pm-1).*Eydinc(2:ssinc-1,1:pm-1)-cbeydinc(2:ssinc-1,1:pm-1).*( (Hzxdinc(2:ssinc-1,1:pm-1)+Hzydinc(2:ssinc-1,1:pm-1) ) - (Hzxdinc(1:ssinc-2,1:pm-1)+Hzydinc(1:ssinc-2,1:pm-1)) ); 

%   update zone Exuinc PML up  & Exdinc PML down
%up
Exuinc(1:ssinc-1,2:pm-1)=caexuinc(1:ssinc-1,2:pm-1).*Exuinc(1:ssinc-1,2:pm-1)+cbexuinc(1:ssinc-1,2:pm-1).*(Hzxuinc(1:ssinc-1,2:pm-1)+Hzyuinc(1:ssinc-1,2:pm-1)-Hzxuinc(1:ssinc-1,1:pm-2)-Hzyuinc(1:ssinc-1,1:pm-2));

Exinc(1:Nxinc-1,Nyinc)=Exinc(1:Nxinc-1,Nyinc) + dt./(epsx(1:Nxinc-1,Nyinc)*ds).*(Hzxuinc(pm:ssinc-pm,1)+Hzyuinc(pm:ssinc-pm,1)-Hzinc(1:Nxinc-1,Nyinc-1));

%down
Exdinc(1:ssinc-1,2:pm-1)=caexdinc(1:ssinc-1,2:pm-1).*Exdinc(1:ssinc-1,2:pm-1)+cbexdinc(1:ssinc-1,2:pm-1).*(Hzxdinc(1:ssinc-1,1:pm-2)+Hzydinc(1:ssinc-1,1:pm-2)-Hzxdinc(1:ssinc-1,2:pm-1)-Hzydinc(1:ssinc-1,2:pm-1));

Exinc(1:Nxinc-1,1)=Exinc(1:Nxinc-1,1) + dt./(epsx(1:Nxinc-1,1)*ds).*(Hzinc(1:Nxinc-1,1)-Hzxdinc(pm:ssinc-pm,1)-Hzydinc(pm:ssinc-pm,1));


%========================update zone Eyrinc PML right & Eylinc PML left
% right
Eyrinc(2:pm-1,1:Nyinc-1)=caeyr(2:pm-1,1:Nyinc-1).*Eyrinc(2:pm-1,1:Nyinc-1)-cbeyr(2:pm-1,1:Nyinc-1).*( (Hzxrinc(2:pm-1,1:Nyinc-1)+Hzyrinc(2:pm-1,1:Nyinc-1) ) - (Hzxrinc(1:pm-2,1:Nyinc-1)+Hzyrinc(1:pm-2,1:Nyinc-1)) ); 

Eyinc(Nxinc,1:Nyinc-1)=Eyinc(Nxinc,1:Nyinc-1)-dt/(epso*ds)*(Hzxrinc(1,1:Nyinc-1)+Hzyrinc(1,1:Nyinc-1)-Hzinc(Nxinc-1,1:Nyinc-1));

Exrinc(1:pm-1,2:Nyinc-1)=caexr(1:pm-1,2:Nyinc-1).*Exrinc(1:pm-1,2:Nyinc-1)+cbexr(1:pm-1,2:Nyinc-1).*(Hzxrinc(1:pm-1,2:Nyinc-1)+Hzyrinc(1:pm-1,2:Nyinc-1)-Hzxrinc(1:pm-1,1:Nyinc-2)-Hzyrinc(1:pm-1,1:Nyinc-2));

%  incident PML up and right
Exrinc(1:pm-1,Nyinc)=caexr(1:pm-1,Nyinc).*Exrinc(1:pm-1,Nyinc)+cbexr(1:pm-1,Nyinc).*(Hzxuinc(ssinc-pm+1:ssinc-1,1)+Hzyuinc(ssinc-pm+1:ssinc-1,1) -Hzxrinc(1:pm-1,Nyinc-1)-Hzyrinc(1:pm-1,Nyinc-1));

%   incident PML down and right
Exrinc(1:pm-1,1)=caexr(1:pm-1,1).*Exrinc(1:pm-1,1)+cbexr(1:pm-1,1).*(Hzxrinc(1:pm-1,1)+Hzyrinc(1:pm-1,1)-Hzxdinc(ssinc-pm+1:ssinc-1,1)-Hzydinc(ssinc-pm+1:ssinc-1,1) );

%   left
Eylinc(2:pm-1,1:Nyinc-1)=caeyl(2:pm-1,1:Nyinc-1).*Eylinc(2:pm-1,1:Nyinc-1)-cbeyl(2:pm-1,1:Nyinc-1).*( (Hzxlinc(2:pm-1,1:Nyinc-1)+Hzylinc(2:pm-1,1:Nyinc-1) ) - (Hzxlinc(1:pm-2,1:Nyinc-1)+Hzylinc(1:pm-2,1:Nyinc-1)) ); 

Eyinc(1,1:Nyinc-1)=Eyinc(1,1:Nyinc-1)-dt/(epso*ds)*(Hzinc(1,1:Nyinc-1) -Hzxlinc(pm-1,1:Nyinc-1)-Hzylinc(pm-1,1:Nyinc-1)  );

Exlinc(1:pm-1,2:Nyinc-1)=caexl(1:pm-1,2:Nyinc-1).*Exlinc(1:pm-1,2:Nyinc-1)+cbexl(1:pm-1,2:Nyinc-1).*(Hzxlinc(1:pm-1,2:Nyinc-1)+Hzylinc(1:pm-1,2:Nyinc-1)-Hzxlinc(1:pm-1,1:Nyinc-2)-Hzylinc(1:pm-1,1:Nyinc-2));

%  incident PML up
Exlinc(1:pm-1,Nyinc)=caexl(1:pm-1,Nyinc).*Exlinc(1:pm-1,Nyinc)+cbexl(1:pm-1,Nyinc).*(Hzxuinc(1:pm-1,1)+Hzyuinc(1:pm-1,1)-Hzxlinc(1:pm-1,Nyinc-1)-Hzylinc(1:pm-1,Nyinc-1));

%  incident PML down
Exlinc(1:pm-1,1)=caexl(1:pm-1,1).*Exlinc(1:pm-1,1)+cbexl(1:pm-1,1).*(Hzxlinc(1:pm-1,1)+Hzylinc(1:pm-1,1)-Hzxdinc(1:pm-1,1)-Hzydinc(1:pm-1,1));

%=====================update eq Ex & Ey=======================

Ex(1:Nx-1,2:Ny-1)=Ex(1:Nx-1,2:Ny-1)+cex(1:Nx-1,2:Ny-1).*(Hz(1:Nx-1,2:Ny-1)-Hz(1:Nx-1,1:Ny-2));

Ey(2:Nx-1,1:Ny-1)=Ey(2:Nx-1,1:Ny-1)-cey(2:Nx-1,1:Ny-1).*(Hz(2:Nx-1,1:Ny-1)-Hz(1:Nx-2,1:Ny-1));

%==interface Ey====Ey(sf,sfd+1:sfu-1) as scattered field
%==interface Ey====Ey(sff,sfd+1:sfu-1) as scattered field
%== Ey(sf,1:Ny-1) and Hz(sf,1:Ny-1)

Ey(sf,sfd+1:sfu-1)=Eycc-cex(sf,sfd+1:sfu-1).*( Hz(sf,sfd+1:sfu-1)-Hzinc(sf,sfd+1:sfu-1)-Hz(sf-1,sfd+1:sfu-1) );   % back

Ey(sff,sfd+1:sfu-1)=Eyff-cex(sff,sfd+1:sfu-1).*( Hz(sff,sfd+1:sfu-1)+Hzinc(sff-1,sfd+1:sfu-1)-Hz(sff-1,sfd+1:sfu-1) ); % front

%===sftf up Ex  among Ex(sf:sff-1,sfu) is total field
Ex(sf:sff-1,sfu)=Exuu+cex(sf:sff-1,sfu).*(Hz(sf:sff-1,sfu)-Hz(sf:sff-1,sfu-1)+Hzinc(sf:sff-1,sfu));

%===sftf down Ex  among Ex(sf:sff-1,sfd+1) is total field
Ex(sf:sff-1,sfd+1)=Exdd+cex(sf:sff-1,sfd+1).*(Hz(sf:sff-1,sfd+1)-Hz(sf:sff-1,sfd)-Hzinc(sf:sff-1,sfd));

%========================Ey PML zone===========================
%============================================================

%========================Ey PML right========================
Eyr(2:pm-1,1:Ny-1)=caeyr(2:pm-1,1:Ny-1).*Eyr(2:pm-1,1:Ny-1)-cbeyr(2:pm-1,1:Ny-1).*( (Hzxr(2:pm-1,1:Ny-1)+Hzyr(2:pm-1,1:Ny-1) ) - (Hzxr(1:pm-2,1:Ny-1)+Hzyr(1:pm-2,1:Ny-1)) ); 

Ey(Nx,1:Ny-1)=Ey(Nx,1:Ny-1)-dt/(epso*ds)*(Hzxr(1,1:Ny-1)+Hzyr(1,1:Ny-1)-Hz(Nx-1,1:Ny-1));
%============================================================

%========================Ey PML left=========================
Eyl(2:pm-1,1:Ny-1)=caeyl(2:pm-1,1:Ny-1).*Eyl(2:pm-1,1:Ny-1)-cbeyl(2:pm-1,1:Ny-1).*( (Hzxl(2:pm-1,1:Ny-1)+Hzyl(2:pm-1,1:Ny-1) ) - (Hzxl(1:pm-2,1:Ny-1)+Hzyl(1:pm-2,1:Ny-1)) ); 

Ey(1,1:Ny-1)=Ey(1,1:Ny-1)-dt/(epso*ds)*(Hz(1,1:Ny-1) -Hzxl(pm-1,1:Ny-1)-Hzyl(pm-1,1:Ny-1)  );

%========================Ey PML up===========================


Eyu(2:ss-1,1:pm-1)=caeyu(2:ss-1,1:pm-1).*Eyu(2:ss-1,1:pm-1)-cbeyu(2:ss-1,1:pm-1).*( (Hzxu(2:ss-1,1:pm-1)+Hzyu(2:ss-1,1:pm-1) ) - (Hzxu(1:ss-2,1:pm-1)+Hzyu(1:ss-2,1:pm-1)) ); 


%========================Ey PML down==========================


Eyd(2:ss-1,1:pm-1)=caeyd(2:ss-1,1:pm-1).*Eyd(2:ss-1,1:pm-1)-cbeyd(2:ss-1,1:pm-1).*( (Hzxd(2:ss-1,1:pm-1)+Hzyd(2:ss-1,1:pm-1) ) - (Hzxd(1:ss-2,1:pm-1)+Hzyd(1:ss-2,1:pm-1)) ); 


%=============================================================
%=======================Ey PML Zone end========================


%========================Ex PML Zone end===========================
%============================================================

%========================Ex PML right========================
Exr(1:pm-1,2:Ny-1)=caexr(1:pm-1,2:Ny-1).*Exr(1:pm-1,2:Ny-1)+cbexr(1:pm-1,2:Ny-1).*(Hzxr(1:pm-1,2:Ny-1)+Hzyr(1:pm-1,2:Ny-1)-Hzxr(1:pm-1,1:Ny-2)-Hzyr(1:pm-1,1:Ny-2));

%  boundary at the upside
Exr(1:pm-1,Ny)=caexr(1:pm-1,Ny).*Exr(1:pm-1,Ny)+cbexr(1:pm-1,Ny).*(Hzxu(ss-pm+1:ss-1,1)+Hzyu(ss-pm+1:ss-1,1) -Hzxr(1:pm-1,Ny-1)-Hzyr(1:pm-1,Ny-1));

%  boundary at the downside
Exr(1:pm-1,1)=caexr(1:pm-1,1).*Exr(1:pm-1,1)+cbexr(1:pm-1,1).*(Hzxr(1:pm-1,1)+Hzyr(1:pm-1,1)-Hzxd(ss-pm+1:ss-1,1)-Hzyd(ss-pm+1:ss-1,1) );
%========================Ex PML left=========================
Exl(1:pm-1,2:Ny-1)=caexl(1:pm-1,2:Ny-1).*Exl(1:pm-1,2:Ny-1)+cbexl(1:pm-1,2:Ny-1).*(Hzxl(1:pm-1,2:Ny-1)+Hzyl(1:pm-1,2:Ny-1)-Hzxl(1:pm-1,1:Ny-2)-Hzyl(1:pm-1,1:Ny-2));

%  boundary at the upside
Exl(1:pm-1,Ny)=caexl(1:pm-1,Ny).*Exl(1:pm-1,Ny)+cbexl(1:pm-1,Ny).*(Hzxu(1:pm-1,1)+Hzyu(1:pm-1,1)-Hzxl(1:pm-1,Ny-1)-Hzyl(1:pm-1,Ny-1));

%  boundary at the downside
Exl(1:pm-1,1)=caexl(1:pm-1,1).*Exl(1:pm-1,1)+cbexl(1:pm-1,1).*(Hzxl(1:pm-1,1)+Hzyl(1:pm-1,1)-Hzxd(1:pm-1,1)-Hzyd(1:pm-1,1));
%========================Ex PML up===========================

Exu(1:ss-1,2:pm-1)=caexu(1:ss-1,2:pm-1).*Exu(1:ss-1,2:pm-1)+cbexu(1:ss-1,2:pm-1).*(Hzxu(1:ss-1,2:pm-1)+Hzyu(1:ss-1,2:pm-1)-Hzxu(1:ss-1,1:pm-2)-Hzyu(1:ss-1,1:pm-2));

Ex(1:Nx-1,Ny)=Ex(1:Nx-1,Ny) + dt./(epsx(1:Nx-1,Ny)*ds).*(Hzxu(pm:ss-pm,1)+Hzyu(pm:ss-pm,1)-Hz(1:Nx-1,Ny-1));

%========================Ex PML down===========================


Exd(1:ss-1,2:pm-1)=caexd(1:ss-1,2:pm-1).*Exd(1:ss-1,2:pm-1)+cbexd(1:ss-1,2:pm-1).*(Hzxd(1:ss-1,1:pm-2)+Hzyd(1:ss-1,1:pm-2)-Hzxd(1:ss-1,2:pm-1)-Hzyd(1:ss-1,2:pm-1));

Ex(1:Nx-1,1)=Ex(1:Nx-1,1) + dt./(epsx(1:Nx-1,1)*ds).*(Hz(1:Nx-1,1)-Hzxd(pm:ss-pm,1)-Hzyd(pm:ss-pm,1));

%=========================Ex PML end========================
%============================================================

%=========================update eq Hz=======================

%=========================SF/TF record matrix========================
Hzcc=Hz(sf,sfd+1:sfu-1);     % back

Hzff=Hz(sff-1,sfd+1:sfu-1);  % front

Hzuu=Hz(sf:sff-1,sfu);       % up

Hzdd=Hz(sf:sff-1,sfd);       % down
%=============================================================

Hz(1:Nx-1,1:Ny-1)=Hz(1:Nx-1,1:Ny-1)+dhz(1:Nx-1,1:Ny-1).*(Ex(1:Nx-1,2:Ny)-Ex(1:Nx-1,1:Ny-1)+Ey(1:Nx-1,1:Ny-1)-Ey(2:Nx,1:Ny-1));

%============================================================
%=========================update SF/TF Hz=====================

Hzinc(1:Nxinc-1,1:Nyinc-1)=Hzinc(1:Nxinc-1,1:Nyinc-1)+dt/(muo*ds)*(Exinc(1:Nxinc-1,2:Nyinc)-Exinc(1:Nxinc-1,1:Nyinc-1)+Eyinc(1:Nxinc-1,1:Nyinc-1)-Eyinc(2:Nxinc,1:Nyinc-1));


%  Hz(sf,sfd+1:sfu-1) as total field
%  Hz(sff-1,sfd+1:sfu-1) as total field

Hz(sf,sfd+1:sfu-1)=Hzcc+dhz(sf,sfd+1:sfu-1).*(Ex(sf,sfd+2:sfu)-Ex(sf,sfd+1:sfu-1)-Ey(sf+1,sfd+1:sfu-1)+Ey(sf,sfd+1:sfu-1)+Eyinc(sf,sfd+1:sfu-1));

Hz(sff-1,sfd+1:sfu-1)=Hzff+dhz(sff-1,sfd+1:sfu-1).*(Ex(sff-1,sfd+2:sfu)-Ex(sff-1,sfd+1:sfu-1)-Ey(sff,sfd+1:sfu-1)+ Ey(sff-1,sfd+1:sfu-1)- Eyinc(sff,sfd+1:sfu-1));

% up side Note: Hz(sf:sff-1,sfu) is scatter field
Hz(sf:sff-1,sfu)=Hzuu+dhz(sf:sff-1,sfu).*(Ex(sf:sff-1,sfu+1)-Ex(sf:sff-1,sfu)+Exinc(sf:sff-1,sfu)+Ey(sf:sff-1,sfu)-Ey(sf+1:sff,sfu));

% down side Note: Hz(sf:sff-1,sfd)is scatter field
Hz(sf:sff-1,sfd)=Hzdd+dhz(sf:sff-1,sfd).*(Ex(sf:sff-1,sfd+1)-Ex(sf:sff-1,sfd)-Exinc(sf:sff-1,sfd+1)+Ey(sf:sff-1,sfd)-Ey(sf+1:sff,sfd));

%=================Hzinc update zone=============================%

% up
Hzyuinc(1:ssinc-1,2:pm-1)=dahzyuinc(1:ssinc-1,2:pm-1).*Hzyuinc(1:ssinc-1,2:pm-1)+dbhzyuinc(1:ssinc-1,2:pm-1).*(Exuinc(1:ssinc-1,3:pm)-Exuinc(1:ssinc-1,2:pm-1));

Hzyuinc(pm:ssinc-pm,1)=dahzyuinc(pm:ssinc-pm,1).*Hzyuinc(pm:ssinc-pm,1)+dbhzyuinc(pm:ssinc-pm,1).*(Exuinc(pm:ssinc-pm,2)-Exinc(1:Nxinc-1,Nyinc));

Hzxuinc(1:ssinc-1,1:pm-1)=dahzxuinc(1:ssinc-1,1:pm-1).*Hzxuinc(1:ssinc-1,1:pm-1)-dbhzxuinc(1:ssinc-1,1:pm-1).*(Eyuinc(2:ssinc,1:pm-1)-Eyuinc(1:ssinc-1,1:pm-1));

%== boundary domain at the right side 
Hzyuinc(ssinc-pm+1:ssinc-1,1)=dahzyuinc(ssinc-pm+1:ssinc-1,1).*Hzyuinc(ssinc-pm+1:ssinc-1,1)+ dbhzyuinc(ssinc-pm+1:ssinc-1,1).*(Exuinc(ssinc-pm+1:ssinc-1,2)- Exrinc(1:pm-1,Nyinc));

%== boundary domain at the left side 
Hzyuinc(1:pm-1,1)=dahzyu(1:pm-1,1).*Hzyuinc(1:pm-1,1)+dbhzyu(1:pm-1,1).*(Exuinc(1:pm-1,2)-Exlinc(1:pm-1,Nyinc));


% down
Hzydinc(1:ssinc-1,2:pm-1)=dahzydinc(1:ssinc-1,2:pm-1).*Hzydinc(1:ssinc-1,2:pm-1)+dbhzydinc(1:ssinc-1,2:pm-1).*(Exdinc(1:ssinc-1,2:pm-1)-Exdinc(1:ssinc-1,3:pm));

Hzydinc(pm:ssinc-pm,1)=dahzydinc(pm:ssinc-pm,1).*Hzydinc(pm:ssinc-pm,1)+dbhzydinc(pm:ssinc-pm,1).*(Exinc(1:Nxinc-1,1)-Exdinc(pm:ssinc-pm,2));

Hzxdinc(1:ssinc-1,1:pm-1)=dahzxdinc(1:ssinc-1,1:pm-1).*Hzxdinc(1:ssinc-1,1:pm-1)-dbhzxdinc(1:ssinc-1,1:pm-1).*(Eydinc(2:ssinc,1:pm-1)-Eydinc(1:ssinc-1,1:pm-1));

%== boundary domain at the right side 
Hzydinc(ssinc-pm+1:ssinc-1,1)=dahzydinc(ssinc-pm+1:ssinc-1,1).*Hzydinc(ssinc-pm+1:ssinc-1,1)+ dbhzydinc(ssinc-pm+1:ssinc-1,1).*(Exrinc(1:pm-1,1)-Exdinc(ssinc-pm+1:ssinc-1,2));

%== boundary domain at the left side 
Hzydinc(1:pm-1,1)=dahzyd(1:pm-1,1).*Hzydinc(1:pm-1,1)+dbhzyd(1:pm-1,1).*(Exlinc(1:pm-1,1)-Exdinc(1:pm-1,2));

% right
Hzyrinc(1:pm-1,1:Nyinc-1)=dahzyr(1:pm-1,1:Nyinc-1).*Hzyrinc(1:pm-1,1:Nyinc-1)+dbhzyr(1:pm-1,1:Nyinc-1).*(Exrinc(1:pm-1,2:Nyinc)-Exrinc(1:pm-1,1:Nyinc-1));

Hzxrinc(2:pm-1,1:Nyinc-1)=dahzxr(2:pm-1,1:Nyinc-1).*Hzxrinc(2:pm-1,1:Nyinc-1)-dbhzxr(2:pm-1,1:Nyinc-1).*(Eyrinc(3:pm,1:Nyinc-1)-Eyrinc(2:pm-1,1:Nyinc-1));

Hzxrinc(1,1:Ny-1)=dahzxr(1,1:Ny-1).*Hzxrinc(1,1:Nyinc-1)-dbhzxr(1,1:Nyinc-1).*(Eyrinc(2,1:Nyinc-1)-Eyinc(Nxinc,1:Nyinc-1));

% left
Hzylinc(1:pm-1,1:Nyinc-1)=dahzyl(1:pm-1,1:Nyinc-1).*Hzylinc(1:pm-1,1:Nyinc-1)+dbhzyl(1:pm-1,1:Nyinc-1).*(Exlinc(1:pm-1,2:Nyinc)-Exlinc(1:pm-1,1:Nyinc-1));

Hzxlinc(1:pm-2,1:Nyinc-1)=dahzxl(1:pm-2,1:Nyinc-1).*Hzxlinc(1:pm-2,1:Nyinc-1)-dbhzxl(1:pm-2,1:Nyinc-1).*(Eylinc(2:pm-1,1:Ny-1)-Eylinc(1:pm-2,1:Nyinc-1));

Hzxlinc(pm-1,1:Nyinc-1)=dahzxl(pm-1,1:Nyinc-1).*Hzxlinc(pm-1,1:Nyinc-1)-dbhzxl(pm-1,1:Nyinc-1).*(Eyinc(1,1:Nyinc-1)- Eylinc(pm-1,1:Nyinc-1));
%========================Hz PML =======================

%========================Hzy PML right========================
Hzyr(1:pm-1,1:Ny-1)=dahzyr(1:pm-1,1:Ny-1).*Hzyr(1:pm-1,1:Ny-1)+dbhzyr(1:pm-1,1:Ny-1).*(Exr(1:pm-1,2:Ny)-Exr(1:pm-1,1:Ny-1));

%========================Hzy PML left========================
Hzyl(1:pm-1,1:Ny-1)=dahzyl(1:pm-1,1:Ny-1).*Hzyl(1:pm-1,1:Ny-1)+dbhzyl(1:pm-1,1:Ny-1).*(Exl(1:pm-1,2:Ny)-Exl(1:pm-1,1:Ny-1));

%========================Hzy PML up==========================
Hzyu(1:ss-1,2:pm-1)=dahzyu(1:ss-1,2:pm-1).*Hzyu(1:ss-1,2:pm-1)+dbhzyu(1:ss-1,2:pm-1).*(Exu(1:ss-1,3:pm)-Exu(1:ss-1,2:pm-1));

Hzyu(pm:ss-pm,1)=dahzyu(pm:ss-pm,1).*Hzyu(pm:ss-pm,1)+dbhzyu(pm:ss-pm,1).*(Exu(pm:ss-pm,2)-Ex(1:Nx-1,Ny));

%== boundary domain at the rightt side 
Hzyu(ss-pm+1:ss-1,1)=dahzyu(ss-pm+1:ss-1,1).*Hzyu(ss-pm+1:ss-1,1)+ dbhzyu(ss-pm+1:ss-1,1).*(Exu(ss-pm+1:ss-1,2)- Exr(1:pm-1,Ny));

%== boundary domain at the left side 
Hzyu(1:pm-1,1)=dahzyu(1:pm-1,1).*Hzyu(1:pm-1,1)+dbhzyu(1:pm-1,1).*(Exu(1:pm-1,2)-Exl(1:pm-1,Ny));

%========================Hzy PML down==========================
Hzyd(1:ss-1,2:pm-1)=dahzyd(1:ss-1,2:pm-1).*Hzyd(1:ss-1,2:pm-1)+dbhzyd(1:ss-1,2:pm-1).*(Exd(1:ss-1,2:pm-1)-Exd(1:ss-1,3:pm));

Hzyd(pm:ss-pm,1)=dahzyd(pm:ss-pm,1).*Hzyd(pm:ss-pm,1)+dbhzyd(pm:ss-pm,1).*(Ex(1:Nx-1,1)-Exd(pm:ss-pm,2));

%== boundary domain at the right side 
Hzyd(ss-pm+1:ss-1,1)=dahzyd(ss-pm+1:ss-1,1).*Hzyd(ss-pm+1:ss-1,1)+ dbhzyu(ss-pm+1:ss-1,1).*(Exr(1:pm-1,1)-Exd(ss-pm+1:ss-1,2));

%== boundary domain at the left side 
Hzyd(1:pm-1,1)=dahzyd(1:pm-1,1).*Hzyd(1:pm-1,1)+dbhzyd(1:pm-1,1).*(Exl(1:pm-1,1)-Exd(1:pm-1,2));
%=============================================================
%========================Hzx PML right========================
Hzxr(2:pm-1,1:Ny-1)=dahzxr(2:pm-1,1:Ny-1).*Hzxr(2:pm-1,1:Ny-1)-dbhzxr(2:pm-1,1:Ny-1).*(Eyr(3:pm,1:Ny-1)-Eyr(2:pm-1,1:Ny-1));
Hzxr(1,1:Ny-1)=dahzxr(1,1:Ny-1).*Hzxr(1,1:Ny-1)-dbhzxr(1,1:Ny-1).*(Eyr(2,1:Ny-1)-Ey(Nx,1:Ny-1));

%========================Hzx PML left=========================
Hzxl(1:pm-2,1:Ny-1)=dahzxl(1:pm-2,1:Ny-1).*Hzxl(1:pm-2,1:Ny-1)-dbhzxl(1:pm-2,1:Ny-1).*(Eyl(2:pm-1,1:Ny-1)-Eyl(1:pm-2,1:Ny-1));
Hzxl(pm-1,1:Ny-1)=dahzxl(pm-1,1:Ny-1).*Hzxl(pm-1,1:Ny-1)-dbhzxl(pm-1,1:Ny-1).*(Ey(1,1:Ny-1)- Eyl(pm-1,1:Ny-1));

%========================Hzx PML up===========================

Hzxu(1:ss-1,1:pm-1)=dahzxu(1:ss-1,1:pm-1).*Hzxu(1:ss-1,1:pm-1)-dbhzxu(1:ss-1,1:pm-1).*(Eyu(2:ss,1:pm-1)-Eyu(1:ss-1,1:pm-1));

Hzxd(1:ss-1,1:pm-1)=dahzxd(1:ss-1,1:pm-1).*Hzxd(1:ss-1,1:pm-1)-dbhzxd(1:ss-1,1:pm-1).*(Eyd(2:ss,1:pm-1)-Eyd(1:ss-1,1:pm-1));

%==========================observation point=================
%==========================NTFF===============================
for nf=1:nfreq                     
%== incident source DFT phase and amplitude      
RD(nf,1)=sin(n*dt*2*pi*f(nf));
RD(nf,2)=cos(n*dt*2*pi*f(nf));
SN(nf,1)=SN(nf,1) + RD(nf,1)*exp(-(4*pi*((n*dt-to)/tau)^2));
SN(nf,2)=SN(nf,2) + RD(nf,2)*exp(-(4*pi*((n*dt-to)/tau)^2));    
%== Hz DFT
Hzefr(nf,ntff:Ny-ntff-1,1)=Hzefr(nf,ntff:Ny-ntff-1,1)+Hz(Nx-ntff,ntff:Ny-ntff-1).*RD(nf,1);
Hzefr(nf,ntff:Ny-ntff-1,2)=Hzefr(nf,ntff:Ny-ntff-1,2)+Hz(Nx-ntff,ntff:Ny-ntff-1).*RD(nf,2);

Hzefl(nf,ntff:Ny-ntff-1,1)=Hzefl(nf,ntff:Ny-ntff-1,1)+Hz(ntff,ntff:Ny-ntff-1).*RD(nf,1);
Hzefl(nf,ntff:Ny-ntff-1,2)=Hzefl(nf,ntff:Ny-ntff-1,2)+Hz(ntff,ntff:Ny-ntff-1).*RD(nf,2);

Hzefu(nf,ntff:Nx-ntff,1)=Hzefu(nf,ntff:Nx-ntff,1)+Hz(ntff:Nx-ntff,Ny-ntff-1)'.*RD(nf,1);
Hzefu(nf,ntff:Nx-ntff,2)=Hzefu(nf,ntff:Nx-ntff,2)+Hz(ntff:Nx-ntff,Ny-ntff-1)'.*RD(nf,2);

Hzefd(nf,ntff:Nx-ntff,1)=Hzefd(nf,ntff:Nx-ntff,1)+Hz(ntff:Nx-ntff,ntff)'.*RD(nf,1);
Hzefd(nf,ntff:Nx-ntff,2)=Hzefd(nf,ntff:Nx-ntff,2)+Hz(ntff:Nx-ntff,ntff)'.*RD(nf,2);

%== Ex DFT
Exefr(nf,ntff+1:Ny-ntff-1,1)=Exefr(nf,ntff+1:Ny-ntff-1,1)+Ex(Nx-ntff,ntff+1:Ny-ntff-1).*RD(nf,1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
Exefr(nf,ntff+1:Ny-ntff-1,2)=Exefr(nf,ntff+1:Ny-ntff-1,2)+Ex(Nx-ntff,ntff+1:Ny-ntff-1).*RD(nf,2);   

Exefl(nf,ntff+1:Ny-ntff-1,1)=Exefl(nf,ntff+1:Ny-ntff-1,1)+Ex(ntff,ntff+1:Ny-ntff-1).*RD(nf,1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
Exefl(nf,ntff+1:Ny-ntff-1,2)=Exefl(nf,ntff+1:Ny-ntff-1,2)+Ex(ntff,ntff+1:Ny-ntff-1).*RD(nf,2);  

Exefu(nf,ntff:Nx-ntff,1)=Exefu(nf,ntff:Nx-ntff,1)+(Ex(ntff:Nx-ntff,Ny-ntff) + Ex(ntff:Nx-ntff,Ny-ntff-1))'.*RD(nf,1)/2;
Exefu(nf,ntff:Nx-ntff,2)=Exefu(nf,ntff:Nx-ntff,2)+(Ex(ntff:Nx-ntff,Ny-ntff-1) + Ex(ntff:Nx-ntff,Ny-ntff))'.*RD(nf,2)/2;

Exefd(nf,ntff:Nx-ntff,1)=Exefd(nf,ntff:Nx-ntff,1)+(Ex(ntff:Nx-ntff,ntff+1) + Ex(ntff:Nx-ntff,ntff))'.*RD(nf,1)/2;
Exefd(nf,ntff:Nx-ntff,2)=Exefd(nf,ntff:Nx-ntff,2)+(Ex(ntff:Nx-ntff,ntff+1) + Ex(ntff:Nx-ntff,ntff))'.*RD(nf,2)/2;

%== Ey DFT
Eyefr(nf,ntff:Ny-ntff-1,1)=Eyefr(nf,ntff:Ny-ntff-1,1)+(Ey(Nx-ntff,ntff:Ny-ntff-1)+Ey(Nx-ntff+1,ntff:Ny-ntff-1)).*RD(nf,1)/2; 
Eyefr(nf,ntff:Ny-ntff-1,2)=Eyefr(nf,ntff:Ny-ntff-1,2)+(Ey(Nx-ntff,ntff:Ny-ntff-1)+Ey(Nx-ntff+1,ntff:Ny-ntff-1)).*RD(nf,2)/2;

Eyefl(nf,ntff:Ny-ntff-1,1)=Eyefl(nf,ntff:Ny-ntff-1,1)+(Ey(ntff,ntff:Ny-ntff-1)+Ey(ntff+1,ntff:Ny-ntff-1)).*RD(nf,1)/2; 
Eyefl(nf,ntff:Ny-ntff-1,2)=Eyefl(nf,ntff:Ny-ntff-1,2)+(Ey(ntff,ntff:Ny-ntff-1)+Ey(ntff+1,ntff:Ny-ntff-1)).*RD(nf,2)/2;

Eyefu(nf,ntff+1:Nx-ntff,1)=Eyefu(nf,ntff+1:Nx-ntff,1)+Ey(ntff+1:Nx-ntff,Ny-ntff-1)'.*RD(nf,1);
Eyefu(nf,ntff+1:Nx-ntff,2)=Eyefu(nf,ntff+1:Nx-ntff,2)+Ey(ntff+1:Nx-ntff,Ny-ntff-1)'.*RD(nf,2);

Eyefd(nf,ntff+1:Nx-ntff,1)=Eyefd(nf,ntff+1:Nx-ntff,1)+Ey(ntff+1:Nx-ntff,ntff)'.*RD(nf,1);
Eyefd(nf,ntff+1:Nx-ntff,2)=Eyefd(nf,ntff+1:Nx-ntff,2)+Ey(ntff+1:Nx-ntff,ntff)'.*RD(nf,2);
end
%===========================plot figure===============================    
if mod(n,10)==0;
grid off
%subplot(2,1,1);
pcolor(-x,y,Hz);shading flat;

title(strcat('Hz at timestep ',num2str(n)));hold off
caxis([-1.0 1.0]);       % set up color axis
view(90,90)  ;           % set up view angle 
hold on
contour(-x,y,cex,1)
hold off
colorbar;                %  add colorbar             
%colormap hot;
%colormap(jet);
axis image;axis off;

%subplot(2,1,2);
%pcolor(-x,y,Hzinc);shading flat;

%title(strcat('Hzinc at timestep ',num2str(n)));hold off
%caxis([-1.0 1.0]);       % set caxis for color bar
%view(90,90)  ;           % set view angle 
%hold on
%contour(-x,y,cex,1)
%hold off
%colorbar;                % put colorbar             
%colormap hot;
%colormap(jet);
%axis image;axis off;

drawnow

end
%========================end of "Plot figure" section =============================

end
%=====================end time loops=========================

%================near to far caculate========================
for nf=1:nfreq                        
store=sqrt(SN(nf,1)^2+SN(nf,2)^2);
SN(nf,2)=-atan2(SN(nf,1),SN(nf,2));
SN(nf,1)=store;
%==Hz
for i=ntff:Ny-ntff-1
store=sqrt(Hzefr(nf,i,1)^2+Hzefr(nf,i,2)^2)/SN(nf,1);
Hzefr(nf,i,2)=-atan2(Hzefr(nf,i,1),Hzefr(nf,i,2))-SN(nf,2);
Hzefr(nf,i,1)=store;
store=sqrt(Hzefl(nf,i,1)^2+Hzefl(nf,i,2)^2)/SN(nf,1);
Hzefl(nf,i,2)=-atan2(Hzefl(nf,i,1),Hzefl(nf,i,2))-SN(nf,2);
Hzefl(nf,i,1)=store;
end

for i=ntff:Nx-ntff
store=sqrt(Hzefu(nf,i,1)^2+Hzefu(nf,i,2)^2)/SN(nf,1);
Hzefu(nf,i,2)=-atan2(Hzefu(nf,i,1),Hzefu(nf,i,2))-SN(nf,2);
Hzefu(nf,i,1)=store;
store=sqrt(Hzefd(nf,i,1)^2+Hzefd(nf,i,2)^2)/SN(nf,1);
Hzefd(nf,i,2)=-atan2(Hzefd(nf,i,1),Hzefd(nf,i,2))-SN(nf,2);
Hzefd(nf,i,1)=store;
end    

%==Ex
for i=ntff+1:Ny-ntff-1
store=sqrt(Exefr(nf,i,1)^2+Exefr(nf,i,2)^2)/(SN(nf,1));
Exefr(nf,i,2)=-atan2(Exefr(nf,i,1),Exefr(nf,i,2))-SN(nf,2);
Exefr(nf,i,1)=store;
store=sqrt(Exefl(nf,i,1)^2+Exefl(nf,i,2)^2)/(SN(nf,1));
Exefl(nf,i,2)=-atan2(Exefl(nf,i,1),Exefl(nf,i,2))-SN(nf,2);
Exefl(nf,i,1)=store;
end

for i=ntff:Nx-ntff
store=sqrt(Exefu(nf,i,1)^2+Exefu(nf,i,2)^2)/(SN(nf,1));
Exefu(nf,i,2)=-atan2(Exefu(nf,i,1),Exefu(nf,i,2))-SN(nf,2);
Exefu(nf,i,1)=store;
store=sqrt(Exefd(nf,i,1)^2+Exefd(nf,i,2)^2)/(SN(nf,1));
Exefd(nf,i,2)=-atan2(Exefd(nf,i,1),Exefd(nf,i,2))-SN(nf,2);
Exefd(nf,i,1)=store;
end

%==Ey
for i=ntff:Ny-ntff-1
store=sqrt(Eyefr(nf,i,1)^2+Eyefr(nf,i,2)^2)/(SN(nf,1));
Eyefr(nf,i,2)=-atan2(Eyefr(nf,i,1),Eyefr(nf,i,2))-SN(nf,2);
Eyefr(nf,i,1)=store;
store=sqrt(Eyefl(nf,i,1)^2+Eyefl(nf,i,2)^2)/(SN(nf,1));
Eyefl(nf,i,2)=-atan2(Eyefl(nf,i,1),Eyefl(nf,i,2))-SN(nf,2);
Eyefl(nf,i,1)=store;
end
    
for i=ntff+1:Nx-ntff
store=sqrt(Eyefu(nf,i,1)^2+Eyefu(nf,i,2)^2)/(SN(nf,1));
Eyefu(nf,i,2)=-atan2(Eyefu(nf,i,1),Eyefu(nf,i,2))-SN(nf,2);
Eyefu(nf,i,1)=store;
store=sqrt(Eyefd(nf,i,1)^2+Eyefd(nf,i,2)^2)/(SN(nf,1));
Eyefd(nf,i,2)=-atan2(Eyefd(nf,i,1),Eyefd(nf,i,2))-SN(nf,2);
Eyefd(nf,i,1)=store;
end

%===============caculated far field scattered field  Hz==================
k=2*pi*f(nf)/c;             % wave vector

j=sqrt(-1);

for phi=0:359;              % integrate the angle from 0 to 360 
jx=0;
jy=0;
jmz=0;   

for p=ntff:Nx-ntff
    r=sqrt((p-Nx/2)^2+(Ny-ntff-1-Ny/2)^2)*ds;
    r2=sqrt((p-Nx/2)^2+(ntff-Ny/2)^2)*ds;
    the=atan2(Ny-ntff-1-Ny/2,p-Nx/2)*180/pi;
    the2=atan2(ntff-Ny/2,p-Nx/2)*180/pi;
    jx=jx + k*sind(phi)*( Hzefu(nf,p,1)*exp(j*Hzefu(nf,p,2))*exp(j*k*r*cosd(phi-the))  -   Hzefd(nf,p,1)*exp(j*Hzefd(nf,p,2))*exp(j*k*r2*cosd(phi-the2) ));
end
for p=ntff:Ny-ntff-1
    r=sqrt((Nx-ntff-Nx/2)^2+(p-Ny/2)^2)*ds;
    r2=sqrt((ntff-Nx/2)^2+(p-Ny/2)^2)*ds;
    the=atan2(p-Ny/2,Nx-ntff-Nx/2)*180/pi;
    the2=atan2(p-Ny/2,ntff-Nx/2)*180/pi;
    jy=jy + k*cosd(phi)*( Hzefr(nf,p,1)*exp(j*Hzefr(nf,p,2))*exp(j*k*r*cosd(phi-the)) - Hzefl(nf,p,1)*exp(j*Hzefl(nf,p,2))*exp(j*k*r2*cosd(phi-the2)) );
end

for p=ntff:Ny-ntff-1
    r=sqrt((Nx-ntff-Nx/2)^2+(p-Ny/2)^2)*ds;;
    the=atan2(p-Ny/2,Nx-ntff-Nx/2)*180/pi;
    the2=atan2(p-Ny/2,ntff-Nx/2)*180/pi;
    jmz=jmz + 2*pi*f(nf)*epso*( Eyefr(nf,p,1)*exp(j*Eyefr(nf,p,2))*exp(j*k*r*cosd(phi-the)) - Eyefl(nf,p,1)*exp(j*Eyefl(nf,p,2))*exp(j*k*r*cosd(phi-the2) ) ); 
end

for p=ntff:Nx-ntff
    r=sqrt((p-Nx/2)^2+(Ny-ntff-1-Ny/2)^2)*ds;
    r2=sqrt((p-Nx/2)^2+(ntff-Ny/2)^2)*ds;
    the=atan2(Ny-ntff-1-Ny/2,p-Nx/2)*180/pi;
    the2=atan2(ntff-Ny/2,p-Nx/2)*180/pi;
    jmz=jmz + 2*pi*f(nf)*epso*(-Exefu(nf,p,1)*exp(j*Exefu(nf,p,2))*exp(j*k*r*cosd(phi-the)) + Exefd(nf,p,1)*exp(j*Exefd(nf,p,2))*exp(j*k*r2*cosd(phi-the2) ));
end

%=============caculate RCS=================================

Fphi=exp(j*pi/4)/sqrt(8*pi*k)*ds*(jx+jy+jmz);  

RCS(nf,phi+1)=2*pi*abs(Fphi)^2;           

tscs(nf,1) = tscs(nf,1) + RCS(nf,phi+1);

end

%=================phi loop=================================

end

%=================nfreq loop===============================

%==========================================================
%subplot(2,1,1);
plot(RCS(1,1:181)/(c/f(1)));
axis([0 180 0 1.21 ]); 
xlabel('angle');
ylabel('RCS/lamda');
tscs=tscs/360;

%subplot(2,1,2);
%plot(10*log10(RCS(1:180)));
%axis([1 180 -40 0 ]); 
%xlabel('角度');
%ylabel('RCS/dBm');
%plot(Hztest);
%tscs
toc;
