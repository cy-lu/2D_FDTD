      !******************************************************************
	! ** 2D_FDTD **            Jin-You, Lu email:chun.lu@ku.ac.ae    
	!******************************************************************
	! Simulatation of a gold-shell silica-core nanocylinder 
	! The dispersion function of gold is modelled with 
	! Three critical point pairs "CP3 model", as reported in
	! Superlattices and Microstructures, 47(1), 60-65. 2010. 
	! By directly running this code, you will get the 
	! extinction/scattering/absorption coefficients of a
	! gold-shell cylinder with outer radius of 50 nm
	! and inner radius of 40 nm in a wide spectral range
	! from 200 nm to 2000 nm.
	! Qext (extinction efficiency) = Cext (extinction cross section) / (A:geometrical project area of the nanostructure) 
	! for a nanosphere: A=pi*radius^2
	! for a nanocylinder: A=2*radius (i.e. a sphere in a 2D FDTD)
	! The definitions of Cext/Cabs/Csca extinction, absorption, and scattering cross section can be found in 
	! Nanoscale Res Lett. 2011; 6(1): 173.
	! **************************************************************************************************
	program  JIYO_2Dfdtd     
	implicit none		   
	integer  Nx			   
	parameter (Nx=500)  
	integer  Ny			   
	parameter (Ny=500)  
	real*8   ds,dt,c,pi,epso,muo,zo,cent_x,cent_y,
	1xdist,ydist
	parameter(ds=1*1D-9)
	integer i,j,k,p,q,ii,kk,nmax,n,nf 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!the variable JJ is used as sqrt(-1) in this code,
	!that's why we did not use it as one of integers
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	parameter (nmax=30000)
	!************field component*********************
	real*8 Ex,Ey,Hz,Eyinc,Hzinc,pulse 
	dimension Ex(Nx,Ny),Ey(Nx,Ny),Hz(Nx,Ny)
	dimension Eyinc(Nx),Hzinc(Nx)
	real*8,dimension(Nx,Ny) :: caex,cbex,caey,cbey,ccex,ccey,dbhz

	complex*8,dimension(Nx,Ny) :: phicp1ex1,phicp2ex1,phicp3ex1,
	1phicp1ey1,phicp2ey1,phicp3ey1

	integer,dimension(Nx,Ny) :: spacex,spacey
	!*******************DFT****NTFF************************
	integer ntff,nfreq
	parameter (ntff=10)
	parameter (nfreq=600)
	real*8 wave,OMEGA,wmax,wmin,wcenter,tau,SN,RD,store,EXT,ABSR,SCA,
	1waves,jr,jl,ju,jd
	dimension wave(nfreq),OMEGA(nfreq),SN(nfreq,2),RD(nfreq,2),
	1EXT(nfreq),ABSR(nfreq),SCA(nfreq),waves(nfreq)
	!======================DFT 1D incident light============================
	real*8 Hzincefr,Hzincefl,Hzincefu,Hzincefd,Eyincefr,Eyincefl
	dimension Hzincefr(nfreq,2),Hzincefl(nfreq,2),
	1Hzincefu(nfreq,Nx-ntff,2),Hzincefd(nfreq,Nx-ntff,2),
	1Eyincefr(nfreq,2),Eyincefl(nfreq,2)
	!======================DFT 1D incident light============================

	real*8 Hzefr,Hzefl,Hzefu,Hzefd,Exefr,Exefl,Exefu,Exefd,Eyefr,
	1Eyefl,Eyefu,Eyefd
	dimension Hzefr(nfreq,Ny-ntff-1,2),Hzefl(nfreq,Ny-ntff-1,2),
	1Hzefu(nfreq,Nx-ntff,2),Hzefd(nfreq,Nx-ntff,2),
	1Exefr(nfreq,Ny-ntff-1,2),Exefl(nfreq,Ny-ntff-1,2),
	1Exefu(nfreq,Nx-ntff,2),Exefd(nfreq,Nx-ntff,2),
	1Eyefr(nfreq,Ny-ntff-1,2),Eyefl(nfreq,Ny-ntff-1,2),
	1Eyefu(nfreq,Nx-ntff,2),Eyefd(nfreq,Nx-ntff,2)

	!=======================pml space parameters===================
	integer ssx,ssy
	real*8 pm,m
	parameter (pm=16)
	parameter (m=4)
	parameter (ssx=Nx+2*(pm-1))
	parameter (ssy=Ny+2*(pm-1))
	real*8 sigmax,sigmar,sigmars,sigmxl,sigmxls,sigmxu,sigmxus,
	1sigmxd,sigmxds
	real*8,dimension(pm,Ny)  :: Exr,Eyr,Hzxr,Hzyr,caexr,cbexr,
	1caeyr,cbeyr,dahzxr,dbhzxr,dahzyr,dbhzyr
	real*8,dimension(pm,Ny)  :: Exl,Eyl,Hzxl,Hzyl,caexl,cbexl,caeyl,
	1cbeyl,dahzxl,dbhzxl,dahzyl,dbhzyl
	real*8,dimension(ssx,pm) :: Exu,Eyu,Hzxu,Hzyu,caexu,cbexu,caeyu,
	1cbeyu,dahzxu,dbhzxu,dahzyu,dbhzyu
	real*8,dimension(ssx,pm) :: Exd,Eyd,Hzxd,Hzyd,caexd,cbexd,caeyd,
	1cbeyd,dahzxd,dbhzxd,dahzyd,dbhzyd
	!======================SFTF parameters=========================
	integer stft,sfl,sfr,sfu,sfd,Nxinc,Nyinc
	parameter (stft=5)
	parameter (sfl=stft)
	parameter (sfr=Nx-stft)
	parameter (sfu=Ny-stft)
	parameter (sfd=stft)
	parameter (Nxinc=Nx)
	parameter (Nyinc=Ny)
	real*8 rbc
	dimension rbc(nmax+2)
	real*8 Eyll,Eyrr,Exuu,Exdd,Hzll,Hzrr,Hzuu,Hzdd
	dimension Eyll(sfu),Eyrr(sfu),Exuu(sfr),Exdd(sfr)
	dimension Hzll(sfu),Hzrr(sfu),Hzuu(sfr),Hzdd(sfr)
	!**********************************************
	real*8 epsr,chi0,A1,phi1,omega1,gamma1,
	1A2,phi2,omega2,gamma2,A3,phi3,omega3,gamma3,sigma,
	1radius,radius2,diameter,distance
	complex*8 jj,deltachicp10,deltachicp20,deltachicp30,
	1cretcp1ex,cretcp2ex,cretcp3ex,cretcp1ey,cretcp2ey,cretcp3ey
	open(unit=23,file="domain.dat")
	open(unit=24,file="EXT.dat")
	open(unit=25,file="SCA.dat")
	open(unit=26,file="ABS.dat")
	open(unit=27,file="Waves.dat")
	epsr=1.0
	chi0=0
	pi=3.141592653589793
	epso=885*1D-14
	muo=4*pi*1D-7
	c=1/sqrt(epso*muo)
	!****************cell************************

	dt=ds/(2*c)
	zo=sqrt(muo/epso)
	!****************initialize******************
	Do i=1,Nx
		Do j=1,Ny
	caex(i,j)=(epsr)/(epsr+chi0)
	cbex(i,j)=(dt/(epso*ds))/(epsr+chi0)
	ccex(i,j)=1/(epsr+chi0)
	caey(i,j)=(epsr)/(epsr+chi0)
	cbey(i,j)=(dt/(epso*ds))/(epsr+chi0)
	ccey(i,j)=1/(epsr+chi0)
	dbhz(i,j)=(dt/(muo*ds))
	spacex(i,j)=0
	spacey(i,j)=0
		end do
	end do
	!material dimension and propertie
	!Taken from the values of gold in Table 1 of 
	!Superlattices and Microstructures, 47(1), 60-65. 2010. 
	!
	epsr=1.1156
	sigma=4.2399257866E16	
	A1=0.5548	
	phi1=2.8463	
	omega1=4.506033811E16
	gamma1=5.0895798919E16
	A2=679.7606
	phi2=-0.0998	
	omega2=3.4587148533E14
	gamma2=3.06427186E13
	A3=3.5244
	phi3=4.6586
	omega3=3.5831868337E15
	gamma3=1.6878397906E15

	
	jj=(0,1)
	chi0=-sigma*dt+
	1real(-jj*2*A1*omega1
     1*exp(-jj*phi1)*(1-exp((-gamma1+jj*omega1)*dt))/(gamma1-jj*omega1)
     1      -jj*2*A2*omega2
     1*exp(-jj*phi2)*(1-exp((-gamma2+jj*omega2)*dt))/(gamma2-jj*omega2)
     1	  -jj*2*A3*omega3
     1*exp(-jj*phi3)*(1-exp((-gamma3+jj*omega3)*dt))/(gamma3-jj*omega3))
      deltachicp10=-jj*2*A1*omega1*exp(-jj*phi1)*
	1(1-exp((-gamma1+jj*omega1)*dt))*(1-exp((-gamma1+jj*omega1)*dt))
     1/(gamma1-jj*omega1)
	deltachicp20=-jj*2*A2*omega2*exp(-jj*phi2)*
	1(1-exp((-gamma2+jj*omega2)*dt))*(1-exp((-gamma2+jj*omega2)*dt))
     1/(gamma2-jj*omega2)
	deltachicp30=-jj*2*A3*omega3*exp(-jj*phi3)*
	1(1-exp((-gamma3+jj*omega3)*dt))*(1-exp((-gamma3+jj*omega3)*dt))
     1/(gamma3-jj*omega3)
	cretcp1ex=exp((-gamma1+jj*omega1)*dt)
	cretcp2ex=exp((-gamma2+jj*omega2)*dt)
	cretcp3ex=exp((-gamma3+jj*omega3)*dt)
	cretcp1ey=exp((-gamma1+jj*omega1)*dt)
	cretcp2ey=exp((-gamma2+jj*omega2)*dt)
	cretcp3ey=exp((-gamma3+jj*omega3)*dt)

	!=========================Shell dimensions============================
	diameter=200*1D-9
	radius=50*1D-9      ! outradius
	radius2=40*1D-9     ! inner radius
	cent_x=(Nx+1)/2.0
	cent_y=(Ny+1)/2.0
	Do p=1,Nx
	  Do q=1,Ny
	   Do ii=-1,1
	    Do kk=-1,1
		   xdist=(p+0.5-(cent_x))-0.333*ii
		   ydist=(q-(cent_y))-0.333*kk
	      if(sqrt(((xdist)*ds)**2+((ydist)*ds)**2).LT.radius) then
				caex(p,q)=(epsr)/(epsr+chi0)
				cbex(p,q)=(dt/(epso*ds))/(epsr+chi0)
				ccex(p,q)=1/(epsr+chi0)
				spacex(p,q)=1
		  end if
	    End Do
	   End Do
		Do ii=-1,1
		 Do kk=-1,1
		     xdist=(p-(cent_x))-0.333*ii
		     ydist=(q+0.5-(cent_y))-0.333*kk
		if(sqrt(((xdist)*ds)**2+((ydist)*ds)**2).LT.radius) then
	
            caey(p,q)=(epsr)/(epsr+chi0)
            cbey(p,q)=(dt/(epso*ds))/(epsr+chi0)
            ccey(p,q)=1/(epsr+chi0)
		  spacey(p,q)=1
	    end if
	     End Do
	   End Do
	  End Do
	End Do	

	!=========================Core dimensions============================
	epsr=2.1
	cent_x=(Nx+1)/2.0
	cent_y=(Ny+1)/2.0
	Do p=1,Nx
	 Do q=1,Ny
        Do ii=-1,1
            Do kk=-1,1
                xdist=(p+1/2-(cent_x))-0.333*ii;
                ydist=(q-(cent_y))-0.333*kk;
			if (sqrt(((xdist)*ds)**2+((ydist)*ds)**2).LT.radius2) then
          
           caex(p,q)=(epsr)/(epsr)
           cbex(p,q)=(dt/(epso*ds))/(epsr)
           ccex(p,q)=1/(epsr)
     
           	spacex(p,q)=0

			end if
            end DO
        end DO
            Do ii=-1,1
             Do kk=-1,1
                xdist=(p-(cent_x))-0.333*ii
                ydist=(q+1/2-(cent_y))-0.333*kk
			if (sqrt(((xdist)*ds)**2+((ydist)*ds)**2).LT.radius2) then

            caey(p,q)=(epsr)/(epsr)
            cbey(p,q)=(dt/(epso*ds))/(epsr)
            ccey(p,q)=1/(epsr)
            
           	spacey(p,q)=0
			end if
            end do
           end do
		end do
	end do


	!=========================PML sigmax============================
	sigmax=0.8*(m+1)/(zo*ds)
		
	!=======================initialize pml parameters============================

	Do i=1,pm
	   
		sigmar=sigmax*((i-1)/(pm-1))**m
		sigmars=(muo/epso)*sigmax*((i-0.5)/(pm-1))**m
		caexr(i,:)=1
		cbexr(i,:)=dt/(epso*ds)
		caeyr(i,:)=(1-sigmar*dt/(2*epso))/(1+sigmar*dt/(2*epso))
		cbeyr(i,:)=(dt/(ds*epso))/(1+sigmar*dt/(2*epso))
		dahzxr(i,:)=(1-sigmars*dt/(2*muo))/(1+sigmars*dt/(2*muo))
		dbhzxr(i,:)=(dt/(ds*muo))/(1+sigmars*dt/(2*muo))
		dahzyr(i,:)=1
		dbhzyr(i,:)=dt/(muo*ds)

		
	END DO	

	Do i=1,pm
		sigmxl=sigmax*((pm-i)/(pm-1))**m
		sigmxls=muo/epso*sigmax*((pm-i-0.5)/(pm-1))**m
		caexl(i,:)=1;
		cbexl(i,:)=dt/(epso*ds);
		caeyl(i,:)=(1-sigmxl*dt/(2*epso))/(1+sigmxl*dt/(2*epso));
		cbeyl(i,:)=(dt/(ds*epso))/(1+sigmxl*dt/(2*epso));
		dahzxl(i,:)=(1-sigmxls*dt/(2*muo))/(1+sigmxls*dt/(2*muo));
		dbhzxl(i,:)=(dt/(ds*muo))/(1+sigmxls*dt/(2*muo));
		dahzyl(i,:)=1;
		dbhzyl(i,:)=dt/(muo*ds);
	ENd DO
	
	Do j=1,pm
	   Do i=1,ssx
	    sigmxu=sigmax*((j-1)/(pm-1))**m;
	    sigmxus=muo/epso*sigmax*((j-0.5)/(pm-1))**m;
	    caexu(i,j)=(1-sigmxu*dt/(2*epso))/(1+sigmxu*dt/(2*epso));
          cbexu(i,j)=(dt/(ds*epso))/(1+sigmxu*dt/(2*epso));
		caeyu(i,j)=1;
		cbeyu(i,j)=dt/(epso*ds);
		dahzxu(i,j)=1;
		dbhzxu(i,j)=dt/(muo*ds);
		dahzyu(i,j)=(1-sigmxus*dt/(2*muo))/(1+sigmxus*dt/(2*muo));
		dbhzyu(i,j)=(dt/(ds*muo))/(1+sigmxus*dt/(2*muo));
	   END DO
	END DO
	
		Do j=1,pm
	      Do i=1,pm
	 sigmxu=sigmax*((pm-i)/(pm-1))**m;
       sigmxus=muo/epso*sigmax*((pm-i-0.5)/(pm-1))**m;
       caeyu(i,j)=(1-sigmxu*dt/(2*epso))/(1+sigmxu*dt/(2*epso));
       cbeyu(i,j)=(dt/(ds*epso))/(1+sigmxu*dt/(2*epso));
       dahzxu(i,j)=(1-sigmxus*dt/(2*muo))/(1+sigmxus*dt/(2*muo));
       dbhzxu(i,j)=(dt/(ds*muo))/(1+sigmxus*dt/(2*muo));
	      END DO
	    END DO

	   DO j=1,pm
		DO i=pm+Nx-1,ssx
	    sigmxu=sigmax*((i-(pm+Nx-1))/(pm-1))**m;
		sigmxus=muo/epso*sigmax*((i-(pm+Nx-1)+0.5)/(pm-1))**m;
		caeyu(i,j)=(1-sigmxu*dt/(2*epso))/(1+sigmxu*dt/(2*epso));
		cbeyu(i,j)=(dt/(ds*epso))/(1+sigmxu*dt/(2*epso));
		dahzxu(i,j)=(1-sigmxus*dt/(2*muo))/(1+sigmxus*dt/(2*muo));
		dbhzxu(i,j)=(dt/(ds*muo))/(1+sigmxus*dt/(2*muo));
	    END DO
	   END DO
	
	      Do j=1,pm
	         Do i=1,ssx
	    sigmxd=sigmax*((pm-j)/(pm-1))**m;
		sigmxds=muo/epso*sigmax*((pm-j-0.5)/(pm-1))**m;
		caexd(i,j)=(1-sigmxd*dt/(2*epso))/(1+sigmxd*dt/(2*epso));
		cbexd(i,j)=(dt/(ds*epso))/(1+sigmxd*dt/(2*epso));
		caeyd(i,j)=1;
		cbeyd(i,j)=dt/(epso*ds);
		dahzxd(i,j)=1;
		dbhzxd(i,j)=dt/(muo*ds);
		dahzyd(i,j)=(1-sigmxds*dt/(2*muo))/(1+sigmxds*dt/(2*muo));
		dbhzyd(i,j)=(dt/(ds*muo))/(1+sigmxds*dt/(2*muo))
	       END DO
	  END DO
	
	    Do j=1,pm
            Do i=1,pm
	    sigmxd=sigmax*((pm-i)/(pm-1))**m;
          sigmxds=muo/epso*sigmax*((pm-i-0.5)/(pm-1))**m;
          caeyd(i,j)=(1-sigmxd*dt/(2*epso))/(1+sigmxd*dt/(2*epso));
          cbeyd(i,j)=(dt/(ds*epso))/(1+sigmxd*dt/(2*epso));
          dahzxd(i,j)=(1-sigmxds*dt/(2*muo))/(1+sigmxds*dt/(2*muo));
          dbhzxd(i,j)=(dt/(ds*muo))/(1+sigmxds*dt/(2*muo))
	     END DO
	  END DO

	
	     Do j=1,pm
              Do i=pm+Nx-1,ssx
			sigmxd=sigmax*((i-(pm+Nx-1))/(pm-1))**m
			sigmxds=muo/epso*sigmax*((i-(pm+Nx-1)+0.5)/(pm-1))**m
			caeyd(i,j)=(1-sigmxd*dt/(2*epso))/(1+sigmxd*dt/(2*epso))
			cbeyd(i,j)=(dt/(ds*epso))/(1+sigmxd*dt/(2*epso))
			dahzxd(i,j)=(1-sigmxds*dt/(2*muo))/(1+sigmxds*dt/(2*muo))
			dbhzxd(i,j)=(dt/(ds*muo))/(1+sigmxds*dt/(2*muo))
	       END DO
	    END DO
     
	!=========================NTFF initialization============================
	wmax=2*3.14159*c/(200*1D-9)
	wmin=2*3.14159*c/(2000*1D-9)
	wcenter=(wmax+wmin)/2
	tau=4/(wmax-wmin)
	   Do i=1,nfreq
      	wave(i)=200*1D-9+(i-1)*1800*1D-9/(nfreq-1);
	    OMEGA(i)=2*pi*c/wave(i);
	   END DO
	 
	!========================field initialization============================
	Ex(:,:)=0.d0
	Ey(:,:)=0.d0
	Hz(:,:)=0.d0
	Eyinc(:)=0.D0
	Hzinc(:)=0.D0
	phicp1ex1(:,:)=(0.d0,0.d0)
	phicp2ex1(:,:)=(0.d0,0.d0)
	phicp3ex1(:,:)=(0.d0,0.d0)
	phicp1ey1(:,:)=(0.d0,0.d0)
	phicp2ey1(:,:)=(0.d0,0.d0)
	phicp3ey1(:,:)=(0.d0,0.d0)
	Eyll(:)=0.d0
	Eyrr(:)=0.d0
	Exuu(:)=0.d0
	Exdd(:)=0.d0

	!=================================FDTD LOOP===============================
	Do n=1,nmax
	pulse=(exp(-1*((n*dt-4*tau)/tau)
	1*((n*dt-4*tau)/tau))*sin(wcenter*(n*dt-4*tau)))
	Hzinc(1)=pulse
	
	Eyll(sfd+1:sfu-1)=Ey(sfl,sfd+1:sfu-1)
	Eyrr(sfd+1:sfu-1)=Ey(sfr,sfd+1:sfu-1)
	Exuu(sfl:sfr-1)=Ex(sfl:sfr-1,sfu)
	Exdd(sfl:sfr-1)=Ex(sfl:sfr-1,sfd+1)

	!=====================Eyinc update zone===============
	rbc(n+1)=Eyinc(Nxinc-1)
	Eyinc(2:Nxinc-1)=Eyinc(2:Nxinc-1)-
	1dt/(epso*ds)*(Hzinc(2:Nxinc-1)-Hzinc(1:Nxinc-2))
	Eyinc(Nxinc)=rbc(n)
	!=================three critical points Loop==========
	Do j=1,Ny
	  Do i=1,Nx
	if(spacex(i,j).EQ.1) then 
	 phicp1ex1(i,j)=deltachicp10*Ex(i,j)+
	1cretcp1ex*phicp1ex1(i,j)
       phicp2ex1(i,j)=deltachicp20*Ex(i,j)+
	1cretcp2ex*phicp2ex1(i,j)
       phicp3ex1(i,j)=deltachicp30*Ex(i,j)+
	1cretcp3ex*phicp3ex1(i,j)
	end if
	if(spacey(i,j).EQ.1) then 
       phicp1ey1(i,j)=deltachicp10*Ey(i,j)+
	1cretcp1ey*phicp1ey1(i,j)
       phicp2ey1(i,j)=deltachicp20*Ey(i,j)+
	1cretcp2ey*phicp2ey1(i,j)
       phicp3ey1(i,j)=deltachicp30*Ey(i,j)+
	1cretcp3ey*phicp3ey1(i,j)
	end if
	   END DO
	END DO
	
	!=====================update eq Ex =======================
	Do j=2,Ny-1
	   Do i=1,Nx-1
	 Ex(i,j)=caex(i,j)*Ex(i,j)+
	1cbex(i,j)*(Hz(i,j)-Hz(i,j-1))+
     1ccex(i,j)*real(phicp1ex1(i,j)+
     1phicp2ex1(i,j)+phicp3ex1(i,j))
	   END DO
	END DO
      !========================SF======Ex=======TF===================
	Ex(sfl:sfr-1,sfu)=Exuu(sfl:sfr-1)+cbex(sfl:sfr-1,sfu)*
	1(Hz(sfl:sfr-1,sfu)-Hz(sfl:sfr-1,sfu-1)+Hzinc(sfl:sfr-1))
      Ex(sfl:sfr-1,sfd+1)=Exdd(sfl:sfr-1)+cbex(sfl:sfr-1,sfd+1)*
	1(Hz(sfl:sfr-1,sfd+1)-Hz(sfl:sfr-1,sfd)-Hzinc(sfl:sfr-1))

	!=========================Ex PML RIGHT==========================
	Exr(1:pm-1,2:Ny-1)=caexr(1:pm-1,2:Ny-1)*Exr(1:pm-1,2:Ny-1)+
	1cbexr(1:pm-1,2:Ny-1)*(Hzxr(1:pm-1,2:Ny-1)-Hzxr(1:pm-1,1:Ny-2)+
     1Hzyr(1:pm-1,2:Ny-1)-Hzyr(1:pm-1,1:Ny-2))
	
       !add UP_PML
	Exr(1:pm-1,Ny)=caexr(1:pm-1,Ny)*Exr(1:pm-1,Ny)+
	1cbexr(1:pm-1,Ny)*(Hzxu(pm+Nx-1:ssx-1,1)+Hzyu(pm+Nx-1:ssx-1,1)-
     1Hzxr(1:pm-1,Ny-1)-Hzyr(1:pm-1,Ny-1))

	
	 !add Down_PML
	Exr(1:pm-1,1)=caexr(1:pm-1,1)*Exr(1:pm-1,1)+
	1cbexr(1:pm-1,1)*(Hzxr(1:pm-1,1)+Hzyr(1:pm-1,1)-
     1Hzxd(pm+Nx-1:ssx-1,pm-1)-Hzyd(pm+Nx-1:ssx-1,pm-1))
	!=========================Ex PML Left==========================
	Exl(1:pm-1,2:Ny-1)=caexl(1:pm-1,2:Ny-1)*Exl(1:pm-1,2:Ny-1)+
	1cbexl(1:pm-1,2:Ny-1)*(Hzxl(1:pm-1,2:Ny-1)-Hzxl(1:pm-1,1:Ny-2)+
     1Hzyl(1:pm-1,2:Ny-1)-Hzyl(1:pm-1,1:Ny-2));
	!add UP_PML
	Exl(1:pm-1,Ny)=caexl(1:pm-1,Ny)*Exl(1:pm-1,Ny)+
	1cbexl(1:pm-1,Ny)*(Hzxu(1:pm-1,1)+Hzyu(1:pm-1,1)-
     1Hzxl(1:pm-1,Ny-1)-Hzyl(1:pm-1,Ny-1));
	!add DOWN_PML
	Exl(1:pm-1,1)=caexl(1:pm-1,1)*Exl(1:pm-1,1)+
	1cbexl(1:pm-1,1)*(Hzxl(1:pm-1,1)+Hzyl(1:pm-1,1)-
     1Hzxd(1:pm-1,pm-1)-Hzyd(1:pm-1,pm-1));
	!=========================Ex PML UP==========================
	Exu(1:ssx-1,2:pm-1)=caexu(1:ssx-1,2:pm-1)*Exu(1:ssx-1,2:pm-1)+
	1cbexu(1:ssx-1,2:pm-1)*(Hzxu(1:ssx-1,2:pm-1)-Hzxu(1:ssx-1,1:pm-2)+
	1Hzyu(1:ssx-1,2:pm-1)-Hzyu(1:ssx-1,1:pm-2));
	Ex(1:Nx-1,Ny)=Ex(1:Nx-1,Ny)+cbex(1:Nx-1,Ny)*(Hzxu(pm:pm+Nx-2,1)+
	1Hzyu(pm:pm+Nx-2,1)-Hz(1:Nx-1,Ny-1));
	!=========================Ex PML UP==========================
	 Exu(1:ssx-1,2:pm-1)=caexu(1:ssx-1,2:pm-1)*Exu(1:ssx-1,2:pm-1)+
	1cbexu(1:ssx-1,2:pm-1)*(Hzxu(1:ssx-1,2:pm-1)-Hzxu(1:ssx-1,1:pm-2)+
     1Hzyu(1:ssx-1,2:pm-1)-Hzyu(1:ssx-1,1:pm-2));
	 Ex(1:Nx-1,Ny)=Ex(1:Nx-1,Ny)+cbex(1:Nx-1,Ny)*(Hzxu(pm:pm+Nx-2,1)+
	1Hzyu(pm:pm+Nx-2,1)-Hz(1:Nx-1,Ny-1))

      !=========================Ey==========================
	Do j=1,Ny-1
	 Do i=2,Nx-1
	Ey(i,j)=caey(i,j)*Ey(i,j)+
	1cbey(i,j)*(-Hz(i,j)+Hz(i-1,j))+
     1ccey(i,j)*real(phicp1ey1(i,j)+
     1phicp2ey1(i,j)+phicp3ey1(i,j))
	 END DO
	END DO
	
	!=======================SF=====Ey=======TF==========================
	Ey(sfl,sfd+1:sfu-1)=Eyll(sfd+1:sfu-1)-
	1cbey(sfl,sfd+1:sfu-1)*(Hz(sfl,sfd+1:sfu-1)-
     1Hz(sfl-1,sfd+1:sfu-1)-Hzinc(sfl));
	Ey(sfr,sfd+1:sfu-1)=Eyrr(sfd+1:sfu-1)-
	1cbey(sfr,sfd+1:sfu-1)*(Hz(sfr,sfd+1:sfu-1)-
     1Hz(sfr-1,sfd+1:sfu-1)+Hzinc(sfr-1))
      !=========================Ey PML RIGHT========================== 
	Eyr(2:pm-1,1:Ny-1)=caeyr(2:pm-1,1:Ny-1)*Eyr(2:pm-1,1:Ny-1)+
	1cbeyr(2:pm-1,1:Ny-1)*(-Hzxr(2:pm-1,1:Ny-1)+Hzxr(1:pm-2,1:Ny-1)-
     1Hzyr(2:pm-1,1:Ny-1)+Hzyr(1:pm-2,1:Ny-1))   
	Ey(Nx,1:Ny-1)=Ey(Nx,1:Ny-1)+cbey(Nx,1:Ny-1)*(-Hzxr(1,1:Ny-1)-
	1Hzyr(1,1:Ny-1)+Hz(Nx-1,1:Ny-1))
	!=========================Ey PML LEFT========================== 
	 Eyl(2:pm-1,1:Ny-1)=caeyl(2:pm-1,1:Ny-1)*Eyl(2:pm-1,1:Ny-1)+
	1cbeyl(2:pm-1,1:Ny-1)*(-Hzxl(2:pm-1,1:Ny-1)+Hzxl(1:pm-2,1:Ny-1)-
     1Hzyl(2:pm-1,1:Ny-1)+Hzyl(1:pm-2,1:Ny-1))
	 Ey(1,1:Ny-1)=Ey(1,1:Ny-1)+cbey(1,1:Ny-1)*(-Hz(1,1:Ny-1)+
	1Hzxl(pm-1,1:Ny-1)+Hzyl(pm-1,1:Ny-1))
	 Ey(Nx,1:Ny-1)=Ey(Nx,1:Ny-1)-(dt/(epso*ds))*(Hzxr(1,1:Ny-1)+
	1Hzyr(1,1:Ny-1)-Hz(Nx-1,1:Ny-1))
	!=========================Ey PML UP========================== 
	Eyu(2:ssx-1,1:pm-1)=caeyu(2:ssx-1,1:pm-1)*Eyu(2:ssx-1,1:pm-1)+
	1cbeyu(2:ssx-1,1:pm-1)*(-Hzxu(2:ssx-1,1:pm-1)+Hzxu(1:ssx-2,1:pm-1)-
     1Hzyu(2:ssx-1,1:pm-1)+Hzyu(1:ssx-2,1:pm-1))
	!=========================Ey PML DoWN========================== 
	Eyd(2:ssx-1,1:pm-1)=caeyd(2:ssx-1,1:pm-1)*Eyd(2:ssx-1,1:pm-1)+
	1cbeyd(2:ssx-1,1:pm-1)*(-Hzxd(2:ssx-1,1:pm-1)+Hzxd(1:ssx-2,1:pm-1)-
     1Hzyd(2:ssx-1,1:pm-1)+Hzyd(1:ssx-2,1:pm-1))
	!=========================Hz  SFTF record matrix=============== 
	   Hzll(sfd+1:sfu-1)=Hz(sfl,sfd+1:sfu-1)
         Hzrr(sfd+1:sfu-1)=Hz(sfr-1,sfd+1:sfu-1)
         Hzuu(sfl:sfr-1)=Hz(sfl:sfr-1,sfu)
         Hzdd(sfl:sfr-1)=Hz(sfl:sfr-1,sfd)
	!=========================Hz========================== 
		Hz(1:Nx-1,1:Ny-1)=Hz(1:Nx-1,1:Ny-1)+
	1dbhz(1:Nx-1,1:Ny-1)*(Ex(1:Nx-1,2:Ny)-Ex(1:Nx-1,1:Ny-1)-
     1Ey(2:Nx,1:Ny-1)+Ey(1:Nx-1,1:Ny-1))
      !=======================incident hz loop==========================  
	Hzinc(1:Nxinc-1)=Hzinc(1:Nxinc-1)+dt/muo/ds*(Eyinc(1:Nxinc-1)-
	1Eyinc(2:Nxinc))
	!=========================Hz SFTF========================== 
	Hz(sfl,sfd+1:sfu-1)=Hzll(sfd+1:sfu-1)+
	1dbhz(sfl,sfd+1:sfu-1)*(Ex(sfl,sfd+2:sfu)-Ex(sfl,sfd+1:sfu-1)-
     1Ey(sfl+1,sfd+1:sfu-1)+Ey(sfl,sfd+1:sfu-1)+Eyinc(sfl))
      Hz(sfr-1,sfd+1:sfu-1)=Hzrr(sfd+1:sfu-1)+
	1dbhz(sfr-1,sfd+1:sfu-1)*(Ex(sfr-1,sfd+2:sfu)-
     1Ex(sfr-1,sfd+1:sfu-1)-Ey(sfr,sfd+1:sfu-1)+
     1Ey(sfr-1,sfd+1:sfu-1)-Eyinc(sfr))
      Hz(sfl:sfr-1,sfu)=Hzuu(sfl:sfr-1)+
	1dbhz(sfl:sfr-1,sfu)*(Ex(sfl:sfr-1,sfu+1)-Ex(sfl:sfr-1,sfu)+
     1Ey(sfl:sfr-1,sfu)-Ey(sfl+1:sfr,sfu))
      Hz(sfl:sfr-1,sfd)=Hzdd(sfl:sfr-1)+
	1dbhz(sfl:sfr-1,sfd)*(Ex(sfl:sfr-1,sfd+1)-Ex(sfl:sfr-1,sfd)+
     1Ey(sfl:sfr-1,sfd)-Ey(sfl+1:sfr,sfd))
      !=========================Hzx PML RIGHT ==========================   
	Hzxr(2:pm-1,1:Ny-1)=dahzxr(2:pm-1,1:Ny-1)*Hzxr(2:pm-1,1:Ny-1)+
	1dbhzxr(2:pm-1,1:Ny-1)*(-Eyr(3:pm,1:Ny-1)+Eyr(2:pm-1,1:Ny-1))
	Hzyr(1:pm-1,1:Ny-1)=dahzyr(1:pm-1,1:Ny-1)*Hzyr(1:pm-1,1:Ny-1)+
	1dbhzyr(1:pm-1,1:Ny-1)*(Exr(1:pm-1,2:Ny)-Exr(1:pm-1,1:Ny-1))
	Hzxr(1,1:Ny-1)=dahzxr(1,1:Ny-1)*Hzxr(1,1:Ny-1)+dbhzxr(1,1:Ny-1)*
	1(-Eyr(2,1:Ny-1)+Ey(Nx,1:Ny-1)) 
      !=========================Hzx PML Left ==========================
	Hzxl(1:pm-2,1:Ny-1)=dahzxl(1:pm-2,1:Ny-1)*Hzxl(1:pm-2,1:Ny-1)+
	1dbhzxl(1:pm-2,1:Ny-1)*(-Eyl(2:pm-1,1:Ny-1)+Eyl(1:pm-2,1:Ny-1)) 
      Hzyl(1:pm-1,1:Ny-1)=dahzyl(1:pm-1,1:Ny-1)*Hzyl(1:pm-1,1:Ny-1)+
	1dbhzyl(1:pm-1,1:Ny-1)*(Exl(1:pm-1,2:Ny)-Exl(1:pm-1,1:Ny-1))  
	Hzxl(pm-1,1:Ny-1)=dahzxl(pm-1,1:Ny-1)*Hzxl(pm-1,1:Ny-1)+
	1dbhzxl(pm-1,1:Ny-1)*(-Ey(1,1:Ny-1)+Eyl(pm-1,1:Ny-1))
	!=========================Hzx PML UP ==========================  
	 Hzxu(1:ssx-1,1:pm-1)=dahzxu(1:ssx-1,1:pm-1)*Hzxu(1:ssx-1,1:pm-1)+
	1dbhzxu(1:ssx-1,1:pm-1)*(-Eyu(2:ssx,1:pm-1)+Eyu(1:ssx-1,1:pm-1))
	 Hzyu(1:ssx-1,2:pm-1)=dahzyu(1:ssx-1,2:pm-1)*Hzyu(1:ssx-1,2:pm-1)+
	1dbhzyu(1:ssx-1,2:pm-1)*(Exu(1:ssx-1,3:pm)-Exu(1:ssx-1,2:pm-1))
	 Hzyu(1:pm-1,1)=dahzyu(1:pm-1,1)*Hzyu(1:pm-1,1)+
	1dbhzyu(1:pm-1,1)*(Exu(1:pm-1,2)-Exl(1:pm-1,Ny))
	 Hzyu(pm:pm+Nx-2,1)=dahzyu(pm:pm+Nx-2,1)*Hzyu(pm:pm+Nx-2,1)+
	1dbhzyu(pm:pm+Nx-2,1)*(Exu(pm:pm+Nx-2,2)-Ex(1:Nx-1,Ny))
	 Hzyu(pm+Nx-1:ssx-1,1)=dahzyu(pm+Nx-1:ssx-1,1)*Hzyu(pm+Nx-1:ssx-1,
	11)+dbhzyu(pm+Nx-1:ssx-1,1)*(Exu(pm+Nx-1:ssx-1,2)-Exr(1:pm-1,Ny))
	!=========================Hz PML DOWN ==========================  
      Hzxd(1:ssx-1,1:pm-1)=dahzxd(1:ssx-1,1:pm-1)*Hzxd(1:ssx-1,1:pm-1)+
	1dbhzxd(1:ssx-1,1:pm-1)*(-Eyd(2:ssx,1:pm-1)+Eyd(1:ssx-1,1:pm-1))
      Hzyd(1:ssx-1,1:pm-2)=dahzyd(1:ssx-1,1:pm-2)*Hzyd(1:ssx-1,1:pm-2)+
	1dbhzyd(1:ssx-1,1:pm-2)*(Exd(1:ssx-1,2:pm-1)-Exd(1:ssx-1,1:pm-2))
	Hzyd(1:pm-1,pm-1)=dahzyd(1:pm-1,pm-1)*Hzyd(1:pm-1,pm-1)+
	1dbhzyd(1:pm-1,pm-1)*(Exl(1:pm-1,1)-Exd(1:pm-1,pm-1))
	Hzyd(pm:pm+Nx-2,pm-1)=dahzyd(pm:pm+Nx-2,pm-1)*
	1Hzyd(pm:pm+Nx-2,pm-1)+dbhzyd(pm:pm+Nx-2,pm-1)*(Ex(1:Nx-1,1)-
     1Exd(pm:pm+Nx-2,pm-1))
	Hzyd(pm+Nx-1:ssx-1,pm-1)=dahzyd(pm+Nx-1:ssx-1,pm-1)*
	1Hzyd(pm+Nx-1:ssx-1,pm-1)+dbhzyd(pm+Nx-1:ssx-1,pm-1)*
     1(Exr(1:pm-1,1)-Exd(pm+Nx-1:ssx-1,pm-1))
	!========================Distrete fourier transfrom =============
	Do nf=1,nfreq
	RD(nf,1)=sin(n*dt*OMEGA(nf))
	RD(nf,2)=cos(n*dt*OMEGA(nf))
      SN(nf,1)=SN(nf,1)+RD(nf,1)*pulse;  !the sine quadrature components  
      SN(nf,2)=SN(nf,2)+RD(nf,2)*pulse;  !the cosine quadrature components
	
	!=======================incident light DFT========================
	Hzincefr(nf,1)=Hzincefr(nf,1)+Hzinc(Nx-ntff)*RD(nf,1)
	Hzincefr(nf,2)=Hzincefr(nf,2)+Hzinc(Nx-ntff)*RD(nf,2)
	Hzincefl(nf,1)=Hzincefl(nf,1)+Hzinc(ntff)*RD(nf,1)
	Hzincefl(nf,2)=Hzincefl(nf,2)+Hzinc(ntff)*RD(nf,2)
	Hzincefu(nf,ntff:Nx-ntff,1)=Hzincefu(nf,ntff:Nx-ntff,1)+
	1RD(nf,1)*Hzinc(ntff:Nx-ntff);
	Hzincefu(nf,ntff:Nx-ntff,2)=Hzincefu(nf,ntff:Nx-ntff,2)+
	1RD(nf,2)*Hzinc(ntff:Nx-ntff);
	Hzincefd(nf,ntff:Nx-ntff,1)=Hzincefd(nf,ntff:Nx-ntff,1)+
	1RD(nf,1)*Hzinc(ntff:Nx-ntff);
	Hzincefd(nf,ntff:Nx-ntff,2)=Hzincefd(nf,ntff:Nx-ntff,2)+
	1RD(nf,2)*Hzinc(ntff:Nx-ntff);
 
	Eyincefr(nf,1)=Eyincefr(nf,1)+Eyinc(Nx-ntff)*RD(nf,1)/2.0+
	1Eyinc(Nx-ntff+1)*RD(nf,1)/2.0
	Eyincefr(nf,2)=Eyincefr(nf,2)+Eyinc(Nx-ntff)*RD(nf,2)/2.0+
	1Eyinc(Nx-ntff+1)*RD(nf,2)/2.0
	Eyincefl(nf,1)=Eyincefl(nf,1)+Eyinc(ntff)*RD(nf,1)/2.0+
	1Eyinc(ntff+1)*RD(nf,1)/2.0
	Eyincefl(nf,2)=Eyincefl(nf,2)+Eyinc(ntff)*RD(nf,2)/2.0+
	1Eyinc(ntff+1)*RD(nf,2)/2.0
	!====================================Hz===================================
  
	Hzefr(nf,ntff:Ny-ntff-1,1)=Hzefr(nf,ntff:Ny-ntff-1,1)+
	1RD(nf,1)*Hz(Nx-ntff,ntff:Ny-ntff-1)
	Hzefr(nf,ntff:Ny-ntff-1,2)=Hzefr(nf,ntff:Ny-ntff-1,2)+
	1RD(nf,2)*Hz(Nx-ntff,ntff:Ny-ntff-1)
 
	Hzefl(nf,ntff:Ny-ntff-1,1)=Hzefl(nf,ntff:Ny-ntff-1,1)+
	1RD(nf,1)*Hz(ntff,ntff:Ny-ntff-1);
	Hzefl(nf,ntff:Ny-ntff-1,2)=Hzefl(nf,ntff:Ny-ntff-1,2)+
	1RD(nf,2)*Hz(ntff,ntff:Ny-ntff-1);
 
	Hzefu(nf,ntff:Nx-ntff,1)=Hzefu(nf,ntff:Nx-ntff,1)+
	1RD(nf,1)*Hz(ntff:Nx-ntff,Ny-ntff-1)
	Hzefu(nf,ntff:Nx-ntff,2)=Hzefu(nf,ntff:Nx-ntff,2)+
	1RD(nf,2)*Hz(ntff:Nx-ntff,Ny-ntff-1)
	 
	Hzefd(nf,ntff:Nx-ntff,1)=Hzefd(nf,ntff:Nx-ntff,1)+
	1RD(nf,1)*Hz(ntff:Nx-ntff,ntff)
	Hzefd(nf,ntff:Nx-ntff,2)=Hzefd(nf,ntff:Nx-ntff,2)+
	1RD(nf,2)*Hz(ntff:Nx-ntff,ntff)
	!====================================Ex============================================
  
	Exefr(nf,ntff:Ny-ntff-1,1)=Exefr(nf,ntff:Ny-ntff-1,1)+
	1RD(nf,1)*Ex(Nx-ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,1)*Ex(Nx-ntff,ntff+1:Ny-ntff)/2.0;
	Exefr(nf,ntff:Ny-ntff-1,2)=Exefr(nf,ntff:Ny-ntff-1,2)+
	1RD(nf,2)*Ex(Nx-ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,2)*Ex(Nx-ntff,ntff+1:Ny-ntff)/2.0
 
	Exefl(nf,ntff:Ny-ntff-1,1)=Exefl(nf,ntff:Ny-ntff-1,1)+
	1RD(nf,1)*Ex(ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,1)*Ex(ntff,ntff+1:Ny-ntff)/2.0
	
	Exefl(nf,ntff:Ny-ntff-1,2)=Exefl(nf,ntff:Ny-ntff-1,2)+
	1RD(nf,2)*Ex(ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,2)*Ex(ntff,ntff+1:Ny-ntff)/2.0
	Do i=ntff,Nx-ntff
	Exefu(nf,i,1)=Exefu(nf,i,1)+
	1RD(nf,1)*Ex(i,Ny-ntff)/2.0+
     1RD(nf,1)*Ex(i,Ny-ntff-1)/2.0
	Exefu(nf,i,2)=Exefu(nf,i,2)+
	1RD(nf,2)*Ex(i,Ny-ntff)/2.0+
     1RD(nf,2)*Ex(i,Ny-ntff-1)/2.0
	END DO
	 
	 Exefd(nf,ntff:Nx-ntff,1)=Exefd(nf,ntff:Nx-ntff,1)+
	1RD(nf,1)*Ex(ntff:Nx-ntff,ntff)/2.0+
     1RD(nf,1)*Ex(ntff:Nx-ntff,ntff+1)/2.0
       Exefd(nf,ntff:Nx-ntff,2)=Exefd(nf,ntff:Nx-ntff,2)+
	1RD(nf,2)*Ex(ntff:Nx-ntff,ntff)/2.0+
     1RD(nf,2)*Ex(ntff:Nx-ntff,ntff+1)/2.0;
	!====================================Ey==================================
	Eyefr(nf,ntff:Ny-ntff-1,1)=Eyefr(nf,ntff:Ny-ntff-1,1)+
	1RD(nf,1)*Ey(Nx-ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,1)*Ey(Nx-ntff+1,ntff:Ny-ntff-1)/2.0
	
	Eyefr(nf,ntff:Ny-ntff-1,2)=Eyefr(nf,ntff:Ny-ntff-1,2)+
	1RD(nf,2)*Ey(Nx-ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,2)*Ey(Nx-ntff+1,ntff:Ny-ntff-1)/2.0
 
	Eyefl(nf,ntff:Ny-ntff-1,1)=Eyefl(nf,ntff:Ny-ntff-1,1)+
	1RD(nf,1)*Ey(ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,1)*Ey(ntff+1,ntff:Ny-ntff-1)/2.0
	Eyefl(nf,ntff:Ny-ntff-1,2)=Eyefl(nf,ntff:Ny-ntff-1,2)+
	1RD(nf,2)*Ey(ntff,ntff:Ny-ntff-1)/2.0+
     1RD(nf,2)*Ey(ntff+1,ntff:Ny-ntff-1)/2.0
 
	Eyefu(nf,ntff:Nx-ntff,1)=Eyefu(nf,ntff:Nx-ntff,1)+
	1RD(nf,1)*Ey(ntff:Nx-ntff,Ny-ntff-1)/2.0+
     1RD(nf,1)*Ey(ntff+1:Nx-ntff+1,Ny-ntff-1)/2.0;
	Eyefu(nf,ntff:Nx-ntff,2)= Eyefu(nf,ntff:Nx-ntff,2)+
	1RD(nf,2)*Ey(ntff:Nx-ntff,Ny-ntff-1)/2.0+
     1RD(nf,2)*Ey(ntff+1:Nx-ntff+1,Ny-ntff-1)/2.0;
 
	Eyefd(nf,ntff:Nx-ntff,1)=Eyefd(nf,ntff:Nx-ntff,1)+
	1RD(nf,1)*Ey(ntff:Nx-ntff,ntff)/2.0+
     1RD(nf,1)*Ey(ntff+1:Nx-ntff+1,ntff)/2.0;
	
	Eyefd(nf,ntff:Nx-ntff,2)= Eyefd(nf,ntff:Nx-ntff,2)+
	1RD(nf,2)*Ey(ntff:Nx-ntff,ntff)/2.0+
     1RD(nf,2)*Ey(ntff+1:Nx-ntff+1,ntff)/2.0;
 	END DO
	!========================FDTD LOOP END=================================
	if(MOD(n,10)==0) then
	write(*,*) n
	end if 
	END DO
	!========================LOOP END=============================
	!========================EXT calculation=============================
	Do nf=1,nfreq
	store=sqrt(SN(nf,1)**2+SN(nf,2)**2)
	SN(nf,2)=-atan2(SN(nf,1),SN(nf,2))
	SN(nf,1)=store
	!=======================incident light=========================
	   store=sqrt(Hzincefr(nf,1)**2+Hzincefr(nf,2)**2)
	   Hzincefr(nf,2)=-atan2(Hzincefr(nf,1),Hzincefr(nf,2))
	   Hzincefr(nf,1)=store
	   store=sqrt(Hzincefl(nf,1)**2+Hzincefl(nf,2)**2)
	   Hzincefl(nf,2)=-atan2(Hzincefl(nf,1),Hzincefl(nf,2))
	   Hzincefl(nf,1)=store
	   
	   store=sqrt(Eyincefr(nf,1)**2+Eyincefr(nf,2)**2)
	   Eyincefr(nf,2)=-atan2(Eyincefr(nf,1),Eyincefr(nf,2))
	   Eyincefr(nf,1)=store
  
	   store=sqrt(Eyincefl(nf,1)**2+Eyincefl(nf,2)**2)
	   Eyincefl(nf,2)=-atan2(Eyincefl(nf,1),Eyincefl(nf,2))
	   Eyincefl(nf,1)=store
		Do i=ntff,Nx-ntff
	store=sqrt(Hzincefu(nf,i,1)**2+Hzincefu(nf,i,2)**2)
	Hzincefu(nf,i,2)=-atan2(Hzincefu(nf,i,1),Hzincefu(nf,i,2))
	Hzincefu(nf,i,1)=store
   
	store=sqrt(Hzincefd(nf,i,1)**2+Hzincefd(nf,i,2)**2)
	Hzincefd(nf,i,2)=-atan2(Hzincefd(nf,i,1),Hzincefd(nf,i,2));
	Hzincefd(nf,i,1)=store;
   
		end DO
	!=============================Hz================================
	Do  i=ntff,Ny-ntff-1
	store=sqrt(Hzefr(nf,i,1)**2+Hzefr(nf,i,2)**2)
	Hzefr(nf,i,2)=-atan2(Hzefr(nf,i,1),Hzefr(nf,i,2));
	Hzefr(nf,i,1)=store
	store=sqrt(Hzefl(nf,i,1)**2+Hzefl(nf,i,2)**2)
	Hzefl(nf,i,2)=-atan2(Hzefl(nf,i,1),Hzefl(nf,i,2))
	Hzefl(nf,i,1)=store
	end DO
	Do i=ntff,Nx-ntff
	store=sqrt(Hzefu(nf,i,1)**2+Hzefu(nf,i,2)**2)
	Hzefu(nf,i,2)=-atan2(Hzefu(nf,i,1),Hzefu(nf,i,2))
	Hzefu(nf,i,1)=store
	store=sqrt(Hzefd(nf,i,1)**2+Hzefd(nf,i,2)**2)
	Hzefd(nf,i,2)=-atan2(Hzefd(nf,i,1),Hzefd(nf,i,2))
	Hzefd(nf,i,1)=store
	end DO
	!=============================Ex================================
	Do i=ntff+1,Ny-ntff-1
	store=sqrt(Exefr(nf,i,1)**2+Exefr(nf,i,2)**2)
	Exefr(nf,i,2)=-atan2(Exefr(nf,i,1),Exefr(nf,i,2))
	Exefr(nf,i,1)=store
	store=sqrt(Exefl(nf,i,1)**2+Exefl(nf,i,2)**2);
	Exefl(nf,i,2)=-atan2(Exefl(nf,i,1),Exefl(nf,i,2))
	Exefl(nf,i,1)=store
	end  Do 
	Do i=ntff,Nx-ntff
	store=sqrt(Exefu(nf,i,1)**2+Exefu(nf,i,2)**2)
	Exefu(nf,i,2)=-atan2(Exefu(nf,i,1),Exefu(nf,i,2))
	Exefu(nf,i,1)=store
      store=sqrt(Exefd(nf,i,1)**2+Exefd(nf,i,2)**2)
      Exefd(nf,i,2)=-atan2(Exefd(nf,i,1),Exefd(nf,i,2))
      Exefd(nf,i,1)=store
      end DO
	!============================Ey==================================
	Do i=ntff,Ny-ntff-1
	 store=sqrt(Eyefr(nf,i,1)**2+Eyefr(nf,i,2)**2)
	 Eyefr(nf,i,2)=-atan2(Eyefr(nf,i,1),Eyefr(nf,i,2))
       Eyefr(nf,i,1)=store
       store=sqrt(Eyefl(nf,i,1)**2+Eyefl(nf,i,2)**2)
       Eyefl(nf,i,2)=-atan2(Eyefl(nf,i,1),Eyefl(nf,i,2))
       Eyefl(nf,i,1)=store
	end  Do
        Do i=ntff+1,Nx-ntff
      store=sqrt(Eyefu(nf,i,1)**2+Eyefu(nf,i,1)**2)
      Eyefu(nf,i,2)=-atan2(Eyefu(nf,i,1),Eyefu(nf,i,2))
      Eyefu(nf,i,1)=store
      store=sqrt(Eyefd(nf,i,1)**2+Eyefd(nf,i,1)**2)
      Eyefd(nf,i,2)=-atan2(Eyefd(nf,i,1),Eyefd(nf,i,2))
      Eyefd(nf,i,1)=store
        end Do
	END DO
	!=========================spectrum calculation=================================
	 Do nf=1,nfreq

	jj=(0,1)
	jr=0
	jl=0
	ju=0
	jd=0
	Do i=ntff,Ny-ntff-1                            
       jr=jr+1/2.0*real(Eyincefr(nf,1)*exp(jj*Eyincefr(nf,2))*
	1Hzefr(nf,i,1)*exp(-jj*Hzefr(nf,i,2))+Eyefr(nf,i,1)*
     1exp(jj*Eyefr(nf,i,2))*Hzincefr(nf,1)*exp(-jj*Hzincefr(nf,2)))
       jl=jl-1/2.0*real(Eyincefl(nf,1)*exp(jj*Eyincefl(nf,2))*
	1Hzefl(nf,i,1)*exp(-jj*Hzefl(nf,i,2))+Eyefl(nf,i,1)*
     1exp(jj*Eyefl(nf,i,2))*Hzincefl(nf,1)*exp(-jj*Hzincefl(nf,2)));    
	end do
      Do i=ntff,Nx-ntff
	ju=ju-1/2.0*real(Exefu(nf,i,1)*exp(jj*Exefu(nf,i,2))*
	1Hzincefu(nf,i,1)*exp(-jj*Hzincefu(nf,i,2)))
	jd=jd+1/2.0*real(Exefd(nf,i,1)*exp(jj*Exefd(nf,i,2))*
	1Hzincefd(nf,i,1)*exp(-jj*Hzincefd(nf,i,2)))
	end DO
	 EXT(nf)=-((jr)+jl+ju+jd)*ds/
	1((SN(nf,1))**2*0.5*sqrt(muo/epso)*2*radius);
	end Do
      !==========================calculate EXT/SCA/ABS coefficients================
	Do nf=1,nfreq
	jr=0;
	jl=0;
	ju=0;
	jd=0;
	 Do i=ntff,Ny-ntff-1                           
       jr=jr+1/2.0*real(Eyefr(nf,i,1)*exp(jj*Eyefr(nf,i,2))*
	1Hzefr(nf,i,1)*exp(-jj*Hzefr(nf,i,2)))
       jl=jl-1/2.0*real(Eyefl(nf,i,1)*exp(jj*Eyefl(nf,i,2))*
	1Hzefl(nf,i,1)*exp(-jj*Hzefl(nf,i,2)))
      end do
	 Do i=ntff,Nx-ntff
       ju=ju-1/2.0*real(Exefu(nf,i,1)*exp(jj*Exefu(nf,i,2))*
	1Hzefu(nf,i,1)*exp(-jj*Hzefu(nf,i,2)))
       jd=jd+1/2.0*real(Exefd(nf,i,1)*exp(jj*Exefd(nf,i,2))*
	1Hzefd(nf,i,1)*exp(-jj*Hzefd(nf,i,2)))
	end do
	ABSR(nf) = -((jr)+jl+ju+jd)*ds/
	1((SN(nf,1))**2*0.5*sqrt(muo/epso)*2*radius);
	end do
	
	Do nf=1,nfreq
	SCA(nf)=EXT(nf)-ABSR(nf)
	end do
	Do nf=1,nfreq
	waves(nf)=1D9*2*pi*c/OMEGA(nf);
	end do
	!============================END=======================

	write(24,42) Ext
42    format(E18.9' '\)
	write(25,43) SCA
43    format(E18.9' '\)
	write(26,44) ABSR
44    format(E18.9' '\)
	write(27,45) waves
45    format(E18.9' '\)
	Do i=1,Nx
	 Do j=1,Ny
	write(23,41) cbex(i,j)
41    format(E18.6' '\)
	end Do
	write(23,' ') 
	end do
	STOP 
		end
