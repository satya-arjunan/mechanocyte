      program mechanocyte 

      use iolibsw
      use molibsw
      use avlibsw
      use cmlibsw

      character ipf*20, opf*20, datafile*20
      character fnum*3, suff*3, pref*3
      real(8) epsl,avdt,cmdt(12),tdump(100),hmax,rmean
      real(8) logInt,area,volume,height,r0,a0
      real(8) PIP2m,PIP3m,PIP3a,PI3Km,PTENm,MessV,MessS,ThetV,ThetS
!
!--read command line to initialize name of input file
!
      call getarg(1,ipf)
      idot=index(ipf,'.')
      if (idot.eq.0) then
         write(*,*)'invalid inputfile name'
         stop
      endif
!--extract numeric field
      fnum=ipf(idot+1:idot+3)
      isuff=(iachar(fnum(1:1))-48)*100+
     1      (iachar(fnum(2:2))-48)*10+
     2       iachar(fnum(3:3))-48
      print *,'isuff=',isuff
!
!--open,read and close the input file
!
      
      open(unit=11,file=ipf,status='old')
      call iordfile(iframe,11)
      close(11)
!--setup data log file
      datafile=ipf
      datafile(idot+1:idot+3)='dat'
      open(31,file=datafile,status='unknown')
!
!--figure out time dumps
!
      logInt=0
      delt=dscp(3,idfrm)
      logInt=dscp(4,idfrm)
!--if delt is zero look for time of dumps in file tdump
      if (delt.eq.zero) then
         open(44,file='tdump',status='old')
         read(44,*) ndump
         do i=1,ndump
            read(44,*) tdump(i)
         enddo
         tstop=tdump(ndump)
c--look for current dump position
         if (time.gt.tstop) then
            print *,'time>tstop',time,tstop
            stop
         else
            do i=1,ndump
               if (tdump(i).gt.time) then
                  idump=i
                  goto 20
               endif
            enddo
         endif
   20    print *,'ndump,idump,tstop',ndump,idump,tstop
         tnext=tdump(idump)
      else
         tstop=time+delt
         tnext=time+logInt
      endif
      print *,'tnext,tstop',tnext,tstop

      idphi=1 !phase rub
      idvis=1 !viscosity
      idvfx=1 !fixed velocity at boundaries
      idgam=1 !surface tension
      idpsi=0 !network contractility
      idsfr=1 !surface forces
      idbfr=0 !body forces
      iddrg=0 !drag at boundaries
      idtrc=0 !boundary tractions
      idhyc=0 !boundary solvent permeability
      !call normalizeCoords
      !call writeInitMeshFile
      !stop
      call initsvec
      call setSurfaceField
      call initMessengerConc
      idebug=0
      call setPhaserub
      call setViscosity
      call getVolume(volume)
      call getArea(area)
      r0=(volume*3.0/(4.0*3.141592))**(1./3.)
      a0=4.0*3.141592*r0**2
      aratio=area/a0
      call setSurfaceTension(aratio)
      call setSlipCondition
      call setSurfaceForce
      epsl=1d-4
      call printFile(isuff, ipf, 0)
  10  call modriver(icyc,epsl,idebug)
      print *,"init icycle:",icyc
      if (icyc.gt.10) goto 10
  100 continue
         call clchm(area)
         call cmstpsiz(cmdt,1d-2)
         tstp=min(cmdt(2),cmdt(4))
!        tstp=cmdt(2)
         if (tstp.ge.(tnext-time)) then
            tstp=tnext-time
            isve=1
         else
            isve=0
         endif
         call dfdriver('eulerian')
         !call avgridmo('lagrangian',0,idebug)
         call avgridmo('lagrangian',idebug)
         call avfield(1,'nw')
         time=time+tstp
         print *,"time:",time
         call setViscosity
         call getVolume(volume)
         call getArea(area)
         r0=(volume*3.0/(4.0*3.141592))**(1./3.)
         a0=4.0*3.141592*r0**2
         aratio=area/a0
         call setSurfaceTension(aratio)
         call setSurfaceForce
         height = volume/area
         open(unit=12,file='test_out')
         do 
           call modriver(icyc,epsl,idebug)
           if (icyc.lt.10) exit
         enddo
         call iowrfile(0,12)
         close(12)
         if (isve.eq.1) then
           call printFile(isuff,ipf,0)
           if (time.ge.tstop*0.999) then
              print *,'time,tstop',time,tstop
              close(31)
              stop
           else
              idump=idump+1
              if (logInt.eq.0) then
                 tnext=tdump(idump)
              else
                 tnext=tnext+logInt
              endif
           endif
         endif
         goto 100
   90 format(i3.3)
      stop
      end


      subroutine printFile(isuff, ipf, isStop)
      use iolibsw
      character ipf*20, opf*20, fnum*3
      idot=index(ipf,'.')
      isuff=isuff+1
      opf=ipf
      write(fnum,90)isuff
      opf(idot+1:idot+3)=fnum
      print *, "file name:", opf
      print *,"time:",time
            hvec(1,lvn,isn) = hvec(1,lvn,isn)*1d-5
      print *,"min x, y, z hvec:",minval(hvec(1,1:3,1:ns)),
     1        minval(hvec(2,1:3,1:ns)),minval(hvec(3,1:3,1:ns))
      print *,"max x, y, z hvec:",maxval(hvec(1,1:3,1:ns)),
     1        maxval(hvec(2,1:3,1:ns)),maxval(hvec(3,1:3,1:ns))
      print *,"Curve min svec(4) l:1,2,3:",minval(svec(4,1,1:ns)),
     1        minval(svec(4,2,1:ns)),minval(svec(4,3,1:ns))
      print *,"Curve max svec(4) l:1,2,3:",maxval(svec(4,1,1:ns)),
     1        maxval(svec(4,2,1:ns)),maxval(svec(4,3,1:ns))
      print *,"MSGR min svec(2) l:1,2,3:",minval(svec(2,1,1:ns)),
     1        minval(svec(2,2,1:ns)),minval(svec(2,3,1:ns))
      print *,"MSGR max svec(2) l:1,2,3:",maxval(svec(2,1,1:ns)),
     1        maxval(svec(2,2,1:ns)),maxval(svec(2,3,1:ns))
      print *,"ThetaN min svec(1) l:1,2,3:",minval(svec(1,1,1:ns)),
     1        minval(svec(1,2,1:ns)),minval(svec(1,3,1:ns))
      print *,"ThetaN max svec(1) l:1,2,3:",maxval(svec(1,1,1:ns)),
     1        maxval(svec(1,2,1:ns)),maxval(svec(1,3,1:ns))
      print *,"viscos min vis() l:1,2,3:",minval(vis(1,1:ns)),
     1        minval(vis(2,1:ns)),minval(vis(3,1:ns))
      print *,"viscos max vis() l:1,2,3:",maxval(vis(1,1:ns)),
     1        maxval(vis(2,1:ns)),maxval(vis(3,1:ns))
      print *,"phaserub min phi() l:1,2,3:",minval(phi(1,1:ns)),
     1        minval(phi(2,1:ns)),minval(phi(3,1:ns))
      print *,"phaserub max phi() l:1,2,3:",maxval(phi(1,1:ns)),
     1        maxval(phi(2,1:ns)),maxval(phi(3,1:ns))
      print *,"surf tension min gam() v,d,e:",minval(gamv(1:nq)),
     1        minval(gamd(1:nq)),minval(game(1:nl))
      print *,"surf tension max gam() v,d,e:",maxval(gamv(1:nq)),
     1        maxval(gamd(1:nq)),maxval(game(1:nl))
      print *,"surf force min sfr v,d,e:",minval(sfrv(1:nq)),
     1        minval(sfrd(1:nq)),minval(sfre(1:nl))
      print *,"surf force max sfr v,d,e:",maxval(sfrv(1:nq)),
     1        maxval(sfrd(1:nq)),maxval(sfre(1:nl))
      open(unit=21,file=opf,status='unknown')
      call iowrfile(0,21)
      close(21)
      if(isStop.eq.1) then
        stop
      endif
   90 format(i3.3)
      end subroutine printFile

      subroutine writeInitMeshFile
      use iolibsw
      open(unit=12,file='mechanocyte.init')
      call iowrfile(0,12)
      close(12)
      end subroutine writeInitMeshFile

      subroutine setViscosity
      use iolibsw
      !increase viscosity to slow down the rounding
      viscosity = 4.39d13
      do isn=1,ns
         do lvn=1,3
            vis(lvn,isn)=max(viscosity*svec(4,lvn,isn), 1238d11)
         enddo
      enddo
      return
      end subroutine setViscosity

      subroutine setPhaserub
      use iolibsw
      phaserub = 1d8
      do isn=1,ns
         do lvn=1,3
            phi(lvn,isn)=max(phaserub*svec(4,lvn,isn), 97414015d1)
         enddo
      enddo
      return
      end subroutine setPhaserub


      subroutine setSurfaceTension(afac)
      use iolibsw 
      !increase surface_tension to increase the rounding
      !decrease surface tension to increase time steps
      surface_tension = 0.0111d-2
      afac3=afac**3
      do iq=1,nq
         gamv(iq)=surface_tension*afac3
         gamd(iq)=surface_tension*afac3
      enddo
      do il=1,nl
         game(il)=surface_tension*afac3
      enddo
      return
      end subroutine setSurfaceTension

      subroutine setSlipCondition
      USE iolibsw
      !vfixv is the ventral element's slip condition
      do iq=1,nq
         vfixv(1,iq)=1d0 !fix (stick) ventral elements in x direction
         vfixv(2,iq)=1d0 !fix (stick) ventral elements in y direction
         vfixv(3,iq)=1d0 !fix (stick) ventral elements in z direction
      enddo

      do il=1,nl
         iq = iqol(il)
         vfixv(1,iq)=0d0 !free (slip) ventral edge elements in x direction
         vfixv(2,iq)=0d0 !free (slip) ventral edge elements in y direction
      enddo
      end subroutine setSlipCondition

      subroutine setSurfaceForce !clsfr
      use iolibsw
      real(8) ThetaEq,Theta0,tau_n
      iThetaN = 12
      Theta0 = 1d-3
      tau_n = 1d1
      iMess = 2
      do iq=1,nq
         sfrv(iq)=0
         dorsal=0d0
         do isn=1,4
            is=isoq(isn,iq)
            !ThetaEq=Theta0*(1d0+svec(iMESS,3,is))
            dorsal=dorsal+svec(iThetaN,3,is)
         enddo
         sfrd(iq)=0.25*dorsal*1d-8
      enddo
      do il=1,nl
         edge=0d0
         do isn=1,2
            is=isol(isn,il)
            do lv=1,3
               !ThetaEq=Theta0*(1d0+svec(iMESS,lv,is))
               edge=edge+svec(iThetaN,lv,is)
            enddo
         enddo
         sfre(il)=edge/6d0*9d-7
      enddo
      return
      end subroutine setSurfaceForce

      subroutine setSurfaceField
      USE iolibsw
      real(8) ventral,dorsal,edge
      do iq=1,nq
         ventral=0d0
         dorsal=0d0
         do isn=1,4
            is=isoq(isn,iq)
            ventral = ventral+snn(0,1,is)/arean(1,is)
            dorsal = dorsal+snn(0,3,is)/arean(3,is) 
         enddo
         vvec(1,iq)=0.25*ventral
         dvec(1,iq)=0.25*dorsal
      enddo
      do il=1,nl
         edge=0d0
         do isn=1,2
            is=isol(isn,il)
            do lv=1,3
               edge=edge+snn(0,lv,is)/arean(lv,is)
            enddo
         enddo
         evec(1,il)=edge/6d0
      enddo
      end subroutine setSurfaceField

      
      subroutine initMessengerConc 
      USE iolibsw
      do isn=1,ns
         do lvn=1,3
            svec(2,lvn,isn) = 0 !messenger
            svec(4,lvn,isn) = 0 !curvature
            if(arean(lvn,isn).gt.0) then
               !set the messenger to be ~20 at the highly curved region:
               svec(4,lvn,isn) = snn(0,lvn,isn)/arean(lvn,isn)*6.6173d-5
               svec(2,lvn,isn) = svec(4,lvn,isn)
            endif 
         enddo
      enddo
      end subroutine initMessengerConc


      subroutine getArea(area)
      USE iolibsw
      real(8) area,surfintv,surfintd,surfinte
      real(8),dimension(3,NSM)::cnode=1d0
      call gosurfintn(cnode,surfintv,surfintd,surfinte)
      area = surfintv+surfintd+surfinte
      end subroutine getArea


      subroutine getVolume(volume)
      USE iolibsw
      real(8) volume
      real(8),dimension(3,NSM)::cnode=1d0
      call govolint(cnode,volume)
      end subroutine getVolume


      subroutine getSurfaceMolecules(kom, cnt)
      USE iolibsw
      real(8) cnt,surfintv,surfintd,surfinte
      real(8) cnode(3,NSM)
      do i=1,ns
         cnode(1,i)=svec(kom,1,i)
         cnode(2,i)=svec(kom,2,i)
         cnode(3,i)=svec(kom,3,i)
      enddo
      call gosurfintn(cnode,surfintv,surfintd,surfinte)
      cnt = surfintv+surfintd+surfinte
      end subroutine getSurfaceMolecules


      subroutine getVolumeMolecules(kom, cnt)
      USE iolibsw
      real(8) cnt
      real(8) cnode(3,NSM)
      do i=1,ns
         cnode(1,i)=svec(kom,1,i)
         cnode(2,i)=svec(kom,2,i)
         cnode(3,i)=svec(kom,3,i)
      enddo
      call govolint(cnode,cnt)
      end subroutine getVolumeMolecules


      subroutine normalizeCoords
      USE iolibsw
      integer in
      real(8) x,y
      do isn=1,ns
         do lvn=1,3
            hvec(1,lvn,isn) = hvec(1,lvn,isn)*0.8
            hvec(2,lvn,isn) = hvec(2,lvn,isn)*0.8
            hvec(3,lvn,isn) = hvec(3,lvn,isn)*0.5
         enddo
      enddo
      end subroutine normalizeCoords


      subroutine initsvec
      USE iolibsw
      integer in
      real(8) x,y
      do isn=1,ns
         do lvn=1,3
            x = (hvec(1,lvn,isn))
            y = (hvec(2,lvn,isn))
            !if(x.lt.5d-6) svec(12,lvn,isn) = 1.444605d12 !PIP2m
            if(y.lt.2.5d-6) then
               svec(12,lvn,isn) = 1.277d13 !PIP2m
               svec(9,lvn,isn) = 0.627d12 !PTENm
            endif
            if(x.lt.-6d-6) then
               svec(8,lvn,isn) = 3.8546d13 !PI3Km
            endif
         enddo
      enddo
      end subroutine initsvec


      subroutine clchm(area)
      use iolibsw
      integer in,iPIP2,iPIP3,iPIP3a,iPI3K,iPTEN,iMESS,iCURV,iThetaN
      real(8) PIPtc,PTENtc,PI3Ktc
      real(8) PIP2m,PIP3m,PIP3a,PI3Km,PTENm,MESS,CURV,ThetaN
      real(8) PIP2,PI3K,PTEN,maxConc,new,area,val
      real(8) k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12
      real(8) ThetaEq,Theta0,tau_n

      iPIP2 = 12
      iPIP3 = 11
      iPIP3a = 10
      iPTEN = 9
      iPI3K = 8
      iCURV = 4
      iMESS = 2
      iThetaN = 1

      k1 = 4d-2
      k2 = 2d-14
      k3 = 2d-13
      k4 = 4d-14
      k5 = 2d-14
      k6 = 5d-14
      k7 = 5d-14
      k8 = 0.09
      k9 = 0.02
      k10 = 0.0001
      k11 = 1d0 !Messenger source rate
      k12 = 1d0 !Messenger decay rate = 1/tau_m
      tau_n = 1d0
      Theta0 = 1d-3


      maxConc = 8.7816d13
      PIPtc = 1.089171d13
      PTENtc = 7.62666d12
      PI3Ktc = 1.44957d13

      call getSurfaceMolecules(iPIP2,val)
      PIP2 = val
      call getSurfaceMolecules(iPIP3,val)
      PIP2 = PIP2 + val
      call getSurfaceMolecules(iPIP3a,val)
      PIP2 = PIP2 + val
      PIP2 = PIPtc-PIP2/area 

      call getSurfaceMolecules(iPI3K,PI3K)
      PI3K = PI3Ktc-PI3K/area

      call getSurfaceMolecules(iPTEN,PTEN)
      PTEN = PTENtc-PTEN/area
c      print *,"PTEN:",int(PTEN*area),"PI3K:",int(PI3K*area),
c     1        "PIP2:",int(PIP2*area)

      do isn=1,ns
         do lvn=1,3
            PIP2m = svec(iPIP2,lvn,isn)
            PIP3m = svec(iPIP3,lvn,isn)
            PIP3a = svec(iPIP3a,lvn,isn)
            PTENm = svec(iPTEN,lvn,isn)
            PI3Km = svec(iPI3K,lvn,isn)
            MESS = svec(iMESS,lvn,isn)
            CURV = svec(iCURV,lvn,isn)
            ThetaN = svec(iThetaN,lvn,isn)

            !Source and decay messenger
            sdot(iMESS,lvn,isn) = k11*CURV-k12*MESS
            sdkr(iMESS,lvn,isn) = -k12

            !Create network
            ThetaEq=Theta0*(1d0+MESS)
            sdot(iThetaN,lvn,isn) = (ThetaEq-ThetaN)*MESS/tau_n
            sdkr(iThetaN,lvn,isn) = -MESS/tau_n



            sdot(iPIP2,lvn,isn) = k1*PIP2-k5*PIP2m*PI3Km+k6*PIP3m*PTENm+
     1                            k7*PIP3a*PTENm-k10*PIP2m
            sdkr(iPIP2,lvn,isn) = -k5*PI3Km-k10
            new = PIP2m + sdot(iPIP2,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPIP2,lvn,isn) = 0
              sdkr(iPIP2,lvn,isn) = 0
            endif

            sdot(iPTEN,lvn,isn) = k2*PTEN*PIP2m-k8*PTENm
            sdkr(iPTEN,lvn,isn) = -k8
            new = PTENm + sdot(iPTEN,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPTEN,lvn,isn) = 0
              sdkr(iPTEN,lvn,isn) = 0
            endif

            sdot(iPIP3a,lvn,isn) = -k3*PIP3a*PI3K+k4*PIP3m*PIP3m-
     1                             k7*PIP3a*PTENm-k9*PIP3a
            sdkr(iPIP3a,lvn,isn) = -k3*PI3k-k7*PTENm-k9
            new = PIP3a + sdot(iPIP3a,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPIP3a,lvn,isn) = 0
              sdkr(iPIP3a,lvn,isn) = 0
            endif

            sdot(iPIP3,lvn,isn) = k3*PIP3a*PI3K-k4*PIP3m*PIP3m+
     1                             k5*PIP2m*PI3Km-k6*PIP3m*PTENm-
     2                             k9*PIP3m
            sdkr(iPIP3,lvn,isn) = -2*k4*PIP3m-k6*PTENm-k9
            new = PIP3m + sdot(iPIP3,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPIP3,lvn,isn) = 0
              sdkr(iPIP3,lvn,isn) = 0
            endif

            sdot(iPI3K,lvn,isn) = k3*PIP3a*PI3K-k9*PI3Km
            sdkr(iPI3K,lvn,isn) = -k9
            new = PI3Km + sdot(iPI3K,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPI3K,lvn,isn) = 0
              sdkr(iPI3K,lvn,isn) = 0
            endif
         enddo
      enddo
      end subroutine clchm
