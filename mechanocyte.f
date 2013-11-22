      program mechanocyte 

      use iolibsw
      use cmlibsw

      character ipf*20, opf*20, datafile*20
      character fnum*3, suff*3, pref*3
      real(8) epsl, avdt, cmdt(12),volint, tdump(100),hmax,rmean
      real(8) logInt,area,volume,height
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
      open(31,file=datafile,status='new')
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

      idphi=0 !phase rub
      idvis=0 !viscosity
      idvfx=0 !fixed velocity at boundaries
      idgam=0 !surface tension
      idpsi=0 !network contractility
      idsfr=0 !surface forces
      idbfr=0 !body forces
      iddrg=0 !drag at boundaries
      idtrc=0 !boundary tractions
      idhyc=0 !boundary solvent permeability
      !call normalizeCoords
      call initsvec
      call setSurfaceField
      call initMessengerConc
      open(unit=12,file='mechanocyte.init')
      call iowrfile(0,12)
      close(12)
      idebug=1
  100 continue
         call getArea(area)
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
         time=time+tstp
         call getVolume(volume)
         call getSurfaceMolecules(12,PIP2m)
         call getSurfaceMolecules(11,PIP3m)
         call getSurfaceMolecules(10,PIP3a)
         call getSurfaceMolecules(9,PTENm)
         call getSurfaceMolecules(8,PI3Km)
         call getVolumeMolecules(2,MessV)
         call getSurfaceMolecules(2,MessS)
         call getVolumeMolecules(1,ThetV)
         call getSurfaceMolecules(1,ThetS)
         height = volume/area
         print *,"time:",time!,"area:",area,"volume:",volume
         print *,"PIP2m:",int(PIP2m),"PI3Km:",int(PI3Km),"PIP3m:",
     1           int(PIP3m),"PIP3a:",int(PIP3a),"PTENm:",int(PTENm),
     2           "MessV:",MessV,"MessStoV:",MessS*height
         print *,"ThetV:",ThetV,"ThetStoV:",ThetS*height
         print *,"maxMess:",maxval(svec(2,1,1:ns)),
     1           maxval(svec(2,2,1:ns)),maxval(svec(2,3,1:ns))
         print *,"minMess:",minval(svec(2,1,1:ns)),
     1           minval(svec(2,2,1:ns)),minval(svec(2,3,1:ns))
         print *,"maxNet:",maxval(svec(1,1,1:ns)),
     1           maxval(svec(1,2,1:ns)),maxval(svec(1,3,1:ns))
         print *,"minNet:",minval(svec(1,1,1:ns)),
     1           minval(svec(1,2,1:ns)),minval(svec(1,3,1:ns))
         open(unit=12,file='test_out')
         call iowrfile(0,12)
         close(12)
         if (isve.eq.1) then
            isuff=isuff+1
            opf=ipf
            write(fnum,90)isuff
            opf(idot+1:idot+3)=fnum
            open(unit=21,file=opf,status='unknown')
            call iowrfile(0,21)
            close(21)
!
!--check if we have reached the end
!
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
            hvec(1,lvn,isn) = hvec(1,lvn,isn)*1d-5
            hvec(2,lvn,isn) = hvec(2,lvn,isn)*1d-5
            hvec(3,lvn,isn) = hvec(3,lvn,isn)*1d-5
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
