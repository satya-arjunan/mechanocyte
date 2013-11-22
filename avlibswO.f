      MODULE AVLIBSW
!
! library of advection/mesh motion routines.
!
      USE IOLIBSW
      USE lupack
!
      real(8),save::xzr(3,3,NSM)!initial node positions
      real(8),save::xnw(3,3,NSM)!node positions after
                                !lagrangian motion by vnw
      real(8),save::xca(3,3,NSM)!node positions after
                                !lagrangian motion by vvol
      real(8),save::xfn(3,3,NSM)!final node positions
      real(8),save::dilfacv(3,NSM)!volume dilation coefficient for node
      real(8),save::dilfaca(3,NSM)!area dilation coefficient for node
      logical,save::avnw
      logical,save::avca
      integer,save::clfix(NSM)
!
      contains!the following subroutines used in advection operations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
      subroutine avgridmo(opt,idebug,iflag)
!
! This subroutine updates the mesh
!
      character*(*) opt !='eul', or 'lag'
      integer idebug
      integer iflag
!
      real(8) a_in, a_out
      real(8) volzr(3,NSM), areazr(3,NSM)
!--move 10% along surface of advected mesh, 90% with rezoned
!--mesh (along tangent volume plane).
      real(8) cfn, cbk
      parameter(cbk=0.05d0) 
      parameter(cfn=1d0-cbk)
      real(8) xback(3,3,NSM)
      integer il,iscl
!
!
!---LAGRANGIAN MOTION OF THE MESH
!
! load initial node positions and update according to velocities
! recall hvec(1:3,:,:) = x-y-z positions
!        hvec(7:9,:,:) = vnx-vny-vnz network velocities
!        hvec(4:6,:,:) = vvrlx-vvrly-vvrlz volume vel. wrt to network
! 
      iflag=0 
      xzr=hvec(1:3,:,:)
      xnw=xzr+tstp*hvec(7:9,:,:)
      xca=xzr+tstp*(hvec(7:9,:,:)+hvec(4:6,:,:))
      avnw=.false.
      avca=.false.
!
!--compute volume and area dilation for network motion
!--(volume motion dilation is zero if there is no flux at bndry)
!
!--calculate and save initial volumes and areas
      call godriver(xzr)
      volzr=voln
      areazr=arean
      call godriver(xnw)
!--compute the dilation factors
      do is=1,ns
         dilfacv(1,is)=volzr(1,is)/voln(1,is)
         dilfacv(2,is)=volzr(2,is)/voln(2,is)
         dilfacv(3,is)=volzr(3,is)/voln(3,is)
         dilfaca(1,is)=areazr(1,is)/arean(1,is)
         dilfaca(3,is)=areazr(3,is)/arean(3,is)
         if (ilos(1,is).ne.0) then !an edge stack
            dilfaca(2,is)=areazr(2,is)/arean(2,is)
         else
            dilfaca(2,is)=0d0
         endif
      enddo
!
      if (index(opt,'eul').gt.0) then !mesh does not move.
         xfn=xzr(:,:,:) !return mesh to original position
      else                            !grid motion is allowed
         xfn=xnw(:,:,:) !mesh follows the network
!
!--BEGIN NON-LAGRANGIAN MOTION OF THE MESH
!
!
         call rezdriver(idebug,xfn,iflag)!rezone the grid
!
!--check that ventral nodes are in z=0 plane
!-- and that dorsal nodes are not too low
         do is=1,ns
            if (abs(xfn(3,1,is)).ge.1d-9) then
               print *,'xfn(3,1,is),is',xfn(3,1,is),is
               print *,'hvec(3,1,is),tstp',hvec(3,1,is),tstp
               stop
            endif
            if (xfn(3,3,is).lt.1d-9) then
              print *,'xfn(3,3,is),is',xfn(3,3,is),is
              print *,'hvec(3,3,is),hvec(9,3,is),tstp',
     1                 hvec(3,3,is),hvec(9,3,is),tstp
              stop
            endif
         enddo
!
!--get all the rezone back surface nodes positions interpolated
!--onto the advected grid
!
         call avback(xnw,xfn,xback)
!
         do is=1,ns
            hvec(1:3,:,is)=cfn*xfn(:,:,is)+cbk*xback(:,:,is)
         enddo
!
      endif
!
      return
      end subroutine avgridmo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezdriver(idebug,xn,iflag)
!
! this subroutine is the driver for the rezoning of the mesh
! 1) the node positions are readjusted vertically so that 
!    all ventral nodes have z=0
! 2) the middle edge nodes are checked for high or low z
!    for low z, the CL is shifted AND
!    for low z or high z, the middle node is shifted (edgadv)
! 3) slide CL to maintain ventral/middle nodes wrt. substratum
!    within limits defined by a_in and a_out. (edgslide)
! 4) the contact line nodes are rezoned to be 
!    equidistant (rezcl)
! 5) middle and dorsal edge nodes are slid so as
!    to form a plane with the ventral (CL) node that
!    is perpendicular to the substratum (rezedg)
! 6) the interior ventral nodes are rezoned
!    to minimize the Winslow functional (rezbot)
! 7) the dorsal nodes and middle edge nodes are rezoned
!    to minimize the Winslow functional (rezback)
! 8) the interior middle nodes are rezoned
!    to minimize the TTM functional (rezone)
!
      implicit none
!
      integer idebug !debug flag
      real(8) xn(3,3,NSM) !
      integer iflag
!
      real(8) cnode(3,NSM)
      real(8) reztol
      parameter(reztol=1d-3)
      integer icmx
      parameter(icmx=200)
      integer icyc,is,lv,nslid,nadv
      real(8) c_newton,dratio,volint,wtot1,wtot2,dwtot,dd
      real(8) a_in, a_out
!
      cnode=1d0
!
! enforce substratum node positions while trying to preserve volume
! this is done because the penalty method for BC leads to very small
! but non-zero velocities.
         do is=1,ns 
            xn(3,3,is)=xn(3,3,is)-xn(3,1,is)
            xn(3,2,is)=xn(3,2,is)-xn(3,1,is)
            xn(3,1,is)=0d0 
         enddo
!
!--deal with the position of middle edge nodes, move CL if needed
         call edgadv(idebug,xn,nadv)
         if (nadv.ne.0) call godriver(xn)
!
!--slide CL to maintain a_in and a_out angles
         a_in=0d0; a_out=0d0
         if (a_in.gt.0d0.or.a_out.gt.0d0) then
            call edgslide(idebug,a_out,a_in,xn,nslid)
            if (nslid.ne.0) call godriver(xn)
         endif
!
!--rezone contact line (ventral edge nodes)
!
      do icyc=1,icmx
         call rezcl(idebug,xn,dd)
         if (dd.lt.0.1d0*reztol.and.icyc.gt.3) exit
         call godriver(xn)
      enddo
      if (idebug.ge.1) print *,'rezcl: dd, icyc',dd,icyc
      if (icyc.ge.icmx) then 
         print *,'rezcl warning',dd
         iflag=1
      endif
!
!--rezone the middle and dorsal edge nodes
!
      call rezedg(idebug,xn)
      call godriver(xn)
!
!--rezone the ventral interior nodes
!
      c_newton=0.5d0 !under-relaxation 
      wtot1=0d0
      do icyc=1,icmx
        call rezbot(idebug,xn,c_newton,wtot2)
        call godriver(xn)
        dwtot=abs(wtot2-wtot1)/wtot2
        if (dwtot.lt.reztol/nq.and.icyc.gt.3) exit
        wtot1=wtot2
      enddo
      if (idebug.ge.1) print *,'rezbot: wtot2, icyc',wtot2,icyc
      if (icyc.ge.icmx) then 
         print *,'rezbot warning',wtot2
         iflag=1
      endif
!
!--rezone the dorsal and middle edge nodes
!
      c_newton=0.5d0 !under-relaxation 
      wtot1=0d0
      do icyc=1,icmx
        call rezback(idebug,xn,c_newton,wtot2)
        call godriver(xn)
        dwtot=abs(wtot2-wtot1)/wtot2
        if (dwtot.lt.reztol/nq.and.icyc.gt.3) exit
        wtot1=wtot2
      enddo
      if (idebug.ge.1) print *,'rezback: wtot2, icyc',wtot2,icyc
      if (icyc.ge.icmx) then
         print *,'rezback warning',wtot2
         iflag=1
      endif
!
!--rezone the middle interior nodes
!
      c_newton=0.5d0 !under-relaxation 
      wtot1=0d0
      do icyc=1,icmx
        call rezone(idebug,xn,c_newton,wtot2)
        call godriver(xn)
        dwtot=abs(wtot2-wtot1)/wtot2
        if (dwtot.lt.reztol/nq.and.icyc.gt.3) exit
        wtot1=wtot2
      enddo
      if (idebug.ge.1) print *,'rezone: wtot2, icyc',wtot2,icyc
      if (icyc.ge.icmx) then 
         print *,'rezone warning',wtot2
         iflag=1
      endif
!
      return
      end subroutine rezdriver 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         BEGIN SUBROUTINES CALLED BY REZDRIVER
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine edgadv(idebug,xn,nadv)
!
! if 1,2,3 are the ventral, middle, and dorsal points of an
! edge stack, and z1=0, then 
!
!        if z2<0.25*z3 then 
!            stack line (1,2,3) goes through the substratum
!            => we move 1 to the other intersection 
!            => we move 2 to z2=0.25*z3  
!
!        if z2>0.75*z3 then 
!           stack line (1,2,3) goes above 3
!           this may lead to problems so we move 2 to z2=0.75*z3
!
      implicit none
!
      integer idebug !debug flag
      real(8) xn(3,3,NSM) !node positions
      integer nadv !number of nodes moved
! 
      real(8) x1,y1,x2,y2,x3,y3,z2,z3,z2new,dx1,dy1,dx2,dy2
      real(8) zeta1, zeta2, H1, H2, H3
      real(8) Lshape,dLdzeta,d2Ldzeta2 !shape function
      integer il, is
!
      nadv=0
      do il=1,nl
         is=isol(1,il)
         z2=xn(3,2,is) 
         z3=xn(3,3,is) 
         if (z3.lt.0d0) then
! we are screwed because the dorsal point is below the substrate
            print *,'edgadv: dorsal point below substrate!',is,z3
            stop
         elseif (z2.lt.0d0) then !we are also screwed
            print *,'edgadv: middle point below substrate!',is,z2
            stop
         elseif (z2.lt.0.25*z3) then  !middle point is low
            x1=xn(1,1,is); y1=xn(2,1,is) 
            x2=xn(1,2,is); y2=xn(2,2,is) 
            x3=xn(1,3,is); y3=xn(2,3,is) 
!--compute the intrinsic coord zeta1 of the 2nd intersection with
!  the substratum: z=(1-zeta**2)*z2+0.5*zeta*(zeta+1)*z3=0
!                    (1+zeta)[z2-zeta(z2-z3/2)=0
            zeta1=1d0/(1d0-0.5d0*z3/z2)
!--get the line shape functions at zeta1
            call getL3(-1,zeta1,H1,dLdzeta,d2Ldzeta2)
            call getL3(0,zeta1,H2,dLdzeta,d2Ldzeta2)
            call getL3(+1,zeta1,H3,dLdzeta,d2Ldzeta2)
!--displacement
            dx1=(H1*x1+H2*x2+H3*x3)-x1
            dy1=(H1*y1+H2*y2+H3*y3)-y1
!--check if the motion is outward or inward by taking the 
!  dot product with the contact line outward normal
            if (idebug.ge.2) then
               if ((dx1*cnn(1,il)+dy1*cnn(2,il)).ge.0d0) then
                  print *,'edgadv: CL moved outward!',is,zeta1
               else
                  print *,'edgadv: CL moved inward!',is,zeta1
               endif
            endif
            xn(1,1,is)=dx1+x1
            xn(2,1,is)=dy1+y1
!-- move middle point such that z2 back to z3/4
            z2new=0.25d0*z3
            call findzeta(0d0,z2,z3,z2new,zeta2)
!--get the line shape functions at zeta2
            call getL3(-1,zeta2,H1,dLdzeta,d2Ldzeta2)
            call getL3(0,zeta2,H2,dLdzeta,d2Ldzeta2)
            call getL3(+1,zeta2,H3,dLdzeta,d2Ldzeta2)
            xn(1,2,is)=H1*x1+H2*x2+H3*x3
            xn(2,2,is)=H1*y1+H2*y2+H3*y3
            xn(3,2,is)=H2*z2+H3*z3
            nadv=nadv+1
         elseif (z2.gt.0.75d0*z3) then !middle point is high
            x1=xn(1,1,is); y1=xn(2,1,is) 
            x2=xn(1,2,is); y2=xn(2,2,is) 
            x3=xn(1,3,is); y3=xn(2,3,is) 
!-- move middle node such that z2 back to (3/4)z3
            z2new=0.75d0*z3
            call findzeta(0d0,z2,z3,z2new,zeta2)
!--get the line shape functions at zeta2
            call getL3(-1,zeta2,H1,dLdzeta,d2Ldzeta2)
            call getL3(0,zeta2,H2,dLdzeta,d2Ldzeta2)
            call getL3(+1,zeta2,H3,dLdzeta,d2Ldzeta2)
            dx2=(H1*x1+H2*x2+H3*x3)-x2
            dy2=(H1*y1+H2*y2+H3*y3)-y2
            xn(1,2,is)=dx2+x2
            xn(2,2,is)=dy2+y2
            xn(3,2,is)=z2new
            if (idebug.ge.2) then
               if ((dx2*cnn(1,il)+dy2*cnn(2,il)).ge.0d0) then
                  print *,'edgadv: middle node moved outward!',is,zeta2
               else
                  print *,'edgadv: middle node moved inward!',is,zeta2
               endif
            endif
            nadv=nadv+1
         endif
      enddo
!
      if (nadv.ne.0.and.idebug.ge.1) 
     1               print *,'edgadv: # of nodes moved:',nadv
!         
      return 
      end subroutine edgadv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine edgslide(idebug,a_out,a_in,xn,nmov)
!
! This subroutine moves nodes at the CL so as to maintain 
! a certain 1->2/substrate  angle.
!
      implicit none
!
      integer idebug !debug flag
      real(8) a_out, a_in !angles for overhang and underhang cases
      real(8) xn(3,3,NSM) !node positions
      integer nmov
!
      real(8) dxy12,dz12,dx12,dy12,dxy13,dz13,dx1,dy1,dx2,dy2,dx3,dy3
      real(8) dx13,dy13,ddxy13
      real(8) dV1,dV2,dV3,dV21,alpha,alpha3,tanalph,a3min,ddotn
      real(8) Bmat(3,3),rhs(3),Bvec(3)
      real(8) d2,d3,clnx,clny,c2,c3,z2,z3,det
      integer il,is
!
      nmov=0
      do il=1,nl
         is=isol(1,il)
         dx12=xn(1,2,is)-xn(1,1,is)
         dy12=xn(2,2,is)-xn(2,1,is)
         ddotn=dx12*cnn(1,il)+dy12*cnn(2,il)
         dxy12=sqrt(dx12**2+dy12**2)
         dz12=xn(3,2,is)-xn(3,1,is)
         alpha=atan(dz12/max(dxy12,vtiny))
         if (alpha.lt.a_out.and.ddotn.gt.0d0) then
!--case where the angle with substrate is smaller than a_out,
!  and the cell has an overhang
!  move node 1 outward to match a_out, 
!  move node 2 inward to conserve volume
            if (abs(dx12).gt.1d-6*dxy12) then
               dx1=0.5d0*(dxy12**2-(dz12/tan(a_out))**2)/
     1                     (dx12+dy12**2/dx12)
               dy1=dx1*dy12/dx12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=-(dV1/(vnn(0,2,is)))/
     1              (vnn(1,2,is)+vnn(2,2,is)*dy12/dx12)
               dy2=dx2*dy12/dx12
            else
               dx1=0d0
               dy1=0.5d0*(dxy12**2-(dz12/tan(a_out))**2)/dy12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=0d0
               dy2=-(dV1/(vnn(0,2,is)))/vnn(2,2,is)
            endif
!           dV2=vnn(0,2,is)*(vnn(1,2,is)*dx2+vnn(2,2,is)*dy2)
            xn(1,1,is)=xn(1,1,is)+dx1
            xn(2,1,is)=xn(2,1,is)+dy1
            xn(1,2,is)=xn(1,2,is)+dx2
            xn(2,2,is)=xn(2,2,is)+dy2
            nmov=nmov+1
            if (idebug.ge.2) then
               print *,'slid bottom edge node is,alpha',is,alpha
            endif
         elseif (alpha.lt.a_in.and.ddotn.lt.0d0) then
!--case where the angle with substrate is smaller than a_in,
!  and the cell has an underhang
!  move node 1 inward to match a_in, 
!  move node 2 outward to conserve volume
            if (abs(dx12).gt.1d-6*dxy12) then
               dx1=0.5d0*(dxy12**2-(dz12/tan(a_in))**2)/
     1                     (dx12+dy12**2/dx12)
               dy1=dx1*dy12/dx12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=-(dV1/(vnn(0,2,is)))/
     1              (vnn(1,2,is)+vnn(2,2,is)*dy12/dx12)
               dy2=dx2*dy12/dx12
            else
               dx1=0d0
               dy1=0.5d0*(dxy12**2-(dz12/tan(a_in))**2)/dy12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=0d0
               dy2=-(dV1/(vnn(0,2,is)))/vnn(2,2,is)
            endif
            dV2=vnn(0,2,is)*(vnn(1,2,is)*dx2+vnn(2,2,is)*dy2)
            xn(1,1,is)=xn(1,1,is)+dx1
            xn(2,1,is)=xn(2,1,is)+dy1
            xn(1,2,is)=xn(1,2,is)+dx2
            xn(2,2,is)=xn(2,2,is)+dy2
            nmov=nmov+1
            if (idebug.ge.2) then
               print *,'retracted bottom edge node is,alpha',is,alpha
            endif
         endif   
      enddo
!
      if (idebug.ge.1.and.nmov.ne.0) then
         print *,nmov,' ventral CL nodes were moved'
      endif
!
      return
      end subroutine edgslide
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezcl(idebug,xn,dd)
!
! This subroutine moves each ventral edge nodes toward the "midpoint"
! between its two neighbors:
! Let the original edge point be 0, the original neighbors A and B
! The target edge point O' is given by the intersection of
! 1) the substratum plane z=0
! 2) the volume conservation plane
!    (x-x0)*vnn_x+(y-y0)*vnn_y+(z-z0)*vnn_z=0
! 1) and 2) (z0=0) give the line (x-x0)*vnn_x+(y-y0)*vnn_y=0; z=0
! 3) the median between A and B
!       [x-(xA+xB)/2](xB-xA)+[y-(yA+yB)/2](yB-yA)=0
!
!  we move a fraction of the way to the target 
!
! User specified:
!   IMPORTANT: at the end of the routine clfix is zeroed out.
!
      implicit none
!
      integer idebug
      real(8) xn(3,3,NSM) !in old, out new node x-y coordinates
      real(8) dd
!
      real(8) dxtg(2,NLM) !displ. of 2nd node of il in tg volume plane
      real(8) dxcl(2,NLM) !displ. of 2nd node of il in contact line
      real(8) cdx,ctg,ccl !coefficients for final displacement
      parameter(cdx=0.25d0)!25% move to target median 
      parameter(ctg=0.9d0)!90% along volume tangent plane
      parameter(ccl=1d0-ctg)!10% along contact line
!
      real(8) x,y,xA,yA,xB,yB,dx,dy,det2,xM,yM,rhs1,rhs2,xt,yt,sint
      real(8) dAB,dxA,dyA,dxB,dyB,dA2,dB2,disp
      integer il,isA,is,isB,lv
!
      do il=1,nl !loop over edges
         isA=isol(1,il)
         is=isol(2,il)
         isB=isol(2,ilol(2,il))
         x=xn(1,1,is)
         y=xn(2,1,is)
         xA=xn(1,1,isA)
         yA=xn(2,1,isA)
         xB=xn(1,1,isB)
         yB=xn(2,1,isB)
         dx=xB-xA
         dy=yB-yA
         dAB=sqrt(dx**2+dy**2)
         sint=sqrt(vnn(1,1,is)**2+vnn(2,1,is)**2)
         xM=0.5*(xA+xB)
         yM=0.5*(yA+yB)
!
!--work on the tangent plane displacement
         if (sint.lt.5d-2) then
! the tangent plane and the substratum essentialy coincide,
! (angle <= 3deg) we solve
! (xnew-x)(yA-yB)+(ynew-y)(xB-xA)=0
! (xnew-xM)(xB-xA)+(ynew-yM)(yB-yA)=0
            det2=-dy**2-dx**2
            rhs1=-x*dy+y*dx
            rhs2=xM*dx+yM*dy
            xt=(rhs1*dy-rhs2*dx)/det2
            yt=(-rhs2*dy-rhs1*dx)/det2
            dxtg(1,il)=xt-x
            dxtg(2,il)=yt-y
         else
! solve:
!       xnew*(xB-xA)+ynew*(yB-yA)=0.5*[(xA+xB)*(xB-xA)+(yA+yB)*(yB-yA)]
! xnew*vnn(1,1,is)+ynew*vnn(2,1,is)=x*vnn(1,1,is)+y*vnn(2,1,is)
            det2=dx*vnn(2,1,is)-dy*vnn(1,1,is)
!check that the plane/substrate intersection and the median are not //
            if (abs(det2).gt.2d-2*dAB) then
               rhs1=xM*dx+yM*dy
               rhs2=x*vnn(1,1,is)+y*vnn(2,1,is)
               xt=(rhs1*vnn(2,1,is)-rhs2*dy)/det2
               yt=(rhs2*dx-rhs1*vnn(1,1,is))/det2
               dxtg(1,il)=xt-x
               dxtg(2,il)=yt-y
            else
               print *,'warning in rezcl: '
               print *,'il,is,dAB,det2',il,is,dAB,det2
               print *,'node number=',(is-1)*3+1
               print *,'vn,vny,vnz',vnn(1:3,1,is)
               print *,'dx,dy',dx,dy
               print *,'sint',sint
               print *,'return to continue'
               read *
            endif
         endif
!
!--work on the displacement along the contact line 
         dxA=x-xA
         dyA=y-yA
         dxB=x-xB
         dyB=y-yB
         dA2=dxA**2+dyA**2
         dB2=dxB**2+dyB**2
         if (dA2.lt.dB2) then !closer to A so we move along the B line
!--compute intersection of B line with median
!   xt*(xB-xA)+yt*(yB-yA)=0.5*[(xA+xB)*(xB-xA)+(yA+yB)*(yB-yA)]       
!   xt*(yB-y)+yt*(x-xB)=x*(yB-y)+y*(x-xB)
            det2=dx*dxB+dy*dyB
            if (abs(det2).lt.1d-16*dB2) then
               print *,'error in rezcl'
               print *,'det2,dB2',det2,dB2
               print *,'is,il',is,il
               stop
            endif
            rhs1=xM*dx+yM*dy
            rhs2=x*(yB-y)+y*(x-xB)
            xt=(rhs1*dxB-rhs2*dy)/det2
            yt=(rhs2*dx+rhs1*dyB)/det2
            dxcl(1,il)=xt-xn(1,1,is)
            dxcl(2,il)=yt-xn(2,1,is)
         elseif(dB2.lt.dA2) then!closer to B so we move along the A line
!--compute intersection of B line with median
!   xt*(xB-xA)+yt*(yB-yA)=0.5*[(xA+xB)*(xB-xA)+(yA+yB)*(yB-yA)]       
!   xt*(yA-y)+yt*(x-xA)=x*(yA-y)+y*(x-xA)
            det2=dx*dxA+dy*dyA
            if (abs(det2).lt.1d-16*dA2) then
               print *,'error in rezcl'
               print *,'det2,dA2',det2,dA2
               print *,'is,il',is,il
               stop
            endif
            rhs1=xM*dx+yM*dy
            rhs2=x*(yA-y)+y*(x-xA)
            xt=(rhs1*dxA-rhs2*dy)/det2
            yt=(rhs2*dx+rhs1*dyA)/det2
            dxcl(1,il)=xt-x
            dxcl(2,il)=yt-y
         else !we are already at the mid point
            dxcl(1,il)=0d0
            dxcl(2,il)=0d0
         endif
      enddo
!
!--check for cases where the node should not be moved
      do il=1,nl
         is=isol(2,il)
         if (clfix(is).gt.0d0) then
            dxtg(1,il)=0d0; dxtg(2,il)=0d0
            dxcl(1,il)=0d0; dxcl(2,il)=0d0
         endif
      enddo
!
!     do is=1,ns
!        hvec(1:3,:,is)=xn(1:3,:,is)
!        hvec(4:6,:,is)=vnn(1:3,:,is)
!        hvec(7:9,:,is)=0d0
!     enddo
!     do il=1,nl
!        is=isol(2,il)
!        hvec(7,1,is)=cdx*(ctg*dxtg(1,il)+ccl*dxcl(1,il))
!        hvec(8,1,is)=cdx*(ctg*dxtg(2,il)+ccl*dxcl(2,il))
!     enddo
!     open(66,file='av.dump')
!     call iowrfile(99,66)
!     close(66)
!
! implement displacements
!
      dd=0d0
      do il=1,nl
         is=isol(2,il)
         isA=isol(1,il)
         isB=isol(2,ilol(2,il))
         dAB=sqrt((xn(1,1,isA)-xn(1,1,isB))**2+
     1            (xn(2,1,isA)-xn(2,1,isB))**2)
         dx=cdx*(ctg*dxtg(1,il)+ccl*dxcl(1,il))
         dy=cdx*(ctg*dxtg(2,il)+ccl*dxcl(2,il))
         dd=max(dd,sqrt(dx**2+dy**2)/dAB)
         xn(1,1,is)=xn(1,1,is)+dx
         xn(2,1,is)=xn(2,1,is)+dy
      enddo
!
      if (idebug.ge.2) print *,'rezcl: dd=',dd
!
      return
      end subroutine rezcl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezedg(idebug,xnew)
!
! this subroutine rezones the middle and dorsal edge nodes
!
!     the dorsal and middle nodes are shifted along their
!     tangent planes so as to lie in the plane perpendicular
!     to the substrate z=0 and containing contact line normal.
!     The displacement keeps z constant
!
      implicit none
!
      integer idebug
      real(8) xnew(3,3,NSM)
!
      real(8) xold(3,3,NSM)
      real(8) det2,rhs1,rhs2
      real(8) xv,yv,xo,yo,zo,zn,tnx,tny,tnz
      real(8) eiznx,eizny,eiznz,teizx,teizy,teizz,eizx,eizy,eizz
      real(8) xA,yA,zA,xB,yB,zB
      real(8) xe,ye,clnx,clny
      integer il,is

!--back up copy of xnew
      xold=xnew
!--loop over edge stacks
      do il=1,nl
         is=isol(1,il)
!--edge stack ventral coord (z=0)
         xe=xold(1,1,is)
         ye=xold(2,1,is)
!--normal to the contact line at isol(1,il)
         clnx=cnn(1,il)
         clny=cnn(2,il)
!--we slide the dorsal edge node along its tangent plane to the
!  plane perpendicular that contains the contact line normal.
!  We require that z remain constant
         xA=xold(1,3,is)
         yA=xold(2,3,is)
         zA=xold(3,3,is)
!--normal to the tangent plane
         tnx=vnn(1,3,is)
         tny=vnn(2,3,is)
         tnz=vnn(3,3,is)
! new position xB, yB, zB given by:
! clny*(xB-xe)-clnx*(yB-ye)=0 (in cl normal perpendicular plane) 
! (xB-xA)*tnx+(yB-yA)*tny+(zB-zA)*tnz=0 (in tangent plane)
! zB=zA
         det2=clny*tny+clnx*tnx
         if (abs(det2).gt.vtiny) then
            rhs1=clny*xe-clnx*ye
            rhs2=xA*tnx+yA*tny
            xnew(1,3,is)=(rhs1*tny+rhs2*clnx)/det2
            xnew(2,3,is)=(rhs2*clny-rhs1*tnx)/det2
            xnew(3,3,is)=zA
         elseif (idebug.ge.1) then
            print *,'rezedg warning:',det2,is
         endif
!
!--do the same for the middle edge node
         xA=xold(1,2,is)
         yA=xold(2,2,is)
         zA=xold(3,2,is)
!--normal to the tangent plane
         tnx=vnn(1,2,is)
         tny=vnn(2,2,is)
         tnz=vnn(3,2,is)
! new position xB, yB, zB given by:
! clny*(xB-xe)-clnx*(yB-ye)=0 (in cl normal perpendicular plane) 
! (xB-xA)*tnx+(yB-yA)*tny+(zB-zA)*tnz=0 (in tangent plane)
! zB=zA
         det2=clny*tny+clnx*tnx
         if (abs(det2).gt.vtiny) then
            rhs1=clny*xe-clnx*ye
            rhs2=xA*tnx+yA*tny
            xnew(1,2,is)=(rhs1*tny+rhs2*clnx)/det2
            xnew(2,2,is)=(rhs2*clny-rhs1*tnx)/det2
            xnew(3,2,is)=zA
         elseif (idebug.ge.1) then
            print *,'rezedg warning:',det2,is
         endif
      enddo
!
      return
      end subroutine rezedg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezbot(idebug,xn,c_newton,wtot)
!
! This subroutine rezones the interior nodes of the 
! ventral surface so as to minimize the Winslow functional
!
      implicit none
!
      integer idebug
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton !Newton method relaxation coefficient
!
      real(8) deps !maximum node motion parameter
      parameter(deps=1d-1)
      real(8) dxn(3,3,NSM)
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23 !metric coeff.
      real(8) dS1,dS2
      real(8) F,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcb, Jcbx,Jcby,Jcbz, dJ
      real(8) Win(NSM),Wtot,Wg
      real(8) Wx(3,NSM),Wy(3,NSM),Wz(3,NSM) !grad of Winslow functional
      real(8) Wxx(3,NSM),Wyy(3,NSM),Wzz(3,NSM) !Hessian of Winslow
      real(8) Wxy(3,NSM),Wxz(3,NSM),Wyz(3,NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
      real(8) HessW(3,3), gradW(3), dvec(3),gradWn
      real(8) HessWtg11,HessWtg12,HessWtg22
      real(8) gradWtg1, gradWtg2, dtg1,dtg2,ddet
      real(8) dmax,dmean
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dscale
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
!--Initialize the W's
      Wtot=0d0
      Wx=0d0; Wy=0d0; Wz=0d0
      Wxx=0d0; Wyy=0d0; Wzz=0d0
      Wxy=0d0; Wxz=0d0; Wyz=0d0
!
!--initialize the displacements
      dxn=0d0
!
      dz1=0d0
      dz2=0d0
      do iq=1,nq !loop over quads
         do isg=1,4 !loop over the ventral gauss points of iq
            dx1=dxSv(1,isg,iq)
            dx2=dxSv(2,isg,iq)
            dy1=dySv(1,isg,iq)
            dy2=dySv(2,isg,iq)
            call Winslow(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F)
            Jcb=dAv(0,isg,iq)
            dJ=1d0/Jcb
            Wg=F*dJ
            Wtot=Wtot+Wg
            do isn=1,4
               istack=isoq(isn,iq)
               dS1=dS4(1,isn,isg)
               dS2=dS4(2,isn,isg)
             call dWin(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F,Jcb,dJ,dS1,dS2,
     3                        Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
               Wx(1,istack)=Wx(1,istack)+Wxi
               Wy(1,istack)=Wy(1,istack)+Wyi
               Wxx(1,istack)=Wxx(1,istack)+Wxxi
               Wyy(1,istack)=Wyy(1,istack)+Wyyi
               Wxy(1,istack)=Wxy(1,istack)+Wxyi
            enddo
         enddo
      enddo
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with c_newton under-relaxation)
!
      do is=1,ns !loop over ventral nodes
         if (ilos(1,is).eq.0) then!not an edge stack
            gradWtg1=-Wx(1,is)
            gradWtg2=-Wy(1,is)
            HessWtg11=Wxx(1,is)
            HessWtg22=Wyy(1,is)
            HessWtg12=Wxy(1,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,1,is)=c_newton*dtg1
            dxn(2,1,is)=c_newton*dtg2
         endif 
      enddo
!
!--write an information dump if required
      if (idebug.ge.2) then      
         if (mod(kount,50).eq.0) then
            do is=1,ns
               hvec(1:3,:,is)=xn(1:3,:,is)
               hvec(4:6,:,is)=vnn(1:3,:,is)
               hvec(7:9,:,is)=dxn(1:3,:,is)
            enddo
            open(66,file='rezbot.dump')
            call iowrfile(99,66)
            close(66)
            print *,'dumped to rezbot.dump',kount
         endif
         kount=kount+1
      endif
! 
!--limiting the displacements
!
      dmean=0d0; dmax=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            dist=sqrt(dxn(1,1,is)**2+dxn(2,1,is)**2)
            dmax=max(dmax,dist)
            dmean=dist+dmean
            dscale=sqrt(abs(arean(1,is)))
            if (dist.gt.deps*dscale) then
               fac=deps*dist/dscale
               dxn(1,1,is)=dxn(1,1,is)*fac
               dxn(2,1,is)=dxn(2,1,is)*fac
            endif
         endif
      enddo
      dmean=dmean/(ns-nl)
      if (idebug.ge.1) print *,'rezbot: Wtot,dmax,dmean',
     1                                  Wtot,dmax,dmean
!
!--move the ventral nodes
      do is=1,ns
         xn(1,1,is)=xn(1,1,is)+dxn(1,1,is)
         xn(2,1,is)=xn(2,1,is)+dxn(2,1,is)
      enddo
!
      return
      end subroutine rezbot
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezback(idebug,xn,c_newton,wtot)
!
! This subroutine rezones the back surface nodes (dorsal
! and middle edge)so as to minimize the Winslow functional
!
      implicit none
!
      integer idebug
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton
!
      real(8) deps !maximum node motion parameter
      parameter(deps=1d-1)
      real(8) dxn(3,3,NSM)
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23 !metric coeff.
      real(8) dS1,dS2
      real(8) F,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcb, Jcbx,Jcby,Jcbz, dJ
      real(8) Wtot,Wg
      real(8) Wx(3,NSM),Wy(3,NSM),Wz(3,NSM) !grad of Winslow functional
      real(8) Wxx(3,NSM),Wyy(3,NSM),Wzz(3,NSM) !Hessian of Winslow
      real(8) Wxy(3,NSM),Wxz(3,NSM),Wyz(3,NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
      real(8) HessW(3,3), gradW(3), dvec(3),gradWn
      real(8) HessWtg11,HessWtg12,HessWtg22
      real(8) gradWtg1, gradWtg2, dtg1,dtg2,ddet
      real(8) u1x,u1y,u1z,u2x,u2y,u2z
      real(8) B(11,11),BLU(11,11),dis(11),vrhs(11)
      integer ipvt(11), isng
      real(8) z2min,z3min
      integer is2min,is3min
      real(8) dmax,dmean
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dist2,dist3,dscale,dscale2,dscale3
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
!--Initialize the W's
      Wtot=0d0
      Wx=0d0; Wy=0d0; Wz=0d0
      Wxx=0d0; Wyy=0d0; Wzz=0d0
      Wxy=0d0; Wxz=0d0; Wyz=0d0
!
!--zero out the displacements
      dxn=0d0
!
!--loop over the dorsal surfaces
      do iq=1,nq 
         do isg=1,4 !loop over the gauss points of iq
            dx1=dxSd(1,isg,iq)
            dx2=dxSd(2,isg,iq)
            dy1=dySd(1,isg,iq)
            dy2=dySd(2,isg,iq)
            dz1=dzSd(1,isg,iq)
            dz2=dzSd(2,isg,iq)
            call Winslow(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F)
            Jcb=dAd(0,isg,iq)
            dJ=1d0/Jcb
            Wg=F*dJ
            Wtot=Wtot+Wg
            do isn=1,4
               istack=isoq(isn,iq)
               dS1=dS4(1,isn,isg)
               dS2=dS4(2,isn,isg)
             call dWin(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F,Jcb,dJ,dS1,dS2,
     3                        Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
               Wx(3,istack)=Wx(3,istack)+Wxi
               Wy(3,istack)=Wy(3,istack)+Wyi
               Wz(3,istack)=Wz(3,istack)+Wzi
               Wxx(3,istack)=Wxx(3,istack)+Wxxi
               Wyy(3,istack)=Wyy(3,istack)+Wyyi
               Wzz(3,istack)=Wzz(3,istack)+Wzzi
               Wxy(3,istack)=Wxy(3,istack)+Wxyi
               Wxz(3,istack)=Wxz(3,istack)+Wxzi
               Wyz(3,istack)=Wyz(3,istack)+Wyzi
            enddo
         enddo
      enddo
!
!--loop over the edge surfaces
      do il=1,nl
         do isg=1,2
            do lvg=1,3
               dx1=dxSe(1,lvg,isg,il)
               dx2=dxSe(2,lvg,isg,il)
               dy1=dySe(1,lvg,isg,il)
               dy2=dySe(2,lvg,isg,il)
               dz1=dzSe(1,lvg,isg,il)
               dz2=dzSe(2,lvg,isg,il)
               call Winslow(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F)
               Jcb=dAe(0,lvg,isg,il)
               dJ=1d0/Jcb
               Wg=F*dJ*wgp(lvg)
               Wtot=Wtot+Wg
               do isn=1,2
                  istack=isol(isn,il)
                  do lvn=1,3
                     dS1=dS6(1,lvn,isn,lvg,isg)
                     dS2=dS6(2,lvn,isn,lvg,isg)
                     call dWin(dx1,dx2,dy1,dy2,dz1,dz2,
     1                         g11,g22,F,Jcb,dJ,dS1,dS2,
     3                         Wxi,Wyi,Wzi,
     3                         Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
                     Wx(lvn,istack)=Wx(lvn,istack)+Wxi*wgp(lvg)
                     Wy(lvn,istack)=Wy(lvn,istack)+Wyi*wgp(lvg)
                     Wz(lvn,istack)=Wz(lvn,istack)+Wzi*wgp(lvg)
                     Wxx(lvn,istack)=Wxx(lvn,istack)+Wxxi*wgp(lvg)
                     Wyy(lvn,istack)=Wyy(lvn,istack)+Wyyi*wgp(lvg)
                     Wzz(lvn,istack)=Wzz(lvn,istack)+Wzzi*wgp(lvg)
                     Wxy(lvn,istack)=Wxy(lvn,istack)+Wxyi*wgp(lvg)
                     Wxz(lvn,istack)=Wxz(lvn,istack)+Wxzi*wgp(lvg)
                     Wyz(lvn,istack)=Wyz(lvn,istack)+Wyzi*wgp(lvg)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with 50% under-relaxation)
!
      z2min=vbig
      z3min=vbig
      do is=1,ns
         if (ilos(1,is).eq.0) then!not an edge stack
!--only dorsal nodes move
!--we project the gradients on vectors tangent to the surface
            u1x=vtg1(1,3,is)
            u1y=vtg1(2,3,is)
            u1z=vtg1(3,3,is)
            u2x=vtg2(1,3,is)
            u2y=vtg2(2,3,is)
            u2z=vtg2(3,3,is)
            gradWtg1=-Wx(3,is)*u1x-Wy(3,is)*u1y-Wz(3,is)*u1z
            gradWtg2=-Wx(3,is)*u2x-Wy(3,is)*u2y-Wz(3,is)*u2z
            HessWtg11=+u1x*u1x*Wxx(3,is)+u1y*u1y*Wyy(3,is)
     1                +u1z*u1z*Wzz(3,is)
     2           +2d0*(u1x*u1y*Wxy(3,is)+u1x*u1z*Wxz(3,is)+
     3                 u1y*u1z*Wyz(3,is))
            HessWtg22=+u2x*u2x*Wxx(3,is)+u2y*u2y*Wyy(3,is)
     1                +u2z*u2z*Wzz(3,is)
     2           +2d0*(u2x*u2y*Wxy(3,is)+u2x*u2z*Wxz(3,is)+
     3                 u2y*u2z*Wyz(3,is))
            HessWtg12=+u1x*u2x*Wxx(3,is)+u1y*u2y*Wyy(3,is)
     1                +u1z*u2z*Wzz(3,is)+
     2                (u1x*u2y+u1y*u2x)*Wxy(3,is)+
     3                (u1x*u2z+u1z*u2x)*Wxz(3,is)+
     4                (u1y*u2z+u1z*u2y)*Wyz(3,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,3,is)=c_newton*(dtg1*u1x+dtg2*u2x)
            dxn(2,3,is)=c_newton*(dtg1*u1y+dtg2*u2y)
            dxn(3,3,is)=c_newton*(dtg1*u1z+dtg2*u2z)
         else !an edge stack
!-middle and dorsal node move in the intersection
! of the tangent plane and the CL normal plane
! with the constraint that dz3/z3=dz2/z2
!
!--Lagrange multiplier method
!
            il=ilos(1,is)
            B=0d0
            B(1,1)=Wxx(2,is); B(1,2)=Wxy(2,is); B(1,3)=Wxz(2,is)
            B(2,1)=Wxy(2,is); B(2,2)=Wyy(2,is); B(2,3)=Wyz(2,is)
            B(3,1)=Wxz(2,is); B(3,2)=Wyz(2,is); B(3,3)=Wzz(2,is)
            B(1,4)=vnn(1,2,is);B(2,4)=vnn(2,2,is);B(3,4)=vnn(3,2,is)
            B(4,1)=vnn(1,2,is);B(4,2)=vnn(2,2,is);B(4,3)=vnn(3,2,is)
            B(1,5)=cnn(2,il);B(2,5)=-cnn(1,il)
            B(5,1)=cnn(2,il);B(5,2)=-cnn(1,il)
            B(6,6)=Wxx(3,is); B(6,7)=Wxy(3,is); B(6,8)=Wxz(3,is)
            B(7,6)=Wxy(3,is); B(7,7)=Wyy(3,is); B(7,8)=Wyz(3,is)
            B(8,6)=Wxz(3,is); B(8,7)=Wyz(3,is); B(8,8)=Wzz(3,is)
            B(6,9)=vnn(1,3,is);B(7,9)=vnn(2,3,is);B(8,9)=vnn(3,3,is)
            B(9,6)=vnn(1,3,is);B(9,7)=vnn(2,3,is);B(9,8)=vnn(3,3,is)
            B(6,10)=cnn(2,il);B(7,10)=-cnn(1,il)
            B(10,6)=cnn(2,il);B(10,7)=-cnn(1,il)
            B(11,3)=xn(3,3,is);B(11,8)=-xn(3,2,is)
            B(3,11)=xn(3,3,is);B(8,11)=-xn(3,2,is)
            call lufa(B,11,BLU,isng,ipvt)
            if (isng.ne.0) print *,
     1         'ill conditioned matrix in rezone, isng=',isng
            vrhs=0d0
            vrhs(1)=-Wx(2,is); vrhs(2)=-Wy(2,is); vrhs(3)=-Wz(2,is)
            vrhs(6)=-Wx(3,is); vrhs(7)=-Wy(3,is); vrhs(8)=-Wz(3,is)
            call lusl(BLU,vrhs,isng,ipvt,11,dis)
            dxn(1,2,is)=c_newton*dis(1)
            dxn(2,2,is)=c_newton*dis(2)
            dxn(3,2,is)=c_newton*dis(3)
            dxn(1,3,is)=c_newton*dis(6)
            dxn(2,3,is)=c_newton*dis(7)
            dxn(3,3,is)=c_newton*dis(8)
!           if (xn(3,3,is).lt.z3min) then
!              is3min=is
!              z3min=xn(3,3,is)
!           endif
!           if (xn(3,2,is).lt.z2min) then
!              is2min=is
!              z2min=xn(3,2,is)
!           endif
         endif 
      enddo
!     print *,'is3min,z,dz',is3min,xn(3,3,is3min),dxn(3,3,is3min)
!     print *,'is3min,z0,dzv',is3min,hvec(3,3,is3min),
!    1                        hvec(9,3,is3min)*tstp
!     print *,'is2min,z,dz',is2min,xn(3,2,is2min),dxn(3,2,is2min)
!
! for debugging might want a dump
      if (idebug.ge.2) then
         if (mod(kount,50).eq.0) then
            do is=1,ns
               hvec(1:3,:,is)=xn(1:3,:,is)
               hvec(4:6,:,is)=vnn(1:3,:,is)
               hvec(7:9,:,is)=dxn(1:3,:,is)
            enddo
            open(66,file='rezback.dump')
            call iowrfile(99,66)
            close(66)
            hvec(7:9,:,1:ns)=0d0
         endif
         kount=kount+1
      endif
!
!--update the node positions while limiting the displacements
!
      dmean=0d0; dmax=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            dist=sqrt(dxn(1,3,is)**2+dxn(2,3,is)**2+dxn(3,3,is)**2)
            dmax=max(dmax,dist)
            dmean=dist+dmean
            dscale=sqrt(abs(arean(3,is)))
            if (dist.gt.deps*dscale) then
               fac=deps*dscale/dist
               dxn(1,3,is)=dxn(1,3,is)*fac
               dxn(2,3,is)=dxn(2,3,is)*fac
               dxn(3,3,is)=dxn(3,3,is)*fac
            endif
         else !an edge stack
            dscale3=min(sqrt(arean(3,is)),xn(3,3,is))
            dist3=sqrt(dxn(1,3,is)**2+dxn(2,3,is)**2+dxn(3,3,is)**2)
            dscale2=min(sqrt(arean(2,is)),xn(3,2,is))
            dist3=sqrt(dxn(1,2,is)**2+dxn(2,2,is)**2+dxn(3,2,is)**2)
            dmax=max(dmax,dist2,dist3)
            dmean=dist3+dist2+dmean
            fac=min(deps*dscale3/dist3,deps*dscale2/dist2)
            if (fac.lt.1d0) then
               dxn(1,3,is)=dxn(1,3,is)*fac
               dxn(2,3,is)=dxn(2,3,is)*fac
               dxn(3,3,is)=dxn(3,3,is)*fac
               dxn(1,2,is)=dxn(1,2,is)*fac
               dxn(2,2,is)=dxn(2,2,is)*fac
               dxn(3,2,is)=dxn(3,2,is)*fac
            endif
         endif
      enddo
      dmean=dmean/(ns+nl)
!
      if (idebug.ge.1) print *,'rezback: Wtot,dmax,dmean',
     1                                   Wtot,dmax,dmean
!
!--displace dorsal nodes
      do is=1,ns
         xn(1,3,is)=xn(1,3,is)+dxn(1,3,is)
         xn(2,3,is)=xn(2,3,is)+dxn(2,3,is)
         xn(3,3,is)=xn(3,3,is)+dxn(3,3,is)
      enddo
!--displace middle edge nodes
      do il=1,nl
         is=isol(1,il)
         xn(1,2,is)=xn(1,2,is)+dxn(1,2,is)
         xn(2,2,is)=xn(2,2,is)+dxn(2,2,is)
         xn(3,2,is)=xn(3,2,is)+dxn(3,2,is)
      enddo
!
      return
      end subroutine rezback
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezone(idebug,xn,c_newton,wtot)
!
! This subroutine rezones the middle interior nodes
! so as to minimize the TTM functional
!
      implicit none
!
      integer idebug
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton
!
      real(8) deps !maximum node motion parameter
      parameter(deps=1d-1)
      real(8) dxn(3,3,NSM)
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23 !metric coeff.
      real(8) dH1,dH2,dH3
      real(8) F,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcb, Jcbx,Jcby,Jcbz, dJ2, dJ3
      real(8) Win(NSM),Wtot,Wg
      real(8) Wx(3,NSM),Wy(3,NSM),Wz(3,NSM) !grad of Winslow functional
      real(8) Wxx(3,NSM),Wyy(3,NSM),Wzz(3,NSM) !Hessian of Winslow
      real(8) Wxy(3,NSM),Wxz(3,NSM),Wyz(3,NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
      real(8) HessW(3,3), gradW(3), dvec(3),gradWn
      real(8) HessWtg11,HessWtg12,HessWtg22
      real(8) gradWtg1, gradWtg2, dtg1,dtg2,ddet
      real(8) u1x,u1y,u1z,u2x,u2y,u2z
      real(8) dmax,dmean
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dscale
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
!--Initialize the W's
      Wtot=0d0
      Win=0d0 
      Wx=0d0; Wy=0d0; Wz=0d0
      Wxx=0d0; Wyy=0d0; Wzz=0d0
      Wxy=0d0; Wxz=0d0; Wyz=0d0
!--zero out displacements
      dxn=0d0
!
      do iq=1,nq !outer loop over iq quadrilateral
         do isg=1,4 !loop over the gauss points of iq
            do lvg=1,3
               dx1=dxV(1,lvg,isg,iq)
               dx2=dxV(2,lvg,isg,iq)
               dx3=dxV(3,lvg,isg,iq)
               dy1=dyV(1,lvg,isg,iq)
               dy2=dyV(2,lvg,isg,iq)
               dy3=dyV(3,lvg,isg,iq)
               dz1=dzV(1,lvg,isg,iq)
               dz2=dzV(2,lvg,isg,iq)
               dz3=dzV(3,lvg,isg,iq)
               call TTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                  g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,Wg)
               Win(iq)=Win(iq)+wgp(lvg)*Wg
               Wtot=Wtot+wgp(lvg)*Wg
               do isn=1,4
                  istack=isoq(isn,iq)
                  dH1=dH(1,2,isn,lvg,isg)
                  dH2=dH(2,2,isn,lvg,isg)
                  dH3=dH(3,2,isn,lvg,isg)
                  call dTTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                      g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,
     2                      dH1,dH2,dH3,
     3                      Wxi,Wyi,Wzi,
     4                      Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
                  Wx(2,istack)=Wx(2,istack)+Wxi*wgp(lvg)
                  Wy(2,istack)=Wy(2,istack)+Wyi*wgp(lvg)
                  Wz(2,istack)=Wz(2,istack)+Wzi*wgp(lvg)
                  Wxx(2,istack)=Wxx(2,istack)+Wxxi*wgp(lvg)
                  Wyy(2,istack)=Wyy(2,istack)+Wyyi*wgp(lvg)
                  Wzz(2,istack)=Wzz(2,istack)+Wzzi*wgp(lvg)
                  Wxy(2,istack)=Wxy(2,istack)+Wxyi*wgp(lvg)
                  Wxz(2,istack)=Wxz(2,istack)+Wxzi*wgp(lvg)
                  Wyz(2,istack)=Wyz(2,istack)+Wyzi*wgp(lvg)
               enddo
            enddo
         enddo
      enddo
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with under-relaxation)
!
      do is=1,ns
         if (ilos(1,is).eq.0) then!not an edge stack
!--middle nodes
            gradW(1)=-Wx(2,is)
            gradW(2)=-Wy(2,is)
            gradW(3)=-Wz(2,is)
            HessW(1,1)=Wxx(2,is)
            HessW(2,1)=Wxy(2,is)
            HessW(3,1)=Wxz(2,is)
            HessW(1,2)=Wxy(2,is)
            HessW(2,2)=Wyy(2,is)
            HessW(3,2)=Wyz(2,is)
            HessW(1,3)=Wxz(2,is)
            HessW(2,3)=Wyz(2,is)
            HessW(3,3)=Wzz(2,is)
            call mat3solv(HessW,gradW,dvec)
            dxn(1,2,is)=c_newton*dvec(1)
            dxn(2,2,is)=c_newton*dvec(2)
            dxn(3,2,is)=c_newton*dvec(3)
         endif 
      enddo
!
! for debugging might want a dump
      if (idebug.ge.2) then
         if (mod(kount,50).eq.0) then
            do is=1,ns
               hvec(1:3,:,is)=xn(1:3,:,is)
               hvec(4:6,:,is)=vnn(1:3,:,is)
               hvec(7:9,:,is)=dxn(1:3,:,is)
            enddo
            open(66,file='rezone.dump')
            call iowrfile(99,66)
            close(66)
            hvec(7:9,:,1:ns)=0d0
         endif
         kount=kount+1
      endif
!
!--update the node positions while limiting the displacements
!
      dmean=0d0
      dmax=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            dist=sqrt(dxn(1,2,is)**2+dxn(2,2,is)**2+dxn(3,2,is)**2)
            dmax=max(dmax,dist)
            dmean=dist+dmean
            dscale=abs(voln(2,is))**(1./3.)
            if (dist.gt.deps*dscale) then
               fac=deps*dscale/dist
               dxn(1,2,is)=dxn(1,2,is)*fac
               dxn(2,2,is)=dxn(2,2,is)*fac
               dxn(3,2,is)=dxn(3,2,is)*fac
            endif
         endif
      enddo
      dmean=dmean/ns
      if (idebug.ge.1) print *,'rezone:  Wtot,dmax,dmean',
     1                                   Wtot,dmax,dmean
!
!--displace the middle nodes
      do is=1,ns
         xn(1,2,is)=xn(1,2,is)+dxn(1,2,is)
         xn(2,2,is)=xn(2,2,is)+dxn(2,2,is)
         xn(3,2,is)=xn(3,2,is)+dxn(3,2,is)
      enddo
!
      return
      end subroutine rezone
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        UTILITIES FOR SUBROUTINES CALLED BY REZDRIVER
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine findzeta(z1,z2,z3,z,zeta)
!
! Given z1 < z2 < z3 node positions (zeta=-1,0,1), find the intrinsic 
! coord zeta that gives position z
!
      implicit none
!
      real(8) z1, z2, z3, z, zeta
!
      real(8) tol
      parameter(tol=1d-6)
      real(8) zetad, zetau, zi
      integer i
!
      if (z.lt.z1.or.z.gt.z3.or.z2.lt.z1.or.z2.gt.z3) then 
         print *,'ERROR in findzeta',z1,z2,z3,z
         stop
      endif
!
      if (z.lt.z2) then
         zetad=-1d0
         zetau=0d0
      elseif (z.gt.z2) then
         zetad=0d0
         zetau=1d0
      else !z=z2 case
         zeta=0d0
         return
      endif
!
      do i=1,25
         zeta=0.5d0*(zetad+zetau)
         zi=0.5d0*zeta*(zeta-1d0)*z1+(1d0-zeta**2)*z2+
     1      0.5d0*zeta*(zeta+1d0)*z3
         if (abs(z-zi).lt.tol) then !close enough
            exit 
         elseif (zi.lt.z) then !zi too small, increase zetad
            zetad=zeta
         elseif (zi.gt.z) then !zi too big, decrease zetau
            zetau=zeta
         endif
      enddo
      if (i.ge.25) then
         print *,'failure to converge in findzeta' 
         print *,'z1,z2,z3,z,zi,zeta',z1,z2,z3,z,zi
         print *,'zeta,abs(z-zi)',zeta,abs(z-zi)
      endif
!
      return
      end subroutine findzeta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Winslow(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F)
!
! This subroutine computes the ingredients for the Winslow 
! functional on a surface (W=F/Jacobian)
!
      implicit none
!
!--derivatives of x,y,z wrt to xi and eta (or zeta) at point on surface
      real(8) dx1,dx2,dy1,dy2,dz1,dz2 
!--metric components and Winslow functional numerator
      real(8) g11,g22,F
!
      g11=dx1*dx1+dy1*dy1+dz1*dz1
      g22=dx2*dx2+dy2*dy2+dz2*dz2
      F=g11+g22
!
      return
      end subroutine Winslow
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dWin(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F,Jcb,dJ,dS1,dS2,
     3                Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
!
! This subroutine computes the xyz derivatives of the Winslow
! functional (W=F/Jcb) on a surface
!
      implicit none
!
!--derivatives of x,y,z wrt to xi and eta (or zeta) at point on surface
      real(8) dx1,dx2,dy1,dy2,dz1,dz2
!--metric components and Winslow function (F) at point
      real(8) g11,g22,F
!--surface Jacobian and its inverse
      real(8) Jcb,dJ
!--derivatives of the surface shape functions wrt to intrinsic coord
      real(8) dS1,dS2
!
!--derivatives of the Winslow functional
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
!
!
!--intermediate derivatives useful for the calculation
      real(8) g11x,g11y,g11z,g22x,g22y,g22z
      real(8) g11xx,g11yy,g11zz,g22xx,g22yy,g22zz
      real(8) Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Nx,Ny,Nz,Nx_y,Nx_z,Ny_x,Ny_z,Nz_x,Nz_y
      real(8) dJ2,dJ3,Jcbx,Jcby,Jcbz,Jcbxx,Jcbyy,Jcbzz,Jcbxy,Jcbxz,Jcbyz
!
!
!               1st derivatives of metric
!--diagonal terms
      g11x=2d0*dS1*dx1; g11y=2d0*dS1*dy1; g11z=2d0*dS1*dz1
      g22x=2d0*dS2*dx2; g22y=2d0*dS2*dy2; g22z=2d0*dS2*dz2
!
!              2nd derivatives of metric
!--diagnonal terms
      g11xx=2d0*dS1*dS1; g11yy=g11xx; g11zz=g11xx
      g22xx=2d0*dS2*dS2; g22yy=g22xx; g22zz=g22xx
!--note: xy cross terms are zero
!
!             gradient of F
!
      Fx=g11x+g22x
      Fy=g11y+g22y
      Fz=g11z+g22z
!
!             Hessian of F
!
      Fxx=g11xx+g22xx
      Fyy=g11yy+g22yy
      Fzz=g11zz+g22zz
!--note: xy cross terms are zero
!
!           gradient of Jacobian
!
      Nx=dy1*dz2-dy2*dz1
      Ny=dz1*dx2-dz2*dx1
      Nz=dx1*dy2-dx2*dy1
!
      Nx_y=dS1*dz2-dS2*dz1
      Nx_z=dy1*dS2-dy2*dS1
      Ny_x=dz1*dS2-dz2*dS1
      Ny_z=dS1*dx2-dS2*dx1
      Nz_x=dS1*dy2-dS2*dy1
      Nz_y=dx1*dS2-dx2*dS1
!
! Nx_yy=Nx_yz=0 etc.
!
      Jcbx=dJ*(Ny*Ny_x+Nz*Nz_x)
      Jcby=dJ*(Nx*Nx_y+Nz*Nz_y)
      Jcbz=dJ*(Nx*Nx_z+Ny*Ny_z)
!
!     Hessian of Jacobian
!
      Jcbxx=dJ*(-Jcbx*Jcbx+Ny_x*Ny_x+Nz_x*Nz_x)
      Jcbyy=dJ*(-Jcby*Jcby+Nx_y*Nx_y+Nz_y*Nz_y)
      Jcbzz=dJ*(-Jcbz*Jcbz+Nx_z*Nx_z+Ny_z*Ny_z)
      Jcbxy=dJ*(-Jcbx*Jcby+Nz_x*Nz_y)
      Jcbxz=dJ*(-Jcbx*Jcbz+Ny_x*Ny_z)
      Jcbyz=dJ*(-Jcby*Jcbz+Nx_y*Nx_z)
!
      dJ2=dJ*dJ
      Wxi=(Fx*Jcb-F*Jcbx)*dJ2
      Wyi=(Fy*Jcb-F*Jcby)*dJ2
      Wzi=(Fz*Jcb-F*Jcbz)*dJ2
      dJ3=dJ*dJ2
      Wxxi=(Fxx*Jcb*Jcb-2d0*Jcbx*(Fx*Jcb-F*Jcbx)-F*Jcb*Jcbxx)*dJ3
      Wyyi=(Fyy*Jcb*Jcb-2d0*Jcby*(Fy*Jcb-F*Jcby)-F*Jcb*Jcbyy)*dJ3
      Wzzi=(Fzz*Jcb*Jcb-2d0*Jcbz*(Fz*Jcb-F*Jcbz)-F*Jcb*Jcbzz)*dJ3
      Wxyi=(-(Fx*Jcby+Fy*Jcbx)*Jcb+2d0*F*Jcbx*Jcby-F*Jcbxy*Jcb)*dJ3
      Wxzi=(-(Fx*Jcbz+Fz*Jcbx)*Jcb+2d0*F*Jcbx*Jcbz-F*Jcbxz*Jcb)*dJ3
      Wyzi=(-(Fy*Jcbz+Fz*Jcby)*Jcb+2d0*F*Jcby*Jcbz-F*Jcbyz*Jcb)*dJ3
!
      return
      end subroutine dWin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine TTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1               g11,g22,g33,g12,g13,g23,
     2               F,Jcb,dJ2,dJ3,Wg)
!
! this subroutine computes the TTM functional (3D Winslow)
!
      implicit none
!
!--derivatives of xyz wrt to xi,eta,zeta 
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
!--metric coefficients
      real(8) g11,g22,g33,g12,g13,g23
!--components of the TTT functional
      real(8) F,Jcb,dJ2,dJ3,Wg
!
      real(8) a(3,3)
!
      g11=dx1*dx1+dy1*dy1+dz1*dz1
      g22=dx2*dx2+dy2*dy2+dz2*dz2
      g33=dx3*dx3+dy3*dy3+dz3*dz3
      g12=dx1*dx2+dy1*dy2+dz1*dz2
      g13=dx1*dx3+dy1*dy3+dz1*dz3
      g23=dx2*dx3+dy2*dy3+dz2*dz3
      F=g11*g22+g11*g33+g22*g33-g12*g12-g13*g13-g23*g23
      a(1,1)=dx1; a(1,2)=dx2; a(1,3)=dx3
      a(2,1)=dy1; a(2,2)=dy2; a(2,3)=dy3
      a(3,1)=dz1; a(3,2)=dz2; a(3,3)=dz3
      Jcb=det3x3(a)
      dJ2=1d0/(Jcb*Jcb)
      dJ3=1d0/(Jcb*Jcb*Jcb)
      Wg=F/Jcb
!
      return
      end subroutine TTM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dTTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,
     2                dH1,dH2,dH3,
     3                Wxi,Wyi,Wzi,
     4                Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
!
!--this subroutine computes the derivatives of the TTM functional
!
      implicit none
!
!--derivatives of xyz wrt to xi,eta,zeta 
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
!--metric and Jacobian coefficients
      real(8) g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3
!--derivatives of the shape function wrt to intrinsic coord.
      real(8) dH1,dH2,dH3
!--derivatives of TTM functional
      real(8) Wxi,Wyi,Wzi
      real(8) Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
!
!--derivatives of the numerator of the TTM functional
      real(8) Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
!--derivatives of the denominator (Jacobian) of the TTM functional
      real(8) Jcbx,Jcby,Jcbz
!--auxillary derivatives of the metric
      real(8) g11x,g11y,g11z,g22x,g22y,g22z,g33x,g33y,g33z
      real(8) g12x,g12y,g12z,g13x,g13y,g13z,g23x,g23y,g23z
      real(8) g11xx,g11yy,g11zz,g22xx,g22yy,g22zz,g33xx,g33yy,g33zz
      real(8) g12xx,g12yy,g12zz,g13xx,g13yy,g13zz,g23xx,g23yy,g23zz
      real(8) a(3,3)
!
!
!               1st derivatives of metric
!--diagonal terms
      g11x=2d0*dH1*dx1; g11y=2d0*dH1*dy1; g11z=2d0*dH1*dz1
      g22x=2d0*dH2*dx2; g22y=2d0*dH2*dy2; g22z=2d0*dH2*dz2
      g33x=2d0*dH3*dx3; g33y=2d0*dH3*dy3; g33z=2d0*dH3*dz3
!
!--off diagonal terms
      g12x=dH1*dx2+dH2*dx1; g12y=dH1*dy2+dH2*dy1; g12z=dH1*dz2+dH2*dz1
      g13x=dH1*dx3+dH3*dx1; g13y=dH1*dy3+dH3*dy1; g13z=dH1*dz3+dH3*dz1
      g23x=dH2*dx3+dH3*dx2; g23y=dH2*dy3+dH3*dy2; g23z=dH2*dz3+dH3*dz2
!
!              2nd derivatives of metric
!--diagnonal terms
      g11xx=2d0*dH1*dH1; g11yy=g11xx; g11zz=g11xx
      g22xx=2d0*dH2*dH2; g22yy=g22xx; g22zz=g22xx
      g33xx=2d0*dH3*dH3; g33yy=g33xx; g33zz=g33xx
!
!--off diagonal terms
      g12xx=2d0*dH1*dH2; g12yy=g12xx; g12zz=g12xx
      g13xx=2d0*dH1*dH3; g13yy=g13xx; g13zz=g13xx
      g23xx=2d0*dH2*dH3; g23yy=g23xx; g23zz=g23xx
!--note: cross terms are zero
!
!             gradient of F
!
      Fx=g11x*g22+g11*g22x+g11x*g33+g11*g33x+g22x*g33+g22*g33x
     1                 -2d0*g12x*g12-2d0*g13x*g13-2d0*g23x*g23
      Fy=g11y*g22+g11*g22y+g11y*g33+g11*g33y+g22y*g33+g22*g33y
     1                 -2d0*g12y*g12-2d0*g13y*g13-2d0*g23y*g23
      Fz=g11z*g22+g11*g22z+g11z*g33+g11*g33z+g22z*g33+g22*g33z
     1                 -2d0*g12z*g12-2d0*g13z*g13-2d0*g23z*g23
!
!             Hessian of F
!--diagonal terms
      Fxx=g11xx*g22+2d0*g11x*g22x+g11*g22xx+
     1    g11xx*g33+2d0*g11x*g33x+g11*g33xx+   
     2    g22xx*g33+2d0*g22x*g33x+g22*g33xx
     3    -2d0*g12x*g12x-2d0*g12xx*g12
     4    -2d0*g13x*g13x-2d0*g13xx*g13
     5    -2d0*g23x*g23x-2d0*g23xx*g23
      Fyy=g11yy*g22+2d0*g11y*g22y+g11*g22yy+
     1    g11yy*g33+2d0*g11y*g33y+g11*g33yy+   
     2    g22yy*g33+2d0*g22y*g33y+g22*g33yy
     3    -2d0*g12y*g12y-2d0*g12yy*g12
     4    -2d0*g13y*g13y-2d0*g13yy*g13
     5    -2d0*g23y*g23y-2d0*g23yy*g23
      Fzz=g11zz*g22+2d0*g11z*g22z+g11*g22zz+
     1    g11zz*g33+2d0*g11z*g33z+g11*g33zz+   
     2    g22zz*g33+2d0*g22z*g33z+g22*g33zz
     3    -2d0*g12z*g12z-2d0*g12zz*g12
     4    -2d0*g13z*g13z-2d0*g13zz*g13
     5    -2d0*g23z*g23z-2d0*g23zz*g23
!
!--off diagonal terms
      Fxy=g11x*g22y+g11y*g22x+g11x*g33y+g11y*g33x+g22x*g33y+g22y*g33x
     1    -2d0*g12x*g12y-2d0*g13x*g13y-2d0*g23x*g23y 
      Fxz=g11x*g22z+g11z*g22x+g11x*g33z+g11z*g33x+g22x*g33z+g22z*g33x
     1    -2d0*g12x*g12z-2d0*g13x*g13z-2d0*g23x*g23z 
      Fyz=g11y*g22z+g11z*g22y+g11y*g33z+g11z*g33y+g22y*g33z+g22z*g33y
     1    -2d0*g12y*g12z-2d0*g13y*g13z-2d0*g23y*g23z 
!
!           gradient of Jacobian
!
!--Jcbx
      a(1,1)=dH1; a(1,2)=dH2; a(1,3)=dH3
      a(2,1)=dy1; a(2,2)=dy2; a(2,3)=dy3
      a(3,1)=dz1; a(3,2)=dz2; a(3,3)=dz3
      Jcbx=det3x3(a)
!--Jcby
      a(1,1)=dx1; a(1,2)=dx2; a(1,3)=dx3
      a(2,1)=dH1; a(2,2)=dH2; a(2,3)=dH3
      a(3,1)=dz1; a(3,2)=dz2; a(3,3)=dz3
      Jcby=det3x3(a)
!--Jcbz
      a(1,1)=dx1; a(1,2)=dx2; a(1,3)=dx3
      a(2,1)=dy1; a(2,2)=dy2; a(2,3)=dy3
      a(3,1)=dH1; a(3,2)=dH2; a(3,3)=dH3
      Jcbz=det3x3(a)
!
!     Hessian of Jacobian=0
!
!  Derivatives of TTM functional
!
      Wxi=(Fx*Jcb-F*Jcbx)*dJ2
      Wyi=(Fy*Jcb-F*Jcby)*dJ2
      Wzi=(Fz*Jcb-F*Jcbz)*dJ2
      Wxxi=(Fxx*Jcb*Jcb-2d0*Jcbx*(Fx*Jcb-F*Jcbx))*dJ3
      Wyyi=(Fyy*Jcb*Jcb-2d0*Jcby*(Fy*Jcb-F*Jcby))*dJ3
      Wzzi=(Fzz*Jcb*Jcb-2d0*Jcbz*(Fz*Jcb-F*Jcbz))*dJ3
      Wxyi=(Fxy*Jcb*Jcb-(Fx*Jcby+Fy*Jcbx)*Jcb+F*2d0*Jcbx*Jcby)*dJ3
      Wxzi=(Fxz*Jcb*Jcb-(Fx*Jcbz+Fz*Jcbx)*Jcb+F*2d0*Jcbx*Jcbz)*dJ3
      Wyzi=(Fyz*Jcb*Jcb-(Fy*Jcbz+Fz*Jcby)*Jcb+F*2d0*Jcby*Jcbz)*dJ3
!
      return
      end subroutine dTTM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!           END OF REZDRIVER RELATED SUBROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avback(xav,xrz,xback)
!
! For every rezoned back node (all dorsal nodes and all middle edge 
! nodes), find the nearest position on the surface of the advected
! mesh and put it in xback
!
      implicit none
!
      real(8) xav(3,3,NSM)  ! coordinates of advected mesh
      real(8) xrz(3,3,NSM)  ! coordinates of rezoned mesh
      real(8) xback(3,3,NSM)! interpolated coordinates of the rezoned
                            ! mesh onto the advected mesh
!
      real(8) xyz(3)
      real(8) xi,eta,zeta
      integer list(NSM)
      integer is,lv,kq,isn
!
      real(8) vnodes(3,3,4)  !v at node (component,level,stack)
      real(8) c(3)  !input intrinsic coordinates
      real(8) v(3)  !interpolated vector
      real(8) dvdc(3,3)  ! vector derivatives wrt to intrinsic coord
!
      xback=xrz !for ventral and interior middle nodes
                !will be the same as the rezoned mesh
!
      do is=1,ns
!--dorsal nodes (level=3)
         xyz=xrz(:,3,is)
         call find_surface(xyz,xav,kq,xi,eta,zeta)
         if (kq.eq.0) then
            print *,'avback: no interpolating element for stack',
     1              is,' level 3'
            stop
         endif
         do isn=1,4
            vnodes(:,:,isn)=xav(:,:,isoq(isn,kq))
         enddo
         c(1)=xi; c(2)=eta; c(3)=zeta
         call getv12(vnodes,c,v,dvdc)
         xback(:,3,is)=v(:)
         if (ilos(1,is).ne.0) then !edge stack, do middle node
            xyz=xrz(:,2,is)
            call find_surface(xyz,xav,kq,xi,eta,zeta)
            if (kq.eq.0) then
               print *,'avback: no interpolating element for stack',
     1                 is,' level 2'
               stop
            endif
            do isn=1,4
               vnodes(:,:,isn)=xav(:,:,isoq(isn,kq))
            enddo
            c(1)=xi; c(2)=eta; c(3)=zeta
            call getv12(vnodes,c,v,dvdc)
            xback(:,2,is)=v(:)
         endif
      enddo
!
!     do is=1,ns
!        hvec(1:3,:,is)=xback(1:3,:,is)
!     enddo
!     open(66,file='avback.dump')
!     call iowrfile(99,66)
!     close(66)
!     print *,'dumped to avback.dump'
!     stop
!
      return
      end subroutine avback
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avrz(xav,xrz,kqrz,xirz,etarz,zetarz)
!
! This subroutine finds the element and intrinsic coordinates
! with which to interpolate the node values of the final
! mesh with respect to the advected mesh 
!
      implicit none
!
      real(8) xav(3,3,NSM) ! coordinates of advected mesh
      real(8) xrz(3,3,NSM) ! coordinates of rezoned mesh
      integer kqrz(3,NSM)  ! interpolating element for rezoned nodes
      real(8) xirz(3,NSM)  ! xi value for rezoned nodes
      real(8) etarz(3,NSM) ! eta value for rezoned nodes
      real(8) zetarz(3,NSM)! zeta value for rezoned nodes
!
      real(8) xyz(3)
      real(8) xi,eta,zeta
      integer list(NSM)
      integer is,lv,kq
!
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
!--middle node (level=2)
!--prepare the list of neighboring elements
            list(1:kqos(is))=lqos(iqos(is):iqos(is)+kqos(is)-1)
            xyz=xrz(:,2,is)
            call find_element(xyz,xav,kqos(is),list,kq,xi,eta,zeta)
            if (kq.eq.0) then ! we look over the entire mesh
               call find_element(xyz,xav,nq,idxarr,kq,xi,eta,zeta)
               if (kq.eq.0) then
                  print *,'no interpolating element found for stack',
     1                    is,' level 2'
                  stop
               endif
            endif
            kqrz(2,is)=kq
            xirz(2,is)=xi
            etarz(2,is)=eta
            zetarz(2,is)=zeta
!
!--dorsal and ventral nodes (level=1 and 3)
            do lv=1,3,2
               xyz=xrz(:,lv,is)
               call find_surface(xyz,xav,kq,xi,eta,zeta)
               if (kq.eq.0) then
                  print *,'no interpolating element found for stack',
     1                    is,' level',lv
                  stop
               endif
               kqrz(lv,is)=kq
               xirz(lv,is)=xi
               etarz(lv,is)=eta
               zetarz(lv,is)=zeta
            enddo
         else !--edge stack case
!--contact line first (ventral node)
            xyz=xrz(:,1,is)
            call find_cl(xyz,xav,kq,xi,eta,zeta)
            kqrz(1,is)=kq
            xirz(1,is)=xi
            etarz(1,is)=eta
            zetarz(1,is)=zeta
!--middle and dorsal node
            do lv=2,3
               xyz=xrz(:,lv,is)
               call find_surface(xyz,xav,kq,xi,eta,zeta)
               kqrz(lv,is)=kq
               xirz(lv,is)=xi
               etarz(lv,is)=eta
               zetarz(lv,is)=zeta
            enddo
         endif
      enddo
!
      return
      end subroutine avrz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avfield(kom,opt)
!
!--this subroutine advects an svec field
!
      implicit none
!
      integer kom !svec field to advect
      character*(*) opt !choices are nw or ca
!
      integer,save::kqnw(3,NSM)  
      real(8),save::xinw(3,NSM), etanw(3,NSM), zetanw(3,NSM)
      integer,save::kqca(3,NSM)  
      real(8),save::xica(3,NSM), etaca(3,NSM), zetaca(3,NSM)
!
      real(8) svtmp(3,NSM), svnew(3,NSM)
      integer nskind, is
!
      nskind=nint(dscp(kom,idspk))
      if (index(opt,'nw').gt.0) then !network advection
!--check if already computed the interpolating coordinates, if not, do it
         if (.not.avnw) then
            call avrz(xnw,hvec(1:3,:,:),kqnw,xinw,etanw,zetanw)
            avnw=.true.
         endif
!--apply dilation coefficient
         if (nskind.eq.1) then !this is a volume variable
            svtmp=svec(kom,:,:)*dilfacv
         elseif (nskind.eq.2) then !this is a surface variable
            svtmp=svec(kom,:,:)*dilfaca
         else 
            print *,'unrecognized variable type'
            print *,'kom=',kom,'   nskind=',nskind
            stop
         endif
!--interpolate
         call avinter(nskind,svtmp,kqnw,xinw,etanw,zetanw,svnew)
      elseif (index(opt,'ca').gt.0) then !bulk advection
         if (.not.avca) then
            call avrz(xca,hvec(1:3,:,:),kqca,xica,etaca,zetaca)
            avca=.true.
         endif
         if (nskind.eq.1.or.nskind.eq.2) then
            svtmp=svec(kom,:,:)
         else 
            print *,'unrecognized variable type'
            print *,'kom=',kom,'   nskind=',nskind
            stop
         endif
!--interpolate
         call avinter(nskind,svtmp,kqca,xica,etaca,zetaca,svnew)
      endif
      svec(kom,:,:)=svnew
!
      return
      end subroutine avfield
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avinter(nskind,svtmp,kqs,xis,etas,zetas,svnew)
!
! This subroutine interpolate an svec field 
!
      implicit none
!
      integer nskind
      real(8) svtmp(3,NSM)
      integer kqs(3,NSM)
      real(8) xis(3,NSM), etas(3,NSM), zetas(3,NSM)
      real(8) svnew(3,NSM)
!
      real(8) xi,eta,zeta
      real(8) dsdxi,dsdeta,dsdzeta
      real(8) snodes(3,4), s
      integer is,lv,i,l,kq,isn,il
!
      if (nskind.eq.1) then !volume variable
         do is=1,ns
            do lv=1,3
               kq=kqs(lv,is)
               xi=xis(lv,is)
               eta=etas(lv,is)
               zeta=zetas(lv,is)
               do i=1,4
                  isn=isoq(i,kq)
                  do l=1,3
                     snodes(l,i)=svtmp(l,isn)
                  enddo
               enddo
               call gets(snodes,xi,eta,zeta,s,dsdxi,dsdeta,dsdzeta)
               svnew(lv,is)=s
            enddo
         enddo
      elseif (nskind.eq.2) then !surface variable
         do is=1,ns  !loop over dorsal and ventral nodes
            do lv=1,3,2
               kq=kqs(lv,is)
               xi=xis(lv,is)
               eta=etas(lv,is)
               zeta=zetas(lv,is)
               do i=1,4
                  isn=isoq(i,kq)
                  do l=1,3
                     snodes(l,i)=svtmp(l,isn)
                  enddo
               enddo
               call gets(snodes,xi,eta,zeta,s,dsdxi,dsdeta,dsdzeta)
               svnew(lv,is)=s
            enddo
         enddo
         do il=1,nl  !loop over mid-point edge nodes
            is=isol(1,il)
            kq=kqs(2,is)
            xi=xis(2,is)
            eta=etas(2,is)
            zeta=zetas(2,is)
            do i=1,4
               isn=isoq(i,kq)
               do l=1,3
                  snodes(l,i)=svtmp(l,isn)
               enddo
            enddo
            call gets(snodes,xi,eta,zeta,s,dsdxi,dsdeta,dsdzeta)
            svnew(2,is)=s
         enddo
      else
         print *,'unable to advect nskind=',nskind
         stop
      endif
!
      return
      end subroutine avinter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avtstep(idebug,avtstp,vmin)
!
! this subroutine computes a recommended advection time step:
! for node x,y,z and neighboring node xn,yn,xn at distance dn
! and joined by unit vector un
!  dt=0.05*dn/abs(v.un)
!
      implicit none
!
      integer idebug !a debug flag
      real(8) avtstp !the recommended time step
      real(8) vmin   !a minimum velocity to set the time step
!
      real(8) ctstp
      parameter(ctstp=0.05d0)
      real(8) ntstp(3,NSM)
      integer tstploc(2)
      integer isp(4), ism(4)
      data isp/2,3,4,1/, ism/4,1,2,3/
      integer iq,is,istack,istackn,istackp
      real(8) x,y,z,vx,vy,vz,xn,yn,zn,dn,unx,uny,unz
      real(8) z2,vz2,z3,vz3
!
      ntstp=vbig
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            istackn=isoq(ism(is),iq)
            istackp=isoq(isp(is),iq)
!
!--ventral node
            x=hvec(1,1,istack)   
            y=hvec(2,1,istack)   
            z=hvec(3,1,istack)   
            vx=hvec(7,1,istack)   
            vy=hvec(8,1,istack)   
            vz=hvec(9,1,istack)   
!--the previous ventral node in the element
            xn=hvec(1,1,istackn)   
            yn=hvec(2,1,istackn)   
            zn=hvec(3,1,istackn)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(1,istack)=min(ntstp(1,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!--the next ventral node in the element
            xn=hvec(1,1,istackp)   
            yn=hvec(2,1,istackp)   
            zn=hvec(3,1,istackp)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(yn-z)/dn
            ntstp(1,istack)=min(ntstp(1,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!--the node above (level 2) in the element
            xn=hvec(1,2,istack)   
            yn=hvec(2,2,istack)   
            zn=hvec(3,2,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(1,istack)=min(ntstp(1,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!
!--middle node
            x=hvec(1,2,istack)   
            y=hvec(2,2,istack)   
            z=hvec(3,2,istack)   
            vx=hvec(7,2,istack)   
            vy=hvec(8,2,istack)   
            vz=hvec(9,2,istack)   
!--the previous middle node in the element
            xn=hvec(1,2,istackn)   
            yn=hvec(2,2,istackn)   
            zn=hvec(3,2,istackn)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(2,istack)=min(ntstp(2,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!--the next middle node in the element
            xn=hvec(1,2,istackp)   
            yn=hvec(2,2,istackp)   
            zn=hvec(3,2,istackp)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(2,istack)=min(ntstp(2,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!--the node below (level 1) in the element
            xn=hvec(1,1,istack)   
            yn=hvec(2,1,istack)   
            zn=hvec(3,1,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(2,istack)=min(ntstp(2,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!--the node above (level 3) in the element
            xn=hvec(1,3,istack)   
            yn=hvec(2,3,istack)   
            zn=hvec(3,3,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(2,istack)=min(ntstp(2,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!
!--dorsal node
            x=hvec(1,3,istack)   
            y=hvec(2,3,istack)   
            z=hvec(3,3,istack)   
            vx=hvec(7,3,istack)   
            vy=hvec(8,3,istack)   
            vz=hvec(9,3,istack)   
!--the previous dorsal node in the element
            xn=hvec(1,3,istackn)   
            yn=hvec(2,3,istackn)   
            zn=hvec(3,3,istackn)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(3,istack)=min(ntstp(3,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!--the next dorsal node in the element
            xn=hvec(1,3,istackp)   
            yn=hvec(2,3,istackp)   
            zn=hvec(3,3,istackp)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(3,istack)=min(ntstp(3,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
!--the node below (level 2) in the element
            xn=hvec(1,2,istack)   
            yn=hvec(2,2,istack)   
            zn=hvec(3,2,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unz=(zn-z)/dn
            ntstp(3,istack)=min(ntstp(3,istack),
     1                          dn/(vmin+abs(vx*unx+vy*uny+vz*unz)))
         enddo
      enddo
!
!--make sure that the middle and top nodes will not descend toward
!  the substratum too fast.
      do istack=1,ns
         z2=hvec(3,2,istack)
         vz2=hvec(9,2,istack)
         if (z2.gt.0d0.and.vz2.lt.0d0) then
            ntstp(2,istack)=min(ntstp(2,istack),-z2/vz2)
         endif
         z3=hvec(3,3,istack)
         vz3=hvec(9,3,istack)
         if (z3.gt.0d0.and.vz3.lt.0d0) then
            ntstp(3,istack)=min(ntstp(3,istack),-z3/vz3)
         endif
      enddo
      avtstp=ctstp*minval(ntstp)
      if (idebug.ge.1) then
         tstploc=minloc(ntstp)
         print *,'avtstp: minimum time step',avtstp
         print *,'avtstp: set by level, stack',tstploc
      endif
!
      return
      end subroutine avtstep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avmidreset(idebug,xn)
!
! this subroutine moves all the middle edge nodes to z2=z3/2
!
      implicit none
!
      integer idebug
      real(8) xn(3,3,NSM)
!
      real(8)  z2, z3, z, zeta, H1, H2, H3, dLdzeta, d2Ldzeta2
      integer il,is
!
      do il=1,nl
         is=isol(1,il)
         z2=xn(3,2,is)
         z3=xn(3,3,is)
         z=z3*0.5d0
         call findzeta(0d0,z2,z3,z,zeta)
!--get the line shape functions at zeta2
         call getL3(-1,zeta,H1,dLdzeta,d2Ldzeta2)
         call getL3(0,zeta,H2,dLdzeta,d2Ldzeta2)
         call getL3(+1,zeta,H3,dLdzeta,d2Ldzeta2)
         xn(1,2,is)=H1*xn(1,1,is)+H2*xn(1,2,is)+H3*xn(1,3,is)
         xn(2,2,is)=H1*xn(2,1,is)+H2*xn(2,2,is)+H3*xn(2,3,is)
         xn(3,2,is)=H2*z2+H3*z3 !should be z3/2
         if (idebug.eq.1) then
            print *,'is,zeta',is,zeta
            print *,'x1,y1,z1',xn(1,1,is),xn(2,1,is),xn(3,1,is)
            print *,'x2,y2,z2',xn(1,2,is),xn(2,2,is),xn(3,2,is)
            print *,'x3,y3,z3',xn(1,3,is),xn(2,3,is),xn(3,3,is)
         endif
      enddo
!
      return
      end subroutine avmidreset
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      RANDOM UTILITIES FOR AVLIB
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function det3x3(a)
! function returns determinant of a 3x3 matrix a
      implicit none
      real(8) det3x3
      real(8) a(3,3)
      det3x3=+a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     1       -a(2,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3))
     2       +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
      end function det3x3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mat3solv(Bmat,rhs,Bvec)
!
! this subroutine solves a 3x3 linear system with
! coefficient Bmat and right hand side rhs and places
! the result in Bvec
!
      implicit none
!
      real(8) Bmat(3,3),rhs(3),Bvec(3)
!
      real(8) cofac(3,3), Inv(3,3)
      real(8) ddet
      integer i, j
!
            
!--compute cofactor matrix
      cofac(1,1)=+Bmat(2,2)*Bmat(3,3)-Bmat(3,2)*Bmat(2,3)
      cofac(2,1)=-Bmat(1,2)*Bmat(3,3)+Bmat(3,2)*Bmat(1,3)
      cofac(3,1)=+Bmat(1,2)*Bmat(2,3)-Bmat(2,2)*Bmat(1,3)
      cofac(1,2)=-Bmat(2,1)*Bmat(3,3)+Bmat(3,1)*Bmat(2,3)
      cofac(2,2)=+Bmat(1,1)*Bmat(3,3)-Bmat(3,1)*Bmat(1,3)
      cofac(3,2)=-Bmat(1,1)*Bmat(2,3)+Bmat(2,1)*Bmat(1,3)
      cofac(1,3)=+Bmat(2,1)*Bmat(3,2)-Bmat(3,1)*Bmat(2,2)
      cofac(2,3)=-Bmat(1,1)*Bmat(3,2)+Bmat(3,1)*Bmat(1,2)
      cofac(3,3)=+Bmat(1,1)*Bmat(2,2)-Bmat(2,1)*Bmat(1,2)
!--determinant
      ddet=1d0/(Bmat(1,1)*cofac(1,1)+Bmat(2,1)*cofac(2,1)+
     1                               Bmat(3,1)*cofac(3,1))
      if (abs(ddet).gt.vbig) then
         print *,'mat3solv: determinant nearly 0!',ddet
      endif
!
      do j=1,3
         do i=1,3
            Inv(i,j)=ddet*cofac(j,i)
         enddo
      enddo
!
      Bvec=0d0
      do j=1,3
         do i=1,3
            Bvec(j)=Bvec(j)+Inv(i,j)*rhs(i)
         enddo
      enddo
! 
      return
      end subroutine mat3solv         
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      END MODULE AVLIBSW!end library of advection routines.
