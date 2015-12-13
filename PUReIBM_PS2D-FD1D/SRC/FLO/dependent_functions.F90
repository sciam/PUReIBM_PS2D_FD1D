!    PUReIBM-PS2D-FD1D is a three-dimensional psudeo-spectral particle-resolved
!    direct numerical simulation solver for detailed analysis of homogeneous
!    fixed and freely evolving fluid-particle suspensions. PUReIBM-PS2D-FD1D
!    is a continuum Navier-Stokes and scalar solvers based on Cartesian grid that utilizes
!    Immeresed Boundary method to represent particle surfuces. The details about the solvers
!    can be found in the below papers in SUBRAMANIAM's group. 
!    Copyright (C) 2015, Shankar Subramaniam, Rahul Garg, Sudheer Tenneti, Bo Sun, Mohammad Mehrabadi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    For acknowledgement, please refer to the following publications:
!     For hydrodynamic solver :
!     (1) TENNETI, S. and SUBRAMANIAM, S., 2014, Particle-resolved direct numerical
!         simulation for gas-solid flow model development. Annu. Rev. Fluid Mech.
!         46 (1) 199-230.
!     (2) M. Mehrabadi, S. Tenneti, R. Garg, and S. Subramaniam, 2015, Pseudo-turbulent 
!         gas-phase velocity fluctuations in homogeneous gas-solid flow: fixed particle
!         assemblies and freely evolving suspensions. J. Fluid Mech. 770 210-246.
!
!     For scalar solver :
!     (3) S. Tenneti, B. Sun, R. Garg, S. Subramaniam, 2013, Role of fluid heating in dense
!         gas-solid flow as revealed by particle-resolved direct numerical simulation.
!         International Journal of Heat and Mass Transfer 58 471-479.

Module dependent_functions
#include "ibm.h"
  USE precision
  USE constants 
  Use global_data
  Use fftw_interface
  USE nlarrays
  implicit none 
Contains
  SUBROUTINE set_interpolation_stencil(pc, ib,ie,jb,je,kb,ke, isch&
       &, ordernew)
    Use global_data, Only :im=>mxf, jm=>my,km=>mz,intx_per, inty_per,intz_per
    IMPLICIT NONE 
    INTEGER, DIMENSION(3), INTENT(in):: pc
    INTEGER, INTENT(out):: ib, ie, jb, je, kb, ke 
    INTEGER, OPTIONAL :: ordernew
    CHARACTER*5, INTENT(in) :: isch 
    INTEGER :: ob2rtmp, ob2ltmp, ordertemp, ordernewtmp
    SELECT CASE(isch)
    CASE('csi')

       ob2rtmp = ob2r
       ob2ltmp = ob2l
       ordertemp = order
!!$       IF(order.NE.3.AND.(pc(1).EQ.1.OR.pc(1).EQ.im&
!!$            &-1.OR.pc(2).EQ.1.OR.pc(2).EQ.jm&
!!$            &-1.OR.pc(3).EQ.1.OR.pc(3).EQ.km-1)) order = 3
       !To make the order at boundary cells for csi equal to 3. 
       ob2l = (order+1)/2
       ob2r = order/2 
       !print*,'ob2l = ',ob2l,ob2r
       
       if(.not.intx_per) then 
          ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
          ie = MIN(im,pc(1) + ob2r)
          IF (ib.EQ.1 ) ie = ib + order - 1
          IF (ie.EQ.im) ib = ie - order + 1
       else 
          if(pc(1).lt.0) then
             ib = pc(1) - ob2r
             ie = pc(1) + (ob2l - 1) !periodic
          else
             ib = pc(1) - (ob2l - 1) !periodic
             ie = pc(1) + ob2r
          end if
       end if

       if(.not.inty_per) then
          
          jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
          je = MIN(jm,pc(2) + ob2r)
          IF (jb.EQ.1 ) je = jb + order - 1
          IF (je.EQ.jm) jb = je - order + 1
       else
          if(pc(2).lt.0) then
             jb = pc(2) - ob2r
             je = pc(2) + (ob2l - 1) !periodic
          else
             jb = pc(2) - (ob2l - 1) !periodic
             je = pc(2) + ob2r
          end if
       end if
       
       
       If(.not.intz_per) then 
          kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
          ke = MIN(km,pc(3) + ob2r)
          IF (kb.EQ.1 ) ke = kb + order - 1
          IF (ke.EQ.km) kb = ke - order + 1
       else
          if(pc(3).lt.0) then
             kb = pc(3) - ob2r
             ke = pc(3) + (ob2l - 1) !periodic
          else
             kb = pc(3) - (ob2l - 1) !periodic
             ke = pc(3) + ob2r
          end if
       end If

       ob2r =  ob2rtmp 
       ob2l = ob2ltmp
       ordernewtmp = order

       !Print*,'ib,ie .... in processing =', ib,ie,jb,je,kb,ke,&
       !     & ordernewtmp
       !PRINT*, 'pc = ',pc(1),pc(2),pc(3)!
       order = ordertemp !reset the order

    CASE('lpi')
       !print*, 'order in set stencil = ', order
       ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
       ie = MIN(im,pc(1) + ob2r)
       if(.not.intx_per) then 
          IF (ib.EQ.1 ) ie = ib + order - 1
          IF (ie.EQ.im) ib = ie - order + 1
       else 
          IF (ib.EQ.1 ) ib = ie - order + 1
          IF (ie.EQ.im) ie = ib + order - 1
       end if

       jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
       je = MIN(jm,pc(2) + ob2r)
       if(.not.inty_per) then
          IF (jb.EQ.1 ) je = jb + order - 1
          IF (je.EQ.jm) jb = je - order + 1
       else
          IF (jb.EQ.1 ) jb = je - order + 1
          IF (je.EQ.jm) je = jb + order - 1
       end if

       kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
       ke = MIN(km,pc(3) + ob2r)
       If(.not.intz_per) then 
          IF (kb.EQ.1 ) ke = kb + order - 1
          IF (ke.EQ.km) kb = ke - order + 1
       else
          IF (kb.EQ.1 ) kb = ke - order + 1
          IF (ke.EQ.km) ke = kb + order - 1
       end If
       ordernewtmp = order
    END SELECT
    IF(PRESENT(ordernew)) ordernew = ordernewtmp
  END SUBROUTINE set_interpolation_stencil


  SUBROUTINE set_interpolation_pstencil(pc, ib,ie,jb,je,kb,ke, isch&
       &, ordernew)
    Use global_data, Only :jm=>my,km=>mz,intx_per, inty_per,intz_per
    IMPLICIT NONE 
    INTEGER, DIMENSION(3), INTENT(in):: pc
    INTEGER, INTENT(out):: ib, ie, jb, je, kb, ke 
    INTEGER, OPTIONAL :: ordernew
    CHARACTER*5, INTENT(in) :: isch 
    INTEGER :: ob2rtmp, ob2ltmp, ordertemp, ordernewtmp, im
    im = mxf - 1
    SELECT CASE(isch)
    CASE('csi')

       ob2rtmp = ob2r
       ob2ltmp = ob2l
       ordertemp = order
!!$       IF(order.NE.3.AND.(pc(1).EQ.1.OR.pc(1).EQ.im&
!!$            &-1.OR.pc(2).EQ.1.OR.pc(2).EQ.jm&
!!$            &-1.OR.pc(3).EQ.1.OR.pc(3).EQ.km-1)) order = 3
       !To make the order at boundary cells for csi equal to 3. 
       ob2l = (order+1)/2
       ob2r = order/2 
       !print*,'ob2l = ',ob2l,ob2r
       
       if(.not.intx_per) then 
          ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
          ie = MIN(im,pc(1) + ob2r)
          IF (ib.EQ.1 ) ie = ib + order - 1
          IF (ie.EQ.im) ib = ie - order + 1
       else 
          if(pc(1).lt.0) then
             ib = pc(1) - ob2r
             ie = pc(1) + (ob2l - 1) !periodic
          else
             ib = pc(1) - (ob2l - 1) !periodic
             ie = pc(1) + ob2r
          end if
       end if

       if(.not.inty_per) then
          
          jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
          je = MIN(jm,pc(2) + ob2r)
          IF (jb.EQ.1 ) je = jb + order - 1
          IF (je.EQ.jm) jb = je - order + 1
       else
          if(pc(2).lt.0) then
             jb = pc(2) - ob2r
             je = pc(2) + (ob2l - 1) !periodic
          else
             jb = pc(2) - (ob2l - 1) !periodic
             je = pc(2) + ob2r
          end if
       end if
       
       
       If(.not.intz_per) then 
          kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
          ke = MIN(km,pc(3) + ob2r)
          IF (kb.EQ.1 ) ke = kb + order - 1
          IF (ke.EQ.km) kb = ke - order + 1
       else
          if(pc(3).lt.0) then
             kb = pc(3) - ob2r
             ke = pc(3) + (ob2l - 1) !periodic
          else
             kb = pc(3) - (ob2l - 1) !periodic
             ke = pc(3) + ob2r
          end if
       end If

       ob2r =  ob2rtmp 
       ob2l = ob2ltmp
       ordernewtmp = order

       !Print*,'ib,ie .... in processing =', ib,ie,jb,je,kb,ke,&
       !     & ordernewtmp
       !PRINT*, 'pc = ',pc(1),pc(2),pc(3)!
       order = ordertemp !reset the order

    CASE('lpi')
       !print*, 'order in set stencil = ', order
       ib = MAX(1 ,pc(1) - (ob2l - 1)) !non-periodic
       ie = MIN(im,pc(1) + ob2r)
       if(.not.intx_per) then 
          IF (ib.EQ.1 ) ie = ib + order - 1
          IF (ie.EQ.im) ib = ie - order + 1
       else 
          IF (ib.EQ.1 ) ib = ie - order + 1
          IF (ie.EQ.im) ie = ib + order - 1
       end if
       
       jb = MAX(1 ,pc(2) - (ob2l - 1)) !non-periodic
       je = MIN(jm,pc(2) + ob2r)
       if(.not.inty_per) then
          IF (jb.EQ.1 ) je = jb + order - 1
          IF (je.EQ.jm) jb = je - order + 1
       else
          IF (jb.EQ.1 ) jb = je - order + 1
          IF (je.EQ.jm) je = jb + order - 1
       end if
       
       kb = MAX(1 ,pc(3) - (ob2l - 1)) !non-periodic
       ke = MIN(km,pc(3) + ob2r)
       If(.not.intz_per) then 
          IF (kb.EQ.1 ) ke = kb + order - 1
          IF (ke.EQ.km) kb = ke - order + 1
       else
          IF (kb.EQ.1 ) kb = ke - order + 1
          IF (ke.EQ.km) ke = kb + order - 1
       end If
       ordernewtmp = order
    END SELECT
    IF(PRESENT(ordernew)) ordernew = ordernewtmp
  END SUBROUTINE set_interpolation_pstencil


  subroutine calc_realdiv(divur)
    Use nlarrays, Only : uf1
    implicit none 

    integer :: i,j,k

    real(prcn), DImension(:,:,:), Intent(out) ::  divur
    do i=1,mx1,1
       do j=1,my2,1
          do k=1,mz,1
             uf1(j,k)=(u(i+1,j,k,1)-u(i,j,k,1))/dx
             uf1(j,k)=uf1(j,k)+half*wy(j)*(u(i,j,k,2)+u(i+1,j,k,2))
             uf1(j,k)=uf1(j,k)+half*wz(k)*(u(i,j,k,3)+u(i+1,j,k,3))
          enddo
       enddo

       call ff2cr(uf1,divur(i,:,:))
    end do
  end subroutine calc_realdiv

  subroutine calc_pressure
    Use nlarrays , Only : pf =>uf1
    Use nlmainarrays, Only : pr => pbcp
!    Use bcsetarrays, Only : ppr
    implicit none 

    integer :: i, j, k, ii, jj, kk, xfac, yfac, zfac
    
!!$    pr(1,1,1) = zero
!!$    
!!$    do k = 1, mz
!!$       do j = 1, my
!!$          do i = 1, mx1
!!$             ii = i-1
!!$             jj = j-1
!!$             kk = k-1
!!$             xfac = 1
!!$             yfac = 1
!!$             zfac = 1
!!$             if(ii.lt.1) then
!!$                ii = 1
!!$                xfac = 0
!!$             end if
!!$             if(jj.lt.1) then
!!$                jj = 1
!!$                yfac = 0
!!$             end if
!!$             if(kk.lt.1) then
!!$                kk = 1
!!$                zfac = 0
!!$             end if
!!$             pr(i,j,k) = pr(ii,jj,kk) + ppr(ii,jj,kk,1)*dx*xfac + ppr(ii,jj,kk,2)*dy*yfac +  ppr(ii,jj,kk,3)*dz*zfac
!!$          end do
!!$       end do
!!$    end do
#if PARALLEL
    do i = 0, nx+1
#else
    do i = 1, nx !mxf
#endif       
!       IF(xperiodic) THEN 
       
       pf(:,:)=p(i,:,:)
       
       call ff2cr(pf(:,:),pr(i,:,:))
!!$       ELSE
!!$          ii = foffset+i
!!$          pf(:,:)=p(ii-1,:,:)
!!$          
!!$          call ff2cr(pf,pr(i,:,:))
!       end IF
       
    end do
  end subroutine calc_pressure

  subroutine calc_forreal(fr)
    implicit none 

    real(prcn), DImension(:,:,:,:), Intent(out) :: fr
    integer :: i,n
    do n=1,ndim
       do i=1,mxf

          call ff2cr(ff(i,:,:,n),fr(foffset+i,:,:,n))
          
       enddo			! loop over i
    end do

  end subroutine calc_forreal

  subroutine calc_velreal_mxf(velr)
    implicit none 

    real(prcn), DImension(:,:,:,:), Intent(out) :: velr
    integer :: i,j,k, n ,ii

    do n=1,ndim

       do i=1,mxf,1
          ii = i+foffset
          call ff2cr(u(ii,:,:,n),velr(i,:,:,n))
          velr(i,:,:,n) = velr(i,:,:,n) + umean(n)
       enddo
    end do

  end subroutine calc_velreal_mxf

  subroutine calc_velreal(velr)
    implicit none 

    real(prcn), DImension(:,:,:,:), Intent(out) :: velr
    integer :: i,j,k, n, idim


#if PARALLEL
	do idim = 1, ndim
		i = 1
		VECSENDRECV(u(i,1,1,idim),1,ucslice,fromproc,1,u(nx+1,1,1,idim),1,toproc,1,decomp_group,status)
		i = nx
		VECSENDRECV(u(i,1,1,idim),1,ucslice,toproc,0,u(0,1,1,idim),1,fromproc,0,decomp_group,status)
	end do
#endif


    do n=1,ndim
#if PARALLEL
       do i=0,nx+1
#else
       do i=1,nx !mx,1       
#endif
          call ff2cr(u(i,:,:,n),velr(i,:,:,n))
          velr(i,:,:,n) = velr(i,:,:,n) + umean(n)
       enddo
    end do

  end subroutine calc_velreal

  subroutine calc_realvort
    Use nlarrays , Only : uf1
    implicit none 

    integer :: i,j,k
  end subroutine calc_realvort
  
  
  subroutine writepres_sph!(pr,ur) 
    USE general_funcs
    USE interpolation

    implicit none 
    INTEGER :: unitno, m, l
    INTEGER, SAVE :: count_routine = 0 
    REAL(prcn) :: stagangle 
    CHARACTER*80 :: FILENAME1
    LOGICAL :: filexist, isopen 
    count_routine  = count_routine + 1
	
    IF(count_routine.eq.3) count_routine = 1

 3020 FORMAT(I1)
    WRITE(FILENAME1, 3020) count_routine	
    FILENAME1 = TRIM(RUN_NAME)//'_cpvstheta_'//TRIM(FILENAME1)//'.dat'
   INQUIRE(FILE=FILENAME1,EXIST=filexist,OPENED=isopen)
       
    unitno = getnewunit(minunitno, maxunitno)
    IF (.NOT.filexist) THEN
       
    OPEN(unitno,FILE=FILENAME1, status='new')
      ELSEIF(filexist.AND..NOT.isopen) THEN
    open(unit=unitno,file=FILENAME1,status='replace')
    ENDIF

    WRITE(unitno, '("# generated at time = ",g12.5)') t
    do m=1,nbody              ! loop over bodies
       write(unitno,*) 'zone'
       do l=1,nbnd
          if(xs(3,l).eq.zero)then
             if(xs(2,l).ge.zero)then
                if(xs(1,l).lt.zero)then
                   stagangle = abs(atan(xs(2,l)/xs(1,l)))
                elseif(xs(1,l).gt.zero)then
                   stagangle = atan(xs(2,l)/xs(1,l))
                   stagangle = pi-stagangle
                elseif(xs(1,l).eq.zero)then
                   stagangle = pi/2.
                endif
                
!                write(unitno,101)stagangle*180./pi,coeffp(m,l)
                !     End loop over all bnd points
             endif
          endif
       enddo
101    FORMAT(10(2x,e20.12))
       
       !     end loop over all bodies
    enddo
       close(unitno, status = "keep")
  end subroutine writepres_sph
  
    
  subroutine deterpressph!(pr,ur) 
    USE general_funcs
    USE interpolation

    Use nlmainarrays, Only : ur=>ubcp, nlr=>nlbcp, onlr=>onlbcp,pr=> pbcp
    Use bcsetarrays, Only : ppr
    !Use boundary_condition
    implicit none 

    !real(prcn), Intent(in) , Dimension(:,:,:) ::  pr
    !real(prcn), Intent(in) , Dimension(:,:,:,:) ::  ur


    real(prcn) ::  xl(ndim), xpb(ndim)
    real(prcn) ::  ul(ndim),unorm(ndim),utang(ndim), ppll(ndim)
    real(prcn) ::  snorm(ndim), stang(ndim)
    real(prcn) ::  unmag,tx,ty,tz, dpdthetamag

    integer :: m,l,n,i,j,k,d,pcell(3), ii,jj,kk

    real(prcn) ::  rad,pl, dfll(ndim), nll(ndim), onll(ndim)
    real(prcn) ::  stagangle
    integer :: isp
    integer :: is(ndim)
    Integer :: ib,ie, jb,je, kb,ke , onew
    open(unit=45,file='presphperip.dat',status='unknown')


   PRINT*,'IN DETER PRESS: NBODY = ', nbody, nbnd 
    write(100,*)'Zone'

    do m=1,nbody              ! loop over bodies
       write(45,*) 'zone'
       do l=1,nbnd

          if(xs(3,l).eq.zero)then
             if(xs(2,l).ge.zero)then
                rad=0.0
                do n=1,ndim

                   xl(n)=xc(m,n)+xs(n,l)*radbdy(m)
                   is(n)=int(xl(n))

                   ul(n)=zero
                   rad=rad+(xs(n,l)*radbdy(m))**2.0
                enddo

                


                rad=dsqrt(rad)
                
                if(xs(1,l).lt.zero)then
                   stagangle = abs(atan(xs(2,l)/xs(1,l)))
                elseif(xs(1,l).gt.zero)then
                   stagangle = atan(xs(2,l)/xs(1,l))
                   stagangle = pi-stagangle
                elseif(xs(1,l).eq.zero)then
                   stagangle = pi/2.
                endif
                do n=1,ndim
                   snorm(n)=(xs(n,l)*radbdy(m))/rad
                enddo
                stang(1) = snorm(2)
                stang(2) = -snorm(1)
                stang(3) = zero

                write(100,'(4(F15.10,1x))')xs(1,l),xs(2,l),xs(3,l)&
                     &,stagangle

                pl = zero
                ppll = zero
                isp=INT(xl(1)-0.5)
                pcell(1) = isp 
                pcell(2:3) = is(2:3)
                xpb(1) = xl(1)-0.5
                xpb(2:3)=xl(2:3)
!                CALL interpolate_pdata(pcell,xpb, ppr,ppll,pl)
!!$                CALL set_interpolation_pstencil(pcell,ib,ie,jb,je,kb,ke&
!!$                     &,interp_scheme, onew) 
!!$                do k = 1, onew
!!$                   do j = 1, onew
!!$                      do i = 1, onew
!!$                         ii = ib+i-1
!!$                         jj = jb+j-1
!!$                         kk = kb+k-1
!!$                         gstencil(i,j,k,1) = ib+(i-1)
!!$                         gstencil(i,j,k,2) = jb+(j-1)
!!$                         gstencil(i,j,k,3) = kb+(k-1)
!!$	     
!!$                         if(ii.lt.1.and.intx_per) THEN
!!$                            !PRINT*,'ii LT 1'
!!$                            ii = mxf+ii-1
!!$                         end if
!!$             
!!$                         if(ii.gt.mxf-1.and.intx_per) THEN 
!!$                            ii = ii-mxf +1
!!$                            !PRINT*,'ii GT mXF-1'
!!$                            !READ(*,*)
!!$                         end if
!!$                         
!!$                         
!!$                         if(jj.lt.1) jj = my+jj
!!$                         if(jj.gt.my) jj = jj-my
!!$                         if(kk.lt.1) kk = mz+kk
!!$                         if(kk.gt.mz) kk = kk-mz
!!$                         prsten(i,j,k) = pr(ii,jj,kk)
!!$                      end do
!!$                   end do
!!$                end do
!!$                CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:ndim)&
!!$                     &,prsten(1:onew,1:onew,1:onew)&
!!$                     &,xpb(1:ndim),pl,onew, interp_scheme,weightp)
                dpdthetamag = zero
                do n = 1, ndim
                   dpdthetamag = dpdthetamag + stang(n)*ppll(ndim)
                end do
                
!!$                write(45,101)stagangle*180./pi,pl/0.5/(upi(1)**2.)
                write(45,101)stagangle*180./pi,dpdthetamag/0.5/(uchar(1)**2.)
             endif
          endif
       enddo

       !     end loop over all bodies
    enddo

101 FORMAT(10(2x,e20.12))
    close(45)


  end subroutine deterpressph


  subroutine interpolate_udata(pc,pos,ib,ie,jb,je,kb,ke,ul&
       &,nll,onll,dfll,flag, ibody,l,onew)

    USE general_funcs
    USE interpolation

    USE bcsetarrays, ONLY :  diffn 
    Use nlmainarrays, Only : ur=>ubcp, nlr=>nlbcp, onlr=>onlbcp!,pr=> pbc 
    implicit none 
    Integer, intent(in) :: pc(3)
    Integer, Intent(in) :: flag, ibody,l
    Real(prcn), Dimension(:), Intent(in) :: pos
    Integer, INtent(out) :: ib, ie, jb,je,kb,ke, onew
    Real(prcn), Intent(out), Dimension(:) :: ul,nll,onll,dfll !, ppll
    Integer :: i, j,k, ii,jj,kk, n
!!$    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in) ::  ppr
    

    CALL set_interpolation_stencil(pc,ib,ie,jb,je,kb,ke&
         &,interp_scheme, onew) 
    if((ib.lt.1.or.ie.gt.mxf).and..not.xperiodic) Print*,'Error in i ....',ib,ie,pc,pos
    do k = 1, onew
       do j = 1, onew
          do i = 1, onew
             ii = ib+i-1
             jj = jb+j-1
             kk = kb+k-1
             gstencil(i,j,k,1) = ib+(i-1)
             gstencil(i,j,k,2) = jb+(j-1)   
             gstencil(i,j,k,3) = kb+(k-1)
#if !PARALLEL
	     if(ii.lt.1) ii = mxf+ii-1
             if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             if(kk.lt.1) kk = mz+kk
             if(kk.gt.mz) kk = kk-mz
!!$          
             LOCAL_INDEX(ii)
#if PARALLEL
             if(flag.eq.1)then
                if(ii.lt.0) then 
                   PRINT*,'udata:ERROR WEST SIDE', myid,frombnd,fromrpr
                   PRINT*,'WEST SIDE: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
                   PRINT*,'WEST SIDE: INDEX: ', l
                   PARALLEL_FINISH()
                   STOP
                else if(ii.gt.nx+1)then
                   PRINT*,'udata:ERROR EAST SIDE', myid,frombnd,fromrpr
                   PRINT*,'EAST SIDE: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
                   PRINT*,'EAST SIDE: INDEX: ', l
                   PARALLEL_FINISH()
                   STOP
               end if
            else
               if(ii.lt.-1) then 
                  PRINT*,'udata:ERROR WEST EXT POINT', myid,frombnd,fromrpr
                  PRINT*,'WEST EXT POINT: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
                  PRINT*,'WEST EXT POINT: INDEX: ', l
                  PARALLEL_FINISH()
                  STOP
                else if(ii.gt.nx+2)then
                   PRINT*,'udata:ERROR EAST EXT POINT', myid,frombnd,fromrpr
                   PRINT*,'EAST EXT POINT: BODY& POSITIONS: ', myid,ibody, pc(1), ii, pos(1) 
                   PRINT*,'EAST EXT POINT: INDEX: ', l
                   PARALLEL_FINISH()
                   STOP
                end if
             end if
#endif
             vsten(i,j,k,1:ndim) = ur(ii,jj,kk,1:ndim)
             if(flag.eq.1) then 
                nlsten(i,j,k,1:ndim) = nlr(ii,jj,kk,1:ndim)
                onlsten(i,j,k,1:ndim) = onlr(ii,jj,kk,1:ndim)
                dfsten(i,j,k,1:ndim) = diffn(ii,jj,kk,1:ndim)
             end if
             
          end do
       end do
    end do
    CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:3),&
         & vsten(1:onew,1:onew,1:onew,1:ndim),pos(1:ndim),ul(1:ndim),onew,&
         & interp_scheme,weightp) 
    if(flag.eq.1) then 
       do n = 1, ndim 
          nll(n) =  array_dot_product(nlsten(1:onew,1:onew,1:onew,n&
               &),weightp(1:onew,1:onew,1:onew)) 
          onll(n)=  array_dot_product(onlsten(1:onew,1:onew,1:onew,n&
               &),weightp(1:onew,1:onew,1:onew))
          dfll(n)=  array_dot_product(dfsten(1:onew,1:onew,1:onew,n&
               &),weightp(1:onew,1:onew,1:onew))
!!$          ppll(n)=  array_dot_product(ppgrsten(1:onew,1:onew,1:onew,n&
!!$               &),weightp(1:onew,1:onew,1:onew))
       end do
    end if
  end subroutine interpolate_udata

  subroutine interpolate_pdata(pc,pos,ppl,pl,l)
    USE general_funcs
    USE interpolation
    Use nlmainarrays, Only : pr=> pbcp 
    Use bcsetarrays, Only : ppr
    implicit none 
    Integer, intent(in) :: pc(3),l
!    REAL(prcn), DIMENSION(:,:,:,:), INTENT(in) ::  ppr
    !REAL(prcn), DIMENSION(:,:,:), INTENT(in) ::  pr
    Real(prcn), Dimension(:), Intent(in) :: pos
    Real(prcn), Intent(out), Dimension(:) :: ppl
    Real(prcn), Intent(out) :: pl
    Integer :: i, j,k, onew, ii,jj,kk, n,ib,ie,jb,je,kb,ke
    

    CALL set_interpolation_pstencil(pc,ib,ie,jb,je,kb,ke&
         &,interp_scheme, onew) 
    do k = 1, onew
       do j = 1, onew
          do i = 1, onew
             ii = ib+i-1
             jj = jb+j-1
             kk = kb+k-1
             gstencil(i,j,k,1) = ib+(i-1)
             gstencil(i,j,k,2) = jb+(j-1)
             gstencil(i,j,k,3) = kb+(k-1)
#if !PARALLEL
	     if(ii.lt.1.and.intx_per) THEN
         !PRINT*,'ii LT 1'
                ii = mxf+ii-1
             end if

             if(ii.gt.mxf-1.and.intx_per) THEN 
                ii = ii-mxf +1
                !PRINT*,'ii GT mXF-1'
                !READ(*,*)
             end if
#endif 

             !The above definition is different compared to that for velocity grid. This is because staggered pressure grid
             !in x- direction, RG 01/19/08
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             if(kk.lt.1) kk = mz+kk
             if(kk.gt.mz) kk = kk-mz
!!$   
             LOCAL_INDEX(ii)
#if PARALLEL
             if(ii.lt.0)then
                PRINT*,'pdata:ERROR WEST SIDE', myid,pos(1),frombnd&
                     &,fromrpr,xstart,xend,l
                PARALLEL_FINISH()
                STOP
             else if(ii.gt.nx+1) then
                PRINT*,'pdata:ERROR EAST SIDE', myid,pos(1),frombnd&
                     &,fromrpr,xstart,xend,l
                PARALLEL_FINISH()
                STOP
             end if
#endif

             ppgrsten(i,j,k,1:ndim) = ppr(ii,jj,kk,1:ndim)
             prsten(i,j,k) = pr(ii,jj,kk)!-DREAL(p(ii,1,1))
          end do
       end do
    end do
    CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:ndim)&
         &,ppgrsten(1:onew,1:onew,1:onew,1:ndim)&
         &,pos(1:ndim),ppl(1:ndim),onew, interp_scheme,weightp)
    pl = array_dot_product(prsten(1:onew,1:onew&
         &,1:onew),weightp(1:onew,1:onew,1:onew)) 
  end subroutine interpolate_pdata

    Subroutine grid_nodes_insphere
    Implicit None 

    Integer :: i,j,k, idim, cor_min(ndim), cor_max(ndim), m, ii, jj,&
         & kk, imin, imax, jmin, jmax, kmin, kmax, &
         & count_solid_tmp, count_fluid_tmp, iphs
    Real(prcn) :: xlr(ndim), xll(ndim), dist, volfracg(nphases)

    
    if(I_AM_NODE_ZERO)Write(*,'(A)')'IN GRID NODES IN SPHERE'
    do iphs = 1, nphases
       phase_array(iphs)%volfracg = zero
       volfracg(iphs) = zero
    end do

    do k = 1, mz
       do j = 1, my
          do i = 0, nx+1 
             fluid_atijk(i,j,k) = .true.
          end do
       end do
    end do
    count_solid  = 0
    count_solid_tmp = 0
    count_fluid_tmp  = nx*my*mz

    do m = 1, nbody
       do idim = 1, ndim 
          if(idim.eq.1) then 
             xlr(idim) = xc(m,idim)  + radbdy(m) + foffset
             xll(idim) = xc(m,idim)  - radbdy(m) + foffset
          else 
             xlr(idim) = xc(m,idim)  + radbdy(m)
             xll(idim) = xc(m,idim)  - radbdy(m) 
          end if
       end do
    
       do idim = 1, ndim 
          cor_min(idim)  = ceiling(xll(idim))
          cor_max(idim) = floor(xlr(idim)) 
       end do
    
       imin = cor_min(1)
       imax = cor_max(1)
       jmin = cor_min(2)
       jmax = cor_max(2)
       kmin = cor_min(3)
       kmax = cor_max(3)
       
!!$       PRINT*,'cor_min(1) = ', cor_min(1), cor_max(1)
!!$       PRINT*,'IMIN, IMAX = ', imin, imax
!!$       PRINT*,'JMIN, JMAX = ', jmin, jmax
!!$       PRINT*,'KMIN, KMAX = ', kmin, kmax
       
       do i = imin, imax 
          ii = i
#if PARALLEL
!!$          if(.not.I_AM_IN_THIS_PROC(ii))then
          if(.not.POINT_IN_PROC(ii))then

             WEST_PERIODIC_POINT_IMAGE(i,ii)
             EAST_PERIODIC_POINT_IMAGE(i,ii)
             if(.not.POINT_IN_PROC(ii)) goto 555
          end if
          
#else
          if(i.lt.1.and.intx_per) ii = mxf+i-1
          if(i.gt.mxf-1.and.intx_per) ii = i-(mxf-1)
#endif
          
          LOCAL_INDEX(ii)

          do j = jmin, jmax 
             jj = j 
             if(j.lt.1.and.inty_per) jj = my+j
             if(j.gt.my.and.inty_per) jj = j-my
             
             do k = kmin, kmax 
                
                kk = k 
                if(k.lt.1.and.intz_per) kk = mz+k
                if(k.gt.mz.and.intz_per) kk = k-mz 
                
                dist = (i - (xc(m,1)+foffset))**two
                dist = dist + (j - xc(m,2))**two + (k - xc(m,3))**two 
                dist = dsqrt(dist)
                IF((dist - radbdy(m)).LE.SMALL_NUMBER) THEN 
                   if(.not.fluid_atijk(ii,jj,kk)) then 
                      PRINT*,'FLUID_ATIJK ALREADY FALSE AT I,J,K=&
                           & ', ii,jj,kk, myid,m
                   ENDIF
                   fluid_atijk(ii,jj,kk)  = .false.
                   !PRINT*,'ijk = ', ii, jj, kk
!!$                   if((1.le.ii).and.(ii.le.nx))then
                   iphs = part_array(m)%iphs
                   volfracg(iphs) = volfracg(iphs) + one
                   count_solid_tmp = count_solid_tmp + 1
                   count_fluid_tmp = count_fluid_tmp - 1
!!$                   end if
                end IF
             end do
          end do
#if PARALLEL
555       continue
          
#endif
       end do
    end do

    VECSENDRECV(fluid_atijk(1,1,1),1,twodlslice,fromproc,1,fluid_atijk(nx+1,1,1),1,toproc,1,decomp_group,status)
    VECSENDRECV(fluid_atijk(nx,1,1),1,twodlslice,toproc,1,fluid_atijk(0,1,1),1,fromproc,1,decomp_group,status)
!!$    PRINT*,'MYID = ', myid, 'count_solid LOCAL = ', count_solid_tmp
!!$    PRINT*,'MYID = ', myid, 'NBODY LOCAL = ', nbody
!!$    PRINT*,'XC(1) = ', myid, XC(1,1), radbdy(1)
    GLOBAL_INT_SUM(count_solid_tmp,count_solid,1,decomp_group)
    GLOBAL_INT_SUM(count_fluid_tmp,count_fluid,1,decomp_group)
    do iphs = 1, nphases
       GLOBAL_DOUBLE_SUM(volfracg(iphs),phase_array(iphs)%volfracg,1,decomp_group)
    end do

    maxvolfrac = REAL(count_solid,prcn)/REAL((mx1*my*mz),prcn)
    do iphs = 1, nphases
       phase_array(iphs)%volfracg = phase_array(iphs)%volfracg&
            &/REAL((mx1*my*mz),prcn)
		if(I_AM_NODE_ZERO) WRite(*,*) 'volfracg = ', (phase_array(iphs)%volfracg)
    end do


#if !PARALLEL
    IF(xperiodic) THEN 
       
       fluid_atijk(mxf,:,:) = fluid_atijk(1,:,:)
       fluid_atijk(0, :, :) = fluid_atijk(mxf-1,:,:)
       !Set the fluid_atijk array correctly for the pseudo inlet and outlet planes correctly for use in scalar non-linear terms calculation
    ELSE	
       !In the case of inflow/outflow B.C.'s, there will be no particles present at the inlet and outlet
       fluid_atijk(0, :, :) = .true.
       fluid_atijk(mxf+1,:,:) = .true.
    ENDIF
#endif

#if 0
#if PARALLEL    
    if(myid.eq.1)then
       OPEN(100,FILE='_fluid.dat')
#else
    if(I_AM_NODE_ZERO)then
       OPEN(100,FILE=TRIM(RUN_NAME)//'_fluid.dat')
#endif
       write(100,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "FX" '
#if PARALLEL
       write(100,*)'ZONE F=POINT, I=', nx+2,  ', J=', my, ', K=', mz
#else
       write(100,*)'ZONE F=POINT, I=', nx/2+2,  ', J=', my, ', K=', mz
#endif

       do k = 1, mz
          do j = 1, my
#if PARALLEL
             do i = 0, nx+1
#else
             do i = nx/2, nx+1
#endif
                if(fluid_atijk(i,j,k))then
                   WRITE(100,*) GLOBAL_INDEX(i), j, k, 1
                else
                   WRITE(100,*) GLOBAL_INDEX(i), j, k, 0
                end if
             end do
          end do
       end do
       close(100,status='keep')
    end if
    if(I_AM_NODE_ZERO) WRITE(*,'(2(2x,A25,i5))')'count_solid=',count_solid, 'count_fluid&
         &=', count_fluid

    PARALLEL_FINISH()
    STOP
#endif

  end Subroutine grid_nodes_insphere


  SUBROUTINE update_nrpr_array(m)
    IMPLICIT NONE
    INTEGER, Intent(in) :: m
    INTEGER :: l,n, iphs
    REAL(prcn) :: xli(ndim), xlo(ndim)
    
    
    iphs = 1!part_array(m)%iphs

    part_array(m)%nrpr_active = phase_array(iphs)%nrpr
    NULLIFY(bndarray)
    bndarray => phase_array(iphs)%bndpts
    nrpr = part_array(m)%nrpr_active

    DO l=1,nrpr
       DO n=1,ndim
          
          xli(n)=xc(m,n)+ bndarray(n,l)*radibdy(m)
          
          xlo(n)=xc(m,n)+ bndarray(n,l)*radobdy(m)
       ENDDO
       CALL check_external_pt(XLO, M, l)
       if(.not.(PART_ARRAY(m)%if_rev(l))) part_array(m)%nrpr_active = part_array(m)%nrpr_active - 1
    END DO
  END SUBROUTINE update_nrpr_array

  subroutine check_external_pt(xloin, LL,irpr)

    USE dem_mod

    implicit none 
    REAL(PRCN), DIMENSION(:), INTENT(IN) :: XLOIN
    INTEGER, INTENT(IN) :: LL, irpr
!!$    LOGICAL, INTENT(OUT) :: REVERSE
    LOGICAL :: PER
    REAL(PRCN), DIMENSION(3) :: XLO
    INTEGER :: IP, PNO, IDIM 
    
    REAL(PRCN) :: DIST
!!$    REVERSE = .TRUE.
    PART_ARRAY(LL)%IF_REV(irpr) = .TRUE.
    PER = .FALSE.
    XLO = XLOIN
    if(XLO(1).lt.one) XLO(1) = mxf+XLO(1)
    if(XLO(1).gt.mxf) XLO(1) = XLO(1)-mxf
    if(XLO(2).lt.one) THEN 
       XLO(2) = my+one+XLO(2)
       !PRINT*, 'XLO2 LT ZERO'
       PER = .TRUE.
    end if

    if(XLO(2).gt.my) THEN 
       XLO(2) = XLO(2)-my-one
       !PRINT*, 'XLO2 GT MY'
       PER = .TRUE.
    end if

    if(XLO(3).lt.one) THEN 
       XLO(3) = mz+one+XLO(3)
       
       PER = .TRUE.
    end if
    
    if(XLO(3).gt.mz) THEN 
       XLO(3) = XLO(3)-mz-one
       
       PER = .TRUE.
    end if
    
    IF(NEIGHBOURS(LL,1).GT.1) THEN 
       DO IP = 2, NEIGHBOURS(LL,1)+1
          PNO = NEIGHBOURS(LL,IP)
          !PRINT*,'PNIO = ', PNO, NEIGHBOURS(LL,NEIGHBOURS(LL,1)+1)
          DIST = ZERO
          DO IDIM = 1, NDIM
             DIST = DIST + (XC(PNO,IDIM)-XLO(IDIM))**2.d0
          end DO
          DIST = SQRT(DIST)
          !PRINT*, 'DIST =', DIST, LL, PNO
          IF(DIST.LT.RADBDY(PNO)) THEN 
             PART_ARRAY(LL)%if_rev(irpr) = .FALSE.
             EXIT
          end IF
       end DO
    end IF
  end subroutine check_external_pt


  SUBROUTINE compute_new_timestep(rks)
    USE nlmainarrays, Only : ubcp, pbcp
    USE scalar_data
    IMPLICIT NONE
    INTEGER, Intent(in) :: rks
    REAL(PRCN) ::  u_max, v_max, w_max, umax_tmp, up_max, uframe_mag
    REAL(PRCN) ::  umax_loc, vmax_loc, wmax_loc  
    REAL(prcn) :: mixmeanslip(ndim),umean_temp(ndim), usmean_temp(ndim), ufmean_temp(ndim),umean_temploc(ndim), usmean_temploc(ndim), ufmean_temploc(ndim), mesh_veltemp(ndim), tempt, utemp
    
    INTEGER :: idim, i, j, k, iphs, partstart, partend, ibody

    DTNEW = LARGE_NUMBER
    u_max = SMALL_NUMBER
    v_max = SMALL_NUMBER
    w_max = SMALL_NUMBER

    umax_loc = SMALL_NUMBER
    vmax_loc = SMALL_NUMBER
    wmax_loc = SMALL_NUMBER

    umax_tmp = LARGE_NUMBER

    UFMEAN_TEMPLOC = zero 
    UMEAN_tempLOC = ZERO
    USMEAN_TeMPLOC = zero 

    mesh_veltemp = mesh_vel
    if(.not.movingcv) mesh_vel = zero
    if(I_AM_NODE_ZERO)WRITE(*,'(A,3(2x,g17.8))')'MESH VEL IN compute_timestep: ', mesh_vel(:) 


#if PARALLEL
    do i = 0, nx+1
#else
       do i = 1,nx
#endif
          DO k = 1,mz
             DO j = 1, my2
                uf1(j,k) = u(i,j,k,1)

                uf2(j,k) = u(i,j,k,2)

                uf3(j,k) = u(i,j,k,3)

             END DO
          END DO
          
          CALL ff2cr(uf1,ubcp(i,:,:,1))
          CALL ff2cr(uf2,ubcp(i,:,:,2))
          CALL ff2cr(uf3,ubcp(i,:,:,3))
          do k = 1, mz
             do j = 1, my 
                ubcp(i,j,k,:) = ubcp(i,j,k,:) + umean(:)
                if((i.gt.0).and.(i.lt.nx+1))then
                   if(FLUID_ATIJK(i,j,k))then
                      do idim = 1, ndim
                         UFMEAN_TEMPLOC(idim)= UFMEAN_TEMPLOC(idim)+ ubcp(i,j,k,idim)
                      end do
                   ELSE
                      USMEAN_TEMPLOC(:)= USMEAN_TEMPLOC(:) + ubcp(i,j,k,:)
                   end if
                   UMEAN_TEMPLOC(:) = UMEAN_TEMPLOC(:) + ubcp(i,j,k,:)

                   umax_loc = MAX(umax_loc, ABS(ubcp(i,j,k,1)-mesh_vel(1)))
                   vmax_loc = MAX(vmax_loc, ABS(ubcp(i,j,k,2)-mesh_vel(2)))
                   wmax_loc = MAX(wmax_loc, ABS(ubcp(i,j,k,3)-mesh_vel(3)))
                end if
             end do
          end do
          
         
          
#if PARALLEL
          if(i.eq.2)then
             do idim=1,ndim
                VECSENDRECV(ubcp(i,1,1,idim),1,urslice,fromproc,1,ubcp(nx+2,1,1,idim),1,toproc,1,decomp_group,status)
!                write (*,*) myid, ': ', 'i: ', i, idim,fromproc ,toproc
             end do
          else if(i.eq.nx-1)then 
             do idim=1,ndim
                VECSENDRECV(ubcp(i,1,1,idim),1,urslice,toproc,0,ubcp(-1,1,1,idim),1,fromproc,0,decomp_group,status)
!                write (*,*) myid, ': ', 'i: ', i, fromproc, toproc, nx
             end do
          end if
#endif
       end do
#if !PARALLEL
       ubcp(mx,:,:,:) = ubcp(1,:,:,:)
#endif
       GLOBAL_DOUBLE_MAX(umax_loc,u_max,1,decomp_group)
       GLOBAL_DOUBLE_MAX(vmax_loc,v_max,1,decomp_group)
       GLOBAL_DOUBLE_MAX(wmax_loc,w_max,1,decomp_group)

       umax_tmp = MAX(u_max,v_max)
       umax_tmp = MAX(umax_tmp,w_max)

       if(umax_tmp*dt/dx.gt.cfl_max) then
          dtnew = cfl_max*dx/umax_tmp
       else
          dtnew = dt 
       endif

		if (move_particles) then
			up_max = SMALL_NUMBER
			do ibody=1, nbody
				up_max = MAX(up_max, sqrt(dot_product(velbdy(ibody,:),velbdy(ibody,:)) ) )
			enddo

			if (up_max*dt/dx > cfl_max) then
				dtnew = cfl_max*dx/umax_tmp
			endif

!			uframe_mag = sqrt(dot_product(frame_vel,frame_vel))
!			if (uframe_mag*dt/dx > cfl_max) then
!				dtnew = cfl_max*dx/uframe_mag
!			endif

		endif

#if PARALLEL       
       GLOBAL_DOUBLE_SUM(UFMEAN_TEMPLOC(1),UFMEAN_TEMP(1),3,decomp_group)
       GLOBAL_DOUBLE_SUM(USMEAN_TEMPLOC(1),USMEAN_TEMP(1),3,decomp_group)
       GLOBAL_DOUBLE_SUM(UMEAN_TEMPLOC(1),UMEAN_TEMP(1),3,decomp_group)
       Write(*,*) myid, UFMEAN_TEMP(1:ndim)
       do idim = 1, ndim
          UFMEAN_TEMP(idim) = UFMEAN_TEMP(idim)/real(count_fluid,prcn)
          USMEAN_TEMP(idim) = USMEAN_TEMP(idim)/real(count_solid,prcn)
          UMEAN_TEMP(idim)= UMEAN_TEMP(idim)/real(mx1*my*mz,prcn)
       end do
       Write(*,"(1i,3d15.7)") myid, UFMEAN_TEMP(1:ndim)

#else
       do idim = 1, ndim
          UFMEAN_TEMP(idim) = UFMEAN_TEMPLOC(idim)/real(count_fluid,prcn)
          USMEAN_TEMP(idim) = USMEAN_TEMPLOC(idim)/real(count_solid,prcn)
          UMEAN_TEMP(idim) = UMEAN_TEMPLOC(idim)/real(mx1*my*mz,prcn)
       end do
#endif
       ufmean_old(:) = ufmean(:)
       ufmean(:) = ufmean_temp(:)

       do idim = 1, ndim 
!!$          usmean_des(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
          if(nbody.gt.0)then
             usmean_act(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
          else
             usmean_act(idim) = zero
          end if
       end do
       
       dufmeandt = -cf*(ufmean(:)-ufmean_old(:))
       
       !dufmeandt = -cf*(ufmean_des(:)-ufmean(:))
       meanslip(:) = (ufmean(:)-usmean_des(:))*(one-maxvolfrac)
       meanslipmod = DSQRT(dot_product(meanslip(1:ndim), meanslip(1:ndim)))
       !usmean can become anything in IBM. SO we dont evaluate mean slip based on the current value of usmean. As long as the particle velocities are correctly attained, usmean_des has been reached and usmean calculated is irrelevant.!
       usmean = zero
       do iphs = 1, nphases
          phase_array(iphs)%mean_spec_vel(1:ndim) = zero
          partstart = phase_array(iphs)%pstart
          partend = phase_array(iphs)%pend
          do idim = 1, ndim
             phase_array(iphs)%mean_spec_vel(idim) = SUM(velbdy(partstart:partend,idim))/real(phase_array(iphs)%npart,prcn)
          end do
!!$          usmean(1:ndim) = usmean(1:ndim) + (phase_array(iphs)&
!!$               &%volfracg)*(phase_array(iphs)%mean_spec_vel(1:ndim))
          usmean(1:ndim) = usmean(1:ndim) + (phase_array(iphs)&
               &%volfrac)*(phase_array(iphs)%mean_spec_vel(1:ndim))
          !WRite(*,*)'volfracg = ', (phase_array(iphs)%volfracg)
       end do

       usmean(:) = usmean(:)/mean_volfrac

       mixmeanslip(:) = (one-maxvolfrac)*(usmean(:)-ufmean(:))
       mixmeanslipmod = DSQRT(DOT_PRODUCT(mixmeanslip(1:ndim),mixmeanslip(1:ndim)))

       if(I_AM_NODE_ZERO)then
          WRITE(*,'(A25,3(2x,g17.8))')'USMEAN DES = ', usmean_des(:)
          WRITE(*,'(A25,3(2x,g17.8))')'USMEAN ACT = ', usmean_act(:)
          WRITE(*,'(A25,3(2x,g17.8))')'USMEAN MIX = ', usmean(:)
          WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN DES = ', ufmean_des(:)
          WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN ACTUAL = ', ufmean_temp(:)
          !WRITE(*,'(A25,3(2x,g17.8))')'UFMEAN CC = ', ufmean_cc(:)
          WRITE(*,'(A25,3(2x,g17.8))')'UMEAN ACTUAL = ', umean(:)
          WRITE(*,'(A25,3(2x,g17.8))')'UMEAN TMP = ', umean_temp(:)
          WRITE(*,'(A25,2(2x,g12.5))') "MAX VEL AND CFL:", umax_tmp,  umax_tmp*dt/dx
          WRITE(*,'(A40,3(2x,g12.5))') "MEANSLIPMOD AND REY(MEANSLIPMOD):", mixmeanslipmod, mixmeanslipmod*char_length/vis
       end if

       !Convert Pressure to real space
       CALL calc_pressure
       
       IF(rks.eq.itrmax)then
          IF(adaptive.and..not.only_dem)then
             DT = DTNEW
          ENDIF
       END IF
       mesh_vel = mesh_veltemp
         
        !! dt = one/300.0d0
        !!  dt = cfl*dx/umax_tmp

       t=t+dt*(coef(rks,1)+coef(rks,2))
       if(iscalon.eq.1)tscal=tscal+dt
       !if(move_particles)S_TIME = t
       utemp = ucharmod 
       tempt = ramp_frac_time*float(mx-1)*dx/utemp
       IF(soft_start) THEN 
          IF (t.LT.tempt) THEN
             cf = -one/dt
             cforig = cf
             !Write(*,*)'I AM HERE : ', tempt, soft_start
             cf=cf*(t/tempt)
             !if(imove.eq.1)move_particles = .FALSE.
          ELSE
             !Write(*,*)'I AM HERE HAHAHAHAHAH', soft_start,tempt,ramp_frac_time
             cf = -one/dt
             cforig=cf
             !if(imove.eq.1)move_particles = .TRUE.
          ENDIF
       ELSE
          cf = -one/dt
          cforig = cf
       end IF
    
       if(I_AM_NODE_ZERO) WRITE(*,'(A, i4,2(2x,g12.5),/,2(A30, g12.5,/))')'IDUMSTEP, t, tend = ', IDUMSTEP, t, tendused ,& 
            'cf original = ', cforig, &
            'cf used = ', cf
#if PARALLEL
       if(iscalon.eq.1)then
          uatminus1(1:my,1:mz) = ubcp(-1,1:my,1:mz,1)
          uatnxp2(1:my,1:mz) = ubcp(nx+2,1:my,1:mz,1)
       end if
#endif
     END SUBROUTINE compute_new_timestep
     
     SUBROUTINE interpolate_fields_to_new_mesh(delta_pos)
       USE nlmainarrays, Only : ubcp, nlbcp, pbcp, onlbcp
       IMPLICIT NONE
       REAL(prcn), Intent(in) :: delta_pos(ndim)
       INTEGER :: i, j, k, vcell(ndim), pcell(ndim), ib,ie,jb,je,kb,ke, idim, onew
       REAL(prcn) :: new_pos(ndim), new_ppos(ndim), p_newmesh, u_newmesh(ndim), onll(ndim), nll(ndim), dfll(ndim), ppll(ndim),pmean

       do k = 1, mz
          new_pos(3) = k + delta_pos(3)
          new_ppos(3) = new_pos(3)
          do j = 1, my
             new_pos(2) = j + delta_pos(2)
             new_ppos(2) = new_pos(2)
             do i = 1, mx1
                new_pos(1) = i + delta_pos(1)
                
                new_ppos(1) = (i + 0.5 + delta_pos(1)) ! Position of the pressure grid point corresponding to vel grid point i
                
                new_ppos(1) = new_ppos(1) - 0.5

                do idim = 1, ndim
                   if(new_pos(idim).lt.zero)then
                      vcell(idim) = int(new_pos(idim)-1)
                   else
                      vcell(idim) = int(new_pos(idim))
                   end if
                   if(new_ppos(idim).lt.zero)then
                      pcell(idim) = int(new_ppos(idim)-1)
                   else
                      pcell(idim) = int(new_ppos(idim))
                   end if
                end do
                p_newmesh = zero
                u_newmesh(1:ndim) = zero
                
                call interpolate_pdata(pcell,new_ppos,ppll,p_newmesh,1)
                call interpolate_udata(vcell,new_pos,ib&
                     &,ie,jb,je,kb,ke,u_newmesh,nll,onll,dfll, 0,1, 1, onew) 
                nlbcp(i,j,k,1:ndim) = u_newmesh(1:ndim)
                onlbcp(i,j,k,1) = p_newmesh
!!$                if(iglobstep.eq.2)then
!!$                   Write(*,*)'pbcp : ', pbcp(i,j,k)
!!$                   Write(*,*)'p on newmesh :', p_newmesh
!!$                   STOP
!!$                end if
             end do
          end do
       end do

       !Store interpolated velocities and pressure fields into ubcp and pbcp
       do idim = 1, ndim
          do k = 1, mz
             do j = 1, my
                do i = 1, mx1
                   ubcp(i,j,k,idim) = nlbcp(i,j,k,idim)
                   if(idim.eq.1)pbcp(i,j,k) = onlbcp(i,j,k,idim)
                enddo
             end do
          end do
       end do
       ubcp(mx,:,:,:) = ubcp(1,:,:,:)
       ! Convert solutions back to complex space
       do idim = 1, ndim
          do i = 1, mx
             ur1(:,:) = ubcp(i,:,:,idim)-umean(idim)
             CALL ff2rc(ur1(:,:),u(i,:,:,idim))
             if(i.le.mx1)CALL ff2rc(pbcp(i,:,:),p(i,:,:))
          end do
       end do
       Write(*,*)'pmean = ', SUM(pbcp(:,:,:))/real(mx1*my*mz,prcn)
     END SUBROUTINE interpolate_fields_to_new_mesh

  SUBROUTINE RUN_TIME_FILE_OPENER(funit, filename, formfile)
		USE general_funcs 

		Implicit NOne 

		Character(LEN=*) :: filename, formfile
		INTEGER, INTENT(out) :: funit 

		LOGICAL :: filexist, isopen


		if (from_post) then
			OPEN(unit=funit,file="tmp_"//trim(filename), form=formfile,status='replace')
		else
			INQUIRE(FILE=FILENAME,EXIST=filexist,OPENED=isopen)
			IF (.NOT.filexist) THEN
				funit = getnewunit(minunitno,maxunitno)
				OPEN(unit=funit,file=FILENAME, form=formfile,status='new')
			ELSEIF(filexist.AND..NOT.isopen) THEN
				IF(irestart.eq.0) THEN 
					funit = getnewunit(minunitno,maxunitno)
					OPEN(unit=funit,file=FILENAME, form=formfile,status='replace')
				ELSEIF(irestart.eq.1) THEN 
					funit = getnewunit(minunitno,maxunitno)
					OPEN(unit=funit,file=FILENAME, form=formfile, POSITION="append")
				ENDIF
			END IF
		endif
  end SUBROUTINE RUN_TIME_FILE_OPENER

  SUBROUTINE SCALAR_RUN_TIME_FILE_OPENER(funit, filename, formfile)
    USE general_funcs 
    
    Implicit NOne 
    
    Character(LEN=*) :: filename, formfile
    INTEGER, INTENT(out) :: funit 
    
    LOGICAL :: filexist, isopen
    
    INQUIRE(FILE=FILENAME,EXIST=filexist,OPENED=isopen)
    
    IF (.NOT.filexist) THEN
       
       funit = getnewunit(minunitno,maxunitno)
       OPEN(unit=funit,file=FILENAME, form=formfile,status='new')
       
    ELSEIF(filexist.AND..NOT.isopen) THEN
       IF(iscal_restart.eq.0) THEN 
          funit = getnewunit(minunitno,maxunitno)
          
          OPEN(unit=funit,file=FILENAME, form=formfile,status='replace')
          
       ELSEIF(iscal_restart.eq.1) THEN 
          funit = getnewunit(minunitno,maxunitno)
          OPEN(unit=funit,file=FILENAME, form=formfile, POSITION="append")
       ENDIF
    END IF
  end SUBROUTINE SCALAR_RUN_TIME_FILE_OPENER

  SUBROUTINE RESTART_FILE_OPENER(funit, filename, formfile)
    USE general_funcs 
    
    Implicit NOne 
    
    Character(LEN=*) :: filename, formfile
    INTEGER, INTENT(out) :: funit 
    
    LOGICAL :: filexist, isopen
    
    INQUIRE(FILE=FILENAME,EXIST=filexist,OPENED=isopen)
    
    IF (.NOT.filexist) THEN
       
       funit = getnewunit(minunitno,maxunitno)
       OPEN(unit=funit,file=FILENAME, form=formfile,status='new')
       
    ELSEIF(filexist.AND..NOT.isopen) THEN
       funit = getnewunit(minunitno,maxunitno)
       
       OPEN(unit=funit,file=FILENAME, form=formfile,status='replace')
       
    END IF
  END SUBROUTINE RESTART_FILE_OPENER


  SUBROUTINE calc_part_statistics(rks) 
    USE general_funcs
    USE randomno
    Implicit None 
    INTEGER :: rks

    REAL(prcn) :: mean_drag(nphases), mean_force(nphases,ndim)
    REAL(prcn) :: fluctv(nbody,ndim),mean_vel(nphases,ndim),&
         & sdev_fluctv(nphases,ndim),sdev_fluctf(nphases,ndim),&
         & flupartslip(ndim), mixmeanvel(ndim), flupartslipmod,&
         & grant(nphases), phasicslip(nphsc2,ndim)
    REAL(prcn) :: mean_contact_force(ndim), fcvel_corr, char_force
    INTEGER :: m,idim,iphs, partstart, partend, npart, phsc2count, jphs
    INTEGER, SAVE :: velinfounit
    LOGICAL, SAVE :: FIRST_TIME = .TRUE.
    CHARACTER(LEN=80) :: FILENAME, formfile

    if(.not.(rks.eq.itrmax))RETURN
    IF(I_AM_NODE_ZERO.and.first_time) THEN     
       FILENAME = TRIM(RUN_NAME)//'_part_info'//'.rst'
       formfile='unformatted'
       CALL  RUN_TIME_FILE_OPENER(unitpartinfo, FILENAME, formfile)

       velinfounit = getnewunit(minunitno,maxunitno)
       FILENAME = TRIM(RUN_NAME)//'_vel_info'//'.dat'
       formfile='formatted'
       CALL  RUN_TIME_FILE_OPENER(velinfounit, FILENAME, formfile)
       
       first_time = .FALSE.
    END IF
    write(unitpartinfo) t-dt
    write(unitpartinfo) xc(1:nbody, 1:ndim)
    write(unitpartinfo) velbdy(1:nbody, 1:ndim)
    write(unitpartinfo) force(1:nbody, 1:ndim)
    write(unitpartinfo) pres(1:nbody, 1:ndim)
    write(unitpartinfo) visc(1:nbody, 1:ndim)
    write(unitpartinfo) contact_force(1:nbody, 1:ndim)
    write(unitpartinfo) frame_vel(1:ndim)
    write(unitpartinfo) frame_accln(1:ndim)
    write(unitpartinfo) ufmean(1:ndim)

    mixmeanvel = zero
    grant = zero
    do iphs = 1, nphases
       phase_array(iphs)%mean_spec_vel(1:ndim) = zero
       partstart = phase_array(iphs)%pstart
       partend = phase_array(iphs)%pend
       do idim = 1, ndim
          phase_array(iphs)%mean_spec_vel(idim) =&
               & SUM(velbdy(partstart:partend,idim))&
               &/real(phase_array(iphs)%npart,prcn)
       end do
       mixmeanvel(1:ndim) = mixmeanvel(1:ndim) + (phase_array(iphs)&
            &%volfracg)*(phase_array(iphs)%mean_spec_vel(1:ndim))
       
       do m = partstart, partend
          do idim = 1, ndim
             grant(iphs) = grant(iphs) + (phase_array(iphs)&
                  &%mean_spec_vel(idim)-velbdy(m,idim))**2.d0
          end do
       end do
       grant(iphs) = grant(iphs)/(three*phase_array(iphs)%npart)
    end do
    mixmeanvel(:) = mixmeanvel(:)/maxvolfrac
    flupartslip(:) = (usmean(:)-ufmean(:))
    flupartslipmod = DOT_PRODUCT(flupartslip(1:ndim),&
         & flupartslip(1:ndim))
    flupartslipmod = DSQRT(flupartslipmod)
    
    phsc2count = 0
    do iphs = 1, nphases-1
       do jphs = iphs+1, nphases
          phsc2count = phsc2count + 1
          phasicslip(phsc2count,1:ndim) = phase_array(jphs)&
               &%mean_spec_vel(1:ndim) - phase_array(iphs)&
               &%mean_spec_vel(1:ndim)
       end do
    end do

    do idim = 1, ndim
       mean_contact_force(idim) = SUM(contact_force(1:nbody,idim))/real(nbody,prcn)
    end do
    char_force = 3.d0*pi*vis*flupartslipmod*(one-maxvolfrac)*char_length

    fcvel_corr = zero
    do m = 1, nbody
       do idim = 1, ndim
          fcvel_corr = fcvel_corr + contact_force(m,idim)*(phase_array(iphs)&
               &%mean_spec_vel(idim)-velbdy(m,idim))
       enddo
    end do
    fcvel_corr = fcvel_corr/real(nbody,prcn)

    Write(velinfounit,'(40(2x,g17.8))')t, t/t_conv/(1-maxvolfrac), t/t_vis,&
         & (one-maxvolfrac)*flupartslipmod*char_length/vis,&
         & mixmeanvel(1:ndim)/flupartslipmod, ((phasicslip(iphs,idim)&
         &/flupartslipmod, idim = 1, ndim),iphs=1,phsc2count),&
         & DSQRT(grant(1:nphases))/flupartslipmod,&
         & mean_contact_force(1:ndim)/char_force, fcvel_corr&
         &/(char_force*flupartslipmod)
    
  END SUBROUTINE calc_part_statistics
  
subroutine gener_lattice_mod(nbody,domlin,xc,dia1, dbdy)

  USE precision 
  USE constants 
  USE global_data, ONLY : RUN_NAME
  
  implicit none
  double precision, intent(in) ::domlin(3),dia1
  double precision, intent(out) :: xc(nbody,3)
  double precision ::dia, dmax, fac , ymin, ymax
  double precision :: rstep, rstep_buff, xtmp, ytmp, ztmp,lngth_spec1, dbdy(nbody), facdmax , dmaxbydmin
  Integer, INTENT(in) :: nbody 
  integer :: nprob1 , i, j, k, ntot, ii, jj, kk, nx, ny, nz, np1, n
  REAL(PRCN) :: ZP, XP, YP
  
  PRINT*, 'IN GENER LATTICE MOD, NBODY = ', NBODY
  
  !doml(:) = DOMLin(:)! - MAXVAL(dbdy(1:nbody))
  dia = 1.2*dia1	! so that particles don't touch each other to begin with
  dmaxbydmin = MAXVAL(DBDY(1:nbody))/MINVAL(DBDY(1:nbody))
  if(ABS(dmaxbydmin-one).lt.small_number) then
     fac = 1.05
  ELSE
     fac = 1.05
  end if
  

  !fac = 1.2
!!$    nx = (domlin(1)-half*dia)/dia
!!$    PRINT*,'nx = ', nx,domlin(1)
!!$    nz = ceiling((real(DOMLIN(3)-half*dia)/dia)) 
!!$    PRINT*,'nx = ', nz    
!!$    np1 = ceiling(real(nbody/nz))+ 1
!!$    ny = ceiling(real(np1/nx)) + 1
    
    !       Specifying the initial distribution of the particles
!!$    PRINT*,'nx = ', nx, 'ny = ', ny, 'nz = ', nz
    n = 1
    dmax =  dia1!MAXVAL(dbdy(1:nbody))
    
    facdmax = 0.01*dmax

    yp = dmax*fac/two
    DO While (n.lt.nbody) 
       !PRINT*,'dmax = ', dmax 
       !nz = ceiling((real(DOML(3)-dmax*fac)/(dmax*fac)))
       !nx = ceiling((doml(1)-dmax*fac)/(dmax*fac))
       nz = ceiling((real(DOMLIN(3)-dmax*fac)/(dmax*fac)))
       nx = ceiling((domlin(1)-dmax*fac)/(dmax*fac))
       

       do k = 1, nz
          zp =  dmax*half +  (k-1)*dmax*fac
          do i = 1, nx
             xp = dmax*half + (i-1)*dmax*fac
             XC(n,1) = xp
             XC(n,2) = yp
             XC(n,3) = zp
             IF(N.Eq.NBODY) GOTO 200
             n = n+1
                       
          end do
          
          !dmax =  MAXVAL(dbdy(n-1:nbody))
          !nx = floor((domlin(1)-dmax*fac)/(dmax*fac))
          !nx = ceiling((domlin(1)-half*dmax*fac)/dmax)

       end do
       
       !dmax =  MAXVAL(dbdy(n-1:nbody))
       yp = yp + dmax*fac

200 CONTINUE       
    end do
    



  open(2222, file = TRIM(RUN_NAME)//"_sphere_center_unfor.inp", form="unformatted", status="replace")

  open(unit=1000,file= TRIM(RUN_NAME)//"_lattice_dist.dat",form='formatted',status='replace')

  write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" '

  do i  = 1, nbody
     write(2222)xc(i,1),xc(i,2),xc(i,3),dbdy(i)*0.5d0
     
     WRITE(1000,'(4(2x,f12.8))') xc(i,1),xc(i,2),xc(i,3),dbdy(i)/2.d0
  end do

  close(2222, status = "keep")
  
  close(1000, status = "keep")
  
  
end subroutine gener_lattice_mod

  
subroutine gener_lattice(nbody1, nbody2,ibidisperse,volfrac1,volfrac2,dia1,dia2&
     &,percent_buff)

  USE precision 
  USE constants 
  USE global_data, ONLY : RUN_NAME

  implicit none

  double precision, intent(in) :: percent_buff, dia1,  dia2, volfrac1, volfrac2
  double precision ::numdens1, numdens2
  double precision :: rstep, rstep_buff, xtmp, ytmp, ztmp,lngth_spec1
  Integer, INTENT(out) :: nbody1, nbody2 
  integer :: nprob1 , i, j, k, ntot, ii, jj, kk

  LOGICAL, INTENT(in) ::ibidisperse
! dia2 = dia1*diaratio

  !    percent_buff = 1.0

  open(2222, file =  TRIM(RUN_NAME)//"_sphere_center_unfor.inp", form="unformatted", status="replace")

  open(unit=1000,file= TRIM(RUN_NAME)//"_lattice_dist.dat",form='formatted',status='replace')

  write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" '
!!!! bidisperse

  numdens1 = 6.*volfrac1*1./(pi*dia1**3)
  numdens2 = 6.*volfrac2*1./(pi*dia2**3)
  PRINT*,'volffracs = ', volfrac1, volfrac2, dia1, dia2
  print*,' Numdens of first species..',numdens1
  print*,' Numdens of second species..',numdens2

!!!! since we assume the domain size to be unity (-0.5,+0.5)

  rstep  = dia1/2.

  !! the hard sphere code does not like
  !! touching spheres; so add a 1% of diameter buffer (0.5% for radius)

  rstep_buff = rstep + percent_buff/2.*rstep/100. 

  nprob1 = int(1./(2.*rstep_buff))
  write(*,*)'number of spheres along one direction...',nprob1

  write(*,*)'Number of rows needed..',numdens1/(nprob1*nprob1)	
  lngth_spec1 = (int(numdens1/(nprob1*nprob1))+1)*2.*rstep_buff

  write(*,*)'Total length along z required...',lngth_spec1

  xtmp = rstep_buff
  ytmp = rstep_buff
  ztmp = rstep_buff

  ntot = 0
  nbody1 = 0 
  nbody2 = 0
  do k = 1, nprob1
     do j = 1, nprob1
        do i = 1, nprob1

           write(2222)xtmp,ytmp,ztmp,dia1/2.d0

           WRITE(1000,'(4(2x,f12.8))') xtmp-0.5,ytmp-0.5d0,ztmp-0.5d0,dia1/2.d0
           ntot = ntot + 1
           nbody1 = nbody1 + 1
           if(ntot.ge.numdens1) then 
              print*,'Total number of particles generated..',ntot
              ii = i
              jj = j
              kk = k
              print*,'final i, j, k...',ii,jj,kk
              goto 40
           endif

           xtmp  = xtmp + 2.*rstep_buff
        end do
        xtmp = rstep_buff
        ytmp = ytmp + 2.*rstep_buff 
     end do
     xtmp = rstep_buff
     ytmp = rstep_buff  
     ztmp = ztmp + 2.*rstep_buff 

     !if(ztmp.gt.volfrac1)then !this works because ztmp \in (0,1) and volfrac \in(0,1)

     if(ztmp.gt.lngth_spec1)then !this works because ztmp \in (0,1) and volfrac \in(0,1)
        print*,'!!!!!!! STOP !!!!!!!!!! : first species has exceeded its domain'
        print*,' total number of particles generated...',ntot
        stop
     end if
  end do

40 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! since we assume the domain size to be unity (-0.5,+0.5)
  if(ibidisperse)then

     rstep  = dia2/2.

     !! the hard sphere code does not like
     !! touching spheres; so add a 1% of diameter buffer (0.5% for radius)

     rstep_buff = rstep + percent_buff/2.*rstep/100. 

     nprob1 = int(1./(2.*rstep_buff))

     xtmp = rstep_buff
     ytmp = rstep_buff
     ztmp = rstep_buff

     ntot = 0

     do k = 1, nprob1
        do j = 1, nprob1
           do i = 1, nprob1
              write(2222)xtmp,ytmp,lngth_spec1+ztmp,dia2/2.d0

              nbody2 = nbody2 + 1

              WRITE(1000,'(4(2x,f12.8))') xtmp-half,ytmp-half,lngth_spec1+ztmp-half,dia2/2.d0

              ntot = ntot + 1
              if(ntot.ge.numdens2) then 
                 print*,'Total number of particles generated..',ntot
                 goto 50
              endif

              xtmp  = xtmp + 2.*rstep_buff
           end do
           xtmp = rstep_buff
           ytmp = ytmp + 2.*rstep_buff 
        end do
        xtmp = rstep_buff
        ytmp = rstep_buff  
        ztmp = ztmp + 2.*rstep_buff 
        !if(volfrac1+ztmp.gt.1.0)then !this works because ztmp \in (0,1) and volfrac \in(0,1)
        if(lngth_spec1+ztmp.gt.1.0)then !this works because ztmp \in (0,1) and volfrac \in(0,1)
           print*,'!!!!!!! STOP !!!!!!!!!! : second species has exceeded its domain'
           print*,' total number of particles generated...',ntot
           stop
	end if
     end do

50   continue
  end if

  close(2222, status = "keep")

  close(1000, status = "keep")


!!!! volfrac of species based on number density

  print*,'Vol frac of first species...',pi*dia1**3/6.*nbody1,' input...',volfrac1
  print*,'Vol frac of second species...',pi*dia2**3/6.*nbody2, ' input...',volfrac2


end subroutine gener_lattice

subroutine scale_to_grid_units(nbody,npart,nphases,my,mbox,xperiodic,percent_buf&
     &,xp_in,rad_in, min_part_sep, toscale)

  USE precision 
  USE constants 
  USE global_data, ONLY : RUN_NAME


  implicit none

  INTEGER, INTENT(in) :: nphases
  INTEGER, INTENT(inout) :: nbody
  INTEGER, INTENT(inout),DIMENSION(nphases):: npart
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nsp
  Logical, Intent(in) ::  xperiodic
  LOGICAL, OPTIONAL :: toscale
  integer :: n,ndim ,i,idim,j,my,mbox,ixcount, partend,partstart, iphs

  real(prcn), intent(in),DIMENSION(nphases) ::  percent_buf
  real(prcn), intent(in) :: min_part_sep
  real(prcn),intent(inout):: xp_in(nbody,3)
  real(prcn),intent(inout) :: rad_in(nbody)
  REAL(prcn) :: dist, overlap, tempx, tempy, tempz
  REAL(prcn) :: rmax(3), L(3), tempr(3), r
  real(prcn), dimension(nbody) :: radbdy
  real(prcn), dimension(nbody,3) :: xp
  integer, dimension(nbody)::  towrite	

  PRINT*,'IN SCALE TO GRID UNITS'
  !PRINT*,'percent_buf = ', percent_buf(1:nphases)
  !nbody =SIZE(xp,1)
  ndim = SIZE(xp,2)
  xp(1:nbody,1:3) = xp_in(1:nbody,1:3)
  radbdy(1:nbody) = rad_in(1:nbody)

  rad_in(1:nbody) = zero 
  xp_in(1:nbody,1:3) = zero 

  ALLOCATE(nsp(nphases))

  n = nbody

  nsp(1:nphases) = npart(1:nphases)

  nbody = 0
  npart(1:nphases) = 0

  do i=1,n
     towrite(i)=1
  enddo

  partstart = 1
  if(.not.present(toscale))toscale = .TRUE.
  
  if(toscale)then
     Write(*,*)' SCALING : ', toscale
     do iphs = 1, nphases
        partend = partstart + nsp(iphs)-1
        do i = partstart, partend
           do idim=1,ndim
              if(idim.eq.1)then
                 xp(i,idim)=xp(i,idim)*(mbox-1)+1.d0
              else
                 xp(i,idim)=xp(i,idim)*(my)+1
              endif
           enddo
           radbdy(i) = radbdy(i)/(one+percent_buf(iphs)*one/100.d0)
           radbdy(i) = radbdy(i)*(my)
        end do
        partstart = partend + 1
     end do
  end if

  Write(*,*)'Max val of XP : ', MAXVAL(XP(1:n,1)), mbox, my
  Write(*,*)'Max val of RADBDY : ', MAXVAL(RADBDY(1:n))
  !check if sphere centers are closer than min_part_sep grid points
  rmax(1) = (mbox-1)/two
  rmax(2:3) = my/two
  L(1) = mbox-1
  L(2:3) = my

  Do i = 1, n
     do j = 1, n
        if(i.ne.j.and.(towrite(j).eq.1))then
           r = zero
           do idim = 1, ndim
              tempr(idim) = xp(i,idim) - xp(j,idim) ! compute the separation
              ! in each dimension
              if((ABS(tempr(idim))).gt.rmax(idim)) then
                 if(tempr(idim).lt.zero) then
                    tempr(idim) = tempr(idim) + L(idim)
                 else
                    tempr(idim) = tempr(idim) - L(idim)
                 end if
              end if
              r = r + tempr(idim)**2.d0
           end do
           r = DSQRT(r)
!!$           if(i.eq.2.and.j.eq.22)then
!!$              write(*,'(A20,3(2x,g12.5))')'Position of i = ', XP(I,1:3)
!!$              write(*,'(A20,3(2x,g12.5))')'Position of j = ', XP(J,1:3)
!!$              Write(*,*) 'r = ', r, tempr(1:ndim), RMAX(1:NDIM)
!!$           end if
           if(r.le.radbdy(i)+radbdy(j)+min_part_sep)then
              towrite(i) = 0
              overlap = r - (radbdy(i)+radbdy(j)+min_part_sep)
              write(*,'(A20,2x,g12.5,A20, 2(2x,I6))')'Spheres closer than ',min_part_sep, '  grid points..', i,j
              Write(*,'(A30,2x, g12.5)')'OVERLAP = ',  overlap/(radbdy(1))
              write(*,'(A20,3(2x,g12.5))')'Position of i = ', XP(I,1:3)
              write(*,'(A20,3(2x,g12.5))')'Position of j = ', XP(J,1:3)
           end if
        end if
     end do
  end Do
  
  partstart = 1
  
  do iphs = 1, nphases
     partend = partstart + nsp(iphs)-1
     do i = partstart, partend
        if(towrite(i).eq.1)then
           nbody = nbody+1
           npart(iphs) = npart(iphs)+1
           rad_in(nbody) = radbdy(i)
           xp_in(nbody,1:3) = xp(i,1:3)
        end if
     end do
     partstart = partend + 1
  end do
  
end subroutine scale_to_grid_units

  end Module dependent_functions
