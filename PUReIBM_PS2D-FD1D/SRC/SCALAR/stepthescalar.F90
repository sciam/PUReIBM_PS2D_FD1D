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

!///////////////////////////////////////////////////////////////////////
!	calculate scalar field at next time level
!       Following the same rk3 stepping procedure.
!-----------------------------------------------------------------------
MODULE stepthescalar
#include "../FLO/ibm.h"
  USE scalar_data
  USE tridiagonal
  USE nlphi
  USE nlmainarrays,  phir=>nlbcp
  USE gauss_seidel
  USE general_funcs
  USE precision 
  IMPLICIT NONE 
CONTAINS 
  SUBROUTINE scalstep(sflag,rks)
    Use nlarrays, Only : uf1,ur1
    IMPLICIT NONE 
    INTEGER, INTENT(in) :: sflag, rks
    INTEGER ::  j,k, isp, i
    Integer, Save :: unit_inn_itn
    Real(prcn) :: dtorig, phimloc, phim
    LOGICAL:: filexist, isopen
    !-----------------------------------------------------------------------
    !	step through wavenumbers 
    if(I_AM_NODE_ZERO)then
       CALL screen_separator(80,'*')
       WRITE(*,'(A25)') 'IN SCALAR ROUTINES'
    end if
    call nlphi_graduphi(rks)
    !PRINT*,'source AVG = ',SUM(sourcephif(1:mx1,1,1,1))/(mx1),(SUM(phif(2:4,1,1,1))+SUM(phif(6:8,1,1,1)))/6.d0

    !DO i = 1, mxf 
    !   PRINT*,'sou = ',i, sourcephif(i,1,1,1)
    !end DO

    !DO i = 1, mxf 
    !   PRINT*,'phif = ',i, phif(i,1,1,1)
    !end DO

    gauss_scal = .true. 
    DO isp = 1,nspmx
       DO j=1,my2
          DO k=1,mz

             CALL phistep(rks,j,k,isp)


             iter_scal = iter_scal+iter_gauss
          ENDDO
       ENDDO
       CALL communicate_scalar(isp)
    ENDDO
    gauss_scal = .false.

    heat_ratio_old(:) = heat_ratio(:)

    heat_ratio(:) = heat_ratio_new(:)

  END SUBROUTINE scalstep

  SUBROUTINE phistep(rks,j,k, isp)

    !	
    INTEGER, INTENT(in) ::  rks,j,k, isp
    INTEGER :: n 

    !-------------------------------------------------------------
    !local variables

    REAL(prcn) ::  clow,chigh
    REAL(prcn) ::  tmp(nx),tmpi(nx)
    REAL(prcn) ::  ut(0:nx+1),uti(0:nx+1)
    REAL(prcn) ::  a(nx),b(nx),c(nx)
    REAL(prcn) ::  d
    REAL(prcn):: auxu(nx), auxv(nx), alpha, beta, auxz(nx), vdotz, vdotyr,vdotyi, tmpphir, tmpphii
    COMPLEX(prcn) ::  r1(nx)
    INTEGER i,il
    REAL(prcn) :: send_temp(3), recv_temp(3)

    alpha = heat_ratio_new(isp) !one
    beta = one/heat_ratio_new(isp) !one

    DO i=1,nx !mx
       r1(i)=czero
       !-----------------------------------------------------------------------
       !set up coefficients of tridiagonal matrix

       a(i)=-dt*gamma(isp)*coef(rks,2)/dx2
       b(i)=one+dt*gamma(isp)*coef(rks,2)*(two/dx2+w2(j,k))
       c(i)=-dt*gamma(isp)*coef(rks,2)/dx2

       auxu(i) = zero
       auxv(i) = zero
       auxz(i) = zero


       !ADD the value of scalar at the previous time step to RHS
       !since CN scheme is being used, add the laplacian of scalar
       !in the y and z directions.
       
       r1(i) = phif(i,j,k,isp)*(one-dt*coef(rks,1)*gamma(isp)*w2(j,k))


       !--------------------------------------------------------------
       !	include convection terms...  

       r1(i)=r1(i) + dt*coef(rks,3)*nlphif(i,j,k,isp)+dt*ffphi(i,j,k,isp)
       r1(i)=r1(i) + dt*coef(rks,4)*onlphif(i,j,k,isp)

       !Add the source term also ....
       if(sourcepresent) r1(i)= r1(i)+ sourcephif(i,j,k,isp)*dt 

    ENDDO
    
    if(xstart.eq.1)a(1) = a(1)*heat_ratio_new(isp)
    if(xend.eq.mx1)c(nx) = c(nx)/heat_ratio_new(isp)
#if PARALLEL
    DO i=1,nx
#else
       DO i=2,nx !mx1
#endif
          !-----------------------------------------------------------------------
          !diffusion in x-direction at the present time...
          !this is done to be consistent with CN scheme. RG-06/01/04

          r1(i) = r1(i) + dt*gamma(isp)*coef(rks,1)*(phif(i-1,j,k,isp)&
               &+phif(i+1,j,k,isp)-two*phif(i,j,k,isp))/dx2 

       ENDDO
       !---------------------------------
       !special cases for i=1,mx (endpoints)

#if !PARALLEL    
       IF(xperiodic) then 
          r1(1) = r1(1) + dt*gamma(isp)*coef(rks,1)*(phif(2,j,k,isp)-two&
               &*phif(1,j,k,isp)+phif(mx1,j,k,isp)*heat_ratio(isp))/dx2 
       end IF
#endif


       !invert tridiagonal matrix to get estimate for u,v,w

       DO i=1,nx !mx
          tmp(i)=dreal(r1(i))
          tmpi(i)=dimag(r1(i))

       ENDDO
#if PARALLEL
       if(xstart.eq.1)then
          auxu(1) = alpha
          auxv(1) = c(nx)!/beta
          b(1) = b(1) - auxu(1)*auxv(1)
       end if
       if(xend.eq.mx1) then
          auxu(nx) = beta
          auxv(nx) = a(1)!/alpha
          b(nx) = b(nx) - auxu(nx)*auxv(nx)
       end if
#else
       auxu(1) = alpha
       auxu(nx) = beta
       auxv(1) = c(nx)/auxu(nx)
       auxv(nx) = a(1)/auxu(1)
       b(1) = b(1) - auxu(1)*auxv(1)
       b(nx) = b(nx) - auxu(nx)*auxv(nx)
#endif
       if(xperiodic) then
#if PARALLEL 
          call mpi_tridag3(a(1:nx),b(1:nx),c(1:nx),tmp(1:nx),auxu(1:nx),tmpi(1:nx),ut(1:nx),auxz(1:nx),uti(1:nx),nx)
          send_temp(1) = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
          send_temp(2) = auxv(1)*ut(1) + auxv(nx)*ut(nx)
          send_temp(3) = auxv(1)*uti(1) + auxv(nx)*uti(nx)
          
          GLOBAL_DOUBLE_SUM(send_temp(1),recv_temp(1),3,decomp_group)
          
          vdotz = recv_temp(1)
          vdotyr = recv_temp(2)
          vdotyi = recv_temp(3)
#else
          call tridag3(a(1:nx),b(1:nx),c(1:nx),tmp(1:nx),auxu(1:nx),tmpi(1:nx),ut(1:nx),auxz(1:nx),uti(1:nx),nx)
          vdotz = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
          vdotyr = auxv(1)*ut(1) + auxv(nx)*ut(nx)
          vdotyi = auxv(1)*uti(1) + auxv(nx)*uti(nx)
#endif
          DO i=1,nx !mx1
             tmpphir = ut(i) - vdotyr*auxz(i)/(one + vdotz)
             tmpphii = uti(i) - vdotyi*auxz(i)/(one + vdotz)

             phif(i,j,k,isp)=dcmplx(tmpphir,tmpphii)
          ENDDO
#if !PARALLEL
          phif(mx,j,k,isp) = phif(1,j,k,isp)/heat_ratio_new(isp)
#endif
       endif
#if 0
       if(k.eq.1.and.j.eq.1)then
          do i = 1, nx
             Write(*,'4(2x,g17.8)')GLOBAL_INDEX(i), dreal(phif(i,j,k,1)), dreal(r1(i))
          end do
!!$          PARALLEL_FINISH()
!!$          STOP
       end if
#endif

     END SUBROUTINE phistep

     SUBROUTINE communicate_scalar(isp)
       Use nlarrays, Only : uf1
       IMPLICIT NONE
       INTEGER, Intent(in) :: isp
       INTEGER :: i,j,k

#if PARALLEL
       i = 1
       VECSENDRECV(phif(i,1,1,isp),1,ucslice,fromproc,1,phif(nx+1,1,1,isp),1,toproc,1,decomp_group,status)
       i = nx
       VECSENDRECV(phif(i,1,1,isp),1,ucslice,toproc,0,phif(0,1,1,isp),1,fromproc,0,decomp_group,status)
       do k = 1, mz
          do j = 1,my2
             if(xstart.eq.1) phif(0,j,k,isp) = heat_ratio_new(isp)*phif(0,j,k,isp)
             if(I_AM_LAST_PROC)phif(nx+1,j,k,isp) = phif(nx+1,j,k,isp)/heat_ratio_new(isp)
          end do
       end do
#endif

     END SUBROUTINE communicate_scalar


   END MODULE stepthescalar
