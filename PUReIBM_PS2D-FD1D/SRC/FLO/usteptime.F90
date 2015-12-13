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

MODULE usteptime
#include "ibm.h"
  USE precision 
  USE constants 
  USE global_data
  USE tridiagonal
  USE gauss_seidel
  implicit none 
Contains
  SUBROUTINE ustep(rks,j,k)
    IMPLICIT NONE 
    INTEGER, INTENT(in) ::  rks,j,k
    INTEGER :: n 

    !-----------------------------------------------------------------------
    !	local variables


    REAL(prcn) ::  clow,chigh
    REAL(prcn) ::  tmp(nx),tmpi(nx)
    REAL(prcn) :: ut(nx),uti(nx)
    REAL(prcn) :: a(nx),b(nx),c(nx)
    REAL(prcn) ::  d
    COMPLEX(prcn) ::  utemp(3), r(nx)


    INTEGER :: i, idim

    REAL(prcn):: auxu(nx), auxv(nx), alpha, beta, auxz(nx), vdotz, vdotyr,vdotyi, tmpur, tmpui

    REAL(prcn) :: send_temp(3), recv_temp(3)

    REAL(prcn) :: rout
    gauss_u = .true.
    !	write(*,*) 'Mean outletvelocity',uout(1,1,1),' before ustep'


    alpha = one
    beta = one


    !-----------------------------------------------------------------------
    !	set up coefficients of tridiagonal matrix
    DO i=1,nx !mx
       !utemp(:) = u(i,j,k,:)

       a(i)=(-g(i)*dt*vis*coef(rks,2)/dx2)
       b(i)=(one+dt*vis*coef(rks,2)*(g(i)*two/dx2+w2(j,k)))
       c(i)=(-g(i)*dt*vis*coef(rks,2)/dx2)


       auxu(i) = zero
       auxv(i) = zero
       auxz(i) = zero
    END DO
#if PARALLEL
    if(xstart.eq.1)then
       auxu(1) = alpha
       auxv(1) = c(nx)/beta
       b(1) = b(1) - auxu(1)*auxv(1)
    end if
    if(xend.eq.mx1) then
       auxu(nx) = beta
       auxv(nx) = c(nx)/alpha
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

    Do idim = 1, ndim
       DO i=1,nx !mx

          !-----------------------------------------------------------------------
          !	diffusion in y- and z- directions
          r(i)=u(i,j,k,idim)*(one-dt*vis*coef(rks,1)*w2(j,k))

          !-----------------------------------------------------------------------
          !	include convection terms

          !       sign change infront of dt done..
          !       to be consistent with Jamals convention
          !       remember nl term is negative..
          !       so negative sign, which should have been in front of dt,
          !       is now absorbed into the nl term
          r(i)=r(i)+dt*coef(rks,3)*nl(i,j,k,idim)+dt*ff(i,j,k,idim)   !*g(i)
          r(i)=r(i)+dt*coef(rks,4)*onl(i,j,k,idim)  !*g(i)

       ENDDO			! end loop over all mx points

#if PARALLEL
       DO i = 1,nx
#else
          DO i=2,nx
#endif

             !-----------------------------------------------------------------------
             !	pressure terms
             if(idim.eq.1)then
                r(i) = r(i)-dt*(coef(rks,1)+coef(rks,2))*(p(i,j,k)-p(i-1,j&
                     &,k))/dx 
             else if(idim.eq.2)then
                r(i)=r(i)-dt*(coef(rks,1)+coef(rks,2))*half*(p(i-1,j,k)&
                     &+p(i,j,k))*wy(j) 
             else if(idim.eq.3)then
                r(i)=r(i)-dt*(coef(rks,1)+coef(rks,2))*half*(p(i-1,j,k)&
                     &+p(i,j,k))*wz(k) 
             end if


             !-----------------------------------------------------------------------
             !	diffusion in x-direction

             r(i)=r(i)+g(i)*dt*vis*coef(rks,1)*(u(i-1,j,k,idim)+u(i+1,j,k&
                  &,idim)-two*u(i,j,k,idim))/dx2 

          ENDDO			!end loop over all 2-mx1 points

          !-----------------------------------------------------------------------
          !	special cases for i=1,mx (endpoints)
          !	no x-direction pressure gradient terms since they are used
          !	in the BCs for the pressure equation


          if(xperiodic) THEN  
#if !PARALLEL
             ! there will be a pressure gradient in x-direction as well
             ! because of the periodic bc    
             if(idim.eq.1)then
                r(1) = r(1)-dt*(coef(rks,1)+coef(rks,2))*(p(1,j,k)-p(mx1,j&
                     &,k))/dx 
             else if(idim.eq.2)then
                r(1) = r(1)-dt*(coef(rks,1)+coef(rks,2))*half*(p(mx1,j,k)&
                     &+p(1,j,k))*wy(j) 
             else if(idim.eq.3)then
                r(1) = r(1)-dt*(coef(rks,1)+coef(rks,2))*half*(p(mx1,j,k)&
                     &+p(1,j,k))*wz(k) 
             end if
#endif
          endif


          if(xperiodic) THEN  
#if !PARALLEL
             r(1)=r(1)+g(1)*dt*vis*coef(rks,1)*(u(2,j,k,idim)-two*u(1,j,k&
                  &,idim)+u(mx1,j,k,idim))/dx2 
#endif
          end if

          !-----------------------------------------------------------------------
          !	implicit timestepping
          !	invert tridiagonal matrix to get estimate for u,v,w

          !solve for u component of velocity
          DO i=1,nx !mx
             tmp(i)=dreal(r(i))
             tmpi(i)=dimag(r(i))
          ENDDO

          if(xperiodic) then
            !gauss_u = .FALSE.
            !if(j.eq.1.and.k.eq.1)gauss_u = .TRUE.
            
#if PARALLEL 
             !PRINT*,' CALLING TRIDAG FROM USTEPTIME FOR j = ', j, ' k = ', k
             call mpi_tridag3(a(1:nx),b(1:nx),c(1:nx),tmp(1:nx),auxu(1:nx),tmpi(1:nx),ut(1:nx),auxz(1:nx),uti(1:nx),nx)
             send_temp(1) = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
             send_temp(2) = auxv(1)*ut(1) + auxv(nx)*ut(nx)
             send_temp(3) = auxv(1)*uti(1) + auxv(nx)*uti(nx)
             GLOBAL_DOUBLE_SUM(send_temp(1),recv_temp(1),3,decomp_group)
             vdotz = recv_temp(1)
             vdotyr = recv_temp(2)
             vdotyi = recv_temp(3)
#else
             !PRINT*,' CALLING TRIDAG FROM USTEPTIME FOR j = ', j, ' k = ', k
             call tridag3(a(1:nx),b(1:nx),c(1:nx),tmp(1:nx),auxu(1:nx),tmpi(1:nx),ut(1:nx),auxz(1:nx),uti(1:nx),nx)
             vdotz = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
             vdotyr = auxv(1)*ut(1) + auxv(nx)*ut(nx)
             vdotyi = auxv(1)*uti(1) + auxv(nx)*uti(nx)
#endif

             do i=1,nx !mx1
                tmpur = ut(i) - vdotyr*auxz(i)/(one + vdotz)
                tmpui = uti(i) - vdotyi*auxz(i)/(one + vdotz)
                
                u(i,j,k,idim) = dcmplx(tmpur,tmpui)
             enddo
#if !PARALLEL
             u(mx,j,k,idim) = u(1,j,k,idim)
#endif
          endif
#if 0
          if(k.gt.1)then
             !             PRINT*,' j = ', j, ' k = ', k, ' idim = ', idim
             !    if(idim.eq.3)then
             do i = 1, nx
!!$                a(i)=-g(i)*dt*vis*coef(rks,2)/dx2
!!$                b(i)=one+dt*vis*coef(rks,2)*(g(i)*two/dx2+w2(j,k))
!!$                c(i)=-g(i)*dt*vis*coef(rks,2)/dx2


!!$               if(i.eq.1)then
!!$                   rout = a(i)*dreal(u(mx1,j,k,idim)) + b(i)*dreal(u(i,j,k,idim)) + c(i)*dreal(u(i+1,j,k,idim))
!!$                else if(i.eq.mx1)then
!!$                   rout = a(i)*dreal(u(i-1,j,k,idim)) + b(i)*dreal(u(i,j,k,idim)) + c(i)*dreal(u(1,j,k,idim))
!!$                else
!!$                   rout = a(i)*dreal(u(i-1,j,k,idim)) + b(i)*dreal(u(i,j,k,idim)) + c(i)*dreal(u(i+1,j,k,idim))
!!$                end if
                Write(*,'4(2x,g17.8)')GLOBAL_INDEX(i), dreal(u(i,j,k,1)), dreal(r(i))
             end do
             PARALLEL_FINISH()
             STOP
          end if
#endif

       END DO
       gauss_u = .false.
     END SUBROUTINE ustep
END MODULE usteptime










