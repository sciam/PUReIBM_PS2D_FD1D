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

module poisson
#include "ibm.h"
  USE precision 
  USE constants 
  USE global_data
  USE tridiagonal
  USE gauss_seidel
  USE bcsetarrays
  IMPLICIT NONE 
  !///////////////////////////////////////////////////////////////////
  !	Calculate the pressure field using AB. The boundary conditions
  !	are such that only the flow field from the previous two
  !	timesteps (at the boundaries) is needed to solve the system.
  !	The pressure is estimated at (t+dt)
  !--------------------------------------------------------------------

  
COntains 

  subroutine pressure(rks,j,k)
    implicit none 
    integer, Intent(in) ::  rks,j,k

    !-----------------------------------------------------------------------
    !	local variables

    real(prcn) ::  a(nx),b(nx),c(nx)
    complex(prcn) :: tmpc, tmpp
    real(prcn) ::  dv(nx),dvi(nx)
    real(prcn) ::  dvt(nx),dvti(nx)
    integer :: i
    REAL(prcn):: auxu(nx), auxv(nx), alpha, beta, auxz(nx),vdotz,vdotyr,vdotyi, tmppr, tmppi, rout
    REAL(prcn) :: send_temp(3), recv_temp(3)
    
    alpha = one
    beta = one
    gauss_p= .true.
    
    do i=1,nx !mx1

       !       should actually be multiplied with negative sign here..
       !       to be consistent with Jamal's convention this is 
       !       not done...since nl's are negative..IOW, negative sign
       !       absorbed into nl term

       a(i)=one/dx2-qrt*w2(j,k)
       b(i)=-two/dx2-half*w2(j,k)
       c(i)=one/dx2-qrt*w2(j,k)
!!$       
!!$
!!$       a(i)=one/dx2
!!$       b(i)=-two/dx2-w2(j,k)
!!$       c(i)=one/dx2

       dv(i)=zero
       dvi(i)=zero
       auxu(i) = zero
       auxv(i) = zero
       auxz(i) = zero
    enddo
    !-----------------------------------------------------------------------
    !	special cases for i=1,mx1 (endpoints)

    !       Sign change done...reverted to Jamals convention

    !-----------------------------------------------------------------------
    !	calculate source terms at interior points

    do i=1,nx !mx1
       tmpc=czero
       tmpc=tmpc+coef(rks,3)*(nl(i+1,j,k,1)-nl(i,j,k,1))/dx + (ff(i+1,j,k,1)-ff(i,j,k,1))/dx
       tmpc=tmpc+coef(rks,3)*wy(j)*(nl(i+1,j,k,2)+nl(i,j,k,2))/2 + wy(j)*(ff(i+1,j,k,2)+ff(i,j,k,2))/2
       tmpc=tmpc+coef(rks,3)*wz(k)*(nl(i+1,j,k,3)+nl(i,j,k,3))/2 + wz(k)*(ff(i+1,j,k,3)+ff(i,j,k,3))/2

       tmpc=tmpc+coef(rks,4)*(onl(i+1,j,k,1)-onl(i,j,k,1))/dx
       tmpc=tmpc+coef(rks,4)*wy(j)*(onl(i+1,j,k,2)+onl(i,j,k,2))/2.d0
       tmpc=tmpc+coef(rks,4)*wz(k)*(onl(i+1,j,k,3)+onl(i,j,k,3))/2.d0

       tmpc=half*(g(i)+g(i+1))*tmpc


       dv(i)=dreal(tmpc)/((coef(rks,1)+coef(rks,2)))
       dvi(i)=dimag(tmpc)/((coef(rks,1)+coef(rks,2)))

    enddo


    !-----------------------------------------------------------------------
    !	integrate pressure field
    
    if(xperiodic) then
       if(j.eq.1.and.k.eq.1)then
#if PARALLEL
          if(xend.eq.mx1)then
             call mpi_tridag(a(1:nx-1),b(1:nx-1),c(1:nx-1),dv(1:nx-1),dvt(1:nx-1),nx-1)
             dvt(nx) = zero
          else
             call mpi_tridag(a(1:nx),b(1:nx),c(1:nx),dv(1:nx),dvt(1:nx),nx)
          end if
#else
          call tridag(a(1:nx-1),b(1:nx-1),c(1:nx-1),dv(1:nx-1),dvt(1:nx-1),nx-1)
          dvt(nx) = zero
#endif
          dvti(1:nx) = zero
          vdotz = zero
          vdotyr = zero
          vdotyi = zero
       else
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
         !gauss_p = .FALSE.
         !if(j.eq.1.and.k.eq.2)gauss_p = .TRUE.
#if PARALLEL          
          !PRINT*,' CALLING TRIDAG FROM POISSON FOR j = ', j, ' k = ', k
          call mpi_tridag3(a(1:nx),b(1:nx),c(1:nx),dv(1:nx),auxu(1:nx),dvi(1:nx),dvt(1:nx),auxz(1:nx),dvti(1:nx),nx)
          send_temp(1) = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
          send_temp(2) = auxv(1)*dvt(1) + auxv(nx)*dvt(nx)
          send_temp(3) = auxv(1)*dvti(1) + auxv(nx)*dvti(nx)
          GLOBAL_DOUBLE_SUM(send_temp(1),recv_temp(1),3,decomp_group)
          vdotz = recv_temp(1)
          vdotyr = recv_temp(2)
          vdotyi = recv_temp(3)
#else
          !PRINT*,' CALLING TRIDAG FROM POISSON FOR j = ', j, ' k = ', k
          call tridag3(a(1:nx),b(1:nx),c(1:nx),dv(1:nx),auxu(1:nx),dvi(1:nx),dvt(1:nx),auxz(1:nx),dvti(1:nx),nx)
          vdotz = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
          vdotyr = auxv(1)*dvt(1) + auxv(nx)*dvt(nx)
          vdotyi = auxv(1)*dvti(1) + auxv(nx)*dvti(nx)
#endif
       end if

       do i=1,nx !mx1
          tmppr = dvt(i) - vdotyr*auxz(i)/(one + vdotz)
          tmppi = dvti(i) - vdotyi*auxz(i)/(one + vdotz)
          tmpp=dcmplx(tmppr,tmppi)
          
          p(i,j,k)=p(i,j,k) + prf*(tmpp - p(i,j,k))
#if 0
          if(k.eq.1.and.j.eq.1) then
            Write(*,'(I,3(2x,g17.8))')GLOBAL_INDEX(i), dreal(p(i,j,k)), dv(i), dreal(ff(i,j,k,1))
          endif
#endif
       enddo
#if 0
       if(k.eq.1.and.j.eq.1) then
         PARALLEL_FINISH()
         STOP
       end if
#endif
    gauss_p= .false.

    RETURN

    if(j.eq.7)then
       do i = 1, nx
       a(i)=one/dx2-qrt*w2(j,k)
       b(i)=-two/dx2-half*w2(j,k)
       c(i)=one/dx2-qrt*w2(j,k)
          if(i.eq.1)then
    !         rout = a(i)*dreal(p(mx1,j,k)) + b(i)*dreal(p(i,j,k)) + c(i)*dreal(p(i+1,j,k))
          else if(i.eq.nx)then
     !        rout = a(i)*dreal(p(i-1,j,k)) + b(i)*dreal(p(i,j,k)) + c(i)*dreal(p(1,j,k))
          else
             rout = a(i)*dreal(p(i-1,j,k)) + b(i)*dreal(p(i,j,k)) + c(i)*dreal(p(i+1,j,k))
          end if
          Write(*,'(4(2x,g17.8))') i,p(i,j,k),dv(i)
       end do
       PARALLEL_FINISH()
       STOP

    end if
    do i = 1, nx
       if(i.eq.1)then
          rout = b(i)*dvt(i) + c(i)*dvt(i+1)
       else if(i.eq.mx1)then
          rout = a(i)*dvt(i-1) + b(i)*dvt(i)
       else
          rout = a(i)*dvt(i-1) + b(i)*dvt(i) + c(i)*dvt(i+1)
       end if
       Write(*,'(4(2x,g17.8))') rout, dv(i)
    end do
    PRINT*,'j = ', j, 'k = ', k
    STOP
    
    end if
    
  end subroutine pressure
  
  subroutine divcorr(rks, j,k)
#if PARALLEL
    USE nlarrays, ONLY: westbuf => uf1, eastbuf => uf2
#endif
    implicit none 
    INTEGER, INTENT(IN) :: rks
    !-----------------------------------------------------------------------
    !	local variables
    real(prcn) ::  a(nx),b(nx),c(nx)

    real(prcn) ::  dv(nx),dvi(nx)
    real(prcn) ::  dv1(nx),dv1i(nx)
    real(prcn) ::  dmax,dmaxi, prftmp, alp_coef_dt
    complex(prcn) ::  tmpc,div(nx), utemp(3), tmpp
    integer ::  i,j,k,n
    REAL(prcn):: auxu(nx), auxv(nx), alpha, beta, auxz(nx), vdotz,vdotyr,vdotyi, tmpdivr, tmpdivi, rout
    REAL(prcn) :: send_temp(3), recv_temp(3)

    !-----------------------------------------------------------------------
    !	calculate divergence of intermediate velocity field
    !	form divergence in Fourier space on pressure grid

    alp_coef_dt = ((coef(rks,1)+coef(rks,2))*dt)
    gauss_phi= .true.

    alpha = one
    beta = one

    do i=1,nx !mx1

       a(i)=one/dx2-qrt*w2(j,k)
       b(i)=-two/dx2-half*w2(j,k)
       c(i)=one/dx2-qrt*w2(j,k)
!!$
!!$       a(i)=one/dx2
!!$       b(i)=-two/dx2-w2(j,k)
!!$       c(i)=one/dx2

       dv(i)=zero
       dvi(i)=zero
       div(i)=czero

       auxu(i) = zero
       auxv(i) = zero
       auxz(i) = zero
    enddo

    do i=1,nx !mx1
       
       tmpc=czero
       tmpc=(u(i+1,j,k,1)-u(i,j,k,1))/dx
       tmpc=tmpc+half*wy(j)*(u(i,j,k,2)+u(i+1,j,k,2))
       tmpc=tmpc+half*wz(k)*(u(i,j,k,3)+u(i+1,j,k,3))
       
       
       !       multiplying with g(i) to remove source term in poisson equation
       
       dv(i)=half*(g(i)+g(i+1))*dreal(tmpc)
       dvi(i)=half*(g(i)+g(i+1))*dimag(tmpc)        

       dv(i) = dv(i)/alp_coef_dt 
       dvi(i) = dvi(i)/alp_coef_dt

    enddo

    !-----------------------------------------------------------------------
    !
    if(xperiodic) then
       if(j.eq.1.and.k.eq.1)then
#if PARALLEL
          if(xend.eq.mx1)then
             call mpi_tridag(a(1:nx-1),b(1:nx-1),c(1:nx-1),dv(1:nx-1),dv1(1:nx-1),nx-1)
             dv1(nx) = zero
          else
             call mpi_tridag(a(1:nx),b(1:nx),c(1:nx),dv(1:nx),dv1(1:nx),nx)
          end if
#else
          call tridag(a(1:nx-1),b(1:nx-1),c(1:nx-1),dv(1:nx-1),dv1(1:nx-1),nx-1)
          dv1(nx) = zero
#endif
      
          dv1i(1:nx) = zero
          vdotz = zero
          vdotyr = zero
          vdotyi = zero
       else
#if PARALLEL
          if(xstart.eq.1)then
             auxu(1) = alpha
             auxv(1) = c(nx)/beta
             b(1) = b(1) - auxu(1)*auxv(1)
          end if
          if(xend.eq.mx1) then
             auxu(nx) = beta
             auxv(nx) = a(1)/alpha
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

#if PARALLEL          
          call mpi_tridag3(a(1:nx),b(1:nx),c(1:nx),dv(1:nx),auxu(1:nx)&
               &,dvi(1:nx),dv1(1:nx),auxz(1:nx),dv1i(1:nx),nx)
          send_temp(1) = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
          send_temp(2) = auxv(1)*dv1(1) + auxv(nx)*dv1(nx)
          send_temp(3) = auxv(1)*dv1i(1) + auxv(nx)*dv1i(nx)
	  GLOBAL_DOUBLE_SUM(send_temp(1),recv_temp(1),3,decomp_group)
          vdotz = recv_temp(1)
          vdotyr = recv_temp(2)
          vdotyi = recv_temp(3)
#else
          call tridag3(a(1:nx),b(1:nx),c(1:nx),dv(1:nx),auxu(1:nx)&
               &,dvi(1:nx),dv1(1:nx),auxz(1:nx),dv1i(1:nx),nx)
          vdotz = auxv(1)*auxz(1) + auxv(nx)*auxz(nx)
          vdotyr = auxv(1)*dv1(1) + auxv(nx)*dv1(nx)
          vdotyi = auxv(1)*dv1i(1) + auxv(nx)*dv1i(nx)
#endif
       end if
    end if

    do i=1,nx !mx1
       tmpdivr = dv1(i) - vdotyr*auxz(i)/(one + vdotz)
       tmpdivi = dv1i(i) - vdotyi*auxz(i)/(one + vdotz)

       div(i)=dcmplx(tmpdivr,tmpdivi)
    enddo

#if PARALLEL
    westbuf(j,k) = div(1)
    eastbuf(j,k) = div(nx)
#endif
    !-----------------------------------------------------------------------
    !	correct intermediate velocity field divergence error
    !	
    !	u(t+dt)=u(int)-grad(phi)
    !
    !-----------------------------------------------------------------------
    !       BEWARE! delta t has been absorbed into phi
    if(xperiodic) then
#if !PARALLEL    
       u(1,j,k,1) = u(1,j,k,1) - alp_coef_dt*(div(1)-div(mx1))/dx
       u(1,j,k,2) = u(1,j,k,2) - alp_coef_dt*wy(j)*half*(div(1)+div(mx1))
       u(1,j,k,3) = u(1,j,k,3) - alp_coef_dt*wz(k)*half*(div(1)+div(mx1))

#else
       u(1,j,k,1) = u(1,j,k,1) - alp_coef_dt*(div(1))/dx!-div(mx1))/dx
       u(1,j,k,2) = u(1,j,k,2) - alp_coef_dt*wy(j)*half*(div(1))!+div(mx1))
       u(1,j,k,3) = u(1,j,k,3) - alp_coef_dt*wz(k)*half*(div(1))!+div(mx1))

#endif
    endif

    do i=2,nx  
       u(i,j,k,1)=u(i,j,k,1)-alp_coef_dt*(div(i)-div(i-1))/dx
       u(i,j,k,2)=u(i,j,k,2)-alp_coef_dt*half*wy(j)*(div(i-1)+div(i))
       u(i,j,k,3)=u(i,j,k,3)-alp_coef_dt*half*wz(k)*(div(i-1)+div(i))

       !u(i,j,k,1)=u(i,j,k,1)-g(i)*(div(i)-div(i-1))/dx
       !u(i,j,k,2)=u(i,j,k,2)-g(i)*half*wy(j)*(div(i-1)+div(i))
       !u(i,j,k,3)=u(i,j,k,3)-g(i)*half*wz(k)*(div(i-1)+div(i))

    enddo

#if !PARALLEL       
    if(xperiodic) u(mx,j,k,:) = u(1,j,k,:)
#endif       
    !-----------------------------------------------------------------------
    !	update pressure field
    !	A positive sign should be here...
    !	P = P* + phi

    do i=1,nx !mx1!=mx-1
       if(i.eq.1)then
#if PARALLEL
!          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((-2*div(i)+div(i+1))/dx2 - w2(j,k)*(two*div(i)+div(i+1))/4.d0)
          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((-2*div(i)+div(i+1))/dx2 - w2(j,k)*div(i))
#else
!          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(mx1)-2*div(i)+div(i+1))/dx2 -w2(j,k)*(div(mx1)+2*div(i)+div(i+1))/4.)
          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(mx1)-2*div(i)+div(i+1))/dx2 -w2(j,k)*div(i))
#endif           
       elseif(i.eq.nx)then
#if PARALLEL
!          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(i-1)-2*div(i))/dx2 - w2(j,k)*(div(i-1)+two*div(i))/4.d0)
          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(i-1)-2*div(i))/dx2 - w2(j,k)*div(i))
#else
!          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(i-1)-2*div(i)+div(1))/dx2 -w2(j,k)*(div(i-1)+2*div(i)+div(1))/4.)
          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(i-1)-2*div(i)+div(1))/dx2 -w2(j,k)*div(i))
#endif
       else
!          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(i-1)-2*div(i)+div(i+1))/dx2 - w2(j,k) * (div(i-1)+2*div(i)+div(i+1))/4.)
          p(i,j,k)=p(i,j,k)+div(i) - (coef(rks,2)*dt)*vis*((div(i-1)-2*div(i)+div(i+1))/dx2 - w2(j,k) * div(i) )
       endif

       !   Changed by Madhu: the positive sign infront of vis/2 was negative 
       !  earlier
       !Reverted the positive sign to negative in front of vis (no division by 2 needed, as it is made consistent with coef's), RG 

       !       divglobal(i,j,k)=div(i)/dt

       !-----------------------------------------------------------------------
       !	calculate new divergence
#if !PARALLEL
       tmpc=czero
       if(i.eq.nx)then
          tmpc=(u(1,j,k,1)-u(i,j,k,1))/dx
          
          tmpc=tmpc+half*wy(j)*(u(i,j,k,2)+u(1,j,k,2))
          tmpc=tmpc+half*wz(k)*(u(i,j,k,3)+u(1,j,k,3))
       else
          tmpc=(u(i+1,j,k,1)-u(i,j,k,1))/dx

          tmpc=tmpc+half*wy(j)*(u(i,j,k,2)+u(i+1,j,k,2))
          tmpc=tmpc+half*wz(k)*(u(i,j,k,3)+u(i+1,j,k,3))
       end if
       divloc = divloc + tmpc
!!$          total_div = total_div + tmpc
       !     	  dv(i)=dreal(g(i)*tmpc)
       !	  dvi(i)=imag(g(i)*tmpc)
#endif

    enddo


    !write(*,*)'new div  = ', dv(8:16),dvi(8:16)

    gauss_phi= .false.

    RETURN

    if(k.gt.1)then
       do i = 1, nx
          a(i)=one/dx2
          b(i)=-two/dx2-w2(j,k)
          c(i)=one/dx2

          if(i.eq.1)then
             rout = a(i)*dreal(div(mx1)) + b(i)*dreal(div(i)) + c(i)*dreal(div(i+1))
          else if(i.eq.mx1)then
             rout = a(i)*dreal(div(i-1)) + b(i)*dreal(div(i)) + c(i)*dreal(div(1))
          else
             rout = a(i)*dreal(div(i-1)) + b(i)*dreal(div(i)) + c(i)*dreal(div(i+1))
          end if
          Write(*,'(4(2x,g17.8))')i, rout, dv(i)
       end do
       STOP
    end if


    do i = 1, nx
       if(i.eq.1)then
          rout = b(i)*auxz(i) + c(i)*auxz(i+1)
       else if(i.eq.mx1)then
          rout = a(i)*auxz(i-1) + b(i)*auxz(i)
       else
          rout = a(i)*auxz(i-1) + b(i)*auxz(i) + c(i)*auxz(i+1)
       end if
       Write(*,'(4(2x,g17.8))') rout, auxu(i)
    end do
    PRINT*,'j = ', j, 'k = ', k
    STOP

  end subroutine divcorr
  
end module poisson








