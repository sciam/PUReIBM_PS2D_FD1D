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

MODULE nl_allflow
!
!Author: Rahul Garg
#include "ibm.h" 
  USE precision 
  USE constants 
  USE global_data
  USE fftw_interface
  Use nlarrays

!  use init_turb

  USE bcsetarrays, fdderiv=>fr  
  
  Use nlmainarrays, Only : ubcp
  IMPLICIT NONE 
  Private
!!$  COMPLEX(prcn) ::  uf1(my2,mz),uf2(my2,mz),uf3(my2,mz)
!!$  COMPLEX(prcn) ::  uf11(my2,mz),uf22(my2,mz),uf33(my2,mz)
!!$  COMPLEX(prcn) ::  uf12(my2,mz),uf13(my2,mz),uf23(my2,mz)
!!$  COMPLEX(prcn) ::  tmp
!!$  REAL(prcn) :: ur1(my,mz),ur2(my,mz),ur3(my,mz)
!!$  REAL(prcn) :: ur11(my,mz),ur22(my,mz),ur33(my,mz)
!!$  REAL(prcn) :: ur12(my,mz),ur13(my,mz),ur23(my,mz)
  INTEGER :: mx2,  mxtmp
  Public :: form_nl
  REAL(PRCN) ::  u_max, v_max, w_max
  
  REAL(PRCN) ::  umax_loc, vmax_loc, wmax_loc  

  REAL(PRCN) ::  umean_tmp(ndim)

  
!  REAL(prcn) :: fplus, fminus
CONTAINS 
  SUBROUTINE form_nl!(ubc,pbc)

    IMPLICIT NONE 
    INTEGER :: i,j,k, im1, ip1, idim,im2,ip2,count, iphs, partstart, partend
    REAL(prcn) :: max_speed, Dplus_vj, Dminus_vj, phi_ri, slope_factor, mesh_veltemp(ndim)
#if !PARALLEL
    REAL(PRCN) :: uplus(mx1), uminus(mx1), U_plus(mx1), U_minus(mx1), slope(mx1)
#else
    REAL(PRCN), ALLOCATABLE, DIMENSION(:) :: uplus, uminus, U_plus, U_minus, slope
#endif    

    slope_factor = one !nl_slope_factor
    phi_ri = one
    if(I_AM_NODE_ZERO)WRITE(*,'(A,2x,g17.8)')'SLOPE FACTOR IN NL = ', slope_factor 

    count = 0
    mesh_veltemp = mesh_vel
    if(.not.movingcv) mesh_vel = zero
    if(I_AM_NODE_ZERO)WRITE(*,'(A,3(2x,g17.8))')'MESH VEL IN NL = ', mesh_vel(:) 

#if !PARALLEL          
    do i = 1,nx
#else
       do i = 0,nx+1
#endif
          do k = 1, mz
             do j = 1, my 
                if(discr_scheme.eq."center")then
                   ur11(j,k)  = (ubcp(i,j,k,1)-mesh_vel(1))*(ubcp(i,j,k,1)-mesh_vel(1))
                end if
                ur22(j,k)  = (ubcp(i,j,k,2)-mesh_vel(2))*(ubcp(i,j,k,2)-mesh_vel(2))
                ur33(j,k)  = (ubcp(i,j,k,3)-mesh_vel(3))*(ubcp(i,j,k,3)-mesh_vel(3))
                ur12(j,k)  = (ubcp(i,j,k,1)-mesh_vel(1))*(ubcp(i,j,k,2)-mesh_vel(2))
                ur13(j,k)  = (ubcp(i,j,k,1)-mesh_vel(1))*(ubcp(i,j,k,3)-mesh_vel(3))
                ur23(j,k)  = (ubcp(i,j,k,2)-mesh_vel(2))*(ubcp(i,j,k,3)-mesh_vel(3))
             end do
          end do
          CALL ff2rc(ur22,uf22)
          CALL ff2rc(ur33,uf33)
          CALL ff2rc(ur12,uf12)
          CALL ff2rc(ur13,uf13)
          CALL ff2rc(ur23,uf23)
          
          do k  = 1, mz 
             do j = 1, my
                if((j.le.my2).and.(i.gt.0).and.(i.le.nx)) then 
                   nl(i,j,k,1) = wy(j)*uf12(j,k) + wz(k)*uf13(j,k)
                   nl(i,j,k,2) = wy(j)*uf22(j,k) + wz(k)*uf23(j,k)
                   nl(i,j,k,3) = wy(j)*uf23(j,k) + wz(k)*uf33(j,k)
                end if
                if(discr_scheme.eq."center")then
                   fdderiv(i,j,k,1) = ur11(j,k)
                   fdderiv(i,j,k,2) = ur12(j,k)
                   fdderiv(i,j,k,3) = ur13(j,k)
                else
                   fdderiv(i,j,k,1:3)= zero
                end if

             end do
          end do

       end do

#if !PARALLEL
       if(discr_scheme.eq."center") fdderiv(mx,:,:,:) = fdderiv(1,:,:,:)
#endif
       !PRINT*,'discretization scheme = ', discr_scheme

       if(discr_scheme.eq."upwind")then
!!$                         xi-1/2    xi+1/2          
!!$                     ......|..........|.......
!!$                     |          |            |
!!$                    xi-1        xi          xi+1   

!!$                   We need to compute (del f/del x)|i, where f = uU, U = {u,v,w}

!!$                      (del f/del x)|i = (f(i+1/2) - f(i-1/2))/dx
!!$ These are like finite volume cells centered around {xi}.
!!$ We cannot compute f(i+1/2) and f(i-1/2) directly since the polynomials are
!!$ discontinuous at these interfaces. 
!!$   For a second order accurate scheme, one must consider piecewise
!!$   linear polynomials.

#if PARALLEL
          if(.not.ALLOCATED(slope))then
             !   BARRIER(decomp_group)
             ALLOCATE(uplus(0:nx), uminus(0:nx), U_plus(0:nx), U_minus(0:nx), slope(0:nx+1))
          end if
#endif

          do k = 1, mz
             do j = 1, my
                do idim = 1,ndim
#if PARALLEL
                   slope(0) = (ubcp(1,j,k,idim)-ubcp(-1,j,k,idim))/(two*dx) 
                   slope(nx+1) = (ubcp(nx+2,j,k,idim)-ubcp(nx,j,k,idim))/(two*dx) 
                   slope(0) = slope(0)*phi_ri*slope_factor
                   slope(nx+1) = slope(nx+1)*phi_ri*slope_factor
#endif
                   do i = 1,nx !mx1
                      im1 = i-1
                      ip1 = i+1
#if !PARALLEL
                      if(im1.lt.1) im1 = mxf+im1-1
                      if(ip1.gt.mxf-1) ip1 = ip1-(mxf-1)
#endif                   
!!$ We have information at the grid locations {xi}. Consider piecewise linear polynomials in the intervals 
!!$   [xi-1/2,xi+1/2].

!!$ Form of the piecewise polynomial:
!!$   P[U](i) = U(i) + Sj*(x-x(i)); Sj is the slope of the  piecewise polynomial. Now how to choose the slope of the polynomial? 
!!$ First Order: Sj = 0

!!$ Second Order: 

                      slope(i) = (ubcp(ip1,j,k,idim)-ubcp(im1,j,k,idim))/(two*dx) 
                      if(limiter.eq."none")then

!!$                Sj = {U(i+1)-U(i-1)}/2dx --> second order slope 
                         phi_ri = one

                      else if(limiter.eq."minmod")then

!!$              Sj = min[(U(i+1)-U(i))/dx, (U(i)-U(i-1))/dx] --> minmod limiter (for very high Re)
#if 0
                         Dplus_vj = (ubcp(ip1,j,k,idim)-ubcp(i,j,k,idim))/dx
                         Dminus_vj = (ubcp(i,j,k,idim)-ubcp(im1,j,k,idim))/dx
                         slope(i) = half*(sgn(Dplus_vj) + sgn(Dminus_vj))
                         slope(i) = slope(i) * MIN(ABS(Dplus_vj),ABS(Dminus_vj))
#endif
                         phi_ri = (ubcp(i,j,k,idim) - ubcp(im1,j,k,idim))/(ubcp(ip1,j,k,idim) - ubcp(i,j,k,idim)+SMALL_NUMBER)
                         if(phi_ri.le.zero)then
                            phi_ri = zero
                         else
                            phi_ri = MIN(one,phi_ri)
                         end if

                      end if
                      slope(i) = slope(i)*phi_ri*slope_factor

!!$ Let us store the fluxes at the interfaces similar to the way we store pressure. That is, flux at the interface (i+1/2) is
!!$ stored at index i.

                      if(idim.eq.1)then
                         
!!$                From the cell i, one can get the following information at the interfaces (i-1/2) and (i+1/2)
!!$                (i-1/2) --> uplus stored at index i-1(im1)
!!$                (i+1/2) --> uminus stored at i
!!$                Uplus(xi+1/2) = P[U](i+1)|(xi+1/2) and
                         
!!$                Uminus(xi+1/2) = P[U](i)|(xi+1/2)

                         uplus(im1) = ubcp(i,j,k,idim) - slope(i)*(dx/two) ! uplus(xi-1/2)
                         uminus(i) =  ubcp(i,j,k,idim) + slope(i)*(dx/two) ! uminus(xi+1/2)

                         U_plus(im1) = uplus(im1)
                         U_minus(i) = uminus(i)
                      else
                         U_plus(im1) = ubcp(i,j,k,idim) - slope(i)*(dx/two)  ! Uplus(i-1/2)
                         U_minus(i) =  ubcp(i,j,k,idim) + slope(i)*(dx/two) ! Uminus(i+1/2)
                      end if
!!$                      if(i.eq.1)then
!!$                         RSENDRECV(slope(i),1,fromproc,1,slope(nx+1),1,toproc,1,decomp_group,status) 
!!$                         
!!$                      else if(i.eq.nx)then 
!!$                         RSENDRECV(slope(i),1,toproc,0,slope(0),1,fromproc,0,decomp_group,status)                    
!!$                      end if

!!$ Evaluate the flux at an interface using the formula :
!!$           f(i+1/2) = 1/2[fplus(i+1/2) + fminus(i+1/2)] -1/2[ max(|uminus|,|uplus|)|(i+1/2){Uplus(i+1/2)-Uminus(i+1/2)}]

!!$           f(i-1/2) = 1/2[fplus(i-1/2) + fminus(i-1/2)] -1/2[ max(|uminus|,|uplus|)|(i-1/2){Uplus(i+1/2)-Uminus(i+1/2)}]
                      
!!$ From the cell i, we can compute the following information:
!!$ At (i-1/2) --> fplus(i-1/2) stored at fdderiv(i-1) or fdderiv(im1)
!!$ At (i+1/2) --> fminus(i+1/2)stored at fdderiv(i) or fdderiv(im1)

                      fdderiv(im1,j,k,idim) = fdderiv(im1,j,k,idim) + half*(uplus(im1)-mesh_vel(1))*(U_plus(im1)-mesh_vel(idim))
                      fdderiv(i,j,k,idim) = fdderiv(i,j,k,idim) + half*(uminus(i)-mesh_vel(1))*(U_minus(i)-mesh_vel(idim))

                      if(i.ne.1)then
!!$ Except for cell 1 (due to PBC), if I am at cell index i, then I am ensured that I now have the complete information 
!!$ for interface (i-1/2).
                         max_speed = MAX(ABS(uminus(im1)-mesh_vel(1)), ABS(uplus(im1)-mesh_vel(1)))
                         fdderiv(im1,j,k,idim) = fdderiv(im1,j,k,idim) - half*max_speed*(U_plus(im1)-U_minus(im1))
                      end if

                   end do

#if !PARALLEL
	           max_speed = MAX(ABS(uminus(mx1)-mesh_vel(1)), ABS(uplus(mx1)-mesh_vel(1)))
                   !max_speed = MAX(ABS(uminus(mx1)), ABS(uplus(mx1)))
                   fdderiv(mx1,j,k,idim) = fdderiv(mx1,j,k,idim) - half*max_speed*(U_plus(mx1)-U_minus(mx1))
                   fdderiv(mx,j,k,idim)  = fdderiv(1,j,k,idim)
#else
                   if(idim.eq.1)then
                      uminus(0) = ubcp(0,j,k,idim) + slope(0)*(dx/two)
                      U_minus(0) = uminus(0)

                      uplus(nx) = ubcp(nx+1,j,k,idim) - slope(nx+1)*(dx/two)
                      U_plus(nx) = uplus(nx)
                   else
                      U_minus(0) = ubcp(0,j,k,idim) + slope(0)*(dx/two)
                      U_plus(nx) = ubcp(nx+1,j,k,idim) - slope(nx+1)*(dx/two)
                   endif
                   fdderiv(0,j,k,idim) = fdderiv(0,j,k,idim) + half*(uminus(0)-mesh_vel(1))*(U_minus(0)-mesh_vel(idim))
                   max_speed = MAX(ABS(uminus(0)-mesh_vel(1)), ABS(uplus(0)-mesh_vel(1)))
                   fdderiv(0,j,k,idim) = fdderiv(0,j,k,idim) - half*max_speed*(U_plus(0)-U_minus(0))
                   fdderiv(nx,j,k,idim) = fdderiv(nx,j,k,idim) + half*uplus(nx)*U_plus(nx)
                   max_speed = MAX(ABS(uminus(nx)-mesh_vel(1)), ABS(uplus(nx)-mesh_vel(1)))
                   fdderiv(nx,j,k,idim) = fdderiv(nx,j,k,idim) - half*max_speed*(U_plus(nx)-U_minus(nx))
#endif
                end do

             end do


          end do

!!$#if PARALLEL
!!$          DEALLOCATE(uplus, uminus, U_plus, U_minus, slope)
!!$#endif

       end if

       DO idim = 1, ndim
          DO i = 1, nx !mx1
             do k = 1, mz
                do j = 1, my 
                   im1 = i-1
                   ip1 = i+1
#if !PARALLEL
                   if(im1.lt.1) im1 = mxf+im1-1
                   if(ip1.gt.mxf-1) ip1 = ip1-(mxf-1)
#endif
                   IF(discr_scheme.eq."center")then
                      ur11(j,k) = (fdderiv(ip1,j,k,idim) - fdderiv(im1,j,k,idim))/(two*dx)
!!$                
                   ELSE IF(discr_scheme.eq."upwind")then
                      ur11(j,k) = (fdderiv(i,j,k,idim) - fdderiv(im1,j,k,idim))/dx

                   END IF
                end do
             end do
             call ff2rc(ur11, uf11)
             
             nl(i,:,:,idim) = nl(i,:,:,idim) + uf11(:,:)
             if(i.eq.1)then
                VECSENDRECV(nl(i,1,1,idim),1,twodcslice,fromproc,1,nl(nx+1,1,1,idim),1,toproc,1,decomp_group,status)
             end if
          end DO
#if !PARALLEL
          nl(nx+1, :,:,idim) = nl(1,:,:,idim)
#endif
       end DO
       mesh_vel = mesh_veltemp
     END SUBROUTINE form_nl
     
END MODULE nl_allflow

