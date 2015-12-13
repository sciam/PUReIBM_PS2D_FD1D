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

Module nlcalc 
#include "ibm.h"
  USE precision 
  USE constants 
  Use nl_allflow
  Use nl_allflow_noncons
  Use boundary_condition
  Use dependent_functions
  Use collision_mod
  implicit none 
  !-----------------------------------------------------------------------
  !	calculate the nonlinear (convective) terms in the NS equations
  !	the velocities are as follows;
  !
  !	u1	u
  !	u2	v
  !	u11	uu
  !	u12	uv
  !	u13	uw
  !	u22	vv
  !	u23	vw
  !	u33	ww	
  !	
  !	use phase-shifting to dealias
  !-----------------------------------------------------------------------
contains   

  subroutine nonlinear(rks)
    Use global_data
    Use fftw_interface
    Use nlmainarrays, Only :  ubc, nlbc, onlbc
    !Use nlarrays, Only : onlf=>uf1
    implicit none 
    integer, Intent(in) ::  rks

    !-----------------------------------------------------------------------
    !	local variables



    !real(prcn) ::  nlr(my,mz)
    !real(prcn) ::  onlr(my,mz)
    real(prcn) ::  eold, tvis, nlmean(ndim), nlmean_loc(ndim) 

    !complex(prcn) ::  onlf(my2,mz)
    complex(prcn) :: nltemp1

    integer i,j,k,l,ii,n,nlflag, unit1
    REAL(prcn) :: utemp, tempt, dt_fluid_orig

	real(prcn) cpu0, cpu1

    !	If nlflag is 1 then use UdotGradU
    !       If nlflag is 2 then use Del.UU

	if (I_AM_NODE_ZERO) CALL CPU_TIME (CPU0)

    nlflag = 1
    !-----------------------------------------------------------------------
    !	reset terms to be overwritten


    !write(*,*)'Beginning of nonlinear inside nl..and idumstep = ', idumstep
    do n=1,ndim
       do i = 1, nx+1
          do j = 1, my2
             do k = 1, mz
                onl(i,j,k,n) = nl(i,j,k,n)
                nl(i,j,k,n) = czero
             end do
          end do
       end do
    end do
    
    if(I_AM_NODE_ZERO)write(*,'(A50)')'CONSERVATIVE: performing del dot uu'
    
    if(.not.only_dem) call form_nl!(ubc,pbc)
    
#if 0        
    t=t+dt*(coef(rks,1)+coef(rks,2))
    if(move_particles)S_TIME = t
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
#endif    
    do n=1,ndim
       
       !Jamals convention of absorbing the negative
       ! sign into the nl term 
       
       if(aliasflag.eq.1)then
          nl(1:nx+1,:,:,n) = -thrd*nl(1:nx+1,:,:,n) !Exposing
          ! the last stride 
       else
          
          nl(1:nx+1,:,:,n) = -nl(1:nx+1,:,:,n)
       endif
    enddo
    
    !-----------------------------------------------------------------------
    !	store nonlinear terms in forcing domain
    
    do n=1,ndim
#if PARALLEL       
       nlbc(0,:,:,n) = zero
       onlbc(0,:,:,n) = zero
#endif
       do i=1,nx+1
          call ff2cr(nl(i,:,:,n), nlbc(i,:,:,n))
          call ff2cr(onl(i,:,:,n),onlbc(i,:,:,n))
#if PARALLEL       
          if(i.eq.nx+1)then
             nlbc(i,:,:,n) = zero
             onlbc(i,:,:,n) = zero
          end if
#endif
       enddo
    enddo
    
    nlmean(1:ndim) = zero
    if(debug_check)then
       do n = 1, ndim 
          nlmean_loc(n) = SUM(nlbc(1:nx,:,:,n))
          GLOBAL_DOUBLE_SUM(nlmean_loc(n),nlmean(n),1,decomp_group)
       end do
       
       nlmean(:) = nlmean(:)/(mx1*my*mz) 
       if(I_AM_NODE_ZERO) WRITE(*,'(A25,3(2x,g12.5))')'NLMEAN = ', nlmean
    end if
    
    ! New dt is generated here so as to resolve the collisions.
    ! This new dt is used to compute the immersed boundary force.
    
    !-----------------------------------------------------------------------
    !	calculate body force terms
    !-----------------------------------------------------------------------

		if (I_AM_NODE_ZERO) then
			CALL CPU_TIME (CPU1) 
			nl_time_current = cpu1-cpu0
			nl_time = nl_time + nl_time_current
		endif

    if(nbody.gt.0)then
       call bcset(rks)
    end if
    
    if(I_AM_NODE_ZERO)write(*,'(A)')'END OF NONLINEAR'
    
  end subroutine nonlinear
end Module nlcalc













