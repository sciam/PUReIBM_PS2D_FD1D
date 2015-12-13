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

module steptheflow
#include "ibm.h"
  !///////////////////////////////////////////////////////////////////////
  !	calculate velocity field at next time level
  !	pressure estimation and divergence correction done in place
  !	(for each wavenumber) to save storage
  !-----------------------------------------------------------------------
  USE precision 
  USE constants 
  Use poisson, Only : pressure, divcorr
  use usteptime 
  Use nlmainarrays, Only : pr => pbcp
  Use nlarrays, Only : uf1,ur1
  Use dependent_functions
  USE dem_mod
Contains
  subroutine velstep(sflag,rks)
    
    Use global_data, Only : my, mz, my2, xend,xstart
    
    SAVE
    
    integer :: j,k, i
    integer, Intent(in) :: sflag,rks
    complex(prcn) :: usumloc, usum
  
    divloc = czero
    total_div = czero


    do j = 1, my2
       do k = 1, mz
          CALL pressure(rks,j,k)
       end do
    end do
    CALL communicate_pressure
    
    do j = 1, my2
       do k = 1, mz
          CALL ustep(rks,j,k)
       end do
    end do
    
    CALL communicate_velocity
    
    do j = 1, my2
       do k = 1, mz
          CALL divcorr(rks,j,k)
       end do
    end do
    CALL correct_bnd_divergence(rks)
    CALL communicate_pressure
    CALL communicate_velocity

#if PARALLEL
    if(debug_check)CALL check_divergence 
#else
    CALL check_divergence 
#endif
    
    if(I_AM_NODE_ZERO) WRITE(*,*) 'EVOLVING UMEAN : '

! COMMENTING THIS, BECAUSE GRAVITY IS CHANGED TO MPG. IF BOTH GRAVITY AND MPG ARE AVAILABLE, THEN UNCOMMENT THIS
    umean(:) = umean(:) + (-mpg(:) + frmean(:) - frame_accln(:))*dt !+ grav(:) 

    if(move_particles)then
       frame_vel(:) = frame_vel(:) + frame_accln(:)*dt
       frame_pos(:) = frame_pos(:) + frame_vel(:)*dt
    end if
    
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A25,3(2x,g17.8))')'FR MEAN @ N  = ', (FRMEAN(j), j = 1, ndim)
       WRITE(*,'(A25,3(2x,g17.8))')'MPG @ N  = ', (MPG(j), j = 1, ndim)
       if(move_particles)then
          WRITE(*,'(A25,3(2x,g17.8))')'FRAME ACCLN @ N  = ', (frame_accln(j), j = 1, ndim)
          WRITE(*,'(A25,3(2x,g17.8))')'FRAME VEL @ N  = ', (frame_vel(j), j = 1, ndim)
          WRITE(*,'(A25,3(2x,g17.8))')'FRAME POS @ N  = ', (frame_pos(j), j = 1, ndim)
       end if
       WRITE(*,'(A25,3(2x,g17.8))')'UMEAN @ n + 1= ', (UMEAN(j), j = 1, ndim)
    end if


#if 0
!!$    CALL calc_pressure(pr(1:nx,1:my,1:mz))
!!$
    if(I_AM_NODE_ZERO)then
       OPEN(1000,FILE=TRIM(RUN_NAME)//'_pressure.dat',status='unknown')
       write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "P" '
       write(1000,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
       !    do i = 1,nx
       !i = 1
       do k=1,mz
          do j=1,my2
!!$             uf1(j,k) = p(i,j,k)
!!$          end do
!!$       end do
!!$       CALL ff2cr(uf1(:,:),ur1(:,:))
!!$       pr(i,:,:) = ur1(:,:)
!!$    end do
!!$    do k=1,mz
!!$       do j=1,my
!!$          !          do i=1,nx !mx
#if !PARALLEL
          do i = 1, nx/2
#else
          do i = 1,nx
#endif
             write(1000,*)(GLOBAL_INDEX(i)),(j),(k),DREAL(p(i,j,k)),DIMAG(p(i,j,k))!,pr(i,j,k)
!!$             !          enddo
         enddo
       enddo
    end do
    close(1000,status='keep')
endif
!!$
PARALLEL_FINISH()
STOP
!    GLOBAL_COMPLEX_SUM(divloc,total_div,1,decomp_group)
!    if(I_AM_NODE_ZERO)  Write(*,'(A25,2(2x,g12.5))') 'MEAN DIV = ', total_div/(mx1*my*mz)
    !    PARALLEL_FINISH()
 !   STOP
#endif


end subroutine velstep

  SUBROUTINE communicate_pressure
    IMPLICIT NONE
    INTEGER :: i
#if PARALLEL
    i = 1
    VECSENDRECV(p(i,1,1),1,ucslice,fromproc,1,p(nx+1,1,1),1,toproc,1,decomp_group,status)
    i = nx
    VECSENDRECV(p(i,1,1),1,ucslice,toproc,0,p(0,1,1),1,fromproc,0,decomp_group,status)
#endif
    RETURN
  END SUBROUTINE communicate_pressure
  
 SUBROUTINE communicate_velocity
   IMPLICIT NONE
    INTEGER :: idim, i
#if PARALLEL
    do idim = 1, ndim
       i = 1
       VECSENDRECV(u(i,1,1,idim),1,ucslice,fromproc,1,u(nx+1,1,1,idim),1,toproc,1,decomp_group,status)
       i = nx
       VECSENDRECV(u(i,1,1,idim),1,ucslice,toproc,0,u(0,1,1,idim),1,fromproc,0,decomp_group,status)
    end do
#endif
    RETURN
  END SUBROUTINE communicate_velocity

  SUBROUTINE correct_bnd_divergence(rks)
    USE nlarrays, ONLY : westsendbuf => uf1, westrecvbuf => uf11, eastsendbuf => uf2, eastrecvbuf => uf22
    IMPLICIT NONE
    INTEGER, Intent(in) :: rks
    INTEGER :: idim, i, j, k
    REAL(prcn) :: alp_coef_dt 
#if PARALLEL
    CSENDRECV(westsendbuf(1,1),my2*mz,fromproc,1,eastrecvbuf(1,1),my2*mz,toproc,1,decomp_group,status)
    CSENDRECV(eastsendbuf(1,1),my2*mz,toproc,1,westrecvbuf(1,1),my2*mz,fromproc,1,decomp_group,status)
    
    alp_coef_dt = ((coef(rks,1)+coef(rks,2))*dt)
    do j = 1, my2
       do k = 1, mz
          p(0,j,k) = westrecvbuf(j,k)
          p(nx+1,j,k) = eastrecvbuf(j,k)
          
          i = 1
          u(i,j,k,1)=u(i,j,k,1)-alp_coef_dt*(-p(i-1,j,k))/dx
          u(i,j,k,2)=u(i,j,k,2)-alp_coef_dt*half*wy(j)*(p(i-1,j,k))
          u(i,j,k,3)=u(i,j,k,3)-alp_coef_dt*half*wz(k)*(p(i-1,j,k))
          
          p(i,j,k)=p(i,j,k) - (coef(rks,2)*dt)*vis*(p(i-1,j,k)/dx2 - w2(j,k)*p(i-1,j,k))

          i = nx
          p(i,j,k)=p(i,j,k) - (coef(rks,2)*dt)*vis*(p(i+1,j,k)/dx2 - w2(j,k)*p(i+1,j,k))
       end do
    end do
#endif
    RETURN
  END SUBROUTINE correct_bnd_divergence
  
  SUBROUTINE check_divergence
    IMPLICIT NONE
    INTEGER :: i, j, k
    COMPLEX(prcn) :: tmpc
    REAL(prcn) :: tmp, divlocmax
    
    divlocmax = zero
    
    do k = 1, mz
       do j = 1, my2
          do i = 1, nx
             tmpc = (u(i+1,j,k,1)-u(i,j,k,1))/dx
             tmpc=tmpc+half*wy(j)*(u(i,j,k,2)+u(i+1,j,k,2))
             tmpc=tmpc+half*wz(k)*(u(i,j,k,3)+u(i+1,j,k,3))
             tmp = DSQRT(dble(tmpc*conjg(tmpc)))
             if(tmp.gt.divlocmax)divlocmax = tmp
             end do 
       end do
    end do
    GLOBAL_DOUBLE_MAX(divlocmax, divmax, 1,decomp_group)
    RETURN
  END SUBROUTINE check_divergence
end module steptheflow
