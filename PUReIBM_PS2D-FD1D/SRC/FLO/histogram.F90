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


SUBROUTINE histogram(u,wt,n,nddu,bldd,brdd,hist)
  !nddu: number of bins for pdf formation
  !bldd: lefh hand side limit .... output
  !brdd: right side limie.... otuput
  USE randomno
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n , nddu
  REAL(prcn), INTENT(in), DIMENSION(:) :: u, wt
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: u_t

  REAL(prcn),  INTENT(out) :: bldd, brdd

  REAL(prcn),  INTENT(out), DIMENSION(:) ::  hist(100)
  REAL(prcn) ::  xdiff

  REAL(prcn) :: vlmt, ave, adev, sdev, var, skew, curt

  INTEGER :: i, ibin
  bldd= 1.e25
  brdd=-1.e25
  ALLOCATE(u_t(n))
  CALL moment1(4, n, n, wt, u, ave,adev,sdev,var,skew&
       &,curt)
  PRINT*,'ave,var,sdev, skew, curt=', ave, var,sdev, skew,curt
  WRITE(*,*)'number of bins in hist..',nddu

  DO i=1,n
     u_t(i)=(u(i)-ave)/sdev
  ENDDO

  CALL moment1(4, n, n, wt, u_t, ave,adev,sdev,var,skew&
       &,curt)
  PRINT*,'ave,var, sdev,skew, curt=', ave, var,sdev, skew,curt
  DO i=1,nddu
     hist(i)=0.0  
  ENDDO
  bldd = MIN(MINVAL(u_t(:)), bldd)
  brdd = MAX(MAXVAL(u_t(:)), brdd)
!!$    DO i=1,n
!!$       bldd = amin1(u(i),bldd)
!!$       brdd = amax1(u(i),brdd)
!!$    ENDDO

  xdiff  = (brdd - bldd)/float(nddu-1)

  DO i=1,n
     ibin = (u_t(i) - bldd)/xdiff + 1
     hist(ibin)=hist(ibin) + wt(i)/xdiff
  ENDDO

  DEALLOCATE(u_t)

END SUBROUTINE histogram

 SUBROUTINE plothist(hist,lb,ub,nbins,iunit,t,tref)
   IMPLICIT NONE
   REAL(prcn), INTENT(in), DIMENSION(:) ::hist

   INTEGER, INTENT(in) ::iunit,nbins
   REAL(prcn) ::  sum_x,lb,ub,dx, t, tref, tmp, tmp2
   INTEGER(prcn) :: i

   sum_x = lb
   dx= (ub-lb)/(float(nbins)-1)
   tmp = one/sqrt(twopi)
   WRITE(iunit,*)'Zone t="',t/tref,'"'
   DO i=1,nbins
      tmp2  = sum_x+dx/2
      WRITE(iunit,*)tmp2,hist(i)!, tmp*exp(-(tmp2*tmp2)/two)
      sum_x = sum_x + dx
   ENDDO

   RETURN
 END SUBROUTINE plothist

