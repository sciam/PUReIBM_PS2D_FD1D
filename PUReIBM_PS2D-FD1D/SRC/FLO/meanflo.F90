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


MODULE meanflo
#include "ibm.h"
  USE precision 
  USE constants 
  USE global_data
  USE tridiagonal
  USE gauss_seidel
  IMPLICIT NONE 
CONTAINS

  SUBROUTINE calc_meanflo(rks)
    
    !       this is a 3d code
    
    IMPLICIT NONE 
    INTEGER, INTENT(in) ::  rks
    
    !-----------------------------------------------------------------------
    !	local variables
    
    INTEGER :: j
    
    umean(1) = umean(1) + (-mpg(1)+frmean(1))*dt
    umean(2) = umean(2) + (-mpg(2)+frmean(2))*dt
    umean(3) = umean(3) + (-mpg(3)+frmean(3))*dt
    
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A25,3(2x,g12.5))')'FR MEAN   = ', (FRMEAN(j), j = 1, ndim)
       WRITE(*,'(A25,3(2x,g12.5))')'MPG   = ', (MPG(j), j = 1, ndim)
       WRITE(*,'(A25,3(2x,g12.5))')'UMEAN = ', (UMEAN(j), j = 1, ndim)
    end if
    RETURN
  end SUBROUTINE calc_meanflo
end MODULE meanflo


    
