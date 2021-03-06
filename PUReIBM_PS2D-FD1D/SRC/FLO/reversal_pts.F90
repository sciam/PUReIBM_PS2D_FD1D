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

module reversal_pts
  USE PRECISION
  
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC:: revpt, parrayrev, revpt_position, update_cell_index
  
  INTEGER, PARAMETER, Save:: maxorder=6
  
  !object reversal point
  Type revpt
     Private
     REAL(prcn), DIMENSION(3):: position
     INTEGER, DIMENSION(3):: cell
     REAL(prcn), DIMENSION(maxorder,maxorder,maxorder) :: wts
  END TYPE revpt
  
  !-------
  ! Have to define a new type (reversal point-array pointer)
  ! in order to declare an array of reversal points-array pointers
  !-------
  TYPE parrayrev
     TYPE(revpt), DIMENSION(:), POINTER:: p
  END TYPE parrayrev
  
  
  !-------
  ! Generic interface definitions
  !-------
  INTERFACE revpt_position
     MODULE PROCEDURE revpt_position_scalar
     MODULE PROCEDURE revpt_position_vector
  END INTERFACE
  
  !-------
  !-------
  ! Access the reversal point position vector
  !-------
  !-------
Contains 
  FUNCTION revpt_position_vector(p)
    REAL(prcn), DIMENSION(3):: revpt_position_vector
    TYPE(revpt), INTENT(IN)::p
    
    revpt_position_vector = p%position
  END FUNCTION revpt_position_vector
  
  !-------
  !-------
  ! Access the individual reversal point  position
  !-------
  !-------
  FUNCTION revpt_position_scalar(p,idirection)
    REAL(prcn):: revpt_position_scalar
    TYPE(revpt), INTENT(IN)::p
    INTEGER:: idirection

    revpt_position_scalar = p%position(idirection)
  END FUNCTION revpt_position_scalar

  !subroutine to update the rev pt cell index. 
  SUBROUTINE update_cell_index(p,newcell)
    TYPE(revpt), INTENT(INOUT):: p
    INTEGER, DIMENSION(3):: newcell
    
    p%cell = newcell
  END SUBROUTINE update_cell_index
end module reversal_pts
