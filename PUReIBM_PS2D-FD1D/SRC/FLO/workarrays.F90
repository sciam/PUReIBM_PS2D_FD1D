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

Module nlarrays 

  Use precision 
  Use constants 
  Use global_data
  implicit none
  !Save 
  Public
  COMPLEX(prcn), DIMENSION(:,:), ALLOCATABLE :: uf1, uf2, uf3, uf11, uf21,uf22, uf33,uf31,uf32
  COMPLEX(prcn), DIMENSION(:,:), ALLOCATABLE :: uf12, uf13, uf23

  REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: ur1,ur2,ur3
#if PARALLEL
  REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: uatminus1, uatnxp2
#endif
  REAL(prcn), DIMENSION(:,:), ALLOCATABLE ::  ur11,ur12, ur13,ur21, ur22,ur23, ur31, ur32,  ur33

end Module nlarrays
!!$
Module bcsetarrays
  Use precision 
  Use constants 
  Use global_data
  implicit none
  !Save 
  
  !Real(prcn), dimension(mxf,my,mz), Target :: dudy,dwdy
  
  Real(prcn),  DIMENSION(:,:,:,:), ALLOCATABLE, target :: omega, fr, ppr, diffn
end Module bcsetarrays

Module nlmainarrays
  Use precision 
  Use constants 
  Use global_data
  implicit none
  !Save 
  Public

  
  real(prcn), DIMENSION(:,:,:,:), ALLOCATABLE,  Target ::  ubc, nlbc, onlbc
  real(prcn), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pbc
  
  real(prcn), Dimension(:,:,:,:), pointer ::  onlbcp,nlbcp, ubcp
  real(prcn), Dimension(:,:,:), pointer ::  pbcp
end Module nlmainarrays

module epsarray
	real(8), allocatable, dimension(:,:,:,:,:) :: ugrad
	real(8), allocatable, dimension(:,:,:)   :: epsij_th_loc,epsij_th, ug_th_loc,ug_th
	real(8), allocatable, dimension(:,:)       :: ugp_th_loc,ugp_th, epsij_surf_loc, epsij_surf, epsij_loc, epsij, epsij_far,epsij_far_loc
	real(8), allocatable, dimension(:)         :: eps_th

	integer, allocatable, dimension(:) :: thbin_loc, thbin

	real(8) :: eps, eps_surf, eps_far
	integer :: nthmax
	

	logical, save :: from_mypost = .false.
end module epsarray



