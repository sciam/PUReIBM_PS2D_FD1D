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

!GLOBAL DATA FOR SCALAR MODULES
!CURRENTLY SCALAR ROUTINES USE A LOT OF ARRAYS FROM THE HYDRODYNAMIC SIDE TO SAVE MEMORY.
!BE CAREFUL!
!AUTHOR: RAHUL GARG 
MODULE scalar_data
#include "../FLO/ibm.h"

  USE precision  
  USE constants 
  USE global_data,  ffphi=>ff
  IMPLICIT NONE 
  SAVE 
  
  INTEGER, PARAMETER :: nscalmax = 1 !maximum number of scalars possible
  REAL(PRCN), dimension(:), ALLOCATABLE :: sourcesink
  COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE ::  phif,nlphif,onlphif,sourcephif
  !COMPLEX(prcn), DIMENSION(:,:,:,:), ALLOCATABLE ::  ffphi
  COMPLEX(prcn), DIMENSION(:,:,:), ALLOCATABLE :: phiin,phiout
  !       Scalar real data 
  !REAL(prcn) :: phireal(mx, my, mz,nspmx)
  REAL(prcn), DIMENSION(:,:,:,:), ALLOCATABLE :: sourcephi
  !REAL(prcn), DIMENSION(:,:,:,:,:), ALLOCATABLE :: gradphi
  !REAL(prcn) :: phi_an(mx,my,mz,nspmx)
  REAL(prcn), DIMENSION(:,:,:), ALLOCATABLE :: surf_scal_value, Nu3, nloutletplane
  REAL(prcn), DIMENSION(:,:), ALLOCATABLE :: phisurfall,&
       & phisurfall_nm1,phisurfall_nm2, flux_body, flux_body2&
       &,bulk_scalar
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: phirmean,fphirmean&
       &,sum_flux_nm1,sum_flux_nm2,sum_flux, GAMMA, phi_fluid_mean,&
       & phimodmean, phi_solid_mean
  
  REAL(PRCN), DIMENSION(:) , ALLOCATABLE :: nu_error_array,&
       & flux_global, flux_global2
  
  
  REAL(prcn) xis,xf,yi,yf,zi,zf,cons,uo,vo, nu_old, nu_error_hist,&
       & nu_error, phimean_des
  
  REAL(prcn), DIMENSION(:) , ALLOCATABLE :: heat_ratio, heat_ratio_new, um

  REAL(prcn), DIMENSION(:) , ALLOCATABLE :: heat_ratio_old

  REAL(prcn), DIMENSION(:,:) , ALLOCATABLE :: phim
  INTEGER, DIMENSION(:) , ALLOCATABLE :: countfl_plane
  
  INTEGER ::  nxi,nxf,nyi,nyf,nzi,nzf,diff,nspmx, iter_scal
  REAL(prcn) :: phisurf,phistream,  Pr_or_Sc_nu, tscal!,gamma(nspmx)
  LOGICAL :: sourcepresent, LUMPED_CAP, zero_flow_in_solid,&
       & dorpr_scal, setphimean
  
  
  !Integer :: inn_itn_step

  NAMELIST/scal_propt/ nspmx,phistream, setphimean, phisurf&
       &,sourcepresent, LUMPED_CAP, Pr_or_Sc_nu, zero_flow_in_solid,&
       & dorpr_scal

END MODULE scalar_data

