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

FUNCTION erfcc_s(x)
  USE nrtype
  USE nrutil, ONLY : poly
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: erfcc_s
  REAL(DP) :: t,z
  REAL(DP), DIMENSION(10) :: coef = (/-1.26551223_dp,1.00002368_dp,&
       0.37409196_dp,0.09678418_dp,-0.18628806_dp,0.27886807_dp,&
       -1.13520398_dp,1.48851587_dp,-0.82215223_dp,0.17087277_dp/)
  z=abs(x)
  t=1.0_dp/(1.0_dp+0.5_dp*z)
  erfcc_s=t*exp(-z*z+poly(t,coef))
  if (x < 0.0) erfcc_s=2.0_dp-erfcc_s
END FUNCTION erfcc_s
   
   
FUNCTION erfcc_v(x)
  USE nrtype
  USE nrutil, ONLY : poly
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(size(x)) :: erfcc_v,t,z
  REAL(DP), DIMENSION(10) :: coef = (/-1.26551223_dp,1.00002368_dp,&
       0.37409196_dp,0.09678418_dp,-0.18628806_dp,0.27886807_dp,&
       -1.13520398_dp,1.48851587_dp,-0.82215223_dp,0.17087277_dp/)
  z=abs(x)
  t=1.0_dp/(1.0_dp+0.5_dp*z)
  erfcc_v=t*exp(-z*z+poly(t,coef))
  where (x < 0.0) erfcc_v=2.0_dp-erfcc_v
END FUNCTION erfcc_v
