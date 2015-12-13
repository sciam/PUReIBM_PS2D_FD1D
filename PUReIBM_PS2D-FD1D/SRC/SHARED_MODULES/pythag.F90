
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
FUNCTION pythag_sp(a,b)
  USE nrtype
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b
  REAL(SP) :: pythag_sp
  REAL(SP) :: absa,absb
  absa=abs(a)
  absb=abs(b)
  if (absa > absb) then
     pythag_sp=absa*sqrt(1.0_sp+(absb/absa)**2)
  else
     if (absb == 0.0) then
        pythag_sp=0.0
     else
        pythag_sp=absb*sqrt(1.0_sp+(absa/absb)**2)
     end if
  end if
END FUNCTION pythag_sp

  FUNCTION pythag_dp(a,b)
    USE nrtype
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: a,b
    REAL(DP) :: pythag_dp
    REAL(DP) :: absa,absb
    absa=abs(a)
    absb=abs(b)
    if (absa > absb) then
       pythag_dp=absa*sqrt(1.0_dp+(absb/absa)**2)
    else
       if (absb == 0.0) then
          pythag_dp=0.0
       else
          pythag_dp=absb*sqrt(1.0_dp+(absa/absb)**2)
       end if
    end if
  END FUNCTION pythag_dp
