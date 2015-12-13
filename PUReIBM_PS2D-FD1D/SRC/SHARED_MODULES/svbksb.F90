
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
SUBROUTINE svbksb_sp(u,w,v,b,x)
  USE nrtype
  USE nrutil, ONLY : assert_eq
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
  REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
  REAL(SP), DIMENSION(:), INTENT(OUT) :: x
  INTEGER(I4B) :: mdum,ndum
  REAL(SP), DIMENSION(size(x)) :: tmp
  mdum=assert_eq(size(u,1),size(b),'svbksb_sp: mdum')
  ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),&
       'svbksb_sp: ndum')
  where (w /= 0.0)
     tmp=matmul(b,u)/w
  elsewhere
     tmp=0.0
  end where
  x=matmul(v,tmp)
END SUBROUTINE svbksb_sp

  SUBROUTINE svbksb_dp(u,w,v,b,x)
    USE nrtype
    USE nrutil, ONLY : assert_eq
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,v
    REAL(DP), DIMENSION(:), INTENT(IN) :: w,b
    REAL(DP), DIMENSION(:), INTENT(OUT) :: x
    INTEGER(I4B) :: mdum,ndum
    REAL(DP), DIMENSION(size(x)) :: tmp
    mdum=assert_eq(size(u,1),size(b),'svbksb_dp: mdum')
    ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),&
         'svbksb_dp: ndum')
    where (w /= 0.0)
       tmp=matmul(b,u)/w
    elsewhere
       tmp=0.0
    end where
    x=matmul(v,tmp)
  END SUBROUTINE svbksb_dp
