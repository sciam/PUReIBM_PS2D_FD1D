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

FUNCTION gammln_s(xx)
	USE nrtype
        USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xx
	REAL(SP) :: gammln_s
	REAL(DP) :: tmp,x
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	call assert(xx > 0.0, 'gammln_s arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	gammln_s=tmp+log(stp*(1.000000000190015_dp+&
		sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
	END FUNCTION gammln_s


	FUNCTION gammln_v(xx)
	USE nrtype
        USE nrutil, ONLY: assert
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), DIMENSION(size(xx)) :: gammln_v
	REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	if (size(xx) == 0) RETURN
	call assert(all(xx > 0.0), 'gammln_v arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	ser=1.000000000190015_dp
	y=x
	do i=1,size(coef)
		y=y+1.0_dp
		ser=ser+coef(i)/y
	end do
	gammln_v=tmp+log(stp*ser/x)
	END FUNCTION gammln_v
