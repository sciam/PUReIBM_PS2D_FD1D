
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

FUNCTION gser_s(a,x,gln)
	USE nrtype
        USE nrutil, ONLY : nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP) :: gser_s
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x)
	INTEGER(I4B) :: n
	REAL(SP) :: ap,del,summ
	if (x == 0.0) then
		gser_s=0.0
		RETURN
	end if
	ap=a
	summ=1.0_sp/a
	del=summ
	do n=1,ITMAX
		ap=ap+1.0_sp
		del=del*x/ap
		summ=summ+del
		if (abs(del) < abs(summ)*EPS) exit
	end do
	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
	if (present(gln)) then
		gln=gammln(a)
		gser_s=summ*exp(-x+a*log(x)-gln)
	else
		gser_s=summ*exp(-x+a*log(x)-gammln(a))
	end if
	END FUNCTION gser_s


	FUNCTION gser_v(a,x,gln)
	USE nrtype
        USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP), DIMENSION(size(a)) :: gser_v
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x)
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(a)) :: ap,del,summ
	LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
	n=assert_eq(size(a),size(x),'gser_v')
	zero=(x == 0.0)
	where (zero) gser_v=0.0
	ap=a
	summ=1.0_sp/a
	del=summ
	converged=zero
	do n=1,ITMAX
		where (.not. converged)
			ap=ap+1.0_sp
			del=del*x/ap
			summ=summ+del
			converged = (abs(del) < abs(summ)*EPS)
		end where
		if (all(converged)) exit
	end do
	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			nrerror('gser: Not enough space for gln')
		gln=gammln(a)
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
	else
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
	end if
	END FUNCTION gser_v
