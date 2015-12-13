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

	FUNCTION gcf_s(a,x,gln)
	USE nrtype
        USE nrutil, ONLY : nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP) :: gcf_s
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: i
	REAL(SP) :: an,b,c,d,del,h
	if (x == 0.0) then
		gcf_s=1.0
		RETURN
	end if
	b=x+1.0_sp-a
	c=1.0_sp/FPMIN
	d=1.0_sp/b
	h=d
	do i=1,ITMAX
		an=-i*(i-a)
		b=b+2.0_sp
		d=an*d+b
		if (abs(d) < FPMIN) d=FPMIN
		c=b+an/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_sp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_sp) <= EPS) exit
	end do
	if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
	if (present(gln)) then
		gln=gammln(a)
		gcf_s=exp(-x+a*log(x)-gln)*h
	else
		gcf_s=exp(-x+a*log(x)-gammln(a))*h
	end if
	END FUNCTION gcf_s


	FUNCTION gcf_v(a,x,gln)
          USE nrtype
          USE nrutil, ONLY : assert_eq,nrerror
          USE nr, ONLY : gammln
          IMPLICIT NONE
          REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
          REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
          REAL(SP), DIMENSION(size(a)) :: gcf_v
          INTEGER(I4B), PARAMETER :: ITMAX=100
          REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
          INTEGER(I4B) :: i
          REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
          LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
          i=assert_eq(size(a),size(x),'gcf_v')
          zero=(x == 0.0)
          where (zero)
             gcf_v=1.0
          elsewhere
             b=x+1.0_sp-a
             c=1.0_sp/FPMIN
             d=1.0_sp/b
             h=d
          end where
          converged=zero
          do i=1,ITMAX
             where (.not. converged)
                an=-i*(i-a)
                b=b+2.0_sp
                d=an*d+b
                d=merge(FPMIN,d, abs(d)<FPMIN )
                c=b+an/c
                c=merge(FPMIN,c, abs(c)<FPMIN )
                d=1.0_sp/d
                del=d*c
                h=h*del
                converged = (abs(del-1.0_sp)<=EPS)
             end where
             if (all(converged)) exit
          end do
          if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
          if (present(gln)) then
             if (size(gln) < size(a)) call &
                  nrerror('gser: Not enough space for gln')
             gln=gammln(a)
             where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
          else
             where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
          end if
	END FUNCTION gcf_v
