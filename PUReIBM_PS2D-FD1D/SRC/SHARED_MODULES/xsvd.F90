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
Module xsvd
  USE nrtype
  Implicit None 
  !private
  Public
Contains
  Subroutine xsvbksb(a,b,sol,N,M)
    !	driver for routine svbksb, which calls routine svdcmp
    USE nrtype
    USE nr
    IMPLICIT NONE
    INTEGER,Intent(in) :: M,n
    REAL(DP), DIMENSION(N,N),Intent(in) :: a
    REAL(DP), DIMENSION(N,M),Intent(in) :: b
    REAL(DP), DIMENSION(N,M),Intent(out) :: sol
    REAL(DP), DIMENSION(N,N) :: u
    INTEGER(I4B) :: i,j,k
    REAL(DP), DIMENSION(N) :: w 
    REAL(DP), DIMENSION(N) :: c,x
    REAL(DP), DIMENSION(N,N) :: v
    REAL(DP) :: wmax,wmin,rcond,tmp
    CHARACTER(3) :: dummy
    !open(7,file='MATRX1.DAT',status='old')
    !do
    !read(7,'(a)') dummy
    !if (dummy == 'END') exit
    !read(7,*)
    !read(7,*) n,m
    !read(7,*)
    !read(7,*) ((a(i,j), j=1,n), i=1,n)
    !read(7,*)
    !read(7,*) ((b(i,j), i=1,n), j=1,m)
    !	copy a into u
    u(1:n,1:n)=a(1:n,1:n)
    !	decompose matrix a
    call svdcmp(u(1:n,1:n),w(1:n),v(1:n,1:n))
    !	find maximum singular value
    wmax=max(maxval(w(1:n)),0.0_dp)
    wmin = (minval(w(1:n)))
    tmp = wmin
    rcond = wmax/wmin
    !	define "small"
    wmin=wmax*(1.0e-5)
    !wmin = 1.0e-6
    !	zero the "small" singular values
    where (w(1:n) < wmin) w(1:n)=0.0
    !	backsubstitute for each right-hand side vector
    do j=1,m
     !  write(*,'(1x,a,i2)') 'Vector number ',j
       c(1:n)=b(1:n,j)
       call svbksb(u(1:n,1:n),w(1:n),v(1:n,1:n),c(1:n),x(1:n))
       !write(*,*) '    Solution vector is:...rcond=',rcond
      ! write(*,'(1x,6f12.6)') (x(i), i=1,n)
      ! write(*,*) '    Original right-hand side vector:'
       !write(*,'(1x,6f12.6)') (c(i), i=1,n)
       !write(*,*) '    Result of (matrix)*(sol''n vector):'
       !c(1:n)=matmul(a(1:n,1:n),x(1:n))
       !write(*,'(1x,6f12.6)') (c(i), i=1,n)
       if(x(1).lt.0.)then
!!$          do k = 1,n
!!$          write(*,'(1x,6e12.5)')(a(k,i),i=1,n)
!!$       end do
!!$          print*,'wmax=',wmax,'  wmin=',tmp,' rcond = ',rcond
!!$          print*,'x(1)=',x(1),x(2)
!!$          write(*,*) '    Solution vector is:...rcond=',rcond
!!$          write(*,'(1x,e12.5)') (x(i), i=1,n)
!!$          write(*,*) '    Original right-hand side vector:'
!!$          write(*,'(1x,6f12.6)') (c(i), i=1,n)
!!$          write(*,*) '    Result of (matrix)*(sol''n vector):'
!!$          c(1:n)=matmul(a(1:n,1:n),x(1:n))
!!$          write(*,'(1x,6f12.6)') (c(i), i=1,n)
         ! read(*,*)
       else 
!!$          do k = 1,n
!!$             write(*,'(1x,6e12.5)')(a(k,i),i=1,n)
!!$          end do
!!$          print*,'not singular','wmax=',wmax,'  wmin=',tmp,' rcond = ',rcond
          !   print*,a
       end if
       
       sol(:,j) = x(:)
    end do
    !write(*,*) '***********************************'
    !write(*,*) 'Press RETURN for next problem'
    !read(*,*)
    !end do
    !close(7)
  END Subroutine xsvbksb
end Module xsvd
 

 

   
