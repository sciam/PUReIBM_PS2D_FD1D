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

module init_turb
#include "ibm.h"

	USE precision 
	USE constants

	use global_data, only : nx, mx1, my, my2, dx, dy, dz, u,&
						&  vis, ndim, run_name, doml, &
						&  maxvolfrac, re, lybyd, &
						&  dia_phys, fluid_atijk, input_type, aliasflag, &
						&  run_name, doml, dia_phys, &
						&  input_type, iglobstep, t, minunitno, &
						&  maxunitno, count_fluid, irestart, rhof, tke_converged
#if PARALLEL
	use global_data , only : err_code, decomp_group, myid, starts, ends
#endif
	use nlmainarrays , only : ubcp, pbcp
	use fftw_interface
	use constants
	use nlarrays
	use randomno
	use general_funcs
	use string_funcs

	type::kappa_type
		real(prcn) :: val
		integer :: num
	end type kappa_type

	real(prcn), allocatable, dimension(:) :: kx, ky, kz

	type(kappa_type),allocatable :: kappa(:)
	real(prcn) :: k0, kmin, kmax, dkappa

	real(prcn),    allocatable, target, dimension(:,:,:,:) :: turb_ur, turb_ur2
	complex(prcn), allocatable, target, dimension(:,:,:,:) :: turb_uf

!	complex(prcn), dimension(:,:,:), allocatable :: uf !,uf2d(:,:)
!	real(prcn),    dimension(:,:,:), allocatable :: ur !,ur2d(:,:)

	real(prcn) :: avr(3)
	real(prcn), save :: epsf, tke_i, epsf_i

	real(prcn) :: eta, l_e, eta_i, l_e_i
	real(prcn) :: c_eta, c_l
	real(prcn) :: divergence, mean_divergence
	character*3  rank_string
	integer :: s_size	

	logical, allocatable :: far_field(:,:,:)

	real(prcn) :: lx,ly,lz
	real(prcn) :: re_lambda, turb_length, turb_intensity, dbyeta, lambda, kmaxeta, u_eta, tau_eta, u_eta_i, tau_eta_i
	integer :: nkappa, ntke, nstat
	logical :: rand_gen, fluc_on, auto_on, kappa_on, pdf_on, detail_stat

	NAMELIST /turb/ re_lambda, turb_length, turb_intensity, dbyeta, &
	              & nkappa, fluc_on, auto_on, kappa_on, pdf_on, rand_gen, nstat, ntke, detail_stat, tke_error
contains

	subroutine allocation
		implicit none
		real(prcn) :: k2, kmode
		integer :: i, j, k, kindex

		lx=DOML(1)
		ly=DOML(2)
		lz=DOML(3)
	
		if (I_AM_NODE_ZERO) write (*,*) "LENGTH OF THE BOX = ", doml

	!	open (unit=50, file="turb.in", action="read")
	!	read (50, NML=turb)
	!	close (50)

		!	initialize wavenumber vectors
		allocate(kx(my2), ky(my), kz(my))

		do i=1,my2
			kx(i)=twopi*dble(i-1)/dble(lx)
		enddo

		do i=1,my/2
			ky(i) =            twopi*dble(i-1)/dble(ly)
			ky(my+1-i) = -1.d0*twopi*dble(i)/dble(ly)
		enddo

		do i=1,my/2
			kz(i) =           twopi*dble(i-1)/dble(lz)
			kz(my+1-i)= -1.d0*twopi*dble(i)/dble(lz)
		enddo

		kmin = kx(2)

		if (aliasflag==0) then
			kmax = my*kmin/2.
			!kmax = sqrt(2.0)/3.*dble(my)*kmin
		else
			kmax = sqrt(2.0)/3.*dble(my)*kmin
		endif
		k0   = kx(2)*(1.1) !sqrt(kx(2)**2+ky(2)**2+kz(2)**2)

#if PARALLEL
		if (.not.I_AM_NODE_ZERO) return
#endif

		! Automatically generates Fourier shells
		!dkappa = k0 !+ky(2)*dconjg(ky(2))+kz(2)*dconjg(kz(2)))
		!nkappa = (kmax-k0)/dkappa+1

		! Numbere of Fourier shells is determined by user
		dkappa=(kmax-k0)/(nkappa-1)
		!dkappa=kmax/nkappa

!		if (dkappa<1.73*kx(2)) then
!			write (*,*) "WARNING: THE ENERGY SPECTRUM SHELLS ARE FINE."
!			write (*,*) "MAKE THE SHELLS COURSER BY DECREASING THE NKAPPA VALUE"
!			STOP
!		endif

		!k0 = dkappa

		allocate(kappa(nkappa))

		! Finding number of nodes on each shell
		kappa(:)%num=0
		do k=1, my
			do j=1, my
				do i=1, my2
					if (i/=1 .or. j/=1 .or. k/=1) then
						k2	  = kx(i)**2+ky(j)**2+kz(k)**2
						kmode = sqrt(k2)
						if (kmode<k0) then
							kindex=1
						else
							kindex = int((kmode-k0)/dkappa)+2
							if (kindex>nkappa) kindex=nkappa
						endif

!						kindex = int(kmode/dkappa)+1
!						if (kindex>nkappa) kindex=nkappa

						if (i==1) then
							kappa(kindex)%num = kappa(kindex)%num + 1
						else
							kappa(kindex)%num = kappa(kindex)%num + 2
						endif
					endif
				end do
			end do
		end do

		write (*,*) "NUMBER OF GRID POINTS = ", sum(kappa(:)%num)

		do i=1,nkappa
			if (kappa(i)%num==0) then
				write (*,*) "WARNING: THEHE IS A FOURIER SHELL IN WHICH THERE IS NO GRID POINT."
				write (*,*) "CHANGE THE 'NKAPPA' VALUE"
				STOP
			endif
		enddo

!		if (ik_on>0) then
!			open(unit=9997,file=trim(run_name)//"_kappa.dat",status="replace",position="append")
!			write (9997,*) "zone"
!			do i=1,nkappa
!				write (9997,"(1D15.7,1I8)") kappa(i)%val,kappa(i)%num
!			enddo
!			close(9997)
!		endif

		c_eta = 0.4
		c_L   = 6.78
	end subroutine allocation

	subroutine specinit
		implicit none
		double precision :: kmode
		complex(prcn) :: theta1, theta2
		complex(prcn) :: alpha, beta
		real(prcn)    :: dk

		! random number generator variables
		integer :: seed, nsample
		real(prcn),allocatable :: randnum0(:)
		real,allocatable    :: rand_num(:), rand_num1(:), rand_num2(:), rand_num3(:)

		! local variable
		integer :: i, j, k, l, m, idim, ierr, unit1, nhbins

		real(prcn) :: k2, umag
		integer :: ikappa,kindex
		real(prcn) :: Re_L, phi, kxy

		real(prcn) :: tmp
		complex(prcn) :: ctmp
		real(prcn), dimension(:), allocatable :: hist, wt
	
		real(prcn), external :: energy_init, integrate_energy, integrate_epsilon
		real(prcn) :: k_x0, k_mag0
		real(prcn) :: k_x1, k_y1, k_z1
		real(prcn) :: phi_s
		real(prcn) :: u_min,u_max,u_prime
		!--------------------------------------
		write (*,*) "GENERATING TURBULENT VELOCITY FLUCTUATIONS"

		lx=DOML(1)
		ly=DOML(2)
		lz=DOML(3)

		call allocation

#if PARALLEL
		if (.not.I_AM_NODE_ZERO) goto 10
#endif

		allocate(turb_ur(my,my,my,ndim))
		allocate(turb_uf(my2,my,my,ndim))
!		allocate(engmag(0:nkappa),epsfmag(0:nkappa))

!		call random_seed()
!		call getnewseed_mod(seed)
!		write(*,'(A,1(I13,1x))') 'seed',seed

		!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
		!IIIII CALCULATING PARAMETERS OF VELOCITY FLUCTUATIONS IIIIIII
		l_e   = turb_length * twopi !ly
		re_l  = re_lambda**2 * 3./20.
		eta   = l_e / ( re_l**0.75 )
		tke   = (re_l * vis / l_e)**2
		epsf  = tke**1.5 / l_e

!		re_lambda = dsqrt(20d0/3d0 * re_l)
		u_prime   = dsqrt(tke * 2d0 / 3d0)
		lambda    = dsqrt(15*vis*u_prime**2/epsf)
		kmaxeta   = kmax*eta
		u_eta	= (epsf*vis)**0.25
		tau_eta	= dsqrt(vis/epsf)

		tke_i  = tke
		epsf_i = epsf
		u_eta_i = u_eta
		tau_eta_i = tau_eta
		eta_i = eta
		l_e_i = l_e

		write (*,*) '-------------------------------'
		write (*,*) 're_lambda = ', re_lambda
		write (*,*) 'Re_L      = ', Re_L
		write (*,*) '-------------------------------'
		write (*,*) 'tke       = ', tke
		write (*,*) 'epsf      = ', epsf_i
		write (*,*) 'u_prime   = ', u_prime
		write (*,*) 'u_eta     = ', u_eta
		write (*,*) 'tau_eta   = ', tau_eta
		write (*,*) '-------------------------------'
		write (*,*) 'Kmin      = ', kmin
		write (*,*) 'Kmax      = ', kmax
		write (*,*) 'dkappa    = ', dkappa
		write (*,*) 'nkappa    = ', nkappa
		write (*,*) '-------------------------------'
		write (*,*) 'eta       = ', eta
		write (*,*) 'eta/dx    = ', eta/dx
		write (*,*) 'lambda    = ', lambda
		write (*,*) 'L_e       = ', l_e
		write (*,*) 'Kmax*eta  = ', kmaxeta
		write (*,*) '-------------------------------'

		write (*,*) "FINDING ENERGY SPECTRUM FUNCTION COEFFICIENTS"
		c_eta = 0.4
		c_L   = 6.78

		call xmnewt(c_eta, c_l, tke, epsf_i, l_e, eta, kmax, vis) !for high reynolds numbers

		write (*,"(1a,2d15.7)") "C_eta, C_l = ", c_eta, c_l
		write (*,*) integrate_energy(kmin,kmax,epsf_i, l_e, eta, c_eta, c_l)

		do i=1, nkappa
			if (i==1) then
				k_x0 = 1e-6
				k_x1 = k0
			else
				k_x0 = k0+(i-2) * dkappa
				k_x1 = k0+(i-1)   * dkappa
			endif

			call xrtnewt(kappa(i)%val,k_x0,k_x1,epsf_i,l_e,eta,c_eta,c_l)
		enddo


!		if (ik_on>0) then
			open(unit=9997,file=trim(run_name)//"_kappa.dat",status="replace",position="append")
			write (9997,*) "zone"
			do i=1,nkappa
				write (9997,"(2D15.7,1I8)") (i-0.5)*dkappa, kappa(i)%val,kappa(i)%num
			enddo
			close(9997)
!		endif


		!##################################################

		turb_uf = czero
		nsample=my2*my*my
!		allocate(rand_num(nsample))
		allocate(randnum0(2))
		allocate(rand_num1(nsample),rand_num2(nsample),rand_num3(nsample))
	
		if (rand_gen) then
!			call uni_dist(randnum0)
			open(unit=999960,file=trim(run_name)//"_seed.d",status="old",action="read")
			read (999960,*) seed
			close(999960)

			call RLuxGo(389, seed, 0, 0)

!			open(unit=999951,file=trim(run_name)//"_rand.dat",status="replace")

			open(unit=999951,file=trim(run_name)//"_rand1.dat",status="replace",form="unformatted")
			open(unit=999952,file=trim(run_name)//"_rand2.dat",status="replace",form="unformatted")
			open(unit=999953,file=trim(run_name)//"_rand3.dat",status="replace",form="unformatted")

!			call uni_dist(rand_num)
!			call uni_dist(rand_num1)
!			call uni_dist(rand_num2)
!			call uni_dist(rand_num3)
!			write (999951,*) rand_num
!			close (999951)

			call RanLux(rand_num1, nsample)
			call RanLux(rand_num2, nsample)
			call RanLux(rand_num3, nsample)

			write (999951) rand_num1
			write (999952) rand_num2
			write (999953) rand_num3

			close (999951)
			close (999952)
			close (999953)

			unit1 = 999955
			open  (unit=unit1, file=trim(run_name)//"_pdf_rand.dat", status="replace")
			write (unit1,*) "variables=V,F(V)"

			if (allocated(wt)) deallocate(wt)
			allocate(wt(nsample))
			do i=1, nsample
				wt(i)=1d0/nsample
			enddo
			nhbins=50
			allocate(hist(nhbins))

			CALL histogram1(dble(rand_num1(:)),wt,nsample,nhbins,u_min,u_max,hist(1:nhbins))
			CALL plothist(hist(1:nhbins),u_min,u_max,nhbins,unit1,real(idim,prcn),1.d0)

			CALL histogram1(dble(rand_num2(:)),wt,nsample,nhbins,u_min,u_max,hist(1:nhbins))
			CALL plothist(hist(1:nhbins),u_min,u_max,nhbins,unit1,real(idim,prcn),1.d0)

			CALL histogram1(dble(rand_num3(:)),wt,nsample,nhbins,u_min,u_max,hist(1:nhbins))
			CALL plothist(hist(1:nhbins),u_min,u_max,nhbins,unit1,real(idim,prcn),1.d0)

			close (unit1)
			deallocate(wt,hist)
		else
!			open(unit=9991,file=trim(run_name)//"_rand.dat",status="old")
			open(unit=999951,file=trim(run_name)//"_rand1.dat",status="old",form="unformatted")
			open(unit=999952,file=trim(run_name)//"_rand2.dat",status="old",form="unformatted")
			open(unit=999953,file=trim(run_name)//"_rand3.dat",status="old",form="unformatted")

			read (999951) rand_num1
			read (999952) rand_num2
			read (999953) rand_num3
			close (999951)
			close (999952)
			close (999953)
		endif

		! generating velocity fluctuations
		l=0
		do k=1, my
!			l=0

			do j=1,my 
				do i=1,my2
!					l=l+3
					if ((i==1) .and. (j==1) .and. (k==1)) then
						turb_uf(1,j,1,1) = cmplx(0.d0,0.d0)
						turb_uf(1,j,1,2) = cmplx(0.d0,0.d0)
						turb_uf(1,j,1,3) = cmplx(0.d0,0.d0)
						kmode = 0.d0
						umag = 0.d0
					else
						l=l+1

						! The real value of wave number for this node is computed
						k2     = kx(i)**2+ky(j)**2+kz(k)**2
						kmode  = sqrt(k2)
						k_mag0 = kmode

						! It is determmined with which shell the energy of this node is associated
						if (kmode<k0) then
!							kindex=0
							kindex=1
						else
							kindex = int((kmode-k0)/dkappa)+2
							if (kindex>nkappa) kindex=nkappa
						endif

!						kindex = int(kmode/dkappa)+1
!						if (kindex>nkappa) kindex=nkappa

						kmode = kappa(kindex)%val

						! The components of the new kappa would be computed to make it consistent with
						! the assumption of the discret energy spectrum on Fourier shells
						k_x1 = kx(i)
						k_y1 = ky(j)
						k_z1 = kz(k)

						phi_s = kmode/k_mag0
						k_x1  = k_x1*phi_s
						k_y1  = k_y1*phi_s
						k_z1  = k_z1*phi_s
						kxy   = dsqrt(k_x1**2+k_y1**2)
						!------------------------


!						theta1	= cmplx(0,twopi*rand_num(l-2))
!						theta2	= cmplx(0,twopi*rand_num(l-1))
!						phi	= twopi*rand_num(l)

						theta1 = cmplx(0,twopi*rand_num1(l))
						theta2 = cmplx(0,twopi*rand_num2(l))
						phi    = twopi*rand_num3(l)

						umag = energy_init(kmode, epsf_i, l_e_i, eta_i, c_eta, c_l) !+energy_init(kmode_plus,epsf,l_e,eta))/2

						! alpha and beta are computed according to the method presented
						! E(k)*dkappa produces the energy of each shell. Then it is distributed
						! uniformly among nodes related to the shell. Acoordingly, the manitude of
						! each velocity fluctuation is computed according to the amount of energy
						! it has.

						dk=dkappa
						if (kindex==1) dk=k0

						alpha = sqrt(2 * umag * dk / kappa(kindex)%num) * exp(theta1) * cos(phi)
						beta  = sqrt(2 * umag * dk / kappa(kindex)%num) * exp(theta2) * sin(phi)

						if (i==1 .and. j==1 .and. k>1) then
							turb_uf(i,j,k,1) = alpha
							turb_uf(i,j,k,2) = beta
							turb_uf(i,j,k,3) = czero
						else
							turb_uf(i,j,k,1) = alpha*kmode*k_y1+beta*k_x1*k_z1
							turb_uf(i,j,k,1) = turb_uf(i,j,k,1)/(kmode*kxy)

							turb_uf(i,j,k,2) = beta*k_y1*k_z1-alpha*kmode*k_x1
							turb_uf(i,j,k,2) = turb_uf(i,j,k,2)/(kmode*kxy)

							turb_uf(i,j,k,3) = -1d0*beta*kxy/kmode
						end if
					end if
				end do
			end do
		end do
		close(9991)
	
		!for i=1,j=1,k<>1 line
		i=1
		j=1
		do k=my/2+1,my
			turb_uf(i,j,k,1) = conjg(turb_uf(i,j,my+2-k,1))
			turb_uf(i,j,k,2) = conjg(turb_uf(i,j,my+2-k,2))
			turb_uf(i,j,k,3) = conjg(turb_uf(i,j,my+2-k,3))
		end do

		!for i=1, j<>1, k=1 line
		i=1
		k=1
		do j=2,my/2
			turb_uf(i,my+2-j,k,1) = dcmplx(dble(turb_uf(i,j,k,1)),(-1.0)*aimag(turb_uf(i,j,k,1)))
			turb_uf(i,my+2-j,k,2) = dcmplx(dble(turb_uf(i,j,k,2)),(-1.0)*aimag(turb_uf(i,j,k,2)))
			turb_uf(i,my+2-j,k,3) = dcmplx(dble(turb_uf(i,j,k,3)),(-1.0)*aimag(turb_uf(i,j,k,3)))
		end do

		!for i=1, J<>1, k<>1 
		 i=1
		 do k=2, my/2
		    do j=2, my/2
		       turb_uf(i,my+2-j,my+2-k,:) = conjg(turb_uf(i,j,k,:))
		       turb_uf(i,my+2-j,k,:)      = conjg(turb_uf(i,j,my+2-k,:))
		    end do
		 end do
    
!^^^^ Mohammad 11-13-2009 ^^^^ existance of this part does not have any effect on the answers
!		turb_uf(1,1,my/2+1,:) = czero
!		turb_uf(1,ny/2+1,1,:) = czero
!		turb_uf(1,ny/2+1,:,1) = czero
!		turb_uf(1,ny/2+1,:,2) = czero
!		turb_uf(1,ny/2+1,:,3) = czero
!
!		turb_uf(1,:,nz/2+1,1) = czero
!		turb_uf(1,:,nz/2+1,2) = czero
!		turb_uf(1,:,nz/2+1,3) = czero
!
!		turb_uf(nx/2+1,:,:,:) = czero     ! i=nxturbby2p1 plane is set to zero
!		if they are added in then, some of the energy would be lost, whi
!-------------------------------------------------

		ubcp = zero

		do idim = 1, ndim
			call ff3cr(turb_uf(:,:,:,idim),turb_ur(:,:,:,idim))
		enddo
	
		call compute_average(0)
		call statistics(0)
!		stop

		if (allocated(rand_num)) deallocate(rand_num)
		if (allocated(rand_num1)) deallocate(rand_num1, rand_num2, rand_num3)
		if (allocated(randnum0)) deallocate(randnum0)

		deallocate(turb_uf)
!		allocate(turb_uf(my,my2,my,ndim))

		do idim=1, ndim
			do i=1, nx
				call ff2rc(turb_ur(i,:,:,idim), u(i,:,:,idim))
			enddo
		enddo
		deallocate(turb_ur)

		write (*,*) "specinit.F90 is done..."
		call screen_separator(80,'I')

10		continue

#if PARALLEL
		BROADCAST_DOUBLE(tke,  1, NODE_ZERO, decomp_group)
		BROADCAST_DOUBLE(epsf, 1, NODE_ZERO, decomp_group)
#endif


		write (*,*) "STOP"
		PARALLEL_FINISH()
		stop

		return
	end subroutine specinit

	subroutine compute_average(call_flag)
		implicit none

		integer :: call_flag
		real(prcn), dimension(:,:,:,:), pointer :: u_r

		real(prcn) :: tmp(3)
		integer :: count, tmp_count
		integer :: i, j, k, idim
		integer :: imax

		if (ntke<1) ntke=20
		if (.not.allocated(kx)) then
			call allocation
			tke_i = tke
			epsf_i = epsf
		endif

		if (call_flag==0) then
			u_r  => turb_ur
!			u_r  => ubcp
			imax =  nx
		else
			if (irestart==0.and.iglobstep==1) then
			elseif (irestart==1.and.iglobstep==1) then
			elseif (mod(iglobstep,ntke)/=0) then
				return
			endif
			u_r => ubcp
			imax = nx
		endif

		avr = 0d0
		count=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					do i=1, imax
						if (fluid_atijk(i,j,k)) then
							avr(idim) = avr(idim) + u_r(i,j,k,idim)
						endif
					end do
				end do
			end do
		enddo
#if PARALLEL
		if (call_flag==1) then
			tmp = 0d0
			GLOBAL_DOUBLE_SUM(avr,tmp,3,decomp_group)
			avr = tmp
		endif
#endif
		avr(:) = avr(:) / count_fluid !((1d0-maxvolfrac)*mx*my*mz)

		write (*,"(1a,3d15.7)") "UAVG = ", avr(:)

	end subroutine compute_average

	subroutine statistics(call_flag)
		implicit none
		integer :: call_flag

		real(prcn) :: kmode, k2, u_prime, k_x0, k_x1, dk
		real(prcn), external :: energy_init, integrate_energy, integrate_epsilon

		real(prcn),    dimension(:),     allocatable :: engmag, epsfmag
		real(prcn),    dimension(:),     allocatable :: hist, wt

		real(prcn),    dimension(:,:,:,:), pointer :: u_r
		complex(prcn), dimension(:,:,:,:), pointer :: u_f

		real(prcn) :: krange, energyf_pre(3)
		real(prcn) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8

		integer :: i, j, k, l, m, idim, ierr, unit1, nhbins, nsample
		integer :: imax, jmax, kkmax, ikrange
		integer	:: ikappa,kindex

		real(prcn) :: l_e, Re_L
		real(prcn) :: umag, phi

		real(prcn) :: var(3), skew, kurt, turb(ndim), bij(ndim), bij_tmp(ndim)
		real(prcn) :: u_min, u_max, tmp, mean_energy
		complex(prcn) :: ctmp
		integer :: unitno
		logical :: filexist, isopen, new_file
		character*30 filename1
		character*7  step_string
		!-------------------------------------------------------------

		if (.not.allocated(kx)) then
			call allocation
			tke_i  = tke
			epsf_i = epsf
		endif

		if (call_flag==0) then
			u_f => turb_uf
			u_r => turb_ur
!			u_r => ubcp
		else
			if (mod(iglobstep,nstat)/=0.or.iglobstep==1) return
!			if (mod(iglobstep, detail_stat)/=0) goto 10
#if !PARALLEL
			if (.not.allocated(turb_ur)) allocate(turb_ur(nx ,my,mz,ndim))
			if (.not.allocated(turb_uf)) allocate(turb_uf(my2,my,mz,ndim))

			u_f=>turb_uf
			u_r=>turb_ur

!			do k=1, mz
!				do j=1, my
!					do i=1, nx
!						u_r(i,j,k,:) = ubcp(i,j,k,:) - avr(:)
!					enddo
!				enddo
!			enddo

!			if (.not.allocated(ur)) allocate(ur(my,my,my),uf(my2,my,my))
!			do idim = 1, ndim
!				do k=1, my
!					do j=1, my
!						do i=1, my
!							ur(i,j,k) = u_r(i,j,k,idim)
!						end do
!					end do
!				end do
!
!				call ff3rc(ur(:,:,:),uf(:,:,:))
!				u_f(:,:,:,idim)=uf(:,:,:)
!			enddo
!			deallocate(ur,uf)
			do idim = 1, ndim
				call ff3rc(u_r(:,:,:,idim),u_f(:,:,:,idim))
			enddo
#endif
		endif

#if PARALLEL
		if (call_flag==1) goto 10
#endif
		tke = 0d0
		epsf = 0d0
		energyf_pre = 0.d0

		do idim=1, ndim
			do k=1, my
				do j=1, my
					do i=1, my2
						k2  = kx(i)**2+ky(j)**2+kz(k)**2
						tmp = dble(u_f(i,j,k,idim) * conjg(u_f(i,j,k,idim)))
						if ( i == 1 ) tmp = tmp/2
				
						tke = tke + tmp
						epsf = epsf + 2 * vis * k2 * tmp
						energyf_pre(idim) = energyf_pre(idim)+tmp
					enddo
				end do
			end do
		enddo

10		if (I_AM_NODE_ZERO) then
			eta  = (vis**3 / epsf)**(0.25d0)
			l_e  = tke**1.5 / epsf
			re_l = l_e * dsqrt(tke) / vis
			re_lambda = dsqrt(20d0/3d0 * re_l)
			u_prime   = dsqrt(tke * 2d0 / 3d0)
			lambda    = dsqrt(15*vis*u_prime**2/epsf)
			kmaxeta   = kmax*eta
			u_eta	= (epsf*vis)**0.25
			tau_eta	= dsqrt(vis/epsf)
		
!			if (call_flag>0.and.(mod(iglobstep,detail_stat)/=0.or.detail_stat==0)) goto 20
	
!			if (input_type=="single-phase") then
!				umeanslip = 0d0
!			endif
!
!			re = umeanslip * (2 * pi) / vis

			if (call_flag==1) then
!				write (*,*) "FINDING ENERGY SPECTRUM FUNCTION COEFFICIENTS"
!				call xmnewt(c_eta,c_l,tke,epsf,l_e,eta,kmax*3,vis)  !FOR HIGH REYNOLDS NUMBERS
				do i=1,nkappa
					if (i==1) then
						k_x0 = 1e-6
						k_x1 = dkappa
					else
						k_x0 = (i-1) * dkappa
						k_x1 = (i)   * dkappa
					endif
					


					call xrtnewt(kappa(i)%val,k_x0,k_x1,epsf_i,l_e_i,eta_i,c_eta,c_l)
				enddo
			endif

!			if (ik_on>0) then
!				open(unit=9997,file=trim(run_name)//"_kappa.dat",status="old",position="append")
!				write (9997,*) "zone"
!				do i=0,nkappa
!					write (9997,"(1D15.7,1I8)") kappa(i)%val,kappa(i)%num
!				enddo
!				close(9997)
!			endif

20			write (*,*) "###### OUTPUT ACCORDING TO SINGLE PHASE FLOW #####"
			write (*,*) '-------------------------------'
			write (*,*) 'C_eta     = ', c_eta
			write (*,*) 'C_L       = ', c_l
			write (*,*) '-------------------------------'
			write (*,*) 'Re_p      = ', re
			write (*,*) 're_lambda = ', re_lambda
			write (*,*) 'Re_L      = ', Re_L
			write (*,*) 'turb_int  = ', turb_intensity
			write (*,*) 'meanslip  = ', umeanslip
			write (*,*) '-------------------------------'
			write (*,*) 'tke       = ', tke
			write (*,*) 'epsf      = ', epsf
			write (*,*) 'vis       = ', vis
			write (*,*) 'u_prime   = ', u_prime
			write (*,*) 'u_eta     = ', u_eta
			write (*,*) 'tau_eta   = ', tau_eta
			write (*,*) '-------------------------------'
			write (*,*) 'Kmin      = ', kmin
			write (*,*) 'Kmax      = ', kmax
			write (*,*) 'dkappa    = ', dkappa
			write (*,*) 'nkappa    = ', nkappa
			write (*,*) '-------------------------------'
			write (*,*) 'eta       = ', eta
			write (*,*) 'dx/eta    = ', dx / eta
			write (*,*) 'eta/L_e   = ', eta / l_e
			write (*,*) 'eta/L     = ', eta / doml(2)
			write (*,*) 'lambda    = ', lambda
			write (*,*) 'lambda/dx = ', lambda / dx
			write (*,*) 'lambda/L_e= ', lambda / l_e
			write (*,*) 'lambda/L  = ', lambda / doml(2)
			write (*,*) 'L_e       = ', l_e
			write (*,*) 'L_e/dx    = ', l_e/dx
			write (*,*) 'L_e/eta   = ', l_e/eta
			write (*,*) "Box Length= ", doml(2)
			write (*,*) 'L_e/Ly    = ', l_e / doml(2)
			write (*,*) 'Kmax*eta  = ', kmaxeta
			write (*,*) '-------------------------------'
			write (*,*) 'eta/eta_i         = ', eta/eta_i
			write (*,*) 'l_e/l_e_i         = ', l_e/l_e_i
			write (*,*) 'tke/tke_i         = ', tke/tke_i
			write (*,*) 'eps/eps_i         = ', epsf/epsf_i
			write (*,*) 'u_eta/u_eta_i     = ', u_eta/u_eta_i
			write (*,*) 'tau_eta/tau_eta_i = ', tau_eta/tau_eta_i

			write (*,*) '-------------------------------'
			if (kmaxeta<1.5) write (*,*) "WARNING: KAMX*ETA<1.5; TE CRITERION IS NOT SATISFIED"

			filename1=trim(run_name)//"_statistics.dat"
			new_file= .false.
			if (irestart==0.and.call_flag==0) then
				new_file = .true.
			else
				inquire(file=trim(filename1),exist=filexist,opened=isopen)
				if (.not.filexist) then
					new_file= .true.
				endif
			endif

			if (new_file) then
				open (unit=9992, file=trim(filename1), status="replace")
				write (9992,"(A)") "variables=t,re_lambda,re_l,tke,u_p,epsf,eta,eta/le,eta/L,lambda,lambda/le, lambda/L,le, &
								& le/L,u_eta,tau_eta,kmaxeta,Ceta,Cl,max_divergence,mean_divergence"
				write (9992,*) "zone"
			else
				open (9992,file=trim(run_name)//"_statistics.dat",form="formatted",status="old",position="append")
			endif

			write (9992,"(21D15.7)") t, re_lambda, re_l, tke, u_prime, epsf, eta, eta/l_e, eta/doml(2), &
						&  lambda, lambda/l_e, lambda/doml(2), l_e, l_e/doml(2), u_eta, &
						&  tau_eta, kmaxeta, c_eta, c_l, divergence, mean_divergence
			close (9992)
		endif

		if (call_flag==1) then
#if PARALLEL
!			BROADCAST_DOUBLE(umeanslip, 1, NODE_ZERO, decomp_group)
			return
#endif
		endif
	
!		if (call_flag>0.and.mod(iglobstep, detail_stat)/=0) return

		allocate(engmag(nkappa),epsfmag(nkappa))
		engmag	 = 0d0
		divergence	 = 0d0
		mean_divergence =0d0
		epsfmag	 = 0d0
		do k=1,my
			do j=1,my
				do i=1,my2
					k2    = kx(i)**2+ky(j)**2+kz(k)**2
					kmode = dsqrt(k2)
					if (kmode<k0) then
!						kindex=0
						kindex=1
					else
						kindex = int((kmode-k0)/dkappa)+2
						if (kindex>nkappa) kindex=nkappa
					endif

!					kindex = int((kmode-k0)/dkappa)+1
!					if (kindex>nkappa) kindex=nkappa
					kmode=kappa(kindex)%val

					tmp =       dble(u_f(i,j,k,1) * conjg(u_f(i,j,k,1)))
					tmp = tmp + dble(u_f(i,j,k,2) * conjg(u_f(i,j,k,2)))
					tmp = tmp + dble(u_f(i,j,k,3) * conjg(u_f(i,j,k,3)))
					if ( i == 1 ) tmp = tmp/2

					engmag(kindex) = engmag(kindex) + tmp
					epsfmag(kindex)= epsfmag(kindex)+ 2*vis*k2*tmp

					ctmp = (cmplx(0,kx(i)) * u_f(i,j,k,1) + cmplx(0,ky(j)) * u_f(i,j,k,2) + cmplx(0,kz(k)) * u_f(i,j,k,3))
					mean_divergence = mean_divergence + sqrt(dble(ctmp*conjg(ctmp)))
					if(sqrt(dble(ctmp*conjg(ctmp))) .gt. divergence) then
						divergence=sqrt(dble(ctmp*conjg(ctmp)))
					end if
				enddo
			end do
		end do

		write (*,*) "energy sum", sum(engmag(:))

		mean_divergence = mean_divergence / count_fluid
		write (*,*) 'Maximum divegence = ' , divergence
		write (*,*) 'Mean divegence = ' , mean_divergence

!		call divergence_3d(call_flag)

		call screen_separator(80,'T')

		!---------Start to output spectra1.dat---------
		call to_string(iglobstep,step_string,s_size)
		do i=1, 7-s_size
			step_string="0"//trim(step_string)
		enddo

		unitno = getnewunit(minunitno,maxunitno)
!		open  (unitno,file=trim(run_name)//"_spectra_"//step_string//".dat",form="formatted",status="replace")
		open  (unitno,file=trim(run_name)//"_spectra.dat",form="formatted",status="replace")
		write (unitno,*) "variables=Index,K,K*eta,K*L,E,E1,E2,eps,eps1,%E"
		write (unitno,*) "zone"

		!integration of discrete energy spectrum function: k0-kmax
		sum1=0d0
!		do i=0,nkappa
		do i=1,nkappa
			kmode = kappa(i)%val
			tmp	 = engmag(i)
			phi  = epsfmag(i)
		
			sum1 = sum1 + tmp
		
			if (i==1) then
				tmp = tmp / k0
				phi = phi / k0
			else
				tmp = tmp / dkappa
				phi = phi / dkappa
			endif

			write(unitno,'(1I,9D15.7)') i,kmode, kmode*eta, kmode*l_e,     &
					& tmp, tmp/(eta*u_eta*u_eta), tmp/(tke*l_e),			&
					& phi, phi/epsf, engmag(i)/tke
		end do

		write (unitno,*) "zone"
		!integration of energy spectrum function: k0-kmax
		krange = kappa(1)%val
		dk = 0.01
		ikrange = int((kmax-krange)/dk)

		do i=1,ikrange
			sum3=integrate_energy(krange+(i-1)*dk,krange+i*dk,epsf,l_e,eta,c_eta,c_l)

			kmode = krange +i*dk
			umag  = energy_init(kmode,epsf,l_e,eta,c_eta,c_l)
			write(unitno,'(1I,9D15.7)') i,kmode, kmode*eta, kmode*l_e,	&
					& umag, umag/(eta*u_eta*u_eta),umag/(tke*l_e),		&
					& 2.0*umag*vis*kmode**2, 2.0*umag*vis*kmode**2/epsf,sum3/tke
		enddo
		close(unitno)

		!integration of energy spectrum function: 0-k0
		sum2 = integrate_energy(1d-6,k0,epsf,l_e,eta,c_eta,c_l)

		!integration of energy spectrum function: k0-kmax
		sum3 = integrate_energy(k0,kmax,epsf,l_e,eta,c_eta,c_l)

		!integration of energy spectrum function: kmax-infinity
		sum4 = integrate_energy(kmax,3*kmax,epsf,l_e,eta,c_eta,c_l)

		!integration for epsilon: 0-k0
		sum5 = integrate_epsilon(1d-6,k0,vis,epsf,l_e,eta,c_eta,c_l)

		!integration for epsilon: k0-kmax
		sum6 = integrate_epsilon(k0,kmax,vis,epsf,l_e,eta,c_eta,c_l)

		!integration for epsilon: kmax-infinity
		sum7 = integrate_epsilon(kmax,3*kmax,vis,epsf,l_e,eta,c_eta,c_l)
	
		write (*,"(A,D)") 'TKE from Fourier space     (1) = ', tke
		write (*,"(A,D)") 'TKE from E.S.F,    0->k0   (2) = ', sum2
		write (*,"(A,D)") 'TKE from E.S.F,   k0->kmax (3) = ', sum3
		write (*,"(A,D)") 'TKE from E.S.F, kmax->inf. (4) = ', sum4
		write (*,"(A,D)") 'TKE from E.S.F,    0->kmax (5) = ', sum2 + sum3
		write (*,"(A,D)") 'TKE INITIAL                    = ', tke_i
		call screen_separator(40,'.')
	
		write (*,"(A,D)") '(1)/(5)                        = ', tke / (sum2 + sum3)
		write (*,"(A,D)") '(1)/TKE_I                      = ', tke / tke_i
		write (*,"(A,D)") '(5)/TKE_I                      = ', (sum2 + sum3) / tke_i

		call screen_separator(50,'*')

		write (*,"(A,D)") 'EPS from Fourier space     (6) = ', epsf
		write (*,"(A,D)") 'EPS from E.S.F,    0->k0   (7) = ', sum5
		write (*,"(A,D)") 'EPS from E.S.F,   k0->kmax (8) = ', sum6
		write (*,"(A,D)") 'EPS from E.S.F, kmax->inf. (9) = ', sum7
		write (*,"(A,D)") 'EPS from E.S.F,    0->kmax (10)= ', sum5 + sum6
		write (*,"(A,D)") 'EPS INITIAL                    = ', epsf_i
		call screen_separator(40,'.')

		write (*,"(A,D)") '(6) /(10)                      = ', epsf / (sum5 + sum6)
		write (*,"(A,D)") '(6) /EPS_I                     = ', epsf / epsf_i
		write (*,"(A,D)") '(10)/EPS_I                     = ', (sum5+sum6) / epsf_i
		call screen_separator(80,'I')
		 !-------------------------------------------------
		 ! checking phase space velocity in fluid space
		 !-------------------------------------------------
	
		if (pdf_on) then
			unit1 = getnewunit(minunitno,maxunitno)
!			open (unitno,file=trim(run_name)//"_pdf_"//step_string//".dat", status="replace")
			open (unitno,file=trim(run_name)//"_pdf.dat", status="replace")
			write (unit1,*) "variables=V,F(V)"
			write (unit1,*) 'Zone t="',0,'"'

			j=50
			do i=-j,j
				phi = 5d0/j * i
				tmp = 1d0/dsqrt(2*pi)*dexp(-phi**2/2)
				write (unit1,"(2D15.7)") phi, tmp
			enddo
		endif

		do idim = 1, ndim
			u_max = maxval(u_r(1:nx,:,:,idim))
			u_min = minval(u_r(1:nx,:,:,idim))
			call gauss3D(u_r(1:nx,1:my,1:mz,idim), avr(idim), var(idim), skew, kurt, nx, my, mz)	

			turb(idim) = var(idim) / 2
			write (*,*) "For direction # ", idim
			write (*,*) "u_min/r.m.s. = " , u_min /dsqrt(var(idim))
			write (*,*) "u_max/r.m.s. = " , u_max /dsqrt(var(idim))
			write (*,*) "Average     = " , avr(idim)
			write (*,*) "Variance    = " , var(idim)
			write (*,*) "Skewness    = " , skew
			write (*,*) "Flatness    = " , kurt
			write (*,'(A,2D15.5)') "0.5<u*u> (physical, Fourier space)= " , turb(idim), energyf_pre(idim)!,turb()-energyf_pre(1) 
			write (*,*) "..............................."

			if (pdf_on) then
				if (.not.allocated(wt)) then
					nsample = my * my * my
					allocate(wt(nsample))
					do i=1, nsample
						wt(i)=1d0/nsample
					enddo
					nhbins=50
					allocate(hist(nhbins))
				endif
			
				CALL histogram1(u_r(:,:,:,idim),wt,nsample,nhbins,u_min,u_max,hist(1:nhbins))
				CALL plothist(hist(1:nhbins),u_min,u_max,nhbins,unit1,real(idim,prcn),1.d0)
				if (idim==3) close(unit1)
			endif
			write (*,*) "-------------------------------"
		enddo
	

		if (allocated(wt)) then
			deallocate(wt)
			deallocate(hist)
		endif
		!-------------------------------------------------

		write (*,"(A,1D)") 'TKE from inverse of the Fourier space = ', sum(turb(1:3))
	 
!		call isotropy_check(u_r,avr,bij,sum(turb(1:3)),my,my,my)
!		write (*,'(A,3D15.7)') 'The Bij is = ',bij
!		write (*,'(A,3D15.7)') 'The diagonal terms are = ',turb/sum(turb(1:3))-1./3.
!		call screen_separator(50,'-')

		if (fluc_on) then
			mean_energy = sum(turb(1:3))

			write (*,*) "WRITING THE VELOCITY FLUCTUATIONS INTO OUTPUT FILE"
			open(unit=9994,file=trim(run_name)//"_velocity.dat",status="replace",form="formatted")
			write (9994,"(a,3(i,a))") "zone i=", nx, "j=", my, "k=", mz, "f=point"

			do k=1, mz
				do j=1, my
					do i=1,nx
						tmp = sqrt(dot_product(u_r(i,j,k,:),u_r(i,j,k,:)))
						write (9994,"(3i6,4d15.7)") i, j, k, tmp/sqrt(mean_energy), u_r(i,j,k,:)/sqrt(mean_energy)
					enddo
				enddo
			enddo
			close (9994)
		endif

!		if (iauto_on==1) then
!			write (*,*) "GENERATING LONGITUDINAL AUTOCORRELATION FUNCTIONS"
!			nhbins=50
!			call autocorrelation(nhbins,0)
!			stop
!		end if
		!-------------------------------------------	
		deallocate(engmag,epsfmag)

		if (call_flag==1) then
			if (allocated(turb_ur)) deallocate(turb_ur)
			if (allocated(turb_uf)) deallocate(turb_uf)
		endif
	end subroutine statistics

	subroutine gauss3D(uin, ave, var, skew, kurt, nx, ny, nz)
	  implicit none
	  integer, intent(in) :: nx,ny,nz
	  double precision, intent(in) :: uin(nx,ny,nz)
	  double precision, intent(inout) :: ave, var, skew, kurt
	  
	  !------------------------
	  ! local variables

	  integer :: i,j,k,ierr
	  double precision :: tmp(4)
	  
	  tmp(:) = 0.d0

	  do k=1,nz            ! decompose along z axis
		  do j=1,ny
		     do i=1,nx
		        tmp(1) = tmp(1) + uin(i,j,k)
		     end do
		  end do
	  end do

	  ave = tmp(1)/dble(nx*ny*nz)
	  
	  tmp(:) = 0.d0
	  skew = 0.d0
	  kurt = 0.d0
	  do k=1,nz            ! decompose along z axis
		  do j=1,ny
		     do i=1,nx
		        tmp(1) = tmp(1) + (uin(i,j,k)-ave)*(uin(i,j,k)-ave)  ! var
		        tmp(2) = tmp(2) + uin(i,j,k)*uin(i,j,k)              ! 
		        tmp(3) = tmp(3) + (uin(i,j,k)-ave)**3                ! skew
		        tmp(4) = tmp(4) + (uin(i,j,k)-ave)**4                ! kurt
		     end do
		  end do
	  end do

	  var  = tmp(1) / dble(nx*ny*nz)
!	  tkein = tmp_r(2)/dble(nx*ny*nz_tot)*0.5
	  skew = tmp(3) / (dble(nx*ny*nz)*var**1.5) !/ var**1.5
	  kurt = tmp(4) / (dble(nx*ny*nz)*var**2) !/ var**2

	  return
	end subroutine gauss3D


#if 0
	subroutine tke_eps(call_flag)
!		use dependent_functions
		implicit none
		integer :: call_flag
		real(prcn), dimension(:,:,:,:), pointer :: ur
		complex(prcn), dimension(:,:,:,:), pointer :: uf

		real(prcn), save, dimension(ndim,ndim) :: uiuj, uiuj_far, eps_ij, R_ij
		real(prcn) :: dt_, tmp_tke, tmp_real, xi, eta1
		real(prcn), save :: tke_old, tke_far
		real(prcn) :: alpha, epsf0, meanf_energy

		real(prcn), dimension(ndim,ndim) :: tmp_tensor
		integer :: dim1, dim2, dim3
		integer :: i, ii, j, jj, k, kk, ip, im, jp, jm, kp, km, idim, iii, jjj, kkk
		logical :: filexist, isopen, new_file
		character*30 filename1, filename2, filename3

		logical :: neighb_insolid

		if (ntke<1) ntke=1
		if (call_flag==0) then
			ur=>turb_ur
			uf=>turb_uf
		else
			if (irestart==0.and.iglobstep==1) then
			elseif (irestart==1.and.iglobstep==1) then
			elseif (mod(iglobstep,ntke)/=0) then
				return
			endif

			ur=>ubcp
			uf=>u
		endif
!		if (near_particle_checked==.false.) then
!			call near_particle_region
!			near_particle_checked = .true.
!		endif


		!^^^^^^^ Computing RSM ^^^^^^^^^^^^^^
		tke_old = tke
		uiuj = 0d0
		uiuj_far = 0d0
		do k=1, mz
			do j=1, my
				do i=1, nx
					if (fluid_atijk(i,j,k)) then
						im = i-1
						ip = i+1
#if !PARALLEL
						if(im<1)  im = nx
						if(ip>nx) ip = 1
#endif               
						do dim1=1, ndim
							do dim2=1, ndim
								uiuj(dim1,dim2) = uiuj(dim1,dim2) + (ur(i,j,k,dim1)-ufmean(dim1)) * (ur(i,j,k,dim2)-ufmean(dim2))
							enddo
						enddo

!						if (far_field(i,j,k)) then
!							do dim1=1, ndim
!								do dim2=1, ndim
!									uiuj_far(dim1,dim2) = uiuj_far(dim1,dim2) + (ur(i,j,k,dim1)-ufmean(dim1)) * (ur(i,j,k,dim2)-ufmean(dim2))
!								enddo
!							enddo
!						endif
					endif
				end do
			end do
		end do

#if PARALLEL
		tmp_tensor = 0d0
		GLOBAL_DOUBLE_SUM(uiuj, tmp_tensor, 9, decomp_group)
		uiuj = tmp_tensor

		tmp_tensor = 0d0
		GLOBAL_DOUBLE_SUM(uiuj_far, tmp_tensor, 9, decomp_group)
		uiuj_far = tmp_tensor
#endif
		uiuj = uiuj / count_fluid
		uiuj_far = uiuj_far / count_fluid

		tke = (uiuj(1,1) + uiuj(2,2) + uiuj(3,3)) / 2
		tke_far = (uiuj_far(1,1) + uiuj_far(2,2) + uiuj_far(3,3)) / 2

		tke_converged = .false.
		if (iglobstep>10.and.abs(tke_old-tke)/tke<tke_error) tke_converged = .true.
		!------------------------------------

		!^^^^^^^ Computing EPSij & Rij^^^^^^^^^^^^^^
		eps_ij = 0d0
#if 0
		call calc_pressure
		do dim1=1, 3
			do dim2=1, 3
				do dim3=1, 3
					do i=1, nx
						im = i-1
#if !PARALLEL
						if (im<1) im = nx
#endif
					
						call derivative2(dim1,dim3,ur1)
						call derivative2(dim2,dim3,ur2)

						do k=1, mz
							do j=1, my
								eps_ij(dim1,dim2) = eps_ij(dim1,dim2) + 2 * vis * (ur1(j,k) * ur2(j,k))
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo

		R_ij   = 0d0
		do dim1=1, 3
			do dim2=1, 3
				do i=1, nx
					im = i-1
#if !PARALLEL
					if (im<1) im = nx
#endif
					call derivative2(dim1,dim2,ur1)
					call derivative2(dim2,dim1,ur2)

					do k=1, mz
						do j=1, my
							R_ij(dim1,dim2) = R_ij(dim1,dim2) + pbcp(i,j,k) * (ur1(j,k)+ur2(j,k))
!							R_ij(dim1,dim2) = R_ij(dim1,dim2) + 0.5*(pbcp(i,j,k)+pbcp(im,j,k)) * (ur1(j,k)+ur2(j,k))
						enddo
					enddo
				enddo
			enddo
		enddo

#if PARALLEL
		tmp_tensor = 0d0
		GLOBAL_DOUBLE_SUM(eps_ij, tmp_tensor, 9, decomp_group)
		eps_ij = tmp_tensor

		tmp_tensor = 0d0
		GLOBAL_DOUBLE_SUM(R_ij, tmp_tensor, 9, decomp_group)
		R_ij = tmp_tensor
#endif
		eps_ij(:,:) = eps_ij(:,:) / count_fluid
		R_ij(:,:)   = R_ij(:,:) / count_fluid

		epsf   = 0.5 * (eps_ij(1,1) + eps_ij(2,2) + eps_ij(3,3))
#endif
	!----------------------------------------------

!		!^^^^^^^ Computing EPSF ^^^^^^^^^^^^^^
!		epsf0   = 0d0
!		do dim1=1, 3
!			do dim2=1, 3
!				do i=1, nx
!					call derivative(dim1,dim2,ur1)
!					call derivative(dim2,dim1,ur2)
!
!					do k=1, mz
!						do j=1, my
!							epsf0 = epsf0 + 2 * vis * (0.5*(ur1(j,k)+ur2(j,k)))**2
!						enddo
!					enddo
!				enddo
!			enddo
!		enddo
!#if PARALLEL
!		tmp_real = 0d0
!		GLOBAL_DOUBLE_SUM(epsf0, tmp_real, 1, decomp_group)
!		epsf0 = tmp_real
!#endif
!		epsf0  = epsf0 / count_fluid
!	!----------------------------------------------

		meanf_energy = half*dot_product(ufmean(:),ufmean(:))
		if(I_AM_NODE_ZERO)then
			call calc_anisotropy(uiuj, xi, eta1) 

			CALL screen_separator(40,'K')
			write (*,'(A,5D15.7)') "T/T_CONV, TKE/MEAN_energy, FAR_FRAC, XI, ETA =", t*umeanslip/dia_phys, tke/meanf_energy, tke_far/tke, xi ,eta1
			write (*,'(A,2D15.7)') "TKE_ERROR: ", abs(tke_old-tke)/tke, tke_error
			CALL screen_separator(40,'K')

			filename1 = trim(run_name)//"_tke_history.dat"
			filename2 = trim(run_name)//"_uiuj_history.dat"
			filename3 = trim(run_name)//"_lumley.dat"

			new_file= .false.
			if (irestart==0.and.iglobstep==1) then
				new_file = .true.
			else
				inquire(file=trim(filename1),exist=filexist,opened=isopen)
				if (.not.filexist) then
					new_file= .true.
				endif
			endif

			if (new_file) then
				open (unit=99900,file=trim(filename1),status="replace")
!				write (99900,"(A)") "variables=t,t_conv, TKE, TKE*, FAR_FRAC, TKE_ERROR"
!				write (99900,"(A)") "variables=t,t_conv, TKE, TKE*, TKE_error"
!									& -(prod-epsf-dTidXi)
				write (99900,*) "zone"

				open (unit=99901,file=trim(filename2),status="replace")
!				write (99901,"(A)") "variables=t,u11,u22,u33,u12,u13,u23,du11,du22,du33,du12,du13, du23"
				write (99901,*) "zone"

				open (unit=99903,file=trim(filename3),status="replace")
!				write (99903,"(A)") "variables=t,xi,eta"
				write (99903,*) "zone"
			else
				open(unit=99900,file=trim(filename1),status='old',position="append")
				open(unit=99901,file=trim(filename2),status='old',position="append")
				open(unit=99903,file=trim(filename3),status='old',position="append")
			endif

			write (99900,'(6D15.7)') t, t*dia_phys/umeanslip, tke, tke/meanf_energy, tke_far/tke, abs(tke_old-tke)/tke

			if (iglobstep>2) then
				write (99901,"(1D15.7)") t*dia_phys/umeanslip
				write (99901,"(6D15.7)") uiuj(1,1)/umeanslip, uiuj(2,2)/umeanslip, uiuj(3,3)/umeanslip, uiuj(1,2)/umeanslip, uiuj(1,3)/umeanslip, uiuj(2,3)/umeanslip
!				write (99901,"(6D15.7)") dui_uj(1,1), dui_uj(2,2), dui_uj(3,3), dui_uj(1,2), dui_uj(1,3), dui_uj(2,3)
!				write (99901,"(6D15.7)") R_ij(1,1), R_ij(2,2), R_ij(3,3), R_ij(1,2), R_ij(1,3), R_ij(2,3)
!				write (99901,"(6D15.7)") -eps_ij(1,1), -eps_ij(2,2), -eps_ij(3,3), -eps_ij(1,2), -eps_ij(1,3), -eps_ij(2,3)
			endif

			write (99903,'(3D15.7)') t*dia_phys/umeanslip, xi, eta1

			close (99900)
			close (99901)
			close (99903)
		endif
	end subroutine tke_eps
#endif
end module init_turb

SUBROUTINE histogram1(u,wt,n,nddu,bldd,brdd,hist)
  !nddu: number of bins for pdf formation
  !bldd: lefh hand side limit .... output
  !brdd: right side limie.... otuput
  USE precision
  USE constants
  USE randomno
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n , nddu
  REAL(prcn), INTENT(in), DIMENSION(n) :: u, wt
  REAL(prcn), DIMENSION(:), ALLOCATABLE :: u_t
  
  REAL(prcn),  INTENT(out) :: bldd, brdd 

  
  REAL(prcn),  INTENT(inout), DIMENSION(nddu) ::  hist
    REAL(prcn) ::  xdiff
    
    REAL(prcn) :: vlmt, ave, adev, sdev, var, skew, curt
    
    INTEGER :: i, ibin
    bldd= 1.e25
    brdd=-1.e25

    ALLOCATE(u_t(n))
    CALL moment1(4, n, n, wt, u, ave,adev,sdev,var,skew&
         &,curt)
!    PRINT*,'ave,var,sdev, skew, curt=', ave, var,sdev, skew,curt
!    WRITE(*,*)'number of bins in hist..',nddu

  DO i=1,n
     u_t(i)=(u(i)-ave)/sdev
  ENDDO

  CALL moment1(4, n, n, wt, u_t, ave,adev,sdev,var,skew&
       &,curt)
!!$  PRINT*,'ave,var, sdev,skew, curt=', ave, var,sdev, skew,curt

    DO i=1,nddu
       hist(i)=0.0  
    ENDDO
    bldd = MIN(MINVAL(u_t(:)), bldd)
    brdd = MAX(MAXVAL(u_t(:)), brdd)
    xdiff = (brdd - bldd)/float(nddu-1)

    DO i=1,n
       ibin = (u_t(i) - bldd)/xdiff + 1
       hist(ibin)=hist(ibin) + wt(i)/xdiff
    ENDDO

    DEALLOCATE(u_t)
END SUBROUTINE histogram1


SUBROUTINE plothist(hist,lb,ub,nbins,iunit,t,tref)
  USE precision
  USE constants
  IMPLICIT NONE
  REAL(prcn), INTENT(in), DIMENSION(nbins) ::hist
  
  INTEGER, INTENT(in) ::iunit,nbins
  REAL(prcn) ::  sum_x,lb,ub,dx, t, tref, tmp, tmp2
  INTEGER(prcn) :: i
  
  sum_x = lb
  dx= (ub-lb)/(float(nbins)-1)
  tmp = one/sqrt(twopi)
  WRITE(iunit,*)'Zone t="',int(t/tref),'"'
  DO i=1,nbins
     tmp2  = sum_x+dx/2
     WRITE(iunit,*)tmp2,hist(i)!, tmp*exp(-(tmp2*tmp2)/two)
     sum_x = sum_x + dx
  ENDDO
  
  RETURN
END SUBROUTINE plothist


real(prcn) function energy_init(kmode,epsf,l_e,eta,c_eta,c_l)
	use precision

	real(prcn), intent(in) :: kmode, epsf, l_e, eta, c_eta, c_l
	real(prcn), parameter  :: C = 1.5, beta = 5.2, p0 = 2.0

!	c_eta = 0.4,c_L = 6.78,

	energy_init = C * epsf**(2./3.) * kmode**(-5./3.)
	energy_init = energy_init*(kmode * l_e/sqrt((kmode*l_e)**2 +&
								& c_L))**(5./3. + p0)
	energy_init = energy_init*exp(-1.d0*beta*(((kmode*eta)**4 + c_eta**4)&
								&**0.25 - c_eta))
	return   
end function energy_init

real(prcn) function integrate_energy(kmin,kmax,eps,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) :: kmin,kmax,eps,l_e,eta,c_eta,c_l
	
	real(prcn), external :: energy_init
	real(prcn) :: dk, umag0, umag, kmode
	integer :: i,ikrange

	dk = 0.001
	ikrange = int((kmax-kmin)/dk)

	integrate_energy  = 0d0
	umag0 = energy_init(kmin,eps,l_e,eta,c_eta,c_l)
	do i=1,ikrange
		kmode            = kmin +i*dk
		umag             = energy_init(kmode,eps,l_e,eta,c_eta,c_l)
		integrate_energy = integrate_energy + (umag+umag0)*dk/2
		umag0            = umag
    end do
end function integrate_energy

real(prcn) function integrate_epsilon(kmin,kmax,vis,epsf,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) :: kmin,kmax,vis,epsf,l_e,eta,c_eta,c_l
	
	real(prcn), external :: energy_init
	real(prcn)	:: dk, umag0, umag, kmode
	integer :: i,ikrange

	dk = 0.001
	ikrange = int((kmax-kmin)/dk)

	integrate_epsilon  = 0.d0
	umag0 = 2 * vis * kmin**2 * energy_init(kmin,epsf,l_e,eta,c_eta,c_l)
	do i=1,ikrange
		kmode             = kmin + i * dk
		umag              = 2 * vis * kmode**2 * energy_init(kmode,epsf,l_e,eta,c_eta,c_l)
		integrate_epsilon = integrate_epsilon + (umag+umag0) * dk / 2
		umag0             = umag
    end do
end function integrate_epsilon




!^^^^^^^^^^^ Subroutines related to the determination of C_eta and C_l
subroutine xmnewt(c_eta,c_l,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	real(prcn),intent(inout) :: c_eta,c_l
	real(prcn),intent(in)    :: tke,eps,l_e,eta,kmax,vis
	INTEGER NTRIAL,N,NP
	real(prcn) TOLX,TOLF
	PARAMETER(NTRIAL=50,TOLX=1.0D-10,N=2,TOLF=1.0D-10,NP=2)
	INTEGER i,j,k,kk
	real(prcn) xx,fjac(NP,NP),fvec(NP),x(NP)

	x(1)=c_eta
	x(2)=c_l

	call usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
!	write(*,'(/1x,t5,a,t14,a,t29,a/)') 'I','X(I)','F'
!	do i=1,N
!		write(*,'(1x,i4,2e15.6)') i,x(i),fvec(i)
!	enddo

	do j=1,NTRIAL
		call mnewt(1,x,N,TOLX,TOLF,tke,eps,l_e,eta,kmax,vis)
		call usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
		write(*,'(/1x,t5,a,t14,a,t29,a/)') 'I','X(I)','F'
		do i=1,N
			write(*,'(1x,i4,2e15.6)') i,x(i),fvec(i)
		enddo
		if (abs(fvec(1))<1D-10.and.abs(fvec(2))<1D-10) exit
	enddo
	c_eta = x(1)
	c_l   = x(2)
end subroutine xmnewt


SUBROUTINE usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER i,n
	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) fjac(n,n),fvec(n),x(n)

	call funcv(n,x,fvec,tke,eps,l_e,eta,kmax,vis)
	call fdjac(n,x,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
END

SUBROUTINE funcv(n,x,f,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER n
	real(prcn) x(n),f(n)

	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) :: c_eta,c_l
	real(prcn) :: krange,kmode,dk,sum1,sum2,umag0,umag
	integer :: i,ikrange
	real(prcn), external  :: energy_init,integrate_energy,integrate_epsilon

	c_eta = x(1)
	c_l   = x(2)
	
	!tke: integral of E(k): 0-kmax
	sum1=integrate_energy(1d-6,kmax,eps,l_e,eta,c_eta,c_l)
	f(1)=tke-sum1

	!eps: integratal of 2*nu*k^2*E(k)*dk (0-kmax)
	sum2=integrate_epsilon(1d-6,kmax,vis,eps,l_e,eta,c_eta,c_l)
	f(2)=eps-sum2
	return
END

SUBROUTINE fdjac(n,x,fvec,df,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER n
	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) df(n,n),fvec(n),x(n)
	real(prcn) small
	PARAMETER (small=1.e-4)
	INTEGER i,j
	real(prcn) h,temp,f(n)

	do j=1,n
		temp=x(j)
		h=small*abs(temp)
		if(h.eq.0.)h=small
		x(j)=temp+h
		h=x(j)-temp
		call funcv(n,x,f,tke,eps,l_e,eta,kmax,vis)
		x(j)=temp
		do i=1,n
			df(i,j)=(f(i)-fvec(i))/h
		enddo
	enddo
	return
END

SUBROUTINE mnewt(ntrial,x,n,tolx,tolf,tke,eps,l_e,eta,kmax,vis)
	use precision
	implicit none
	INTEGER n,ntrial
	real(prcn),intent(in) :: tke,eps,l_e,eta,kmax,vis
	real(prcn) tolf,tolx,x(n)
	INTEGER i,k,indx(n)
	real(prcn) d,errf,errx,fjac(n,n),fvec(n),p(n)

	do k=1,ntrial
		call usrfun(x,n,fvec,fjac,tke,eps,l_e,eta,kmax,vis)
		errf=0.
		do i=1,n
			errf=errf+abs(fvec(i))
		enddo
		if(errf.le.tolf)return
		do i=1,n
			p(i)=-fvec(i)
		enddo
		call ludcmp(fjac,n,indx,d)
		call lubksb(fjac,n,indx,p)
		errx=0.
		do i=1,n
			errx=errx+abs(p(i))
			if (x(i)+p(i)>0) then
				x(i)=x(i)+p(i)
			else
				x(i)=x(i)*0.5
			endif
		enddo
		if(errx.le.tolx)return
	enddo
	return
END


SUBROUTINE lubksb(a,n,indx,b)
	use precision
	implicit none
	INTEGER n
	integer indx(n)
	real(prcn) a(n,n),b(n)
	INTEGER i,ii,j,ll
	real(prcn) sum
	ii=0
	do i=1,n
		ll=indx(i)
		sum=b(ll)
		b(ll)=b(i)
		if (ii.ne.0)then
			do j=ii,i-1
				sum=sum-a(i,j)*b(j)
			enddo
		else if (sum.ne.0.) then
			ii=i
		endif
		b(i)=sum
	enddo
	do i=n,1,-1
		sum=b(i)
		do j=i+1,n
			sum=sum-a(i,j)*b(j)
		enddo
		b(i)=sum/a(i,i)
	enddo
	return
END

SUBROUTINE ludcmp(a,n,indx,d)
	use precision
	implicit none
	INTEGER n
	integer indx(n),NMAX
	real(prcn) d,a(n,n),TINY
	PARAMETER (NMAX=500,TINY=1.0e-20)
	INTEGER i,imax,j,k
	real(prcn) aamax,dum,sum,vv(n)
	d=1.
	do i=1,n
		aamax=0.
		do j=1,n
			if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
		enddo
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
	enddo
	do j=1,n
		do i=1,j-1
			sum=a(i,j)
			do k=1,i-1
				sum=sum-a(i,k)*a(k,j)
			enddo
			a(i,j)=sum
		enddo
		aamax=0.
		do i=j,n
			sum=a(i,j)
			do k=1,j-1
				sum=sum-a(i,k)*a(k,j)
			enddo
			a(i,j)=sum
			dum=vv(i)*abs(sum)
			if (dum.ge.aamax) then
				imax=i
				aamax=dum
			endif
		enddo
		if (j.ne.imax)then
			do k=1,n
				dum=a(imax,k)
				a(imax,k)=a(j,k)
				a(j,k)=dum
			enddo
			d=-d
			vv(imax)=vv(j)
		endif
		indx(j)=imax
		if(a(j,j).eq.0.)a(j,j)=TINY
		if(j.ne.n)then
			dum=1./a(j,j)
			do i=j+1,n
				a(i,j)=a(i,j)*dum
			enddo
		endif
	enddo
	return
END
!--------------------------------------------------

!^^^^^^^^^^^ Subroutines related to the mean value theorem
subroutine xrtnewt(root,k0,kmax,epsf,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) epsf,l_e,eta,c_eta,c_l
	INTEGER N,NBMAX
	real(prcn):: X1,X2
	PARAMETER(N=100,NBMAX=1)
	INTEGER i,nb
	real(prcn):: integrate_energy,f,rtnewt,root,xacc,xb1(NBMAX),xb2(NBMAX)
	EXTERNAL funcd,f,integrate_energy

	real(prcn) :: dk,sum,umag0,umag,kmode,k0,kmax,const
	integer :: ikrange

	const = integrate_energy(k0,kmax,epsf,l_e,eta,c_eta,c_l)
	const = const / (kmax-k0)

	X1=k0
	X2=kmax
	nb=NBMAX
	call zbrak(f,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,X1,X2,N,xb1,xb2,nb)
!	write(*,'(/1x,a)') 'Roots of BESSJ0:'
!	write(*,'(/1x,t19,a,t31,a/)') 'x','F(x)'
	do i=1,nb
		xacc=(1.0e-6)
		root=rtnewt(funcd,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,xb1(i),xb2(i),xacc)
!		write(*,'(1x,a,i2,2x,f12.6,e16.4)') 'Root ',i,root,f(root,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	enddo
END

function f(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	use precision
	implicit none
	real(prcn) :: x,f
	real(prcn) :: k0,kmax,const,tke,epsf,l_e,eta,c_eta,c_l
	real(prcn),external :: energy_init

	f=energy_init(x,epsf,l_e,eta,c_eta,c_l)-const
end

function df(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,f)
	use precision
	implicit none
	real(prcn)  :: x,f,df
	real(prcn)  :: k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn),external  :: fdjac1
	
!	df=2*x
	df=fdjac1(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,f)
end

function fdjac1(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,fvec)
	use precision
	implicit none
	INTEGER n
	real(prcn) fdjac1
	real(prcn) x,fvec,df
	real(prcn) k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn) small
	PARAMETER (small=1.e-4)
	INTEGER i,j
	real(prcn) h,temp
	real(prcn),external :: f

	temp=x
	h=small*abs(temp)
	if(h.eq.0.)h=small
	x=temp+h
	h=x-temp

	fdjac1=(f(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)-fvec)/h
	x=temp
	return
END

SUBROUTINE funcd(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,fn,dfn)
	use precision
	implicit none
	real(prcn)  :: fn,dfn,x
	real(prcn)  :: k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn),external :: f,df

	fn  =  f(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	dfn = df(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,fn)
	return
END

FUNCTION rtnewt(funcd,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,x1,x2,xacc)
	use precision
	implicit none
	INTEGER JMAX
	real(prcn):: rtnewt,x1,x2,xacc
	real(prcn)::k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	EXTERNAL funcd
	PARAMETER (JMAX=200)
	INTEGER j
	real(prcn):: df,dx,f
	rtnewt=.5*(x1+x2)
	do j=1,JMAX
		call funcd(rtnewt,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,f,df)
		dx=f/df
		rtnewt=rtnewt-dx
!		if((x1-rtnewt)*(rtnewt-x2).lt.0.) pause
!		if((x1-rtnewt)*(rtnewt-x2).lt.0.) write (*,*) 'rtnewt jumped out of brackets'
		if(abs(dx).lt.xacc) return
	enddo
!	write (*,*) 'rtnewt exceeded maximum iterations'
!	pause 
END

SUBROUTINE zbrak(fx,k0,kmax,const,epsf,l_e,eta,c_eta,c_l,x1,x2,n,xb1,xb2,nb)
	use precision
	implicit none
	INTEGER n,nb
	real(prcn):: x1,x2,xb1(nb),xb2(nb)
	real(prcn) :: k0,kmax,const,epsf,l_e,eta,c_eta,c_l
	real(prcn),EXTERNAL :: fx
	INTEGER i,nbb
	real(prcn):: dx,fc,fp,x
	nbb=0
	x=x1
	dx=(x2-x1)/n
	fp=fx(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
	do i=1,n
		x=x+dx
		fc=fx(x,k0,kmax,const,epsf,l_e,eta,c_eta,c_l)
		if(fc*fp.le.0.) then
			nbb=nbb+1
			xb1(nbb)=x-dx
			xb2(nbb)=x
			if(nbb.eq.nb) exit
		endif
		fp=fc
	enddo
	nb=nbb
	return
END



