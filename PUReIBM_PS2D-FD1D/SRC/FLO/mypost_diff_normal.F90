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

module mypost_process
#include "ibm.h"
	use global_data, only : nx, mx1, my, my2, dx, dy, dz, u,&
						&  vis, ndim, run_name, doml, &
						&  maxvolfrac, re, lybyd, &
						&  dia_phys, fluid_atijk, input_type, aliasflag, &
						&  run_name, doml, dia_phys, &
						&  input_type, iglobstep, t, minunitno, &
						&  maxunitno, count_fluid, irestart, rhof, wy, wz, &
						&  usmean, ufmean, pres_avg, visc_avg, bndarray
#if PARALLEL
	use global_data, only : err_code, devisc_comp_group, myid, starts, ends
#endif
	use nlmainarrays , only : ubcp, pbcp
	use fftw_interface
	use constants
	use nlarrays
	use randomno
	use general_funcs
	use string_funcs
	use init_turb , only : avr, tke, iauto_on, autocorrelation
	use boundary_condition
	use dependent_functions
	use bcsetarrays

contains

subroutine post
	implicit none

	type :: phi_data
		real(8) :: ph
		real(8) :: eps_visc, eps_pres
		real(8), dimension(:,:), pointer :: visc_comp
		real(8), dimension(:), pointer:: pres_comp

		real(8), dimension(3) :: pgrad,vel
		real(8), dimension(3) :: tdiv
		real(8) :: pres
	endtype phi_data

	type :: theta_data
		real(8) :: th
		integer :: nphi
		type (phi_data), dimension(:), allocatable :: phi
	endtype theta_data
		
	type :: surface_stress
		integer :: nth
		type (theta_data), dimension(:), pointer :: theta
	endtype surface_stress

	type(surface_stress) :: surf, out1, out2

	type(theta_data), pointer :: layer_m, layer_c, layer_p

	real(8), allocatable :: ugrad(:,:,:,:,:)
	integer, allocatable, dimension(:) :: thbin_loc, thbin
	real(8), allocatable, dimension(:) :: eps_th_loc, eps_th
	real(8), allocatable, dimension(:) :: visc_th_loc, visc_th, pres_th_loc, pres_th
	real(8), allocatable, dimension(:,:) :: ugp_th_loc, ugp_th, eps_visc_phi, eps_pres_phi

	real(8), pointer :: surf_tens(:,:), out1visc_comp(:,:), out2visc_comp(:,:), visc_comp(:,:)
	real(8), pointer :: out1pres_comp, out2pres_comp, pres_comp
!	real(8), dimension(:,:), pointer :: bndarray

	real(8) :: tmp1, tmp2, rad, r, ro, ro2, x, y, z, dth, th, dphi, phi, dist
	integer :: ith, nth, thp, thm, iphi, nphi, phim, phip
	integer :: iphs, nbnd, idim
	integer :: i, j, k, l, m, n, dim1, dim2

	integer :: ib, ie, jb, je, kb, ke, onew, ii, jj, kk
	
	integer, dimension(ndim) :: is, io, io2
	integer, dimension(ndim) :: vcellb, vcello, vcello2
	integer, dimension(ndim) :: pcellb, pcello, pcello2
	real(8), dimension(ndim) :: xl, xlo, xlo2
	real(8), dimension(ndim) :: xlp, xlop, xlo2p
	real(8), dimension(ndim) :: ul, ulo, ulo2, ul_p, ulo_p, ulo2_p
	real(8), dimension(ndim) :: pgl, pglo, pglo2
	real(8) :: pl, plo, plo2
	real(8) :: rdummy(ndim), umagnorm, umagnormo, umagnormo2, viscm(ndim,ndim), viscp(ndim,ndim), presp, presm, div(ndim), gradp(ndim)
	real(8) :: ddr, ddth, ddphi, phi_star

	real(8), pointer :: pr(:,:,:)

	logical :: velterm


	real(8), allocatable :: test(:,:,:)


	allocate(ugrad(nx,my,mz,ndim,ndim))
	

	do idim=1, ndim
		do i=1, nx
			call ff2cr(u(i,:,:,idim),ubcp(i,:,:,idim))
			ubcp(i,:,:,idim) = ubcp(i,:,:,idim) + umean(idim)
		enddo
	enddo
	call transfer_velocity

	!GENERATING VELOCITY GRADIENT
	if (I_AM_NODE_ZERO) write (*,*) "CALCULATING UGRAD"
	call calc_ugrad
	call calc_pressure
	pr => pbcp

	tmp1 = 0.5*phase_array(nphases)%dia/doml(2)*my
	nth  = nint(twopi*tmp1*f2/2)

	allocate(ugp_th_loc(0:nth,ndim),ugp_th(0:nth,ndim))
	allocate(thbin_loc(0:nth),thbin(0:nth))
!	allocate(eps_th_loc(0:nth),eps_th(0:nth))
	allocate(visc_th_loc(0:nth),visc_th(0:nth))
	allocate(pres_th_loc(0:nth),pres_th(0:nth))
	allocate(eps_visc_phi(0:nth,0:7), eps_pres_phi(0:nth,0:7))
	
	thbin_loc=0
	thbin=0

	ugp_th = 0d0
	ugp_th_loc = 0d0

!	eps_th_loc=0d0
!	eps_th=0d0

	visc_th = 0d0
	visc_th_loc = 0d0

	pres_th = 0d0
	pres_th_loc = 0d0

	eps_visc_phi = 0d0
	eps_pres_phi = 0d0

	do m=1, nbody
		iphs = 1
		nbnd = phase_array(iphs)%nbnd
		bndarray => phase_array(iphs)%bndpts

		dr = dr * dx		
		r   = radbdy(m) * dx
		ro  = r+dr	! * dx
		ro2 = r+two*dr !* dx

		write (*,*) r, ro, ro2
		write (*,*) radbdy(m), radobdy(m), rado2bdy(m)
		read (*,*)
	
		do l=1,nbnd
			tmp1 = 0.5*phase_array(nphases)%dia/doml(2)*my
			nth  = nint(twopi*tmp1*f2/4)*2

			x = bndarray(1,l)
			y = bndarray(2,l)
			z = bndarray(3,l)

			rad = sqrt(x**2+y**2+z**2)

			th  = acos(abs(x)/rad)
			if (x>0) th = pi-th
!			th  = acos(abs(z)/rad)
!			if (z<0) th = pi-th

			dth = twopi/2/nth
			ith = nint(th/dth)
			th  = dth*ith

			if (th<0.or.th>pi) write (*,*) "ERROR IN THETA"
			if (ith>nth)       write (*,*) "ERROR IN NTHETA"

			nphi=nint(twopi*tmp1*sin(th)*f2/4)*4

			if (abs(sin(th))>1d-4) then
				phi = acos(abs(y/rad/sin(th)))
				if (abs(y/rad/sin(th))>1d0) phi=0d0
				if (y>=0d0.and.z>=0d0) then
				elseif (y<0d0.and.z>=0d0) then
					phi = pi - phi
				elseif (y<0d0.and.z<0d0) then
					phi = pi + phi
				elseif (y>=0d0.and.z<0) then
					phi = twopi - phi
				endif

!				phi = acos(abs(x/rad/sin(th)))
!				if (abs(x/rad/sin(th))>1d0) phi=0d0
!
!				if (x>=0d0.and.y>=0d0) then
!				elseif (x<0d0.and.y>=0d0) then
!					phi = pi - phi
!				elseif (x<0d0.and.y<0d0) then
!					phi = pi + phi
!				elseif (x>=0d0.and.y<0) then
!					phi = twopi - phi
!				endif

				dphi = twopi/nphi
				iphi = nint(phi/dphi)
!				if (iphi>=nphi) iphi = 0
				phi  = dphi*iphi
			else
				phi  = 0
				nphi = 1
				iphi = 0
			endif


!			if (abs(th * 180d0/pi-176.4)<1) then
			if (ith==49) then
				write (*,"(4D20.12,4I6)") th * 180d0/pi, dth * 180d0/pi, phi * 180d0/pi, dphi * 180d0/pi, iphi, nphi, ith, nth
			endif

			if (phi<0.or.phi>twopi) write (*,*) "ERROR IN PHI"
			if (iphi>=nphi) then
				write (*,*) "ERROR IN NPHI"
				write (*,*) ith, iphi, nphi
				stop
			endif

!			goto 10


			if (.not.associated(surf%theta)) then
!				write (*,*) "ALLOCATING THETA ARRAY", nth
				surf%nth = nth
!				surf%dth = dth
				allocate(surf%theta(0:nth))

				out1%nth = nth
!				out1%dth = dth
				allocate(out1%theta(0:nth))

				out2%nth = nth
!				out2%dth = dth
				allocate(out2%theta(0:nth))
			endif 

			if (.not.allocated(surf%theta(ith)%phi)) then
!				write (*,*) "ALLOCATING PHI ARRAY", ith, nphi
				surf%theta(ith)%nphi = nphi
!				surf%theta(ith)%dphi = dphi
				allocate(surf%theta(ith)%phi(0:nphi-1))

				out1%theta(ith)%nphi = nphi
!				out1%theta(ith)%dphi = dphi
				allocate(out1%theta(ith)%phi(0:nphi-1))

				out2%theta(ith)%nphi = nphi
!				out2%theta(ith)%dphi = dphi
				allocate(out2%theta(ith)%phi(0:nphi-1))

				do i=0, nphi-1
					allocate(surf%theta(ith)%phi(i)%visc_comp(3,3), surf%theta(ith)%phi(i)%pres_comp(1))
					allocate(out1%theta(ith)%phi(i)%visc_comp(3,3), out1%theta(ith)%phi(i)%pres_comp(1))
					allocate(out2%theta(ith)%phi(i)%visc_comp(3,3), out2%theta(ith)%phi(i)%pres_comp(1))
				enddo
			endif

			do n=1,ndim
				!     location of surface points
				xl(n)=xc(m,n)+ bndarray(n,l)*radbdy(m)
				is(n)=INT(xl(n))
				ul(n)=zero

				!     location of r+dr points
				xlo(n)=xc(m,n)+ bndarray(n,l)*radobdy(m)
				io(n)=INT(xlo(n))
				ulo(n)=zero

				!     location of r+2*dr points
				xlo2(n)=xc(m,n)+ bndarray(n,l)*rado2bdy(m)
				io2(n)=INT(xlo2(n))
				ulo2(n)=zero

				rdummy(n) = zero
!				onll(n)=zero
!				nll(n)=zero
!				ppll(n)=zero
!				dfll(n)=zero
			enddo

			xlp(1)   = xl(1)-0.5
			xlp(2:3) = xl(2:3)

			xlop(1)   = xlo(1)-0.5
			xlop(2:3) = xlo(2:3)

			xlo2p(1)   = xlo2(1)-0.5
			xlo2p(2:3) = xlo2(2:3)


			do n=1, ndim
				if (xl(n).lt.zero) then 
					vcellb(n) = int(xl(n)-1)
				else 
					vcellb(n) = int(xl(n))
				endif
				if (xlp(n).lt.zero) then 
					pcellb(n) = int(xlp(n)-1)
				else 
					pcellb(n) = int(xlp(n))
				end if

				if(xlo(n).lt.zero) then 
					vcello(n) = int(xlo(n)-1)
				else 
					vcello(n) = int(xlo(n))
				endif
				if (xlop(n).lt.zero) then 
					pcello(n) = int(xlop(n)-1)
				else 
					pcello(n) = int(xlop(n))
				end if

				if(xlo2(n).lt.zero) then 
					vcello2(n) = int(xlo2(n)-1)
				else 
					vcello2(n) = int(xlo2(n))
				endif
				if(xlo2p(n).lt.zero) then 
					pcello2(n) = int(xlo2p(n)-1)
				else 
					pcello2(n) = int(xlo2p(n))
				end if
			enddo

#if PARALLEL
			! R POINTS
			xltemp = xl(1)
			vcelltemp = vcellb(1)
			if(l.eq.FOCUS_POINT)then
				PRINT*,' xl = ', myid, xltemp
			end if
			if(.not.CELL_IN_VEL_GRID(vcelltemp))then
				WEST_PERIODIC_IMAGE(vcellb(1),vcelltemp,xl(1),xltemp)
				EAST_PERIODIC_IMAGE(vcellb(1),vcelltemp,xl(1),xltemp)
				if(l.eq.FOCUS_POINT)then
					PRINT*,' xl IMAGES = ', myid,xltemp
				end if

				if(.not.CELL_IN_VEL_GRID(vcelltemp)) goto 777
			end if
			vcellb(1) = vcelltemp
			xl(1) = xltemp

			! R+DR POINTS
			xltemp = xlo(1)
			vcelltemp = vcello(1)
			if(l.eq.FOCUS_POINT)then
				PRINT*,' xl = ', myid, xltemp
			endif
			if(.not.CELL_IN_VEL_GRID(vcelltemp))then
				WEST_PERIODIC_IMAGE(vcello(1),vcelltemp,xlo(1),xltemp)
				EAST_PERIODIC_IMAGE(vcello(1),vcelltemp,xlo(1),xltemp)
				if(l.eq.FOCUS_POINT)then
					PRINT*,' xl IMAGES = ', myid,xltemp
				end if

				if(.not.CELL_IN_VEL_GRID(vcelltemp)) goto 777
			end if
			vcello(1) = vcelltemp
			xo(1) = xltemp

			! R+2 * DR POINTS
			xltemp = xlo2(1)
			vcelltemp = vcello2(1)
			if(l.eq.FOCUS_POINT)then
				PRINT*,' xl = ', myid, xltemp
				end if
			if(.not.CELL_IN_VEL_GRID(vcelltemp))then
				WEST_PERIODIC_IMAGE(vcello2(1),vcelltemp,xlo2(1),xltemp)
				EAST_PERIODIC_IMAGE(vcello2(1),vcelltemp,xlo2(1),xltemp)
				if(l.eq.FOCUS_POINT)then
					PRINT*,' xl IMAGES = ', myid,xltemp
				end if

				if(.not.CELL_IN_VEL_GRID(vcelltemp)) goto 777
			end if
			vcello2(1) = vcelltemp
			xo2(1) = xltemp

!		! INNER POINTS
!		vcelltemp = vcelli(1)
!		xltemp  = xli(1)
!		if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
!			PRINT*,' INNER REV PT = ', xli(1),vcelltemp, myid, m
!			PRINT*,' EXTERNAL REV PT = ', xlo(1), myid,m
!		end if
!
!		if(.not.CELL_IN_PROC(vcelltemp))then
!			WEST_PERIODIC_IMAGE(vcelli(1),vcelltemp,xli(1),xltemp)
!			EAST_PERIODIC_IMAGE_MOD(vcelli(1),vcelltemp, xli(1), xltemp)
!			if(.not.CELL_IN_PROC(vcelltemp))then
!				if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
!					PRINT*,' INNER REVERSAL PT = ', xli(1),vcelltemp, myid, m
!				end if
!				goto 777
!			end if
!		end if
!		if(EAST_NO_MANS_LAND(vcelli(1)).or.EAST_NO_MANS_LAND(vcelltemp)) then 
!			velterm = .not.CONCAVE(xli,1,m)
!		else if(WEST_NO_MANS_LAND(vcelli(1)).or.WEST_NO_MANS_LAND(vcelltemp))then
!			velterm = CONCAVE(xli,1,m)
!		else
!			velterm = .TRUE.
!		end if
!		vcelli(1) = vcelltemp
!		xli(1) = xltemp
!
!		! OUTER POINTS
!!		if(velterm)then
!			vcelltemp = vcello(1)
!			xltemp = xlo(1)
!			if(.not.RPR_CELL_IN_PROC(vcelltemp))then
!				WEST_PERIODIC_IMAGE(vcello(1),vcelltemp, xlo(1),xltemp)
!				EAST_PERIODIC_IMAGE_MOD(vcello(1),vcelltemp,xlo(1),xltemp)
!				if(.not.RPR_CELL_IN_PROC(vcelltemp))then
!					if(I_AM_NODE_ZERO)then
!						if(vcelltemp.eq.mxf-3)then
!							vcelltemp = vcelltemp-(mxf-1)+1
!							xltemp = xltemp-(mxf-1)
!						endif
!					else
!						PRINT*,' ERROR WITH EXTERNAL POINT IN THIS PROCESSOR : ', myid, m, l, xlo(1), vcelltemp,vcello(1),xli(1)
!					end if
!				end if
!			end if
!			vcello(1) = vcelltemp
!			xlo(1) = xltemp
!!		endif
!
!		if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
!			PRINT*,' VELTERM = ', l,velterm, myid, m
!			PARALLEL_FINISH()
!			STOP
!		endif
#else
			velterm = .TRUE.
#endif
		!---------------------------------------

			!ON THE OUTTER SPHERE
			call interpolate_pdata(pcello,xlop,pglo,plo,l)
			out1%theta(ith)%phi(iphi)%pres_comp(1) = plo

			call interpolate_udata(vcello,xlo,ib,ie,jb,je,kb,ke,ulo,rdummy,rdummy,rdummy,0,m,l,onew)

			umagnormo  = dot_product(ulo(:),cd(:,l))
			do idim=1, ndim
				ulo_p(idim)  = ulo(idim)  - umagnormo  * cd(idim,l)
			enddo

			out1visc_comp => out1%theta(ith)%phi(iphi)%visc_comp
			out1visc_comp = 0d0
			do k = 1, onew 
				do j = 1, onew
					do i = 1, onew
						ii = ib+i-1
						jj = jb+j-1
						kk = kb+k-1
#if !PARALLEL
						if(ii.lt.1) ii = mxf+ii-1
						if(ii.gt.mxf-1) ii = ii-(mxf-1)
#endif
						if(jj.lt.1) jj = my+jj
						if(jj.gt.my) jj = jj-my
						if(kk.lt.1) kk = mz+kk
						if(kk.gt.mz) kk = kk-mz 
						LOCAL_INDEX(ii)                

						do dim1=1,ndim
							do dim2=1, ndim
								out1visc_comp(dim1,dim2) = out1visc_comp(dim1,dim2) + weightp(i,j,k) * &
								& vis * (ugrad(ii,jj,kk,dim1,dim2)+ugrad(ii,jj,kk,dim2,dim1))
							enddo
						enddo
					enddo
				enddo
			enddo

			!ON THE 2nd OUTTER SPHERE
			call interpolate_pdata(pcello2,xlo2p,pglo2,plo2,l)
			out2%theta(ith)%phi(iphi)%pres_comp(1) = plo2

			call interpolate_udata(vcello2,xlo2,ib,ie,jb,je,kb,ke,ulo2,rdummy,rdummy,rdummy,0,m,l,onew)  

			umagnormo2 = dot_product(ulo2(:),cd(:,l))
			do idim=1, ndim
				ulo2_p(idim) = ulo2(idim) - umagnormo2 * cd(idim,l)
			enddo

			out2visc_comp => out2%theta(ith)%phi(iphi)%visc_comp
			out2visc_comp = 0d0
			do k = 1, onew 
				do j = 1, onew
					do i = 1, onew
						ii = ib+i-1
						jj = jb+j-1
						kk = kb+k-1
#if !PARALLEL
						if(ii.lt.1) ii = mxf+ii-1
						if(ii.gt.mxf-1) ii = ii-(mxf-1)
#endif
						if(jj.lt.1) jj = my+jj
						if(jj.gt.my) jj = jj-my
						if(kk.lt.1) kk = mz+kk
						if(kk.gt.mz) kk = kk-mz 
						LOCAL_INDEX(ii)                

						do dim1=1,ndim
							do dim2=1, ndim
								out2visc_comp(dim1,dim2) = out2visc_comp(dim1,dim2) + weightp(i,j,k) * &
								& vis * (ugrad(ii,jj,kk,dim1,dim2)+ugrad(ii,jj,kk,dim2,dim1))
							enddo
						enddo
					enddo
				enddo
			enddo

			!visc_compUTING THE VISCOUS PART OF THE STRESS ON THE SURFACE OF THE PARTICLE TENSOR BASE ON CARTISION COORDINATES
			!ON THE SURFACE
			call interpolate_pdata(pcellb,xlp,pgl,pl,l)
			surf%theta(ith)%phi(iphi)%pres_comp(1) = pl

			call interpolate_udata(vcellb,xl,ib,ie,jb,je,kb,ke,ul,rdummy,rdummy,rdummy,0,m,l,onew) 

			umagnorm   = dot_product(ul(:),cd(:,l))
			do idim=1, ndim
				ul_p(idim)   = ul(idim)   - umagnorm   * cd(idim,l)
			enddo

			do idim=1, ndim
				tmp1 = (-3*ul_p(idim)+4*ulo_p(idim)-ulo2_p(idim))/(2*dr)

				ugp_th_loc(ith,idim) = ugp_th_loc(ith,idim) + tmp1
!				ugp_loc(nth,nphi,idim) = tmp1
			enddo

			thbin_loc(ith) = thbin_loc(ith)+1

!			if (ith==0) then
!				write (*,*) th*180/pi, phi*180/pi
!				write (*,"(6D15.7)") ul,ul_p
!				write (*,"(6D15.7)") ulo,ulo_p
!				write (*,"(6D15.7)") ulo2,ulo2_p
!				write (*,"(3D15.7)") ugp_th_loc(ith,:)
!				write (*,"(1I6)") thbin_loc(ith)
!				write (*,*)
!				read(*,*)
!			endif

			surf%theta(ith)%th = th
			surf%theta(ith)%phi(iphi)%ph = phi
			visc_comp => surf%theta(ith)%phi(iphi)%visc_comp
			visc_comp = 0d0
			do k = 1, onew 
				do j = 1, onew
					do i = 1, onew
						ii = ib+i-1
						jj = jb+j-1
						kk = kb+k-1
#if !PARALLEL
						if(ii.lt.1) ii = mxf+ii-1
						if(ii.gt.mxf-1) ii = ii-(mxf-1)
#endif
						if(jj.lt.1) jj = my+jj
						if(jj.gt.my) jj = jj-my
						if(kk.lt.1) kk = mz+kk
						if(kk.gt.mz) kk = kk-mz 
						LOCAL_INDEX(ii)                

						do dim1=1,ndim
							do dim2=1, ndim
								visc_comp(dim1,dim2) = visc_comp(dim1,dim2) + weightp(i,j,k) * vis * (ugrad(ii,jj,kk,dim1,dim2)+ugrad(ii,jj,kk,dim2,dim1))
							enddo
						enddo
					enddo
				enddo
			enddo

			call rotate_tensor(out1visc_comp,th,phi)
			call rotate_tensor(out2visc_comp,th,phi)
			call rotate_tensor(visc_comp,th,phi)

			nullify(visc_comp, out1visc_comp, out2visc_comp)
10		enddo

!		goto 20

		do ith=1, surf%nth-1
			thp = ith+1
			thm = ith-1

			layer_m => surf%theta(thm)
			layer_c => surf%theta(ith)
			layer_p => surf%theta(thp)

			do iphi=0, layer_c%nphi-1
				visc_comp     => surf%theta(ith)%phi(iphi)%visc_comp
				out1visc_comp => out1%theta(ith)%phi(iphi)%visc_comp
				out2visc_comp => out2%theta(ith)%phi(iphi)%visc_comp

				pres_comp     => surf%theta(ith)%phi(iphi)%pres_comp(1)
				out1pres_comp => out1%theta(ith)%phi(iphi)%pres_comp(1)
				out2pres_comp => out2%theta(ith)%phi(iphi)%pres_comp(1)
				
				!SEARCHING IN THE PREVIOUS LAYER
				do i=0, layer_m%nphi-1
					if (i/=layer_m%nphi-1) then
						if (layer_c%phi(iphi)%ph>=layer_m%phi(i)%ph.and.layer_c%phi(iphi)%ph<=layer_m%phi(i+1)%ph) then
							dist = (layer_c%phi(iphi)%ph-layer_m%phi(i)%ph) / (layer_m%phi(i+1)%ph-layer_m%phi(i)%ph)
							viscm(:,:)= (1d0-dist) * layer_m%phi(i)%visc_comp(:,:) + dist * layer_m%phi(i+1)%visc_comp(:,:)
							presm     = (1d0-dist) * layer_m%phi(i)%pres_comp(1)   + dist * layer_m%phi(i+1)%pres_comp(1)
							exit
						endif
					else
						if (layer_c%phi(iphi)%ph>=layer_m%phi(i)%ph.and.layer_c%phi(iphi)%ph<=twopi) then
							dist = (layer_c%phi(iphi)%ph-layer_m%phi(i)%ph) / (twopi-layer_m%phi(i)%ph)
							viscm(:,:)= (1d0-dist) * layer_m%phi(i)%visc_comp(:,:) + dist * layer_m%phi(0)%visc_comp(:,:)
							presm     = (1d0-dist) * layer_m%phi(i)%pres_comp(1)   + dist * layer_m%phi(0)%pres_comp(1)
							exit
						endif
					endif
				enddo
				if (i==layer_m%nphi) then
					write (*,*) "INTERPOLATION FAILED FOR PHI ON THE PREVIOUS LAYER"
					stop
				endif

				!SEARCHING IN THE NEXT LAYER
				do i=0, layer_p%nphi-1
					if (i/=layer_p%nphi-1) then
						if (layer_c%phi(iphi)%ph>=layer_p%phi(i)%ph.and.layer_c%phi(iphi)%ph<=layer_p%phi(i+1)%ph) then
							dist = (layer_c%phi(iphi)%ph-layer_p%phi(i)%ph) / (layer_p%phi(i+1)%ph-layer_p%phi(i)%ph)
							viscp(:,:)= (1d0-dist) * layer_p%phi(i)%visc_comp(:,:) + dist * layer_p%phi(i+1)%visc_comp(:,:)
							presp     = (1d0-dist) * layer_p%phi(i)%pres_comp(1)   + dist * layer_p%phi(i+1)%pres_comp(1)
							exit
						endif
					else
						if (layer_c%phi(iphi)%ph>=layer_p%phi(i)%ph.and.layer_c%phi(iphi)%ph<=twopi) then
							dist = (layer_c%phi(iphi)%ph-layer_p%phi(i)%ph) / (twopi-layer_p%phi(i)%ph)
							viscp(:,:)= (1d0-dist) * layer_p%phi(i)%visc_comp(:,:) + dist * layer_p%phi(0)%visc_comp(:,:)
							presp     = (1d0-dist) * layer_p%phi(i)%pres_comp(1)   + dist * layer_p%phi(0)%pres_comp(1)
							exit
						endif
					endif
				enddo
				if (i==layer_p%nphi) then
					write (*,*) "INTERPOLATION FAILED FOR PHI ON THE NEXT LAYER"
					stop
				endif

				phip = iphi+1
				phim = iphi-1
				if (phip>=layer_c%nphi-1) phip=0
				if (phim<0)     phim=layer_c%nphi-1

				do dim1=1, ndim
					ddr = (-3 * r**2 * visc_comp(dim1,1) + 4 * ro**2 * out1visc_comp(dim1,1) - ro2**2 * out2visc_comp(dim1,1)) / (2*dr)

					ddth = (sin(layer_p%th) * viscp(dim1,2) - sin(layer_m%th) * viscm(dim1,2)) / (layer_p%th-layer_m%th)

					ddphi = (layer_c%phi(phip)%visc_comp(dim1,3)-layer_c%phi(phim)%visc_comp(dim1,3)) / (2*layer_c%phi(1)%ph)

!					if (iph==0.or.iph==layer_c%nphi-1) then
!						ddphi = (layer_c%phi(phip)%visc_comp(dim1,3)-layer_c%phi(phim)%visc_comp(dim1,3)) / (2*layer_c%phi(1)%ph)
!					else
!						ddphi = (layer_c%phi(phip)%visc_comp(dim1,3)-layer_c%phi(phim)%visc_comp(dim1,3)) / (layer_c%phi(phip)%ph-layer_c%phi(phim)%ph)
!					endif

					div(dim1) = 1d0/r**2 * ddr + 1d0/(r*sin(layer_c%th)) * ddth + 1d0/(r*sin(layer_c%th)) * ddphi
				enddo

				layer_c%phi(iphi)%tdiv = div

				ddr = (-3 * pres_comp + 4 * out1pres_comp - out2pres_comp) / (2*dr)
				ddth = (presp - presm) / (layer_p%th-layer_m%th)
				ddphi = (layer_c%phi(phip)%pres_comp(1)-layer_c%phi(phim)%pres_comp(1)) / (2*layer_c%phi(1)%ph)
!				ddphi = (layer_c%phi(phip)%pres_comp(1)-layer_c%phi(phim)%pres_comp(1)) / (layer_c%phi(phip)%ph-layer_c%phi(phim)%ph)

				gradp(1) = -ddr
				gradp(2) = -1d0/r * ddth
				gradp(3) = -1d0/(r*sin(layer_c%th)) * ddphi

				layer_c%phi(iphi)%pgrad = gradp

				ul(:)=usmean(:)-ufmean(:)
				call rotate_vector(ul,layer_c%th,layer_c%phi(iphi)%ph)
				layer_c%phi(iphi)%vel=ul

				tmp1 = dot_product(div,ul)
				layer_c%phi(iphi)%eps_visc = tmp1 / (tke*umeanslip/dia_phys)

				tmp1 = dot_product(gradp,ul)
				layer_c%phi(iphi)%eps_pres = tmp1 / (tke*umeanslip/dia_phys)

!				layer_c%phi(iphi)%eps_pres(1) = ul(1) !ddr
!				layer_c%phi(iphi)%eps_pres(2) = ul(2) !1d0/r * ddth
!				layer_c%phi(iphi)%eps_pres(3) = ul(3) !1d0/(r*sin(layer_c%th)) * ddphi

				visc_th_loc(ith) = visc_th_loc(ith) + layer_c%phi(iphi)%eps_visc
				pres_th_loc(ith) = pres_th_loc(ith) + layer_c%phi(iphi)%eps_pres

				nullify(visc_comp, out1visc_comp, out2visc_comp)
				nullify(pres_comp, out1pres_comp, out2pres_comp)
			enddo
		
			do j=0, 7
				phi_star = j*pi/4
				do iphi=0, layer_c%nphi-1
					if (iphi/=layer_c%nphi-1) then
						if (layer_c%phi(iphi)%ph<=phi_star.and.phi_star<=layer_c%phi(iphi+1)%ph) then
							dist = (phi_star-layer_c%phi(iphi)%ph) / (layer_c%phi(iphi+1)%ph-layer_c%phi(iphi)%ph)

							tmp1 = (1d0-dist) * layer_c%phi(iphi)%eps_visc + dist * layer_c%phi(iphi+1)%eps_visc
							tmp2 = (1d0-dist) * layer_c%phi(iphi)%eps_pres + dist * layer_c%phi(iphi+1)%eps_pres

							eps_visc_phi(ith,j) = tmp1 !/ (tke*umeanslip/dia_phys)
							eps_pres_phi(ith,j) = tmp2 !/ (tke*umeanslip/dia_phys)
							exit
						endif
					else
						if (layer_c%phi(iphi)%ph<=phi_star.and.phi_star<=twopi) then
							dist = (phi_star-layer_c%phi(iphi)%ph) / (twopi-layer_c%phi(iphi)%ph)

							tmp1 = (1d0-dist) * layer_c%phi(iphi)%eps_visc + dist * layer_c%phi(0)%eps_visc
							tmp2 = (1d0-dist) * layer_c%phi(iphi)%eps_pres + dist * layer_c%phi(0)%eps_pres

							eps_visc_phi(ith,j) = tmp1 !/ (tke*umeanslip/dia_phys)
							eps_pres_phi(ith,j) = tmp2 !/ (tke*umeanslip/dia_phys)
							exit
						endif
					endif
				enddo
				if (iphi==layer_c%nphi) then
					write (*,*) "INTERPOLATION PROBLEM IN SPECIFIED PHI: ", ith, j
					stop
				endif
			enddo

			nullify(layer_m, layer_c, layer_p)
		enddo
20		nullify(bndarray)
	enddo

	write (*,"(3D15.7)") usmean(:)-ufmean(:)

	do ith=0, nth
		write (*,*) ith, surf%theta(ith)%th*180/pi
!		do iphi=0, surf%theta(ith)%nphi-1
!			do i=1, 3
!				write (*,"(3D15.7)") surf%theta(ith)%phi(iphi)%visc_comp(i,1), surf%theta(ith)%phi(iphi)%visc_comp(i,2), surf%theta(ith)%phi(iphi)%visc_comp(i,3)
!			enddo
!			write (*,*)
!		enddo



		do iphi=0, surf%theta(ith)%nphi-1
			write (*,"(3D15.7)") surf%theta(ith)%phi(iphi)%eps_pres, surf%theta(ith)%phi(iphi)%eps_visc, surf%theta(ith)%phi(iphi)%eps_pres + surf%theta(ith)%phi(iphi)%eps_visc
!			write (*,"(9D15.7)") surf%theta(ith)%phi(iphi)%tdiv, 0d0, surf%theta(ith)%phi(iphi)%vel, 0d0, surf%theta(ith)%phi(iphi)%eps_visc
		enddo

!			write (*,"(2D15.7)") surf%theta(ith)%phi(iphi)%eps_visc, surf%theta(ith)%phi(iphi)%eps_pres
		read (*,*)
	enddo


#if PARALLEL

#else
	ugp_th(:,:) = ugp_th_loc(:,:)
	thbin(:) = thbin_loc(:)
	visc_th(:) = visc_th_loc(:)
	pres_th(:) = pres_th_loc(:)
#endif

	do ith=0, nth
!		write (*,"(6D15.7)") ugp_th_loc(ith,:),ugp_th_loc(ith,:)
!		read (*,*)
		ugp_th(ith,:) = ugp_th(ith,:) / thbin(ith) / (umeanslip/dia_phys)
		visc_th(ith)  = visc_th(ith)  / thbin(ith) !/ (tke*umeanslip/dia_phys)
		pres_th(ith)  = pres_th(ith)  / thbin(ith) !/ (tke*umeanslip/dia_phys)
	enddo

	if (I_AM_NODE_ZERO) then
		open (unit=1, file="ug_eps_th.dat", status="replace")
		write (1,"(A)") "variables=ith,nth,theta,u1,u2,u3,umag,eps_visc,eps_pres,eps"
		write (1,"(A)") "zone"
		do ith=0, nth
			tmp1 = dsqrt(ugp_th(ith,1)**2+ugp_th(ith,2)**2+ugp_th(ith,3)**2)
			write (1,"(2I6,8D15.7)") ith, thbin(ith), surf%theta(ith)%th*180/pi, ugp_th(ith,1), ugp_th(ith,2), ugp_th(ith,3), tmp1, visc_th(ith), pres_th(ith), visc_th(ith)+pres_th(ith)
		enddo
		close(1)

		open (unit=2, file="eps_phi.dat", status="replace")
		write (2,"(A)") "variables=ith,theta,eps_pres,eps_visc,eps"
		do j=0, 7
			write (2,"(A)") "zone"
			do ith=0, nth
				write (2,"(1I6,4D15.7)") ith, surf%theta(ith)%th*180/pi, eps_visc_phi(ith,j), eps_pres_phi(ith,j), eps_visc_phi(ith,j)+eps_pres_phi(ith,j)
			enddo
		enddo
		close(2)
	endif

contains

	subroutine rotate_tensor(tensor,th,phi)
		implicit none
		real(8), intent(in) :: th,phi
		real(8), intent(inout) :: tensor(ndim,ndim)

		real(8) :: tmp(3,3)
		real(8) :: rotate(3,3)

		integer :: dim1,dim2,dim3,dim4

		!ROTATION OF THE STRESS TENSOR ALIGNED WITH THE SPHERICAL COORDINATE
		rotate(1,1) =-cos(th)
		rotate(2,1) = sin(th) * cos(phi)
		rotate(3,1) = sin(th) * sin(phi)

		rotate(1,2) = sin(th)
		rotate(2,2) = cos(th) * cos(phi)
		rotate(3,2) = cos(th) * sin(phi)

		rotate(1,3) = 0d0
		rotate(2,3) =-sin(phi)
		rotate(3,3) = cos(phi)

		tmp = 0d0
		do dim1=1, ndim
			do dim2=1, ndim
				do dim3=1, ndim
					do dim4=1, ndim
						tmp(dim1,dim2) = tmp(dim1,dim2) + rotate(dim3,dim1) * rotate(dim4,dim2) * tensor(dim3,dim4)
					enddo
				enddo
			enddo
		enddo
		tensor(:,:) = tmp(:,:)
	end subroutine rotate_tensor

	subroutine rotate_vector(vec,th,phi)
		implicit none
		real(8), intent(in) :: th,phi
		real(8), intent(inout) :: vec(ndim)

		real(8) :: tmp(3)
		real(8) :: rotate(3,3)

		integer :: dim1,dim2,dim3,dim4

		!ROTATION OF THE STRESS TENSOR ALIGNED WITH THE SPHERICAL COORDINATE
		rotate(1,1) =-cos(th)
		rotate(2,1) = sin(th) * cos(phi)
		rotate(3,1) = sin(th) * sin(phi)

		rotate(1,2) = sin(th)
		rotate(2,2) = cos(th) * cos(phi)
		rotate(3,2) = cos(th) * sin(phi)

		rotate(1,3) = 0d0
		rotate(2,3) =-sin(phi)
		rotate(3,3) = cos(phi)

		tmp = 0d0
		do dim1=1, ndim
			do dim2=1, ndim
				tmp(dim1) = tmp(dim1) + rotate(dim2,dim1) * vec(dim2)
			enddo
		enddo
		vec(:) = tmp(:)
	end subroutine rotate_vector

	subroutine calc_ugrad
		implicit none
		integer :: ip, im, dir

		do i=1, nx
			if (I_AM_NODE_ZERO.and.mod(i,10)==0) write (*,*) " PLANE # = ", i
			ip=i+1
			im=i-1
#if !PARALLEL
			if (ip>nx) ip=1
			if (im<1)  im=nx
#endif
			do dir=1, ndim
				do idim=1, ndim
					do k=1, mz
						do j=1, my2
							if (dir==1) then
								uf1(j,k) = (u(ip,j,k,idim)-u(im,j,k,idim)) / (two*dx)
							elseif (dir==2) then
								uf1(j,k) = u(i,j,k,idim) * wy(j)
							elseif (dir==3) then
								uf1(j,k) = u(i,j,k,idim) * wz(k)
							endif
						enddo
					enddo
					call ff2cr(uf1,ur1)
					ugrad(i,:,:,idim,dir) = ur1(:,:)
				enddo
			enddo
		enddo

	
!		open (unit=1, file=trim(run_name)//"_ugrad.dat", status="replace")
!		write (1,*) "variables=i,j,k,u11,u12,u13" !,u21,u22,u23,u31,u32,u33"
!		write (1,*) "zone  i=",nx," j=",my," k=",mz," f=point"
!
!		do k=1,mz		
!			do j=1, my
!				do i=1, nx
!!					do idim=1, ndim
!						write (1,"(3I6,3D15.7)") i,j,k,ugrad(i,j,k,1,:)
!!					enddo
!				enddo
!			enddo
!		enddo
!		close(1)


	end subroutine calc_ugrad

	subroutine transfer_velocity
		implicit none
		integer :: left, right, l
		real(8), allocatable :: trans_s(:), trans_r(:)
#if PARALLEL
		left = myid-1
		right = myid+1

		if (left<0) left=left+nproc
		if (right>nproc-1) right=right-nproc

		allocate(trans_s(my*mz*3), trans_r(my*mz*3))

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					trans_s(l)=ubcp(1,j,k,idim)
				enddo
			enddo
		enddo

	
		call MPI_SENDRECV(trans_s,my*mz*3,MPI_DOUBLE_PRECISION,left,0,trans_r,my*mz*3,MPI_DOUBLE_PRECISION,right,0,decomp_group,status,err_code)

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					ubcp(nx+1,j,k,idim)=trans_r(l)
				enddo
			enddo
		enddo

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					trans_s(l)=ubcp(nx,j,k,idim)
				enddo
			enddo
		enddo
	
		call MPI_SENDRECV(trans_s,my*mz*3,MPI_DOUBLE_PRECISION,right,1,trans_r,my*mz*3,MPI_DOUBLE_PRECISION,left,1,decomp_group,status,err_code)

		l=0
		do idim=1, ndim
			do k=1, mz
				do j=1, my
					l=l+1
					ubcp(0,j,k,idim)=trans_r(l)
				enddo
			enddo
		enddo
#else
		ubcp(nx+1,:,:,:) = ubcp(1,:,:,:)
#endif	
	end subroutine transfer_velocity

end subroutine post
end module mypost_process
