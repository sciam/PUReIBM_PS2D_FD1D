
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
						&  usmean, ufmean, pres_avg, visc_avg
#if PARALLEL
	use global_data, only : err_code, decomp_group, myid, starts, ends
#endif
	use nlmainarrays , only : ubcp, pbcp
	use fftw_interface
	use constants
	use nlarrays
	use randomno
	use general_funcs
	use string_funcs
	use dependent_functions
	use bcsetarrays

	logical, allocatable :: far_field(:,:,:)
	logical, save :: near_particle_checked=.false.
contains

subroutine post
	implicit none
	
	real(8), allocatable :: dissip(:,:,:), fluc(:,:,:)
	real(8), allocatable, dimension(:) :: pres_part, visc_part, int_tke
	real(8) :: sijsij, sijsij_loc, sijsij_far, sijsij_far_loc

	character*40 filename0
	
	real(8) :: tmp1, tmp2
	integer :: i, j, k, m, idim, iph
	integer :: im, ip, jm, jp, km, kp, ii, jj, kk, iii, jjj, kkk
	integer :: dim1, dim2, dim3
	logical :: neighb_insolid

	allocate (dissip(nx,my,mz),fluc(nx,my,mz))

	if (near_particle_checked==.false.) then
		call near_particle_region
		near_particle_checked = .true.
	endif

	dissip = 0d0
	fluc   = 0d0

	sijsij_loc = 0d0
	sijsij_far_loc = 0d0
!	generating thetaij
	if (I_AM_NODE_ZERO) write (*,*) "GENERATING SijSij"
	do i=1, nx
		if (I_AM_NODE_ZERO.and.mod(i,20)==0) write (*,*) " PLANE # = ", i
		do dim1=1, ndim
			do dim2=1, ndim
				call derivative_1st(i,dim1,dim2,ur1)
				call derivative_1st(i,dim2,dim1,ur2)

				do k=1, mz
					do j=1, my
						if (fluid_atijk(i,j,k)) then
							tmp1 = 2*vis*((ur1(j,k)+ur2(j,k))*half)**2
!							tmp1 = vis*(ur1 (j,k))**2
							sijsij_loc = sijsij_loc + tmp1
			
							dissip(i,j,k) = dissip(i,j,k) + tmp1
							fluc(i,j,k) = dot_product(ubcp(i,j,k,:)-ufmean(:),ubcp(i,j,k,:)-ufmean(:))/2

							if (far_field(i,j,k)) then
								sijsij_far_loc = sijsij_far_loc + tmp1
							endif
						else
							dissip(i,j,k) = -(vis*(umeanslip/dia_phys)**2)
							fluc(i,j,k)   = -(dot_product(ufmean(:),ufmean(:))/2)
						endif
					enddo
				enddo
			enddo
		enddo
	enddo
	dissip = dissip / (vis*(umeanslip/dia_phys)**2)
	fluc   = fluc   / (dot_product(ufmean(:),ufmean(:))/2)

#if PARALLEL
	sijsij = 0d0
	GLOBAL_DOUBLE_SUM(sijsij_loc,sijsij,1,decomp_group)

	sijsij_far = 0d0
	GLOBAL_DOUBLE_SUM(sijsij_far_loc,sijsij_far,1,decomp_group)
#else
	sijsij     = sijsij_loc
	sijsij_far = sijsij_far_loc
#endif
	sijsij     = sijsij     / count_fluid / (vis*(umeanslip/dia_phys)**2)
	sijsij_far = sijsij_far / count_fluid / (vis*(umeanslip/dia_phys)**2)

	allocate(pres_part(nphases+1), visc_part(nphases+1), int_tke(nphases+1))

	int_tke = 0d0
	pres_part=0d0
	visc_part=0d0

	do m=1, nbody
		if (m>phase_array(1)%npart) then
			iph = 2
		else
			iph = 1
		endif

		do idim=1, ndim
			pres_part(iph) = pres_part(iph) + pres(m,idim) * (ufmean(idim)-usmean(idim)) !/ phase_array(iph)%npart
			visc_part(iph) = visc_part(iph) + visc(m,idim) * (ufmean(idim)-usmean(idim)) !/ phase_array(iph)%npart
			int_tke(iph)   = int_tke(iph)   + (visc(m,idim)+pres(m,idim)) * (ufmean(idim)-usmean(idim)) !/ phase_array(iph)%npart

			pres_part(nphases+1) = pres_part(nphases+1) + pres(m,idim) * (ufmean(idim)-usmean(idim))
			visc_part(nphases+1) = visc_part(nphases+1) + visc(m,idim) * (ufmean(idim)-usmean(idim))
			int_tke(nphases+1)   = int_tke(nphases+1)   + (visc(m,idim)+pres(m,idim)) * (ufmean(idim)-usmean(idim))
		enddo
	enddo
	
	do i=1, nphases+1
		pres_part(i)=pres_part(i) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)**3) / (vis*(umeanslip/dia_phys)**2)
		visc_part(i)=visc_part(i) / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)**3) / (vis*(umeanslip/dia_phys)**2)

		int_tke(i) = int_tke(i)   / (doml(1)*doml(2)*doml(3)*(1-maxvolfrac)**3) / (vis*(umeanslip/dia_phys)**2)
	enddo

	if (I_AM_NODE_ZERO) then
		open (unit=1, file=trim(run_name)//"_sijsij.dat",status="replace")
		write (1,"(A,1D15.7)") "SijSij          = ", sijsij
		write (1,"(A,1D15.7)") "SijSij_FAR_FRAC = ", sijsij_far/sijsij
		Write (1,"(A,1D15.7)") "INT_TKE         = ", int_tke(nphases+1)
		Write (1,"(A,1D15.7)") "VISC_PART_FRAC  = ", visc_part(nphases+1)/int_tke(nphases+1)
		Write (1,"(A,1D15.7)") "PRES_PART_FRAC  = ", pres_part(nphases+1)/int_tke(nphases+1)
		close (1)

		open (unit=2, file=trim(run_name)//"_dissip_field.dat",status="replace")
		write (2,*) "variables=x,y,z,tke,eps"
		write (2,*) "zone I=,", nx, " J=", my, " K=", mz
		do k=1, mz
			do j=1, my
				do i=1 ,nx
					write (2,"(3I6,2D15.7)") i, j, k, fluc(i,j,k), dissip(i,j,k)
				enddo
			enddo
		enddo
		close (2)
	endif


	contains

	subroutine derivative_1st(pos,dim1,dim2,duidxj)
		implicit none
		integer, intent(in)  :: pos, dim1, dim2
		real(8), intent(out) :: duidxj(my,mz)
		
		real(8) :: tmp1, tmp2
		integer :: i,j,k
		
		i = pos
		ip=i+1
		im=i-1
		
#if !PARALLEL
		if(im<1)  im = nx
		if(ip>nx) ip = 1
#endif       
		if (dim2==1) then
			do k=1, mz
				do j=1, my
					duidxj(j,k)=(1./2/dx)* (ubcp(ip,j,k,dim1)-ubcp(im,j,k,dim1))
				enddo
			enddo
		else
			do k=1, mz
				do j=1, my2
					if (dim2==2) then
						uf3(j,k) = u(i,j,k,dim1) * wy(j)
					else
						uf3(j,k) = u(i,j,k,dim1) * wz(k)
					endif
				enddo
			enddo
			call ff2cr(uf3,ur3)
			duidxj = ur3
		endif
	end subroutine derivative_1st

end subroutine post

subroutine near_particle_region
	implicit none
	
	integer :: i, j, k, ii, jj, kk, m, idim, imin, imax, jmin, jmax, kmin, kmax
	real(8), dimension(ndim) :: xlr, xll
	integer, dimension(ndim) :: cor_min, cor_max
	real(8) :: near_particle_rad, dist

	near_particle_rad = 1.05
	allocate(far_field(0:nx+1,my,mz))

	do k = 1, mz
		do j = 1, my
			do i = 0, nx+1 
				far_field(i,j,k) = .true.
			end do
		end do
	end do

    do m = 1, nbody
       do idim = 1, ndim 
          if(idim.eq.1) then 
             xlr(idim) = xc(m,idim)  + radbdy(m) + foffset
             xll(idim) = xc(m,idim)  - radbdy(m) + foffset
          else 
             xlr(idim) = xc(m,idim)  + radbdy(m)
             xll(idim) = xc(m,idim)  - radbdy(m) 
          end if
       end do
    
       do idim = 1, ndim 
          cor_min(idim)  = ceiling(xll(idim))
          cor_max(idim) = floor(xlr(idim)) 
       end do
    
       imin = cor_min(1)
       imax = cor_max(1)
       jmin = cor_min(2)
       jmax = cor_max(2)
       kmin = cor_min(3)
       kmax = cor_max(3)
       
       
       do i = imin, imax 
          ii = i
#if PARALLEL
          if(.not.POINT_IN_PROC(ii))then

             WEST_PERIODIC_POINT_IMAGE(i,ii)
             EAST_PERIODIC_POINT_IMAGE(i,ii)
             if(.not.POINT_IN_PROC(ii)) goto 555
          end if
          
#else
          if(i.lt.1.and.intx_per) ii = mxf+i-1
          if(i.gt.mxf-1.and.intx_per) ii = i-(mxf-1)
#endif
          
          LOCAL_INDEX(ii)

          do j = jmin, jmax 
             jj = j 
             if(j.lt.1.and.inty_per) jj = my+j
             if(j.gt.my.and.inty_per) jj = j-my
             
             do k = kmin, kmax 
                
                kk = k 
                if(k.lt.1.and.intz_per) kk = mz+k
                if(k.gt.mz.and.intz_per) kk = k-mz 
                
                dist = (i - (xc(m,1)+foffset))**two
                dist = dist + (j - xc(m,2))**two + (k - xc(m,3))**two 
                dist = dsqrt(dist)
                IF((dist - radbdy(m)*near_particle_rad).LE.SMALL_NUMBER) THEN 
                   if(far_field(ii,jj,kk).eq..false.) then 
                      PRINT*,'FLUID_ATIJK ALREADY FALSE AT I,J,K=&
                           & ', ii,jj,kk, myid,m
                   ENDIF
                   far_field(ii,jj,kk)  = .false.
                end IF
             end do
          end do
#if PARALLEL
555       continue
          
#endif
       end do
    end do
end subroutine near_particle_region

end module mypost_process

