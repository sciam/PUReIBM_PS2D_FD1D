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


MODULE initialize_flo
#include "ibm.h"
  USE precision 
  USE constants 
  USE scalar_init
  USE general_funcs
  USE randomno
  USE hard_sphere
  USE postproc_funcs
  USE initialize
  !use init_turb

!  use mypost_process

  IMPLICIT NONE 
  !Read input files for initial velocity field and simulation
  !parameters. Check initial flow field for continuity and 
  !correct for any divergence. 
  
  INTEGER, PARAMETER, PRIVATE :: nbin_max = 200, mis_gof=120
  INTEGER, PRIVATE :: gofunit, gofavgunit
  REAL(PRCN), PRIVATE ::  gofr(nbin_max), rho_est(nbin_max)
  LOGICAL, PRIVATE :: rescaling
CONTAINS
  SUBROUTINE initflo
    USE dependent_functions
    USE restart_funcs
    USE fftw_interface
    USE global_data
    USE geom_init
    USE nlarrays, ONLY  : uf1, ur1
    Use nlmainarrays, Only : ubcp, nlbcp, onlbcp, pbcp , ubc, nlbc, onlbc, pbc
    USE maternmod
    USE collision_mod

    !^^^ Mohammad 10-16-2009 ^^^
!	use mypost_process

!    use string_funcs
!    use init_turb, only : avr, tke, turb_uf
!    use mypost_process
!   use steptheflow
    !---------------------------

use bcsetarrays
use fftw_interface


    IMPLICIT NONE  

    !-----------------------------------------------------------------------
    !local variables


    REAL(prcn) ::  jt(ndim,ndim), flagu,dum,&
         & gammaloc, ldomtmp(3), upi_int(ndim), uchar_tmp

    INTEGER :: i,j,k,l,n,m, mytmp, mbox,idim, partstart, partend, iphs
    
    REAL(prcn) ::   umax_tmp, lchar, xli(ndim), xlo(ndim), vfrac_mean !,umeanslip

    REAL(prcn) :: vbox

!^^^ Mohammad 12-30-2009 ^^^ to initializa_flo, is used in dtibm
!	real(8) :: dt_turb
!	integer :: ii
!	real(8) :: tke,tmp
!	real(8), allocatable:: turb_ur2(:,:,:,:)

	character*80 filename1, filename2
#if PARALLEL
!	INTEGER :: starts(nproc), ends(nproc)
	complex(8),allocatable :: trans_buf(:)

!	real(8),external :: tke_phase
!	real(8):: avr(3)
!	integer :: count,tmp_count
	integer :: node_num,iproc
	integer :: strlen

	allocate(starts(0:nproc-1), ends(0:nproc-1))
#endif


integer :: count

!---------------------------
    FIRST_PASS = .TRUE.
    
    soft_start = .true.
    move_particles = .FALSE.
    !if(imove.eq.1)soft_start = .FALSE.
    base_state_saved = .FALSE.
    compute_auto_corr = .FALSE.
    rhoofs0 = one
    stepbody=1
    gauss_u = .false.
    gauss_p = .false.
    gauss_phi = .false.
    gauss_p = .false.
    if(only_dem) move_particles = .TRUE.
    if(.not.mean_vel_to_particles)movingcv = .FALSE.
    if(imove.eq.1) movingcv = .FALSE.
    
    CALL set_interpolation_scheme(2) !DO it before allocating memory
    ! for gstencil and rest of the arrays required for interpolation.
    ! This iniliazes order variable required for allocating other
    ! arrays. Do it for the largest shcmes that is gonna be used. 
    
    if(irestart.eq.0)then
       count_restart = 0
    end if
    ! In this routine, if MY is undefined, it will be set based on LYBYD and DBYDX. If input_type is random, then LYBYD is based on the average diameter and DBYDX is based on the smallest diameter.
    
    CALL generate_configuration

    if(POST_NO_FLOW_MEM_ALLOC)RETURN

    CALL DOMAIN_1D_DECOMP(mx-1,nproc,myid,xstart,xend)
    
    nx = xend-xstart+1 ! For serial case nx = mx-1
    
    
    if(I_AM_NODE_ZERO) CALL screen_separator(80,'#')

#if PARALLEL
    starts(myid) = xstart
    ends(myid) = xend
    CALL MPI_ALLGATHER(xstart, 1, MPI_INTEGER, starts, 1, MPI_INTEGER, decomp_group, err_code)
    CALL MPI_ALLGATHER(xend,   1, MPI_INTEGER, ends,   1, MPI_INTEGER, decomp_group, err_code)
    
    if(I_AM_NODE_ZERO)then
       do i = 0, nproc-1
          Write(*,'(4(A8,2x,I6))')'RANK = ', i, 'XSTART=', starts(i), 'XEND = ', ends(i), 'NX = ', ends(i)-starts(i)+1
       end do
    end if
#else
    Write(*,'(4(A8,2x,I6))')'RANK = ', 0, 'XSTART=', xstart, 'XEND = ', xend, 'NX = ', nx
#endif

    if(I_AM_NODE_ZERO) CALL screen_separator(80,'#')
    BARRIER(decomp_group)
    
    vbox = doml(1)*doml(2)*doml(3)
    
    nphsc2 = nphases*(nphases-1)/2
    partstart = 1
    DO iphs = 1, nphases
       phase_array(iphs)%pstart = partstart
       if(irestart.eq.0)phase_array(iphs)%dia = radbdy(partstart)*two*dx
       !PRINT*,'radius = ', iphs, phase_array(iphs)%dia, dx
       phase_array(iphs)%volfrac = (phase_array(iphs)%npart)*pi*(phase_array(iphs)%dia**3.d0)/(6.d0*vbox)
       
       partend = partstart + phase_array(iphs)%npart-1 
       phase_array(iphs)%pend = partend
       do m = PARTSTART,PARTEND
          part_array(m)%iphs = iphs
       ENDDO
       partstart = partend + 1
    END DO
    if(I_AM_NODE_ZERO)then
       if((TRIM(input_type).eq.'random').and.(TRIM(psd_type).eq.'bidisp'))then
          Write(*,'(A30,I6, 2x,g17.8)')'VOLFRAC RATIO FOR BIDISP CASE = ', myid, phase_array(2)%volfrac/phase_array(1)%volfrac
       end if
    end if

    mean_volfrac = zero
    
    do iphs = 1, nphases
       mean_volfrac = mean_volfrac + phase_array(iphs)%volfrac
    end do

    if(.not.ALLOCATED(u))then
       CALL alloc_mem !allocate memory for main data 
       
       !ALSO ALLOCATE SCALAR 
       IF(ISCALON.EQ.1) call alloc_scalar_mem
    end if
    
    
    CALL grid_nodes_insphere
    
    if(I_AM_NODE_ZERO)Write(*,'(A)')'OUT OF GRID_NODES_INSPHERE'
    if(I_AM_NODE_ZERO) WRITE(*,'(2(2x,A25,i10))')'count_solid=',count_solid, 'count_fluid&
         &=', count_fluid
    
    if(I_AM_NODE_ZERO) Write(*,'(A10,2(2x,g17.8))')'MAXVOLFRAC = ',&
         & maxvolfrac, count_solid
    GRAV(:) = zero


#if 0
		if (iturbon==1) then
			turbvelunit = getnewunit(minunitno,maxunitno)
			call gener_filename(filename1,"/data003/mech_eng3/mohammad/TEST_ISOTURB/MIS1_velocity.dat")

			open(turbvelunit, filename1, action="read", format="unformatted")
			read re_lambda, Re_L, tke_i, epsf_i, vis, u_eta_i, tau_eta_i, kmin, kmax, eta_i, lambda_i, l_e_i, kmaxeta

			do k=1, mz
				do j=1, my
					do i=1,mx1
						write (9994) tmpi, tmpj, tmpk, tmpvec(:)
						if (starts(myid)<=i.and.i<=ends(myid)) 
					enddo
				enddo
			enddo
			close (9994)

			do i=1, nx
#endif


		if (impose_grav) then
			set_mpg = .TRUE.

!			mpg(:) = (-0.2638240D-02) * cos(pi*flo_ang(1:3)/180.d0)
			mpg(:) = (9.81) * cos(pi*flo_ang(1:3)/180.d0)
			grav(:) = - mpg(:) / (rhof*(one-maxvolfrac) + rhos*maxvolfrac)

!			GRAV(1:3) = -0.2638240D-02*cos(pi*flo_ang(1:3)/180.d0)
!			mpg(:) = (rhof*(one-maxvolfrac) + rhos*maxvolfrac)*GRAV(:)

			vis = DSQRT(ABS(rhos/rhof-one)* sqrt(dot_product(grav(:),grav(:))) * char_length**3.d0/archno)
			if(I_AM_NODE_ZERO)Write(*,'(2(A25,2x,g17.8))')'To attain Arch. no = ', archno, 'fluid viscosity is set to ', vis

			uchar_tmp = sqrt( sqrt(dot_product(grav(:),grav(:))) * abs(rhos/rhof-1) * char_length) * (1-maxvolfrac)
			uchar(1:3) = uchar_tmp*cos(pi*flo_ang(1:3)/180.d0)
			ucharmod = SQRT(uchar(1)**2.d0 + uchar(2)**2.d0+uchar(3)**2.d0)    
			fsslip(1:ndim) = uchar(1:ndim)/(one-maxvolfrac)

			fsslipmod = sqrt(dot_product(fsslip,fsslip))
		else
			!uchar_tmp = ((Re*vis)/char_length)!*cos(pi*flo_ang(1:3)/180.d0)
			! New formulation to make IBM GI

!			fsslipmod = one
!			uchar_tmp = (one-maxvolfrac)*fsslipmod
!			vis = uchar_tmp*char_length/Re 
			if(I_AM_NODE_ZERO) Write(*,*) 'Changing viscosity to get uchar = 1. New vis : ', vis
        
			uchar_tmp = ((Re*vis)/char_length)
			uchar(1:3) = uchar_tmp*cos(pi*flo_ang(1:3)/180.d0)
			ucharmod = SQRT(uchar(1)**2.d0 + uchar(2)**2.d0+uchar(3)**2.d0)    
			fsslip(1:ndim) = uchar(1:ndim)/(one-maxvolfrac)

			!dchar = dia_phys
		endif

	if(I_AM_NODE_ZERO)then
		write (*,"(1a,10d15.7)") 'CHAR LENGTH, UCHAR_TEMP = ', char_length, uchar_tmp
		write (*,"(1a,10d15.7)") 'UCHAR(1:3), UCHARMOD    = ', uchar, ucharmod
		write (*,"(1a,10d15.7)") 'FSSLIP(1:3), FSSLIPMOD  = ', fsslip(1:ndim), fsslipmod

		write (*,"(1a,3d15.7)") "GRAVITY = ", grav(:)
		write (*,"(1a,3d15.7)") "MPG     = ", mpg(:)
	end if
    
    ! Assign mean and fluctuating velocities to particles
    if(ReT.gt.SMALL_NUMBER) gran_temp = (ReT*vis/char_length)**2.d0
    if(irestart.eq.0)then
       if(TRIM(input_type).eq."random")then
          CALL generate_particle_velocities
       else if(TRIM(input_type).eq."default")then
!          if (I_AM_NODE_ZERO) Write(*,*)'DEFAULT CASE. VELOCITIES ARE READ FROM sphr_center.inp'
          CALL generate_particle_velocities
       else if(TRIM(input_type).eq."lubtest")then
			if (I_AM_NODE_ZERO) then
				Write(*,*)'LUBRICATION TEST CASE CASE. CHARACTERISTIC VEL SCALE IN RE IS REL. VELOCITY.'
				Write(*,*)'OPPOSITE VELOCIIES ARE ASSIGNED.'
			endif
          velbdy(1,1:ndim) = uchar(1:ndim)/two
          velbdy(2,1:ndim) = -uchar(1:ndim)/two
       else
          if(mean_vel_to_particles)then
             do idim=1,ndim
                do m = 1, nbody
                   velbdy(m,idim) =  -uchar(idim)/(one-maxvolfrac)
                end do
             end do
          else
             do m = 1, nbody
                velbdy(m,1:ndim) = zero
             end do
          Endif
       endif
       
       do idim = 1, ndim
          if(nbody.gt.0)then
             usmean_des(idim) = SUM(velbdy(1:nbody, idim))/real(nbody, prcn)
          else
             usmean_des(idim) = zero
          end if
       end do
       
       ! Assign mean velocity to the fluid.             
       if(mean_vel_to_particles)then
          do idim = 1, ndim
             ufmean_des(idim) = zero
          end do
       else
          do idim = 1, ndim
             ufmean_des(idim) = uchar(idim)/(one-maxvolfrac)
          end do
       end if
       if(TRIM(input_type).eq.'lubtest')ufmean_des(1:ndim) = zero
       
       ! Transform velocities to get homogeneous boundary conditions on particle surface
       do idim = 1, ndim
!!$          do m = 1, nbody
!!$             velbdy(m,idim) = (velbdy(m,idim) - solid_vel(idim))/fsslipmod
!!$          end do
!!$          usmean_des(idim) = (usmean_des(idim)-solid_vel(idim))/fsslipmod
!!$          ufmean_des(idim) = (ufmean_des(idim)-solid_vel(idim))/fsslipmod
          
          ! Compute the volumetric mean velocity

			if (set_mpg.or.initialize_zero) then
				umean(idim) = zero
			else
				umean(idim) = (one-maxvolfrac)*ufmean_des(idim) + maxvolfrac*usmean_des(idim)
			endif
       enddo
       
       frame_accln(1:ndim) = zero
       frame_vel(1:ndim) = zero
       frame_pos(1:ndim) = one
       CALL INPUT_CHECK

!		if (I_AM_NODE_ZERO) then
!			write (*,"(3d15.7)") 'USMEAN DESIRED : ', usmean_des(1:ndim)
!			write (*,"(3d15.7)") 'UfMEAN DESIRED : ', ufmean_des(1:ndim)
!		endif

		!^^^ Mohammad 10-16-2009 ^^^^^^^^^^^^^^^^^^^
		if (iturbon==1) then
!			call specinit
!			call transfer_vel
		endif
		!-------------------------------------------

    end if


#if 0
		do i=1, 1
			do j=1, 1l
				do k=1, 1
					do m=1, nbody
						write (1,"(4d15.7,1i)") xc(m,1)/dbydx+(i-1)*lybyd, &
					&									xc(m,2)/dbydx+(j-1)*lybyd, &
					&									xc(m,3)/dbydx+(k-1)*lybyd, radbdy(m)/dbydx, 1
					enddo
				enddo
			enddo
		enddo

		write (1,*) "zone"
		do i=1, 1
			do j=1, 1
				do k=1, 1
					do m=1, nbody
						if (xc(m,2)/dbydx<=lybyd/4.or.xc(m,2)/dbydx>=lybyd*3/4) &
					&	write (1,"(4d15.7,1i)") xc(m,1)/dbydx+(i-1)*lybyd, &
					&									xc(m,2)/dbydx+(j-1)*lybyd, &
					&									xc(m,3)/dbydx+(k-1)*lybyd, radbdy(m)/dbydx, 1
					enddo
				enddo
			enddo
		enddo

		close(1)
	endif
   PARALLEL_FINISH()
   STOP
#endif

#if 0
	if (I_AM_NODE_ZERO) then
		WRITE(*,*) 'WRITING IMMERSED BODY DATA'
		OPEN(unit=1,file=TRIM(RUN_NAME)//'_sphr_center.inp',status='replace', Action='write')
!		write (1,"(1a,1d15.7)") "VOLFRAC=", maxvolfrac
!		write (1,"(1i,1a)") nbody, " atoms"
!		write (1,"(1i,1a)") nphases, " atom types"
!		write (1,"(1a)") "atom_style = sphere"
!		write (1,"(2f8.2,1A)") zero, lybyd, " xlo xhi"
!		write (1,"(2f8.2,1A)") zero, lybyd, " ylo yhi"
!		write (1,"(2f8.2,1A)") zero, lybyd, " zlo zhi"
!		write (1,"(1a)") "atoms"
!		write (1,"(1a)") " "
		do i=1, nbody
!			write (1,"(2i,5d15.7)") i, 1, dia_phys, 1000.0, xc(i,:)/dbydx
			write (1,"(4d15.7)") xc(i,:)/dbydx, radbdy(i)/dbydx
		enddo
!		write (1,"(1a)") "Velocities"
!		write (1,"(1a)") " "
!		do i=1, nbody
!			write (1,"(1i,6d15.7)") i, velbdy(i,:), 0d0, 0d0, 0d0
!		enddo
		close(1)
	endif
   PARALLEL_FINISH()
   STOP
#endif





    
    dtorig = dt
    
    urf = 1.0d0 ! momentum equation under-relaxation factor
    prf = 1.0d0!0.8d0
    !r = 0.5*dia_phys/doml(2)*my
    upi = zero 
    if(I_AM_NODE_ZERO.AND.nbody.gt.0)then
       open(unit=2001,file=TRIM(RUN_NAME)//'_sphr_vel_out.dat',form='formatted',status='unknown')
       write(2001,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',   &
            &    ' "UX" '
       DO m=1,nbody
          WRITE(2001,'(6(2x,f12.8))') xc(m,1), xc(m,2), xc(m,3), velbdy(m,1), velbdy(m,2), velbdy(m,3)
       enddo
       close(2001, status='keep')
    end if
    if(nphases.gt.0)then
       r = 0.5*phase_array(nphases)%dia/doml(2)*my ! Resolve the largest sphere and store in the array of first phase
       
       do iphs = 1, 1!nphases
          !r = half*(phase_array(iphs)%dia)/doml(2)*my
          if(I_AM_NODE_ZERO)Write(*,*)'RADIUS IN QUADGENER = ', r
          NULLIFY(phase_array(iphs)%bndpts)
          if(ALLOCATED(xs))DEALLOCATE(xs)
          if(ALLOCATED(cd))DEALLOCATE(cd)
          CALL quadgener(r,dr,f1,f2)

          nrpr = nbnd
          if(I_AM_NODE_ZERO)then
             WRITE(*,'(A,2(2x,I5))')'No. of boundary points for phase ',iphs, nbnd           
             WRITE(ounit,'(A,2(2x,I5))')'No. of boundary points for phase ',iphs, nbnd           
          end if
       
          IF (ASSOCIATED(phase_array(iphs)%bndpts)) THEN
             DEALLOCATE(phase_array(iphs)%bndpts)
             IF (nbnd.GT.0) then 
                ALLOCATE(phase_array(iphs)%bndpts(ndim,nbnd))
             end IF
          ELSE
             IF(nbnd.GT.0) ALLOCATE(phase_array(iphs)%bndpts(ndim,nbnd))
          ENDIF
          
          phase_array(iphs)%nbnd = nbnd
          phase_array(iphs)%nrpr = nrpr
          
          do i = 1, nbnd
             do idim = 1, ndim
                phase_array(iphs)%bndpts(idim,i) = xs(idim,i)
             end do
          end do
          !READ(*,*)
       enddo
       if(ALLOCATED(xs))DEALLOCATE(xs)
    end if
    
    do m = 1, nbody
       radibdy(m) = radbdy(m)-dr
       radobdy(m) = radbdy(m)+dr
       rado2bdy(m) = radbdy(m)+two*dr
    end do
    
    do m= 1, nbody
       iphs = 1!part_array(m)%iphs
       nrpr = phase_array(iphs)%nrpr
       NULLIFY(part_array(m)%if_rev)
       IF (ASSOCIATED(part_array(m)%if_rev)) THEN
          DEALLOCATE(part_array(m)%if_rev)
          IF (nrpr.GT.0) then 
             ALLOCATE(part_array(m)%if_rev(nrpr))
          end IF
       ELSE
          IF(nrpr.GT.0) ALLOCATE(part_array(m)%if_rev(nrpr))
       ENDIF
    end do

if (igeometry==1) return
    
    IF(irestart.EQ.1) then
       CALL read_restart
        write(*,*)'u',u(1:150,1,1,1)

       CALL INPUT_CHECK

#if 0
	write (*,*) "YYYYYYYYYOOOOOOOOOOOHHHHHHHHOOOOOOOO"
	write (*,*) from_post

		!^^^ Mohammad 10-16-2009 ^^^^^^^^^^^^^^^^^^^
		if (iturbon==1) then
			call specinit
			call transfer_vel
		endif
		!-------------------------------------------

!		if (from_post) then
			write (*,*) "UMEAN=",umean(:)
			call calc_velreal(ubcp)
			
			call post
!		endif

#endif
    END IF
    
    voldom = doml(1)*doml(2)*doml(3)
    
    umeanslip = SQRT(uchar(1)**2.d0 + uchar(2)**2.d0+uchar(3)**2.d0)/(one-maxvolfrac)
    if(I_AM_NODE_ZERO)WRITE(*,'(A,2(2x,g12.6))') 'UMEANSLIP AND MAXVOLFRAC = ', UMEANSLIP, maxvolfrac
    
    !This mean slip will have same expression regardless of velocity given to particles or fluid phase.
    t_conv = (char_length)/(umeanslip*(one-maxvolfrac))
    !min = MINVAL(
    tendused = tend*doml(1)/(umeanslip*(one-maxvolfrac))
    
    t_min = t_conv
    
    
    umax_tmp = MAX(umeanslip, MAXVAL(ABS(velbdy(1:nbody,1:ndim))))
    
    if(I_AM_NODE_ZERO) WRITE(*,'(A25,g12.6)') 'MAXIMUM VELOCITY  = ', umax_tmp
    
    dt_tmp_conv = (dx)/(umax_tmp)!/100
    
    
    dt_tmp = dt_tmp_conv 
    
    !lchar = dchar!*(one-maxvolfrac)
    lchar  = char_length
    t_vis =(lchar*lchar)/vis
    
    diff_ts = .false.
    

    
    dt_tmp_vis = large_number
    dt_tmp_diff = large_number
    
    if(t_vis-t_min.LT.small_number) then 
       if(I_AM_NODE_ZERO) print*,'T_VIS-T_MIN = ', T_VIS-T_MIN
       t_min = t_vis
       diff_ts = .true.
       lchar = ((one-maxvolfrac)*doml(1)*doml(2)*doml(3))**(one/three)
       !lchar = lchar*doml
       tendused = 0.2*tend*lchar*lchar/vis
       
       
       dt_tmp = dt_tmp_vis
    end if
    
    lchar = dx 

    if(Re.lt.one)then
       
       !dt_tmp_vis = (lchar*lchar*(one-maxvolfrac))/vis
       
       !dt_tmp_vis = t_vis/50.d0 !(lchar*lchar*(one-maxvolfrac))/vis
       
       dt_tmp_vis = Re*dbydx*dt_tmp_conv/(cfl*vfl) 
       
       !vfl is the number
       ! of time steps used to resolve the viscous time scale
       
    end if
    
    if(I_AM_NODE_ZERO) WRITE(*,'(A25,g12.6)') 'RE = tvis/tcon = ', t_vis/t_conv
    
    !t_grav = DSQRT(two*dchar/9.81)
    t_grav = DSQRT(two*char_length/9.81)
    dt_tmp_grav = LARGE_NUMBER
    if(impose_grav)dt_tmp_grav = DSQRT(two*dx/9.81)
    
    !dt_tmp_vis = (vfl*cfl*lchar*lchar)/vis
    !t_vis =( 4.d0*r*r)/vis
    
    !end if
    
    IF (iscalon.EQ.1)THEN
       CALL INITSCAL ! To Intialize the Scalar Field
       !dt_tmp_diff = dt_tmp_diff*(1.d0-maxvolfrac) 
    ELSE
       dt_tmp_diff = LARGE_NUMBER
       t_diff = LARGE_NUMBER 
    ENDIF

    dtcoll = LARGE_NUMBER
    if(nbody.gt.0)then
       XLENGTH = DOML(1)
       ZLENGTH = DOML(3)
       YLENGTH = DOML(2)
       DES_EN_INPUT(:) = coeff_rest
       DES_ET_INPUT(:) = zero
       
       CALL des_time_march(.TRUE.)
    end if
    DO m=1,nbody            ! loop over bodies
       CALL update_nrpr_array(m)
    END DO
    
    IF(irestart.eq.0) then
!!$       IF (imove.eq.1.and.(.not.(TRIM(collision_type).eq."none"))) then
!!$          dtcoll = dtsolid_orig
!!$       end IF
       
       S_TIME = ZERO 
       AUTO_CORR_SEP_TIME = ZERO
       if(I_AM_NODE_ZERO)then
          PRINT*, 'MAXIMUM OVERLAP = ', OVERLAP_MAX
          !write(*,'(A25,2x,i5)')'Default # of Rev pts',nrpr
       end if
       
       if(I_AM_NODE_ZERO)then
          WRITE(*,'(A25,2x,g12.5)')'DT_CONV = ', DT_TMP_CONV
          WRITE(*,'(A25,2x,g12.5)')'DT_VIS = ', DT_TMP_VIS
          WRITE(*,'(A25,2x,g12.5)')'DT_GRAV = ', DT_TMP_GRAV
          WRITE(*,'(A25,2x,g12.5)')'DT_DIFF = ', DT_TMP_DIFF
          WRITE(*,'(A25,2x,g12.5)')'DT_COLLISIONAL = ', DTCOLL
       end if
       dt = MIN(DT_TMP_CONV, DT_TMP_VIS, DT_TMP_GRAV,DT_TMP_DIFF)
       dt = cfl*dt
       
       dt = MIN(dt, DTCOLL)
       
       if(I_AM_NODE_ZERO) WRITE(*,'(A25,2x,g12.5)')'DT CHOSEN = ', DT
       
       cf=-1.d0/dt
       cforig = cf
       if(I_AM_NODE_ZERO) WRITE(*,*) 'cf=',cf
       if(I_AM_NODE_ZERO) WRITE(ounit,'(2x,A,g12.5)') 'cf=',cf

       
    

       
       ferror = one
       fold = zero
       
       do iphs = 1, nphases
          phase_array(iphs)%ferror = one
          phase_array(iphs)%fold = zero
          phase_array(iphs)%ferror_array(1:nerr_steps) = one
          phase_array(iphs)%ferror_hist = one
          if(imove.eq.1)then
             phase_array(iphs)%grant_error = one
             phase_array(iphs)%grant_old = gran_temp
             phase_array(iphs)%grant_array(1:nerr_steps) = one
             phase_array(iphs)%gran_error_hist = one
          end if
       end do
       
       ferror_array(:) = one
       ferror_hist = one
       source_hydro(:) = zero 
    end IF
    
    ramp_frac_time = ramp_frac_steps*ucharmod*dt/(float(mx-1)*dx)
    if(I_AM_NODE_ZERO) WRITE(*,'(A25,2x,g12.5)')'ramping time steps ', ramp_frac_steps
    if(I_AM_NODE_ZERO) WRITE(*,'(A25,2x,g12.5)')'ramp_frac_time ', ramp_frac_time

    mesh_vel(1:ndim) = usmean_des(1:ndim)
    

    IF(.NOT.XPERIODIC) THEN 
       
       DO i=1,mx-mbuffer
          g(i)=one
       ENDDO

       DO i=mx-mbuffer,mx
          IF(XPERIODIC) PRINT*,'IF IAM PRINTING FOR PERIODIC CASE, THEN C&
               &OME TO ME IN g(x) CALCULATION' 
          j=mx-i
          g(i)=(one-dcos(twopi*j/(mbuffer*2.d0)))/two
          !g(i)=one
       ENDDO
       g(mx)=zero
    ELSE 
       g(1:nx+1) = one
    end IF
    
    !-----------------------------------------------------------------------

    dx2=dx*dx

    !-----------------------------------------------------------------------
    !	define wavenumbers: Introduced by Shankar Subramaniam, 2001

    !WRITE(*,*)'Using new definition of wavenumbers'
    !WRITE(*,*)'(consistent with FFTW)'

    DO i=1,my2
       wy(i)=dcmplx(zero,twopi*DBLE(i-1)/(DBLE(my)*dy))
    ENDDO

    DO i=1,mz/2
       wz(i)=dcmplx(zero,twopi*DBLE(i-1)/(DBLE(mz)*dz))
       wz(mz+1-i)=dcmplx(zero,-twopi*DBLE(i)/(DBLE(mz)*dz))
    ENDDO

    !WRITE(*,*)'wy: ',wy
    !WRITE(*,*)'wz: ',wz
    !-----------------------------------------------------------------------
    !	calculate phase-shift coefficients
    IF(aliasflag.eq.1) THEN 
       DO k=1,mz
          DO j=1,my2
             shiftyz(j,k)=EXP((wy(j)*dy+wz(k)*dz)/three)
          ENDDO
       ENDDO
    end IF

    !-----------------------------------------------------------------------
    !	calculate diffusion coefficient at each grid point

    DO k=1,mz
       DO j=1,my2
          w2(j,k)=-dreal(wy(j)*wy(j)+wz(k)*wz(k))
       ENDDO
    ENDDO








#if 0


do k=1, mz
	do j=1, my
		fr(1,j,k,2) = sin((j-1)*dy)
		fr(1,j,k,3) = sin((k-1)*dz)
	enddo
enddo

write (*,*) "HERE..."
do count=1, skip_num
	if (mod(count,100)==0) write (*,*) count

	do idim=2, 3
		CALL ff2rc(fr(1,:,:,idim),ff(1,:,:,idim))
	enddo

	do k=1, mz
		do j=1, my2
			ff(1,j,k,2) = ff(1,j,k,2) * wy(j)
			ff(1,j,k,3) = ff(1,j,k,3) * wz(k)
		enddo
	enddo

	do idim=2, 3
		CALL ff2cr(ff(1,:,:,idim), fr(1,:,:,idim))
	enddo
enddo
#endif
#if 0
open(unit=1,file="force.dat",status="replace")
write (1,"(4(1a,1i))"), "zone j=", my, " k=", mz, " f=point"
do k=1, mz
	do j=1, my
		write (1,"(2i,2d15.7)") j,k,fr(1,j,k,2:3)
	enddo
enddo
close(1)
!stop
#endif





















    !-----------------------------------------------------------------------
    !	define maximum wavenumber for spectrum truncation

    wmax2=w2(my2,1)
    !	wy(my2)=czero
    !	wz(mz/2+1)=czero

111 FORMAT(e12.4)
112 FORMAT(e12.4,5x,e12.4)

    !-----------------------------------------------------------------------
    !	initialize boundary conditions in Fourier space

!!$    DO n=1,ndim
!!$       DO k=1,mz
!!$          DO j=1,my2
!!$
!!$             uin(j,k,n)=u(1,j,k,n)
!!$             uout(j,k,n)=u(mx,j,k,n)
!!$
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO




    if(rk3) then
       !FOR RK3 TIME STEPPING
       coef(1,1)=4.d0/15.d0
       coef(1,2)=4.d0/15.d0 !For Alpha, Diff term
       coef(1,3)=8.d0/15.d0 !for Gamma, NL term at n
       coef(1,4)=zero       !For Zeta, ONL TERM

       coef(2,1)=1.d0/15.d0
       coef(2,2)=1.d0/15.d0
       coef(2,3)=5.d0/12.d0
       coef(2,4)=-17.d0/60.d0

       coef(3,1)=1.d0/6.d0
       coef(3,2)=1.d0/6.d0
       coef(3,3)=3.d0/4.d0
       coef(3,4)=-5.d0/12.d0
       itrmax = 3
    else
       !FOR EULER TIME STEPPING 
       coef(1,1)=half
       coef(1,2)=half
       coef(1,3)=three/two
       coef(1,4)=-half
       itrmax = 1
    end if

#if PARALLEL
    vel_converted = .FALSE.
#endif
    
    if(I_AM_NODE_ZERO)then
       write(*,'(/,3(2x,A25,i4,/),5(2x,A25, g12.5,/),3(2x,A25, i4,/),2x,A25,3(g12.5),/,A25,2x,g12.5)') &
            '# of I grid pts = ', mx, &
            '# of J grid pts = ', my, &
            '# of K grid pts = ', mz, &
            'dx = ', dx, &
            'dy = ', dy, &
            'dz = ', dz, &
            'Radius(#of grid pts.) = ', R, &
            'Re = ', Re, &
            'reversl points (yes=1)', dorpr, &
            'mxf = ', mxf, &
            'foffset = ', foffset, &
            'upi = ', uchar(1:3), &
            'dia_phys = ', dia_phys
       
       
       Write(*,*) 'TEND = ', tendused
       Write(*,'(A,i10)') 'MINIMUM NO. STEPS BASED ON CURRENT TIME STEP =  ', NINT(tendused/dt)
       write(*,'(6(2x,A25,g12.5,/))') &
            'Convective Time scale =', t_conv, &
            'Viscous Time scale =', t_vis, &
            'Gravity Time Scale =', t_grav,&
            'Diffusive Time scale =', t_diff, &
            'Min. Time scale =', t_min, &
            'Time steP Size =', dt
    !   CALL write_input(ounit)
       
       write(ounit,'(3(2x,A25,i4,/),5(2x,A25, g12.5,/),3(2x,A25, i4,/),2x,A25,3(g12.5))') &
            '# of I grid pts = ', mx, &
            '# of J grid pts = ', my, &
            '# of K grid pts = ', mz, &
            'dx = ', dx, &
            'dy = ', dy, &
            'dz = ', dz, &
            'Radius(#of grid pts.) = ', R, &
            'Re = ', Re, &
            'reversl points (yes=1)', dorpr, &
            'mxf = ', mxf, &
            'foffset = ', foffset, &
            'upi = ', uchar(1:3)
       write(ounit,'(5(2x,A25,g12.5,/))') &
            'Convective Time scale =', t_conv, &
            'Viscous Time scale =', t_vis, &
            'Diffusive Time scale =', t_diff, &
            'Min. Time scale =', t_min, &
            'Time step size =', dt
    end if

    if(TRIM(input_type).eq."single-phase")then
       frmean(1:ndim) = zero
       mpg(1:ndim) = zero
    endif
    !saveitns = INT(saveitns/100.d0*dia_phys/(ucharmod*dt))
    !if(I_AM_NODE_ZERO)Write(*,*)' NEW SAVE ITNS = ', saveitns
#if 0
  contains
	subroutine transfer_vel
		implicit none

#if PARALLEL
		if (I_AM_NODE_ZERO) then
			! velocity fluctuations for node zero
			u(1:nx,:,:,:) = u(1:nx,:,:,:) + turb_uf(1:nx,:,:,:)
			
			! sending velocity fluctuations to other processes
			do iproc=1,nproc-1
				node_num=my*my2*(ends(iproc)-starts(iproc)+1)
				allocate(trans_buf(node_num))

				do idim=1,ndim
					l=0
					do k=1,my
						do j=1,my2
							do i=starts(iproc),ends(iproc)
								l=l+1
								trans_buf(l)=turb_uf(i,j,k,idim)
							enddo
						enddo
					enddo
					call mpi_send(trans_buf(1),node_num,mpi_double_complex,iproc,iproc,comm_group,err_code)	
!						write (*,*) "0 sending to ",iproc
				enddo
				deallocate(trans_buf)
			enddo
		else
			! recieving velocity fluctuations from node zero
			node_num=my*my2*nx
			allocate(trans_buf(node_num))

			do idim=1,ndim
				call mpi_recv(trans_buf(1),node_num,mpi_double_complex,node_zero,myid,comm_group,status,err_code)
!					write (*,*) myid," receiving from 0"
				l=0
				do k=1,my
					do j=1,my2
						do i=1,nx
							l=l+1
							u(i,j,k,idim) = u(i,j,k,idim) + trans_buf(l)
						enddo
					enddo
				enddo
			enddo
			deallocate(trans_buf)
		endif

!			call communicate_velocity

		do idim = 1, ndim
			i = 1
			VECSENDRECV(u(i,1,1,idim),1,ucslice,fromproc,1,u(nx+1,1,1,idim),1,toproc,1,decomp_group,status)
			i = nx
			VECSENDRECV(u(i,1,1,idim),1,ucslice,toproc,0,u(0,1,1,idim),1,fromproc,0,decomp_group,status)
		end do
#else
		u(1:nx,:,:,:) = u(1:nx,:,:,:) + turb_uf(1:nx,:,:,:)
		u(nx+1,:,:,:) = u(1,:,:,:)
#endif
		if (allocated(turb_uf)) deallocate(turb_uf)
	
	end subroutine transfer_vel
#endif
  END SUBROUTINE initflo


  SUBROUTINE generate_configuration
    USE global_data
    USE postproc_funcs
    USE hard_sphere
    USE maternmod
    USE collision_mod
    USE randomno
    USE dependent_functions
    IMPLICIT NONE

    REAL(prcn) ::  rad, temp, dia_gcg1,dia_gcg2,vol_frac1_gcg, vol_frac2_gcg, vol_tmp

    REAL(prcn), ALLOCATABLE, DIMENSION(:,:,:,:) ::  velr
    REAL(prcn) :: ldomtmp(3), x1,x2, umf0(ndim), rsf(ndim,ndim)

    INTEGER :: i,j,k,n,m,  mytmp, mbox, nbdymax, idim, tmp_unit, iphs, nsim
    REAL(PRCN) :: lambda_p,  rad_tmp, final_vol_frac, pvel_var1, pvel_var2

    REAL(prcn), dimension(:,:), ALLOCATABLE :: xc_temp, vtemp
    CHARACTER*80 :: TEMPCHAR
    CHARACTER*80 :: filename, FILENAME1
    
    REAL(prcn) :: part_sep, dummyx, dummyy, dummyz, dummyr
	integer :: dummyi

    INTEGER, dimension(:), ALLOCATABLE :: npart

	real(prcn) :: conf, confint, int_dist_avg, int_dist_var, int_dist_sd, int_dist2

	logical, allocatable :: contact(:,:)
	real(prcn) :: max_overlap

    mytmp = my
    mbox = mxf


!    nbins = 200 

    rescaling = .true.

    if(I_AM_NODE_ZERO)then
       gofunit = getnewunit(minunitno,maxunitno)

       OPEN(gofunit,file=TRIM(RUN_NAME)//'_gof_used.dat', form = 'formatted')

       gofavgunit = getnewunit(minunitno,maxunitno)
       
       OPEN(gofavgunit,file= TRIM(RUN_NAME)//'_gof_avg.dat', form = 'formatted')
    end if
    IF(nbins.GT.nbin_max) THEN
       WRITE(*,*) "INCREASE nbin_max in initialize flo "
       STOP
    end IF
    
    !radbdy(1:nbody) = 0.0625
    !allocate the radbdy array to the size of the number of bodies
    IF(irestart.eq.0) THEN 
       IF(input_type.eq."mat") THEN 

		CALL screen_separator(80,'I')
		write (*,*) "IN MATERN CARD-CORE PROCESS"

		nsim = mis_mat
		if(my.eq.undefined_I)then
			my = lybyd * dbydx
			if(MOD(my,2).ne.0) then
				my = my+1
				PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
			endif
			!if(my/2.ne.0)my = my + 1
			mz = my
			mx = my+1
			mxf = mx 
			my2 =  my/2+1
			mx1 = mx - 1 
			IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
		end if

		if (.not.allocated(gofr_avg)) then
			allocate(gofr_avg(nbins), gofr_mis(nsim, nbins))
			allocate(rad_bin(nbins))
			allocate(int_dist(nsim))

			gofr_avg = zero
			gofr_mis = zero
			rad_bin  = zero
			int_dist = zero
		endif

          doml(2) = lybyd*dia_phys
          dx = doml(2)/(my)
          dy = dx
          dz = dx 

          doml(1) = (mx-1)*dx
          doml(3) = doml(2)
          
2002      continue
          rad_tmp = half*(real(my))/lybyd
          !rad_tmp = half*(one)/lybyd
          vol_tmp =  phiavg
          lambda_p = -one/(fourthirdpi*(hhat*rad_tmp)**3)

				write (*,"(1A,1D15.7)") "RAD_TMP = ", rad_tmp
				write (*,"(1A,1D15.7)") "LAMBDA_P = ", lambda_p

          if(one.gt.(hhat**3.0)*phiavg) then 
!             lambda_p = lambda_p * log (one-fourthirdpi*(hhat*rad_tmp*2)**3*phiavg)
             lambda_p = lambda_p * log (one-(hhat)**3*phiavg)

			write (*,"(1A,1D15.7)") "LAMBDA_P = ", lambda_p
          else
             write(*,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'
             write(ounit,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'

			write (*,"(1A,1D15.7)") "HHAT = ", hhat
			write (*,"(1A,1D15.7)") "PHIAVG = ", hhat
			write (*,"(1A,1D15.7)") "(HHAT**3.0)*PHIAVG = ", (hhat**3.0)*phiavg
			write (*,"(1A,1D15.7)") "LAMBDA_P = ", lambda_p

             PARALLEL_FINISH()
             STOP
          end if

          
          if(I_AM_NODE_ZERO) then
             ldomtmp(1) = mx1 !-1
             ldomtmp(2:3) = my

             nbdymax = 2*NINT(lambda_p*ldomtmp(1)*ldomtmp(2)*ldomtmp(3))
             Write(*,*)'MAX NUMBER OF BODIES IN MATERN : ', nbdymax

			if (allocated(xc_temp)) deallocate(xc_temp)
			allocate(xc_temp(nbdymax,3))

             CALL matern(3, ibordel, Ldomtmp, lambda_p, rad_tmp, hhat, nsim, nbdymax, xc_temp, nbody)!0.006 
          end if
          BROADCAST_INT(nbody,1,NODE_ZERO,decomp_group)

          CALL ALLOC_BND_RELATED

          if(I_AM_NODE_ZERO)then
             do n = 1, ndim
                xc(1:nbody,n) = xc(1:nbody,n) + one
             enddo
             vol_frac1 = vol_tmp
             PRINT*,'NUMBER OF BODIES AFTER MATERN = ', nbody
             WRITE(ounit,*) 'NUMBER OF BODIES AFTER MATERN = ', nbody
             !STOP

             xc(1:nbody,1:ndim)  = xc_temp(1:nbody, 1:ndim)
             radbdy(1:nbody)  = rad_tmp

             open(unit=2000,file=TRIM(RUN_NAME)//'_sphr_center_mat.dat&
                  &',form='formatted',status='unknown')
             !open(unit=2000,file='sphr_center.inp',form='formatted',status='unknown')
             write(2000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',   &
                  &    ' "UX" '
             DO m=1,nbody
                WRITE(2000,'(4(2x,f12.8))') xc(m,1), xc(m,2), xc(m,3), radbdy(m)
             enddo

             CLOSE(2000, status="keep")
          end if
	     CALL ALLOC_PHASE_RELATED
          phase_array(1)%npart = nbody
          
          final_vol_frac = nbody*fourthirdpi*(rad_tmp**3)/(ldomtmp(1)*ldomtmp(2)*ldomtmp(3))

#if 0
          IF(final_vol_frac.lt.matern_treshold*phiavg) THEN 
             CALL DEALLOC_BND_RELATED
             CALL DEALLOC_PHASE_RELATED

             PRINT*,'REDOING MATERN'
             goto 2002
          end IF
#endif


          if(I_AM_NODE_ZERO)then
             PRINT*,'number of bodies after scaling= ', nbody
             Write(*,*)'final volume fraction after Matern: ', final_vol_frac
             WRITE(ounit,*)'number of bodies after scaling = ', nbody
             
             !NOw calculate the gof 
             !CALL calc_gofr(nbody, xc(1:nbody,1:3), radbdy(1:nbody), mytmp, mbox,&
             !     & xperiodic,nrbins, 3, rescaling, gofr(1:nrbins),&
             !     & rho_est(1:nrbins), rad_bin(1:nrbins)) 


			gofr_avg = zero 
			do j=1,nbins
				gofr_avg(j) = gofr_avg(j) + one/dble(nsim)*sum(gofr_mis(1:nsim,j))
			end do

			do j=1,nbins
				IF(nsim.ge.2) THEN 
					conf=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* sum((gofr_mis(1:nsim,j)- gofr_avg(j))**2))
				end IF

				write(gofavgunit,'(10(E20.10,1x))') rad_bin(j), rad_bin(j)*lybyd, gofr_avg(j), conf
			end do
			close(gofavgunit, status = "keep") 

			int_dist_avg = sum(int_dist(1:nsim))/nsim

			if (nsim>2) then
				call get_confin(nsim, confint)
				int_dist_var = sum( (int_dist(1:nsim)- int_dist_avg) **2) / nsim
				int_dist_sd  = sqrt(int_dist_var)

				conf = int_dist_sd / sqrt(float(nsim)) * confint
			else
				conf = zero
			endif

!			call int_dist_gofr(nbins, rad_bin(1:nbins), gofr_avg(1:nbins), int_dist2)

			filename = "NMIS_interparticle.dat"
			open (unit=1,file=trim(filename),status="replace",action="write")
			write (1,"(4D15.7)") final_vol_frac, int_dist_avg, conf, int_dist2
			close (1)

			IF(GOF_AVG) then 
				if(I_AM_NODE_ZERO)WRITE(*,*) 'GOF AVG IS TRUE, SO STOPPING THE SIMULATION'
				PARALLEL_FINISH()
				STOP 
			end IF
          end if
          !STOP

       ELSE IF(TRIM(input_type).eq."risers") THEN 
          
          doml(2) = widthbyd*dia_phys

          Write(*,*)'DOING RISER SIMULATION'

          if(my.eq.undefined_I)then
             my = widthbyd*dbydx
             if(MOD(my,2).ne.0) then
                my = my+1
                PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
             endif
             
             mz = my
             my2 =  my/2+1
             IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
          end if
          
          if(mx.eq.undefined_I)then
             mx = aspect_ratio*my + 1
             mxf = mx 
             mx1 = mx - 1 
             IF(I_AM_NODE_ZERO)PRINT*,'MX UNDEFINED: SETTING TO:', MX
          end if
          
          dx = doml(2)/(my)
          dy = dx
          dz = dx 
          
          doml(1) = (mx-1)*dx
          doml(3) = doml(2)
          
          nphases = 1
          CALL ALLOC_PHASE_RELATED
          
          phase_array(1)%dia = dia_phys
          phase_array(1)%volfrac = volume_fraction
          char_length = dia_phys
          

          rad_tmp = half*(my)/lybyd
          vol_tmp =  volume_fraction
          lambda_p = (-one/(fourthirdpi*(hhat**3)*(rad_tmp**3)))
          if(one.gt.(hhat**3.0)*volume_fraction) then 
             lambda_p = lambda_p*log(one-(hhat**3.0)*volume_fraction)
          else
             write(*,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'
             write(ounit,*) 'Exiting the simulation: because matern distribution not possible at these values of hhat and vol fraction'
             PARALLEL_FINISH()
             STOP
          end if
          
          if(I_AM_NODE_ZERO) then
             ldomtmp(1) = mx1 !-1
             ldomtmp(2:3) = my
!!$             ldomtmp(1) = doml(1)
!!$             ldomtmp(2:3) = doml(2:3)
             
             nbdymax = 2*NINT(lambda_p*ldomtmp(1)*ldomtmp(2)*ldomtmp(3))
             Write(*,*)'MAX Number of bodies in Matern : ', nbdymax
             
             if(.not.ALLOCATED(xc_gener))ALLOCATE(xc_gener(nbdymax,3)&
                  &, rad_gener(nbdymax))
             
20020        continue             
             CALL matern(3, ibordel,Ldomtmp,lambda_p, rad_tmp,hhat,mis_mat ,nbdymax,xc_gener,nbody)!0.006 
             final_vol_frac = nbody*fourthirdpi*(rad_tmp**3)/(ldomtmp(1)*ldomtmp(2)*ldomtmp(3))
             IF(final_vol_frac.lt.0.9*vol_frac1) THEN 
                PRINT*,'REDOING MATERN'
                goto 20020
             end IF
             Write(*,*)'number of bodies after Matern= ', nbody
             Write(*,*)'final volume fraction after Matern: ', final_vol_frac
             
#if 0             
             percent_buf(1:nphases) = zero
             min_part_sep = zero
             ALLOCATE(npart(nphases))
             npart(1) = nbody
             write(*,*)'nbody: ', nbody
             Write(*,*)'DOML(:) ', DOML(:)
             
             read(*,*)

             DO IDIM = 1, 3
                XC_GENER(1:NBODY,IDIM) = XC_GENER(1:NBODY, IDIM)&
                     &/DOML(IDIM)
             end DO
             
             rad_gener(1:NBODY)=rad_tmp/DOML(2)

             CALL scale_to_grid_units(nbody,npart(1:nphases),nphases,my,mxf,xperiodic&
                  &,percent_buf(1:nphases),xc_gener(1:nbody,1:3),rad_gener(1:nbody),&
                  & min_part_sep, toscale=.TRUE.)
             Write(*,*),'nbody: ', nbody
#endif                  
             CALL ALLOC_BND_RELATED
             do n = 1, ndim
                xc(1:nbody,n) = xc_gener(1:nbody,n) + one
             enddo
             radbdy(1:nbody) = rad_tmp
          end if

          do iphs = 1, nphases
             phase_array(iphs)%npart = nbody
          end do
          
          BROADCAST_INT(nbody,1,NODE_ZERO,decomp_group)
          
          if(.not.I_AM_NODE_ZERO)CALL ALLOC_BND_RELATED
          
          do idim = 1, ndim
             BROADCAST_DOUBLE(xc(1,idim),nbody,NODE_ZERO,decomp_group)
          end do
          BROADCAST_DOUBLE(radbdy(1),nbody,NODE_ZERO,decomp_group)

       ELSEIF(input_type.eq."random") then 
          CALL generate_psd_config
          
          do idim = 1, ndim
             BROADCAST_DOUBLE(xc(1,idim),nbody,NODE_ZERO,decomp_group)
          end do
          BROADCAST_DOUBLE(radbdy(1),nbody,NODE_ZERO,decomp_group)

       ELSEIF(input_type.eq."simple") then
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'S')
             WRITE(*,*) 'IN SIMPLE CUBIC'
          end if
          nbody = 1
          nphases = 1
          CALL ALLOC_PHASE_RELATED
          phase_array(1)%npart = nbody
          phase_array(1)%dia = dia_phys
          phase_array(1)%volfrac = vol_frac1
          char_length = dia_phys
          CALL ALLOC_BND_RELATED
          
          lybyd = (pi/(6.d0*vol_frac1))**(one/three)
          doml(2) = lybyd*dia_phys
          
          if(my.eq.undefined_I)then
             my = lybyd * dbydx
             if(MOD(my,2).ne.0) then
               my = my+1
               PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
             endif
             !if(my/2.ne.0)my = my + 1
             mz = my
             mx = my+1
             mxf = mx 
             my2 =  my/2+1
             mx1 = mx - 1 
             IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
          end if

          xc(1:nbody, 1) = mx1/two+one
          xc(1:nbody,2) =  my/2.d0+one
          xc(1:nbody,3) =  mz/2.d0+one
          !xc(1,1:3) = my/2.d0+one
          xperiodic = .true.
          
          
          dx = doml(2)/(my)
          dy = dx
          dz = dx 
          doml(1) = (mx-1)*dx
          doml(3) = doml(2)
          
          radbdy(1:nbody) = half*(my)/lybyd
          !xc(1:nbody, 1) =  radbdy(1:nbody)!mx1/two+one
          if(I_AM_NODE_ZERO)then
             PRINT*,'VOL_FRAC1: ', VOL_FRAC1
             !PRINT*,'VOL_FRAC1 in SIMPLE ACTUAL: ', nbody*
             PRINT*,'LYBYD AND RADBDY',lybyd, radbdy(1)
             CALL screen_separator(80,'S')
          end if
       ELSEIF(input_type.eq."fcc") then
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'FCC')
             WRITE(*,*) 'IN FACE CENTERED CUBIC'
          end if
          nbody = 4
          nphases = 1
          CALL ALLOC_PHASE_RELATED
          phase_array(1)%npart = nbody
          phase_array(1)%dia = dia_phys
          phase_array(1)%volfrac = vol_frac1
          char_length = dia_phys
          CALL ALLOC_BND_RELATED
          
          lybyd= (4.d0*pi/(6.d0*vol_frac1))**(one/three)
          
          if(my.eq.undefined_I)then
             my = lybyd * dbydx
             if(MOD(my,2).ne.0) then
               my = my+1
               PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
             endif
             !if(my/2.ne.0)my = my + 1
             mz = my
             mx = my+1
             mxf = mx 
             my2 =  my/2+1
             mx1 = mx - 1 
             IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
          end if

          xc(1,1:3) = 1.d0
          xc(2,1) = my/2.d0+one
          xc(2,2) = my/2.d0+one
          xc(2,3) = 1.d0
          xc(3,1) = my/2.d0+one
          xc(3,2) = 1.d0
          xc(3,3) = my/2.d0+one
          xc(4,1) = 1.d0
          xc(4,2) = my/2.d0+one
          xc(4,3) = my/2.d0+one

          xperiodic = .true.
          
          doml(2) = lybyd*dia_phys
          dx = doml(2)/(my)
          dy = dx
          dz = dx 
          doml(1) = (mx-1)*dx
          doml(3) = doml(2)
          
          radbdy(1:nbody) = half*(my)/lybyd
          if(I_AM_NODE_ZERO)then
             PRINT*,'VOL_FRAC1: ', VOL_FRAC1
             !PRINT*,'VOL_FRAC1 in SIMPLE ACTUAL: ', nbody*
             PRINT*,'LYBYD AND RADBDY',lybyd, radbdy(1)

             CALL screen_separator(80,'FCC')
          end if

       ELSEIF(input_type.eq."default") then
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'D')
             WRITE(*,*) 'IN DEFAULT FRESH START'
             
             tmp_unit =  getnewunit(minunitno,maxunitno)
             OPEN(tmp_unit,file=TRIM(RUN_NAME)//'_sphr_center.inp',status='old', Action='read')

             vol_frac1 = zero
             nbody = 0
            READ(tmp_unit,*)
125          continue
             READ(tmp_unit,*,END=225) dummyx, dummyy, dummyz, dummyr !, dummyi !M ,dummyx,dummyy,dummyz,dummyx

!if (dummyy<=lybyd/8.or.dummyy>=lybyd*7/8) then
!				if (dummyx<lybyd.and.dummyy<lybyd.and.dummyz<lybyd) nbody = nbody + 1
!endif
             goto 125
225          continue
             CLOSE(tmp_unit, status='keep')
             Write(*,*)'Number of bodies in the DEFAULT CASE are : ', nbody
          end if
          BROADCAST_INT(nbody,1,NODE_ZERO,decomp_group)
          
          IF (nbody.NE.0) THEN
             CALL ALLOC_BND_RELATED
             if(I_AM_NODE_ZERO)then
                OPEN(tmp_unit,file=TRIM(RUN_NAME)//'_sphr_center.inp',status='old', Action='read')
                WRITE(*,*) 'Reading immersed surface data'
               READ(tmp_unit,*)
					i=0
					do
						READ(tmp_unit,*) dummyx, dummyy, dummyz, dummyr !, dummyi

!if (dummyy<=lybyd/8.or.dummyy>=lybyd*7/8) then
!						if (dummyx<lybyd.and.dummyy<lybyd.and.dummyz<lybyd) then 
							i=i+1
							xc(i,1) = dummyx*dbydx
							xc(i,2) = dummyy*dbydx
							xc(i,3) = dummyz*dbydx
!							color(i) = dble(dummyi)

							velbdy(i,1:ndim) = zero
!M                   READ(tmp_unit,*) xc(i,1),xc(i,2),xc(i,3), radbdy(i), &
!M                        &velbdy(i,1), velbdy(i,2), velbdy(i,3)
!M                   vol_frac1 = vol_frac1 + (pi*(dia_phys**3.d0))/(6.d0)
!						endif
!endif
						if (i==nbody) exit
					end do
                
                CLOSE(tmp_unit,STATUS='keep')
!                dummyx = MINVAL(xc(1:nbody,1))
!                do m = 1, nbody
!                   xc(m, 1) = xc(m, 1) - dummyx + 1
!                end do
             end if
             nphases = 1
             CALL ALLOC_PHASE_RELATED
             phase_array(1)%npart = nbody
             phase_array(1)%dia = dia_phys
             phase_array(1)%volfrac = vol_frac1
             char_length = dia_phys
          end if
          
          do idim = 1, ndim
             BROADCAST_DOUBLE(xc(1,idim),nbody,NODE_ZERO,decomp_group)
             BROADCAST_DOUBLE(velbdy(1,idim),nbody,NODE_ZERO,decomp_group)
          end do
#if 1
          ALLOCATE(npart(nphases))
          do iphs = 1, nphases
             npart(iphs) = nbody
             percent_buf(iphs) = zero
          end do
#endif             
          if(I_AM_NODE_ZERO) PRINT*,'PART VOLUME = ', VOL_FRAC1
          doml(2) = lybyd*dia_phys
          
          if(my.eq.undefined_I)then
             my = lybyd * dbydx
             if(MOD(my,2).ne.0) then
                my = my+1
                PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
             endif
             !if(my/2.ne.0)my = my + 1
             mz = my
             mx = my+1
             mxf = mx 
             my2 =  my/2+1
             mx1 = mx - 1 
             IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
          end if

         !Convert the centers to lie between 0 and 1
          
          do m = 1, nbody
!             xc(m,1:ndim) = xc(m,1:ndim) * real(my,prcn)/(lybyd*dbydx)
             radbdy(m) = half*real(my,prcn)/lybyd !M half*one/lybyd
          end do

          dx = doml(2)/(my)
          dy = dx
          dz = dx 
          doml(1) = (mx-1)*dx
          doml(3) = doml(2)
#if 0
          CALL scale_to_grid_units(nbody,npart(1:nphases),nphases,my,mxf,xperiodic&
               &,percent_buf(1:nphases),xc(1:nbody,1:3),radbdy(1:nbody),&
               & min_part_sep, toscale=.TRUE.) 
          
          do iphs = 1, nphases
             phase_array(iphs)%npart = npart(iphs)
          end do
#endif             !STOP
          
          if(I_AM_NODE_ZERO)then
             WRITE(*,*) 'OUT OF DEFAULT SPHR CONFIG'
             CALL screen_separator(80,'D')
          end if
          
       ELSEIF(input_type.eq."agglomerate") then
          
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'A')
             WRITE(*,*) 'INPUT TYPE AGGLOMERATE FRESH START'
             
             tmp_unit =  getnewunit(minunitno,maxunitno)
             OPEN(tmp_unit,file=TRIM(RUN_NAME)//'_sphr_center.inp',status='old', Action='read')
             
             vol_frac1 = zero
             nbody = 0
             !READ(tmp_unit,*)
1250         continue
             READ(tmp_unit,*,END=2250) dummyx, dummyy, dummyz
             nbody = nbody + 1
             goto 1250
2250         continue
             CLOSE(tmp_unit, status='keep')
             Write(*,*)'Number of bodies in the DEFAULT CASE are : ', nbody
          end if
          BROADCAST_INT(nbody,1,NODE_ZERO,decomp_group)
          
          doml(2) = lybyd*dia_phys
          if(my.eq.undefined_I)then
             my = lybyd * dbydx
             if(MOD(my,2).ne.0) then
                my = my+1
                PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
             endif
             !if(my/2.ne.0)my = my + 1
             mz = my
             mx = my+1
             mxf = mx 
             my2 =  my/2+1
             mx1 = mx - 1 
             IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
          end if
          
          dx = doml(2)/(my)
          dy = dx
          dz = dx 
          doml(1) = (mx-1)*dx
          doml(3) = doml(2)
          
          IF (nbody.NE.0) THEN
             CALL ALLOC_BND_RELATED
             if(I_AM_NODE_ZERO)then
                OPEN(tmp_unit,file=TRIM(RUN_NAME)//'_sphr_center.inp',status='old', Action='read')
                WRITE(*,*) 'Reading immersed surface data'
                !READ(tmp_unit,*)
                do i = 1, nbody
                   READ(tmp_unit,*) xc(i,1),xc(i,2),xc(i,3)
                   velbdy(i,1:ndim) = zero
                end do
                
                CLOSE(tmp_unit,STATUS='keep')
                
                do m = 1, nbody
                   do idim = 1, ndim
                      xc(m, idim) = xc(m, idim)/DOML(idim) + half
                   end do
                end do
             end if
             nphases = 1
             CALL ALLOC_PHASE_RELATED
             phase_array(1)%npart = nbody
             phase_array(1)%dia = dia_phys
             phase_array(1)%volfrac = vol_frac1
             char_length = dia_phys
          end if
          
          do idim = 1, ndim
             BROADCAST_DOUBLE(xc(1,idim),nbody,NODE_ZERO,decomp_group)
             BROADCAST_DOUBLE(velbdy(1,idim),nbody,NODE_ZERO,decomp_group)
          end do
          
          ALLOCATE(npart(nphases))
          do iphs = 1, nphases
             npart(iphs) = nbody
             percent_buf(iphs) = zero
          end do
          
          do m = 1, nbody
             radbdy(m) = half*(one/lybyd)
          end do
          
          
          CALL scale_to_grid_units(nbody,npart(1:nphases),nphases,my,mxf,xperiodic&
               &,percent_buf(1:nphases),xc(1:nbody,1:3),radbdy(1:nbody),&
               & min_part_sep) 
          
          do iphs = 1, nphases
             phase_array(iphs)%npart = npart(iphs)
          end do
          
          if(I_AM_NODE_ZERO)then
             WRITE(*,*) 'OUT OF AGGLOMERATE FRESH START'
             CALL screen_separator(80,'A')
          end if
       else if(TRIM(input_type).eq.'lubtest')then
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'L')
             WRITE(*,*) 'IN LUBRICATION THEORY TEST --> FRESH START'
!!$             tmp_unit =  getnewunit(minunitno,maxunitno)
!!$             OPEN(tmp_unit,file=TRIM(RUN_NAME)//'_lubtest.inp',status='old', Action='read')
!!$             READ(tmp_unit,*)nphases
          end if
          nphases = 2
          CALL ALLOC_PHASE_RELATED
          ALLOCATE(npart(nphases))
          do iphs = 1, nphases
             phase_array(iphs)%npart = 1
             npart(iphs) = 1
             percent_buf(iphs) = zero
          end do
          phase_array(1)%dia = dia_phys
          phase_array(2)%dia = phase_array(1)%dia*dia_ratio
          
          char_length = two*(phase_array(1)%dia*phase_array(2)%dia)/(phase_array(1)%dia+phase_array(2)%dia)

          DOML(2) = lybyd*phase_array(2)%dia
          
          DOML(1) = DOML(2)
          DOML(3) = DOML(2)
          
          nbody = 2
          CALL ALLOC_BND_RELATED
          
          do m = 1, nbody
             radbdy(m) = phase_array(m)%dia/two
             xc(m,2) = DOML(2)/two
             xc(m,3) = DOML(2)/two
          end do
          part_sep = hbyd*phase_array(2)%dia
          
          xc(1,1) = (DOML(2)/two - part_sep/two - radbdy(1))
          xc(2,1) = (DOML(2)/two + part_sep/two + radbdy(2))
          
          do m = 1, nbody
             radbdy(m) = radbdy(m)/DOML(2)
             do idim = 1, ndim
                xc(m,idim) = xc(m,idim)/DOML(idim)
             end do
          end do
          
          if(I_AM_NODE_ZERO)then
             Write(*,'(A25,2x,g17.8)')'LENGTH OF THE BOX IS',DOML(2)
             Write(*,'(A25,2(2x,g17.8))')'RADII ',radbdy(1), radbdy(2)
             Write(*,'(A25,2(2x,g17.8))')'PHYSICAL LOCATIONS ',xc(1,1), xc(2,1)
          end if
          
          if(my.eq.undefined_I)then
             my = doml(2)/phase_array(1)%dia * dbydx
             if(MOD(my,2).ne.0) then
                my = my+1
                PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
             endif
             !if(my/2.ne.0)my = my + 1
             mz = my
             mx = my+1
             mxf = mx 
             my2 =  my/2+1
             mx1 = mx - 1 
             IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
          end if

          
          dx = doml(2)/real(my,prcn)
          
          dy = dx
          dz = dx
          
          !IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
          
          min_part_sep = zero
          hbydx = part_sep/dx
          Write(*,'(A25,2(2x,g17.8))')'SEPARATION IN GRID UNITS : ',hbydx
          CALL scale_to_grid_units(nbody,npart(1:nphases),nphases,my,mxf,xperiodic&
               &,percent_buf(1:nphases),xc(1:nbody,1:3),radbdy(1:nbody),&
               & min_part_sep) 
          
          if(I_AM_NODE_ZERO)then
             PRINT*,'RADII IN GRID UNITS: ',radbdy(1), radbdy(2)
             CALL screen_separator(80,'L')
          end if
       else if(TRIM(input_type).eq.'single-phase')then
          DOML(1:ndim) = two*pi
!          CALL screen_separator(80,'T')
!          WRITE(*,*) 'IN SINGLE PHASE TURBULENCE RUN --> FRESH START'
!          if(my.eq.undefined_I)then
!             Write(*,'(A)')'MY IS UNDEFINED FOR SINGLE-PHASE RUN. DEFINE IT In floparam.in AND TRY AGAIN.'
!             PARALLEL_FINISH()
!             STOP
!          end if
!          mz = my
!          mx = my+1
!          mxf = mx 
!          my2 =  my/2+1
!          mx1 = mx - 1 

          dx = doml(2)/real(my,prcn)
          dy = dx
          dz = dx
          nbody = 0
          nphases = 0
          char_length = DOML(2)
       end IF!ihardsphere 
       
    ELSE
       if(I_AM_NODE_ZERO)then
          open(unit=1050,file=TRIM(RUN_NAME)//'_RESTART',form="formatted"      &
               ,status="unknown")    
          read(1050, *) count_restart
          close(1050, status='keep')

          if(imove==1 .and. .not.moving_from_fixed_bed) then
             WRITE(FILENAME1,'(I1)')count_restart
             
             open(1050,file=TRIM(RUN_NAME)//'_sphere_config_'//TRIM(FILENAME1)//'.rst',form=rstsave, status="unknown")
          else  
             open(1050, file = TRIM(RUN_NAME)//'_sphere_config.rst', form="unformatted", status="old")
          end if
          
          read(1050) nbody
          if(TRIM(input_type).eq."random".or.TRIM(input_type).eq.'risers')then

             READ(1050)nphases

             CALL ALLOC_PHASE_RELATED
             do iphs = 1, nphases
                READ(1050) phase_array(iphs)%npart
             end do
             do iphs = 1, nphases
                READ(1050) phase_array(iphs)%dia
             end do

             READ(1050)DOML(1:ndim)
		

             WRITE(*,*) 'FROM RESTART: DOML = ', DOML(1:3)
             if(imove.eq.1)then

					!^^^^^^^ 06-24-2010 Mohammad: RESTARTING A FIXED BED FOR MOVING CASE ^^^^          
					if (moving_from_fixed_bed) then
						MOVE_PARTICLES = .true.
					else
						read(1050)MOVE_PARTICLES
					endif
					!------------------------------------------------------------------------

             end if
          else if(TRIM(input_type).eq."lubtest")then
             nphases = 2
             CALL ALLOC_PHASE_RELATED
          else if(TRIM(input_type).eq."single-phase")then
             nphases = 0
          else
             nphases = 1
             CALL ALLOC_PHASE_RELATED
          end IF
          READ(1050) char_length
          WRITE(*,*) 'FROM RESTART: NBODY = ', nbody
          WRITE(ounit,*) 'FROM RESTART: NBODY = ', nbody
       end if
       BROADCAST_INT(count_restart,1,NODE_ZERO,decomp_group)
       BROADCAST_INT(nbody,1,NODE_ZERO,decomp_group)
       BROADCAST_INT(nphases,1,NODE_ZERO,decomp_group)
       BROADCAST_DOUBLE(char_length,1,NODE_ZERO,decomp_group)
#if PARALLEL
       if(TRIM(input_type).eq.'random'.or.TRIM(input_type).eq.'risers') BROADCAST_DOUBLE(DOML(1),ndim,NODE_ZERO,decomp_group)
#endif
       if(.not.I_AM_NODE_ZERO)then
          if(nphases.gt.0)CALL ALLOC_PHASE_RELATED
       end if
       IF (nbody.NE.0) THEN
          CALL ALLOC_BND_RELATED
          if(I_AM_NODE_ZERO)then
             READ(1050) xc(1:nbody,1:3)
             
             READ(1050) radbdy(1:nbody)
             READ(1050) velbdy(1:nbody,1:3)

		     !^^^^^^^ 06-24-2010 Mohammad: RESTARTING A FIXED BED FOR MOVING CASE ^^^^          
!            if(imove.eq.1) READ(1050) frame_vel(1:ndim)
!            if(imove.eq.1) READ(1050) frame_pos(1:ndim)

		     if (imove.eq.1) then
				if (moving_from_fixed_bed) then
					frame_vel(1:ndim) = 0d0
					frame_pos(1:ndim) = 0d0
				else
					read (1050) frame_vel(1:ndim)
					read (1050) frame_pos(1:ndim)
				endif
			endif
			!------------------------------------------------------------------------


          end if
       end IF
       if(I_AM_NODE_ZERO) close(1050, status = "keep")
       IF (nbody.NE.0) THEN
          do idim = 1, ndim
             BROADCAST_DOUBLE(xc(1,idim),nbody,NODE_ZERO,decomp_group)
             BROADCAST_DOUBLE(velbdy(1,idim),nbody,NODE_ZERO,decomp_group)
          end do
          BROADCAST_DOUBLE(radbdy(1),nbody,NODE_ZERO,decomp_group)
       END IF

       IF(input_type.eq."random".or.TRIM(input_type).eq.'risers') then
          do iphs = 1, nphases
             BROADCAST_INT(phase_array(iphs)%npart,1,NODE_ZERO,decomp_group)
             BROADCAST_DOUBLE(phase_array(iphs)%dia,1,NODE_ZERO,decomp_group)
          end do
          
       ELSEIF(input_type.eq."simple") then
          xperiodic = .true.
          phase_array(1)%npart = nbody
          phase_array(1)%dia = dia_phys
          phase_array(1)%volfrac = vol_frac1
          lybyd= (pi/(6.d0*vol_frac1))**(one/three)
          if(I_AM_NODE_ZERO)PRINT*,'LYBYD IN SIMPLE',lybyd
          doml(2) = lybyd*dia_phys
          
       ELSEIF(input_type.eq."fcc") then
          xperiodic = .true.
          lybyd= (4.d0*pi/(6.d0*vol_frac1))**(one/three)
          phase_array(1)%npart = nbody
          phase_array(1)%dia = dia_phys
          phase_array(1)%volfrac = vol_frac1
          if(I_AM_NODE_ZERO) PRINT*,'LYBYD IN FCC',lybyd
          doml(2) = lybyd*dia_phys
       ELSEIF(input_type.eq."default") then 
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'D')
             WRITE(*,*) 'IN DEFAULT SPHR CONFIG RESTART'
          end if
          vol_frac1 = zero 
          do i = 1, nbody
             vol_frac1 = vol_frac1 + (pi*(dia_phys**3.d0))/(6.d0)
#if 0
             if(I_AM_NODE_ZERO)then
                WRITE(*,*) 'NBODY = ', NBODY
                WRITE(*,*) 'XC = ', XC(1,:)
                write(*,*) 'RADBDY and VELBDY =', RADBDY(I), VELBDY(I,:)
             end if
#endif
          end do
          phase_array(1)%npart = nbody
          phase_array(1)%dia = dia_phys
          phase_array(1)%volfrac = vol_frac1
          
          doml(2) = lybyd*dia_phys
          
          PRINT*,'PART VOL =  ', VOL_FRAC1
          
          if(I_AM_NODE_ZERO)then
             WRITE(*,*) 'OUT OF DEFAULT SPHR CONFIG RESTART'
             CALL screen_separator(80,'D')
          end if
       ELSEIF(input_type.eq."agglomerate") then 
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'A')
             WRITE(*,*) 'IN AGGLOMERATE RESTART'
          end if
          vol_frac1 = zero 
          do i = 1, nbody
             vol_frac1 = vol_frac1 + (pi*(dia_phys**3.d0))/(6.d0)
#if 0
             if(I_AM_NODE_ZERO)then
                WRITE(*,*) 'NBODY = ', NBODY
                WRITE(*,*) 'XC = ', XC(1,:)
                write(*,*) 'RADBDY and VELBDY =', RADBDY(I), VELBDY(I,:)
             end if
#endif
          end do
          phase_array(1)%npart = nbody
          phase_array(1)%dia = dia_phys
          phase_array(1)%volfrac = vol_frac1
          
          doml(2) = lybyd*dia_phys
          
          PRINT*,'PART VOL =  ', VOL_FRAC1
          
          if(I_AM_NODE_ZERO)then
             WRITE(*,*) 'OUT OF AGGLOMERATE RESTART'
             CALL screen_separator(80,'A')
          end if
       ELSEIF(input_type.eq."lubtest") then 
          if(I_AM_NODE_ZERO)then
             CALL screen_separator(80,'L')
             WRITE(*,*) 'IN LUBRICATION TEST RESTART'
          end if
          
          do i = 1, nbody
             if(I_AM_NODE_ZERO)then
                WRITE(*,*) 'NBODY = ', NBODY
                WRITE(*,*) 'XC = ', XC(1,:)
                write(*,*) 'RADBDY and VELBDY =', RADBDY(I), VELBDY(I,:)
             end if
          end do
          do iphs = 1, nphases
             phase_array(iphs)%npart = 1
          end do
          phase_array(1)%dia = dia_phys
          phase_array(2)%dia = phase_array(1)%dia*dia_ratio
          
          DOML(2) = lybyd*phase_array(2)%dia
          
          if(I_AM_NODE_ZERO)then
             WRITE(*,*) 'OUT OF LUBRICATION TEST RESTART'
             CALL screen_separator(80,'L')
          end if
       else if(TRIM(input_type).eq.'single-phase')then
          DOML(1:ndim) = two*pi
          CALL screen_separator(80,'T')
          WRITE(*,*) 'IN SINGLE PHASE TURBULENCE RUN --> RESTART'
          if(my.eq.undefined_I)then
             Write(*,'(A)')'MY IS UNDEFINED FOR SINGLE-PHASE RUN. DEFINE IT In floparam.in AND TRY AGAIN.'
             PARALLEL_FINISH()
             STOP
          end if
          mz = my
          mx = my+1
          mxf = mx 
          my2 =  my/2+1
          mx1 = mx - 1 
       end IF
       
       if(my.eq.undefined_I)then
          !if((TRIM(input_type).eq.'random').and.(TRIM(psd_type).eq.'bidisp'))then
          if((TRIM(input_type).eq.'random').and.(TRIM(psd_type).eq.'csd'))then !Need to be done correctly for CSD
             my = DOML(2)/phase_array(1)%dia * dbydx
          else if(TRIM(input_type).eq.'risers')then
             my = widthbyd*dbydx
          else if(TRIM(input_type).eq.'lubtest')then
             my = DOML(2)/phase_array(1)%dia * dbydx
          else
             my = lybyd * dbydx
          end if
          if(MOD(my,2).ne.0) then
             my = my+1
             PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
          endif
          !if(my/2.ne.0)my = my + 1
          mz = my
          my2 =  my/2+1
          
          IF(I_AM_NODE_ZERO)PRINT*,'MY UNDEFINED: SETTING TO:', MY
       end if
       
       if(mx.eq.undefined_I)then
          if(TRIM(input_type).eq.'risers')then
             mx = aspect_ratio*my + 1
          else
             mx = my+1
          end if
          mxf = mx 
          mx1 = mx - 1
       end if
       
       dx = doml(2)/my
       doml(1) = (mx-1)*dx
       doml(3) = doml(2)
       dy = dx
       dz = dx
       if(FROM_POST)Write(*,*)'DONE WITH RESTART'
#if 0
	if (nbody>0) then
!		nbins = 200
		if (.not.allocated(gofr_avg)) then
			nsim = 1
			allocate(gofr_avg(nbins), gofr_mis(nsim, nbins))
			allocate(rad_bin(nbins))
			allocate(int_dist(nsim))

			gofr_avg = zero
			gofr_mis = zero
			rad_bin  = zero
			int_dist = zero
		endif

		WRITE(*,*) 'GENERATING GOFR' 
		if (.not.allocated(contact)) allocate(contact(nbody,nbody))
		CALL calculate_gofr_homog(nbody,xc(1:nbody,1:3), contact, my, mxf, nbins, .true., gofr(1:nbins), rad_bin(1:nbins), max_overlap)

		do j = 1, nsim
			gofr_mis(j,:) = gofr(:)
		end do

		gofr_avg = zero 
		do j=1,nbins
			gofr_avg(j) = gofr_avg(j) + one/dble(nsim)*sum(gofr_mis(1:nsim,j))
		end do
     
		gofavgunit = getnewunit(minunitno,maxunitno)
		open(gofavgunit,file= trim(run_name)//'_gof_avg.dat', form = 'formatted')


		do j=1,nbins
			write(gofavgunit,'(3d15.7)') rad_bin(j), rad_bin(j)*lybyd, gofr_avg(j)
		end do
		close(gofavgunit, status = "keep") 
	endif

#endif       
    ENDIF!irestart
    
    if(I_AM_NODE_ZERO.and.nbody.ne.0)then
       open(unit=2000,file=TRIM(RUN_NAME)//'_sphr_center_out.dat',form='formatted',status='unknown')
       write(2000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',   &
            &    ' "UX" ', ' "IBDY" '
!!$       write(2000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',   &
!!$            &    ' "UX" ', ' "UY" ', ' "UZ" ', ' "PHI" '
       
       DO m=1,nbody
          !WRITE(2000,'(4(2x,g17.8),2x,i3)') xc(m,1), xc(m,2), xc(m,3), radbdy(m),m
          WRITE(2000,'(6(2x,g17.8))') xc(m,1), xc(m,2), xc(m,3), radbdy(m)
       enddo
       
       CLOSE(2000, status="keep")
    end if
  end subroutine generate_configuration

  SUBROUTINE generate_particle_velocities
    USE global_data
    IMPLICIT NONE
    
    REAL(PRCN) :: lambda_p, gran_energy_out,&
         & diff_vel_spec1(ndim), diff_vel_spec2(ndim), temp_spec1,&
         & temp_spec2, mean_gt, umean_test1(ndim), umean_test2(ndim),&
         & mean_vel(ndim),gran_energy_out1,gran_energy_out2,&
         & var_test(ndim), hist(20), fmin, fmax, umf0(ndim),&
         & rsf(ndim,ndim), rsf1(ndim,ndim), rsf2(ndim,ndim), x1, x2,&
         & vfrac_mean,variance1,variance2
    REAL(prcn) :: gran_temp_spec(nphases),diff_vel_spec(nphases,ndim)&
         &, constant1(ndim), constant2, test_mix_var,&
         & test_mixmean_vmag, test_mean_spec_vel(nphases,ndim),&
         & test_mixmean_vel(ndim), test_spec_var(nphases)

    REAL(prcn), dimension(:,:), ALLOCATABLE ::  vtemp,vtemp1,vtemp2&
         &,veltemp

    REAL(prcn), dimension(:), ALLOCATABLE ::  ftemp, wt 

    INTEGER ::  unit1, ncount(ndim), count_temp(ndim),idim, i, j, iphs, pstart, pend, m

    velbdy(1:nbody,:) = zero

    !if(nphases.gt.1)mean_vel_to_particles = .TRUE.
    !Assign Mean velocity to the particles
    if(mean_vel_to_particles)then
       if(TRIM(psd_type).eq."discrete".or.TRIM(psd_type).eq."bidisp")then
          
          if(equal_momentum)then
             constant2 = zero
             do iphs = 1, nphases
                constant2 = constant2 + phase_array(iphs)%npart&
                     &/phase_array(iphs)%volfracg
             end do
             constant2 = constant2*maxvolfrac/((one-maxvolfrac)&
                  &*real(nbody,prcn)) + real(nphases,prcn)/maxvolfrac
             do idim = 1, ndim
!!$                constant1(idim) = maxvolfrac*uchar(idim)/((one&
!!$                     &-maxvolfrac)*real(nphases,prcn))
                constant1(idim) = uchar(idim)/((one-maxvolfrac)&
                     &*(constant2))                   
                
             end do
             do iphs = 1, nphases
                do idim = 1, ndim
                   phase_array(iphs)%mean_spec_vel(idim) = constant1(idim)/phase_array(iphs)%volfracg
                end do
             end do
          else
             do iphs = 1, nphases
                do idim = 1, ndim
                   phase_array(iphs)%mean_spec_vel(idim) =&
                        & -uchar(idim)/(one-maxvolfrac) ! Equal
                   ! velocities
!!$                   phase_array(iphs)%mean_spec_vel(idim) =&
!!$                        & uchar(idim) ! Equal velocities
                   
                end do
             end do
          end if
          pstart = 1
          do iphs = 1, nphases
             do idim=1,ndim
                diff_vel_spec(iphs,idim) = phase_array(iphs)%mean_spec_vel(idim)- uchar(idim)
             end do
          end do
          
          do m = 1, nbody
             iphs = part_array(m)%iphs
             velbdy(m,1:ndim) = phase_array(iphs)%mean_spec_vel(1:ndim)
          end do

       else if(psd_type.eq."csd")then
          
          do idim = 1, ndim  
             velbdy(1:nbody,idim) = uchar(idim)!/(one-maxvolfrac) ! Assigning equal velocities to all the spheres for the time being. Need to decide on Re_m(r).
          end do
       end if
    else
       do iphs = 1, nphases
          do idim = 1, ndim
             phase_array(iphs)%mean_spec_vel(idim) = zero
          end do
       end do
       velbdy(1:nbody,1:ndim) = zero
    end if

    ! Generate velocity distribution to the particles

    if(ReT.gt.SMALL_NUMBER)then
       gran_temp = (ReT*vis/dia_phys)**2.d0
       if(I_AM_NODE_ZERO)then
          if(psd_type.eq."csd")then
             do j=1,ndim
                umf0(j)=0.0
                do i=1,ndim
                   if(i.eq.j)then
                      rsf(i,j)= gran_temp
                   else
                      rsf(i,j)=0.0
                   endif
                enddo
             enddo
             ALLOCATE(vtemp(nbody,ndim))
             CALL jn_dist(vtemp,SIZE(vtemp,1),ndim,umf0,rsf) ! Maxwellian
             ! distribution for the particle velocity fluctuations
             do m = 1, nbody
                velbdy(m,1:ndim) = velbdy(m,1:ndim) + vtemp(m,1:ndim)
             end do
             DEALLOCATE(vtemp)
          else if(TRIM(psd_type).eq."discrete".or.TRIM(psd_type).eq."bidisp")then
!!$	    constant2 = three*SUM(volfrac(1:nphases))*gran_temp
             constant2 = mean_volfrac*gran_temp
             do iphs = 1, nphases
                 constant2 = constant2 !-  DOT_PRODUCT(diff_vel_spec(iphs,1:ndim),diff_vel_spec(iphs,1:ndim))*volfrac(iphs)
             end do
             constant2 = constant2/real(nphases,prcn)
             pstart = 1
             do iphs = 1, nphases
                gran_temp_spec(iphs) = constant2/(phase_array(iphs)%volfrac)
!!$	       gran_temp_spec(iphs) = constant2/(three*volfrac(iphs))
                do j=1,ndim
                   umf0(j) = phase_array(iphs)%mean_spec_vel(j)
                   do i=1,ndim
                      if(i.eq.j)then
                         rsf(i,j)= gran_temp_spec(iphs)
                      else
                         rsf(i,j)=0.0
                      endif
                   enddo
                enddo
                ALLOCATE(vtemp(phase_array(iphs)%npart,ndim))
                CALL jn_dist(vtemp,SIZE(vtemp,1),ndim,umf0,rsf) ! Maxwellian
                pend = pstart + phase_array(iphs)%npart - 1
                i = 1
                do m = pstart, pend
                   velbdy(m,1:ndim) = vtemp(i,1:ndim)
                   i = i+1
                end do
                pstart = pend + 1
                DEALLOCATE(vtemp)
             end do
 	       endif
       endif
    endif

    do idim = 1, ndim
       BROADCAST_DOUBLE(velbdy(1,idim),nbody,NODE_ZERO,decomp_group)
    end do
  end SUBROUTINE generate_particle_velocities

  SUBROUTINE INPUT_CHECK
    IMPLICIT NONE
    REAL(prcn) :: gran_temp_spec(nphases),diff_vel_spec(nphases,ndim)&
         &, constant1(ndim), constant2, test_mix_var,&
         & test_mixmean_vmag, test_mean_spec_vel(nphases,ndim),&
         & test_mixmean_vel(ndim), test_spec_var(nphases),&
         & volfrac(nphases), vfracmean, vbox,&
         & dia(nphases), diag(nphases)
    INTEGER :: pstart, pend, iphs, m, idim, npart(nphases)
    
    vbox = doml(1)*doml(2)*doml(3)
    volfrac = zero
    npart = 0
    do m = 1, nbody
       iphs = part_array(m)%iphs
       npart(iphs) = npart(iphs)+1
       dia(iphs) = two*radbdy(m)*dx
       diag(iphs) = two*radbdy(m)
       volfrac(iphs) = volfrac(iphs) + pi*(dia(iphs))**3.d0/6.d0
    end do

    do iphs = 1, nphases
       volfrac(iphs) = volfrac(iphs)/vbox
    end do

    vfracmean = SUM(volfrac(1:nphases))
    if(I_AM_NODE_ZERO)then
       CALL screen_separator(80,'I')
       if(.not.TRIM(psd_type).eq.'csd')then
          Write(*,'((A25,2x,I6))')'Number of PHASES = ', nphases
          do iphs = 1, nphases
             Write(*,'((A25,2x,I6, A5,2(2x,I8)))')'No. of Particles in phase ', iphs,&
                  & ' = ', npart(iphs),(phase_array(iphs)%npart)
          end do
          
          do iphs = 1, nphases
             Write(*,'((A25,2x,I6, 2x, A5,2(2x,g17.8)))')'Diameter of phase ', iphs, ' = ', dia(iphs), diag(iphs)
          end do

          do iphs = 1, nphases
             Write(*,'((A25,2x,I6, 2x, A5,2(2x,g17.8)))')'Volume fraction &
                  &of phase ', iphs, ' = ', volfrac(iphs)&
                  &,(phase_array(iphs)%volfracg)
          end do
       end if
       Write(*,'((A25,2(2x,g17.8)))')'Mean Volume fraction = ', vfracmean, maxvolfrac
    end if
    
    if(TRIM(psd_type).eq."discrete".or.TRIM(psd_type).eq."bidisp")then
       test_mean_spec_vel(1:nphases,1:ndim) = zero
       test_mixmean_vel = zero
       pstart = 1
       do iphs = 1, nphases
          pend = pstart + phase_array(iphs)%npart- 1
          do m = pstart, pend
             test_mean_spec_vel(iphs,1:ndim) = test_mean_spec_vel(iphs,1:ndim) + velbdy(m,1:ndim)
          end do
          pstart = pend + 1
          test_mean_spec_vel(iphs,1:ndim) = test_mean_spec_vel(iphs,1:ndim)/real(phase_array(iphs)%npart, prcn)
          test_mixmean_vel(1:ndim) =  test_mixmean_vel(1:ndim) + test_mean_spec_vel(iphs,1:ndim)*phase_array(iphs)%volfracg/maxvolfrac
       end do
       do idim = 1, ndim
          test_mixmean_vel(idim) = (one-maxvolfrac)*(test_mixmean_vel(idim)-&
               & ufmean_des(idim))
       end do

       test_mixmean_vmag = DSQRT(DOT_PRODUCT(test_mixmean_vel(1:ndim)&
            &,test_mixmean_vel(1:ndim)))

!!$       PRINT*,'test_mean_mag = ', test_mixmean_vmag
!!$       PRINT*,'test_mean_mag = ', phase_array(1)%volfracg,&
!!$            & phase_array(2)%volfracg, phase_array(1)%volfracg&
!!$            &+phase_array(2)%volfracg, maxvolfrac
       if(I_AM_NODE_ZERO)then
          Write(*,'(2(A25,2x,g17.8,/))')&
               'Re_m Desired = ', Re, &
               'Re_m Output = ', test_mixmean_vmag*char_length/vis 
       end if
       
       test_mix_var = zero
       test_spec_var = zero
       pstart = 1
       do iphs = 1, nphases
          pend = pstart + phase_array(iphs)%npart - 1
          do m = pstart, pend
             do idim = 1, ndim
                test_spec_var(iphs) = test_spec_var(iphs) + (velbdy(m,idim)-test_mean_spec_vel(iphs,idim))**2.d0
             end do
          end do
          test_spec_var(iphs) = test_spec_var(iphs)/real(phase_array(iphs)%npart,prcn)
          pstart = pend + 1
       end do
       do iphs = 1, nphases
          test_mix_var = test_mix_var + test_spec_var(iphs)*phase_array(iphs)%volfracg/maxvolfrac
       end do
       if(I_AM_NODE_ZERO)then
          Write(*,'(2(A25,2x,g17.8,/))')&
               'T Desired = ', gran_temp, &
               'T Output = ',  test_mix_var/three
          
          Write(*,'(2(A25,2x,g17.8,/))')&
               'Re_T Desired = ', ReT, &
               'Re_T Output = ', DSQRT(test_mix_var/three)*char_length/vis 
       end if

       test_mix_var = zero
       test_spec_var = zero
       pstart = 1
       do iphs = 1, nphases
          pend = pstart + phase_array(iphs)%npart- 1
          do m = pstart, pend
             do idim = 1, ndim
                test_spec_var(iphs) = test_spec_var(iphs) + (velbdy(m,idim)-test_mixmean_vel(idim))**2.d0
             end do
          end do
          test_spec_var(iphs) = test_spec_var(iphs)/real(phase_array(iphs)%npart,prcn)
          pstart = pend + 1
       end do
       do iphs = 1, nphases
          test_mix_var = test_mix_var + test_spec_var(iphs)*phase_array(iphs)%volfrac/mean_volfrac
       end do
       
!!$       Write(*,'(2(A25,2x,g17.8,/))')&
!!$            'Re_T Desired = ', ReT, &
!!$            'Re_T Output HRENYA DEF= ', DSQRT(test_mix_var/three)*dia_phys/vis 
    ELSE IF(TRIM(psd_type).eq."csd")then
       test_mixmean_vel = zero
       do idim = 1, ndim
          test_mixmean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
       end do
       
       test_mixmean_vmag = DSQRT(DOT_PRODUCT(test_mixmean_vel(1:ndim),test_mixmean_vel(1:ndim)))
       if(I_AM_NODE_ZERO)then
          Write(*,'(2(A25,2x,g17.8,/))')&
               'Re_m Desired = ', Re, &
               'Re_m Output = ', test_mixmean_vmag*char_length/vis 
       end if
       test_mix_var = zero
       do m = 1, nbody
          do idim = 1, ndim
             test_mix_var = test_mix_var + (velbdy(m,idim)- test_mixmean_vel(idim))**2.d0
          end do
       end do
       test_mix_var = test_mix_var/real(nbody,prcn)
       if(I_AM_NODE_ZERO)then       
          Write(*,'(2(A25,2x,g17.8,/))')&
               'Re_T Desired = ', ReT, &
               'Re_T Output = ', DSQRT(test_mix_var/three)*char_length/vis 
       end if
    END if

    if(I_AM_NODE_ZERO)  CALL screen_separator(80,'I')
  END SUBROUTINE INPUT_CHECK
         
  
  SUBROUTINE alloc_mem

    USE nlarrays
    USE  bcsetarrays
    USE  nlmainarrays
    USE global_data 
    
    IMPLICIT NONE 
    INTEGER :: rg, npt, i,j,k, iphs
!!$    revgroup = 2 
!!$    ALLOCATE(revp(revgroup))
!!$    ALLOCATE(pgrev  (revgroup))
!!$
!!$    revp(1:2) = nrpr
!!$    DO rg = 1,revgroup
!!$       npt = revp(rg)
!!$       ALLOCATE(pgrev  (rg)%p(npt))
!!$    ENDDO

    ALLOCATE(gstencil  (order,order,order,3), vstencil(order,order&
         &,order,3), prsten(order,order,order))

    ALLOCATE(vsten(order,order,order,3), nlsten(order,order,order,3)&
         &,onlsten(order,order,order,3), dfsten(order,order,order,3),&
         & ppgrsten(order,order,order,3))
     
    ALLOCATE( gradphisten(order, order, order, 3, 3))
    ALLOCATE(sourcesten(order,order,order,3))
    !Print*,'In alloc_mem_initflo ', size(gstencil,1), size(gstencil,2)
    !Read(*,*)
    !nlarrays 

    !global_data

if (igeometry==0) then 
    
#if !PARALLEL
    ALLOCATE(u(nx+1,my2,mz,ndim),p(nx,my2,mz))!,uin(my2,mz,ndim),uout(my2,mz,ndim))
#else
    ALLOCATE(u(0:nx+1,my2,mz,ndim), p(0:nx+1,my2,mz))
#endif
    ! Make sure to create the appropriate strides at the end of the subroutine below.
    ALLOCATE(nl(nx+1,my2,mz,ndim),onl(nx+1,my2,mz,ndim))
    ALLOCATE(ff(nx+1,my2,mz,ndim))

    ALLOCATE(wy(my2),wz(mz))

    IF(aliasflag.eq.1) THEN 
       ALLOCATE(shifty(my2),shiftz(mz),shiftyz(my2,mz))
    end IF
    
!!$    ALLOCATE(wtbar(mxf,my,mz,2),wtbar_tmp(mxf,my,mz))
    ALLOCATE(w2(my2,mz), g(nx+1))
    !nlarrays 
    ALLOCATE(uf1(my2,mz),uf2(my2,mz),uf3(my2,mz))
    
    ALLOCATE(uf11(my2,mz),uf22(my2,mz),uf33(my2,mz))
    
    ALLOCATE(uf12(my2,mz),uf13(my2,mz),uf23(my2,mz))

    ALLOCATE(uf21(my2,mz),uf31(my2,mz),uf32(my2,mz))
    
    ALLOCATE(ur1(my,mz),ur2(my,mz),ur3(my,mz))
    
    ALLOCATE(ur11(my,mz),ur22(my,mz),ur33(my,mz))
    
    ALLOCATE(ur12(my,mz),ur13(my,mz),ur23(my,mz))

    ALLOCATE(ur21(my,mz),ur31(my,mz),ur32(my,mz))

#if !PARALLEL
    ALLOCATE(fr(nx+1,my,mz, ndim), ppr(nx,my,mz,ndim))
    
    ALLOCATE(diffn(nx+1,my,mz,ndim))
    
    !ALLOCATE(omega(mxf,my,mz, ndim))
    
    !nlmain arrays 
    ALLOCATE(ubc(nx+1,my,mz,ndim), nlbc(nx+1,my,mz,ndim), onlbc(nx+1,my,mz,ndim), pbc(nx,my,mz))
#else
    ALLOCATE(fr(0:nx+1,my,mz, ndim), ppr(0:nx+1,my,mz,ndim))
    
    ALLOCATE(diffn(0:nx+1,my,mz,ndim))
    
    !ALLOCATE(omega(mxf,my,mz, ndim))
    
    !nlmain arrays 
    ALLOCATE(ubc(-1:nx+2,my,mz,ndim), nlbc(0:nx+1,my,mz,ndim), onlbc(0:nx+1,my,mz,ndim), pbc(0:nx+1,my,mz))
    !extra buffers needed to store the velocity for interpolating the velocity at external points
#endif

    do k = 1, mz
       do j = 1, my2
#if !PARALLEL          
          do i = 1,nx+1
#else
          do i = 0,nx+1
#endif
             u(i,j,k,:) = czero
#if !PARALLEL
             if(i.le.mxf-1)  p(i,j,k)=  czero
#else             
             p(i,j,k)=  czero
#endif
             if(i.gt.0)then
                nl(i,j,k,:)=  czero
                onl(i,j,k,:)=  czero
                ff(i,j,k,:)=  czero
             end if
             if(i.eq.1) then 
                
                !uin(j,k,:) = czero 
                !uout(j,k,:) = czero 
                uf1(j,k) = czero  
                uf2(j,k) = czero  
                uf3(j,k) = czero  
                uf11(j,k) = czero 
                uf12(j,k) = czero 
                uf13(j,k) = czero 
                uf22(j,k) = czero 
                uf23(j,k) = czero 
                uf33(j,k) = czero 
                
             end if
          end do
       end do
    end do

    do k = 1, mz
       do j = 1, my
#if !PARALLEL          
          do i = 1,nx+1
#else
          do i = 0,nx+1
#endif
             !omega(i,j,k,:) = zero
             fr(i,j,k,:)=  zero
#if !PARALLEL
             if(i.le.nx) then
#endif
                ppr(i,j,k,:)=  zero
                pbc(i,j,k)=  zero
#if !PARALLEL
             endif
#endif


             diffn(i,j,k,:) = zero 
             onlbc(i,j,k,:)=  zero
             nlbc(i,j,k,:)=  zero
             ubc(i,j,k,:)=  zero
             
             
             if(i.eq.1) then 
                
                ur1(j,k) = zero  
                ur2(j,k) = zero  
                ur3(j,k) = zero  
                ur11(j,k) = zero 
                ur12(j,k) = zero 
                ur13(j,k) = zero 
                ur22(j,k) = zero 
                ur23(j,k) = zero 
                ur33(j,k) = zero 
                
             end if
          end do
#if PARALLEL
          ubc(-1,j,k,:) = zero
          ubc(nx+2,j,k,:) = zero
#endif
       end do
    end do
! Make sure to create the vectors in an appropriate manner. Create the vectors with the correct strides.
    ! With incorrect strides,there wont be any memory violations but the results will be completely awry.
    ! Remember this is FORTRAN :(

    ubcp=>ubc 
    nlbcp=>nlbc
    onlbcp =>onlbc
    pbcp => pbc 
    

endif

    ALLOCATE(fluid_atijk(0:nx+1,my,mz))

    CREATE_2D_LSLICE(my*mz,1,nx+2,twodlslice)
    COMMIT(twodlslice)
    
    CREATE_2D_RSLICE(my*mz,1,nx+2,twodrslice)
    COMMIT(twodrslice)
    CREATE_2D_RSLICE(my*mz,1,nx+4,urslice)
    COMMIT(urslice)
    CREATE_2D_CSLICE(my2*mz,1,nx+1,twodcslice)
    COMMIT(twodcslice)          
    CREATE_2D_CSLICE(my2*mz,1,nx+2,ucslice)
    COMMIT(ucslice)          
    ALLOCATE(ferror_array(nerr_steps))


    do iphs = 1, nphases
       ALLOCATE(phase_array(iphs)%ferror_array(nerr_steps))
       if(imove.eq.1)ALLOCATE(phase_array(iphs)%grant_array(nerr_steps))
    end do

    
  END SUBROUTINE alloc_mem

  SUBROUTINE set_interpolation_scheme(choice)
    USE global_data, ONLY : scheme, interp_scheme, order,ob2l,ob2r&
         &,gstencil,intx_per, inty_per,intz_per, vstencil, vsten,&
         & prsten, ppgrsten, nlsten, onlsten, dfsten
    IMPLICIT NONE 
    INTEGER, INTENT(in) :: choice
    INTEGER :: order_orig
    order_orig = order
    IF(choice.EQ.1) THEN 
       interp_scheme = 'lpi'
       scheme = '4-order'
    ELSE IF(choice.EQ.2) THEN 
       interp_scheme = 'lpi'
       scheme = '2-order'
    ELSE IF(choice.EQ.3) THEN 
       interp_scheme = 'csi'
       scheme = '4-order'
    ENDIF
    SELECT CASE(scheme)
    CASE("2-order")
       order = 2
    CASE("3-order")
       order = 3
    CASE("4-order")
       order = 4
    CASE("5-order")
       order = 5
    CASE("6-order")
       order = 6
    END SELECT

    IF(ALLOCATED(gstencil).AND.order_orig.NE.order) THEN 
       DEALLOCATE(gstencil) 
       ALLOCATE(gstencil  (order,order,order,3))
    END IF

    IF(ALLOCATED(vstencil).AND.order_orig.NE.order) THEN 
       DEALLOCATE(vstencil) 
       ALLOCATE(vstencil  (order,order,order,3))
    END IF

    IF(ALLOCATED(prsten).AND.order_orig.NE.order) THEN 
       DEALLOCATE(prsten) 
       ALLOCATE(prsten  (order,order,order))
    END IF

    IF(ALLOCATED(vsten).AND.order_orig.NE.order) THEN 
       DEALLOCATE(vsten) 
       ALLOCATE(vsten  (order,order,order,3))
    END IF

    IF(ALLOCATED(nlsten).AND.order_orig.NE.order) THEN 
       DEALLOCATE(nlsten) 
       ALLOCATE(nlsten  (order,order,order,3))
    END IF

    IF(ALLOCATED(onlsten).AND.order_orig.NE.order) THEN 
       DEALLOCATE(onlsten) 
       ALLOCATE(onlsten  (order,order,order,3))
    END IF

    IF(ALLOCATED(dfsten).AND.order_orig.NE.order) THEN 
       DEALLOCATE(dfsten) 
       ALLOCATE(dfsten  (order,order,order,3))
    END IF

    IF(ALLOCATED(ppgrsten).AND.order_orig.NE.order) THEN 
       DEALLOCATE(ppgrsten) 
       ALLOCATE(ppgrsten  (order,order,order,3))
    END IF
    IF(ALLOCATED(gradphisten).AND.order_orig.NE.order) THEN 
       DEALLOCATE(gradphisten) 
       ALLOCATE(gradphisten  (order,order,order,3,3))
    END IF
    intx_per = .FALSE.
    IF(xperiodic)     intx_per = .TRUE.
    inty_per = .True.
    intz_per = .True.
    ob2l = (order+1)/2
    ob2r = order/2 
  END SUBROUTINE set_interpolation_scheme

  SUBROUTINE generate_psd_config
    USE functions
    IMPLICIT NONE
    REAL(prcn) :: davg,sigmad,drms,lbydrms,phi_calc,drms_calc,min_dia,max_dia, davg_calc,sigma_calc
    REAL(prcn), ALLOCATABLE, DIMENSION(:) :: FRANDN, DBDY
    REAL(prcn), ALLOCATABLE, DIMENSION(:) :: FTEMP,Wt
    REAL(prcn), ALLOCATABLE, DIMENSION(:) :: dia_inc,vfrac_inc, radii
    REAL(prcn), ALLOCATABLE, DIMENSION(:) :: moments,weights,abscissa

    REAL(prcn) :: fmin,fmax,hist(15), sum_volfrac
    REAL(prcn) :: vbox, x1, x2, sigd3, sigd2, conf
    INTEGER :: SIZE_FRANDN, IDIM, NCOUNT, M, unit1,iphs, partstart, partend, sum_part
    INTEGER :: nsim, imis_gof, j

	character*50 filename
	real(prcn) :: confint, int_dist_avg, int_dist_var, int_dist_sd, int_dist2

    
    davg = dia_phys
    sigmad = sigmabydavg*davg

    if(TRIM(psd_type).eq.'csd')nphases = UNDEFINED_I

    if(nphases.eq.1)sigmad = zero
    
    drms = davg*(DSQRT(1 + SQR(sigmad*davg)))
    
    
    !lbydrms = DOML(2)/drms
    !PRINT*,'DRMS = ', DRMS
    
    if(I_AM_NODE_ZERO)then
       if(TRIM(psd_type).eq.'csd')then
          char_length = davg
          DOML(2) = lybyd*davg
          doml(1) = doml(2)
          doml(3) = doml(2)
          
          PRINT*,'LYBYD = ', LYBYD, 'DAVG = ', davg, 'DOML(2) = ', DOML(2)
          vbox = doml(1)*doml(2)*doml(3)

          nbody = NINT(phiavg*6.d0*(lybyd**3.d0)/pi)
          if(I_AM_NODE_ZERO)PRINT*,'NBODY IN CSD CASE  = ', NBODY
          SIZE_FRANDN = 4*NBODY
          ALLOCATE(FRANDN(SIZE_FRANDN),DBDY(nbody))
          
599       CONTINUE
          CALL norm_dist(FRANDN(1:SIZE_FRANDN))
          
          ncount = 0 
          
          Do m = 1, SIZE_FRANDN
             IF(ABS(FRANDN(m)).lt.3.d0)then
                ncount = ncount + 1
                DBDY(NCOUNT) = FRANDN(M)*sigmad + davg
             END IF
             IF(NCOUNT.eq.NBODY)EXIT
          End Do
          IF(NCOUNT.lt.NBODY) then
             WRITE(*,*)'NCOUNT = ', 'ncount, NBODY = ', NBODY
             WRITE(*,*)'REDOING NORMAL DISTRIBUTION'
             GOTO 599
          END IF
          MIN_DIA = MINVAL(DBDY)
          MAX_DIA = MAXVAL(DBDY)
          DAVG_CALC = ZERO
          DRMS_CALC = ZERO
          PHI_CALC = ZERO
          Do m = 1, nbody
             davg_calc = davg_calc + dbdy(m)/REAL(nbody,prcn)
             drms_calc = drms_calc+ SQR(dbdy(m))
             phi_calc = phi_calc + pi*dbdy(m)*SQR(dbdy(m))/6.d0
          End Do
          drms_calc = DSQRT(drms_calc/real(nbody,prcn))
          phi_calc = phi_calc/(DOML(1)*DOML(2)*DOML(3))
          sigma_calc = DSQRT(SQR(drms_calc)- SQR(davg_calc))
          if(ABS(sigma_calc-sigmabydavg*davg).gt.0.01*sigmabydavg*davg)then
             Write(*,'(A80)')'REQUIRED SIGMA/<D> NOT REACHED. DOING NORM DIST AGAIN'
             goto 599
          end if
          Write(*,'(3(2x,A25,g17.8,/))')&
               'MIN DIA = ', min_dia, &
               'MAX DIA =', max_dia,&
               'RATIO = ', max_dia/min_dia
          Write(*,'(2(2x,A25,g17.8,/))')&
               'DIA AVG CALC = ', davg_calc, &
               'DIA AVG INPUT =', davg

          Write(*,'(2(2x,A25,g17.8,/))')&
               'VOL FRAC CALC = ', phi_calc, &
               'VOL FRAC INPUT =', phiavg
          Write(*,'(2(2x,A25,g17.8,/)))')&
               'DRMS CALC = ', drms_calc, &
               'DRMS INPUT =', drms
          Write(*,'(2(2x,A25,g17.8,/)))')&
               'SIGMA/<D> CALC = ', DSQRT(SQR(drms_calc/davg_calc)-one), &
               'SIGMA/<D> INPUT =', sigmabydavg
!!$             unit1=getnewunit(minunitno,maxunitno)
!!$             OPEN(unit=unit1,FILE='psd.dat',status='unknown')
!!$             ALLOCATE(ftemp(nbody),wt(nbody))
!!$             do m=1,nbody
!!$                wt(m) = 1/real(nbody,prcn)
!!$             end do
!!$             ftemp(1:nbody) = dbdy(1:nbody)
!!$             CALL histogram(ftemp,wt,nbody,15,fmin,fmax,hist)
!!$             CALL plothist(hist(1:15),fmin,fmax,15,unit1,1.d0,1.d0)
!!$             close(unit1,status='keep')
          nphases = nbody
          CALL ALLOC_PHASE_RELATED
          do iphs = 1, nphases
             phase_array(iphs)%npart = 1
             phase_array(iphs)%dia = dbdy(iphs)
             phase_array(iphs)%volfrac = (phase_array(iphs)%npart)*pi*(phase_array(iphs)%dia**3.d0)/(6.d0*vbox)
          end do

       else if(psd_type.eq.'discrete')then
          Write(*,'(A,2x,I6)')'GENERATING A REPRESENTATION OF GAUSSIAN CSD WITH PHASES = ', nphases
          ALLOCATE(moments(2*nphases),weights(nphases),abscissa(nphases))

          char_length = davg
          DOML(2) = lybyd*davg
          doml(1) = doml(2)
          doml(3) = doml(2)
          
          PRINT*,'LYBYD = ', LYBYD, 'DAVG = ', davg, 'DOML(2) = ', DOML(2)
          vbox = doml(1)*doml(2)*doml(3)
          
          nbody = INT(phiavg*6.d0*(lybyd**3.d0)/pi)
          char_length = davg
          if(I_AM_NODE_ZERO)PRINT*,'NBODY IN DISCRETE CASE  = ', NBODY
          moments(:) = zero
          moments(1) = 1
          if(nphases.gt.1)then
             moments(3) = 1
             do iphs = 3,nphases
                moments(2*iphs-1) = moments(2*iphs-3)*(2*iphs-3)
             end do
          end if
          PRINT*,'MOMENTS =', moments(:)
          CALL hermite_polynomial(nphases,moments,weights,abscissa) 
          CALL ALLOC_PHASE_RELATED
          ALLOCATE(DBDY(nbody))
          sum_part = 0
          do iphs = 2, nphases
             phase_array(iphs)%npart = INT(weights(iphs)*nbody)
             sum_part = sum_part+phase_array(iphs)%npart
          end do
!!$             npart(nphases) = nbody - SUM(npart(1:nphases-1))
          if(nphases.gt.1)then
             phase_array(1)%npart = nbody - sum_part
          else
             phase_array(1)%npart = nbody
          end if
          sum_volfrac = zero
          do iphs = 1, nphases
             phase_array(iphs)%dia = davg + abscissa(iphs)*sigmad
             phase_array(iphs)%volfrac = phase_array(iphs)%npart*pi*(phase_array(iphs)%dia**3.d0)/(6.d0*vbox)
             sum_volfrac = sum_volfrac + phase_array(iphs)%volfrac
          end do
          
          Write(*,'(2(2x,A25,g17.8,/))')&
               'VOL FRAC CALC = ', sum_volfrac, &
               'VOL FRAC INPUT =', phiavg
          Write(*,'(2(2x,A25,I6,/))') 'NBODY CALC = ', sum_part+phase_array(1)%npart, 'NBODY INPUT =', nbody
          
       else if(TRIM(psd_type).eq.'bidisp') then
          
          nphases = 2
          ibidisperse = .TRUE.
          CALL ALLOC_PHASE_RELATED
          dia_ratio = yalpha(2)/yalpha(1)
          volfrac_rat = (dia_ratio-yalpha(2))/(yalpha(2)-1)

          Write(*,'(A)')'GENERATING A BIDISPERSE SUSPENSION WITH : '
          Write(*,'(A, 2x, g17.8)')'VOLUME FRACTION RATIO = ', volfrac_rat
          Write(*,'(A, 2x, g17.8)')'DIAMETER RATIO = ', dia_ratio
          
          
          phase_array(1)%volfrac = phiavg/(one+volfrac_rat)
          phase_array(2)%volfrac = volfrac_rat*(phase_array(1)%volfrac)
          
          x1 = phase_array(1)%volfrac/phiavg
          x2 = phase_array(2)%volfrac/phiavg
          
          phase_array(1)%dia = dia_phys*(x1 + x2/dia_ratio)
          phase_array(2)%dia = dia_ratio*(phase_array(1)%dia)

          !DOML(2) = lybyd*phase_array(2)%dia
          DOML(2) = lybyd*dia_phys
          doml(1) = doml(2)
          doml(3) = doml(2)
          PRINT*,'LYBYD = ', LYBYD, 'DIA2 = ',phase_array(2)%dia , 'DOML(2) = ', DOML(2)
          vbox = doml(1)*doml(2)*doml(3)


          phase_array(1)%npart = INT(6.d0*vbox*phase_array(1)%volfrac/(pi*phase_array(1)%dia**3.d0))
          phase_array(2)%npart = INT(6.d0*vbox*phase_array(2)%volfrac/(pi*phase_array(2)%dia**3.d0))

          nbody = phase_array(1)%npart + phase_array(2)%npart

!^^^^^^SINGLE PARTiCL MOVING ^^^^^^^^
!          phase_array(1)%dia = dia_phys/3
!         phase_array(2)%dia = dia_phys
!          phase_array(1)%volfrac = pi/6 * (phase_array(1)%dia/doml(2))**3
!          phase_array(2)%volfrac = phiavg
!          phase_array(1)%npart = 1
!          phase_array(2)%npart = INT(6.d0*vbox*phase_array(2)%volfrac/(pi*phase_array(2)%dia**3.d0))
!          nbody = phase_array(1)%npart + phase_array(2)%npart
!		  phiavg = phase_array(1)%volfrac + phase_array(2)%volfrac
!------------------------------------
          
          WRITE(*, '(A25,2(2x,g17.8))')'Diameters (ORIGINAL)=', phase_array(1)%dia, phase_array(2)%dia
          WRITE(*, '(A25,2(2x,g17.8))')'Volume fractions (ORIGINAL)=', phase_array(1)%volfrac, phase_array(2)%volfrac
          WRITE(*, '(A25,2(2x,I6))') 'No. of Particles =', phase_array(1)%npart, phase_array(2)%npart
          Write(*,'(2(2x,A25,I6,/))') 'NBODY CALC = ', phase_array(1)%npart+phase_array(2)%npart
          
          ALLOCATE(DBDY(nbody))
       end if

       
       ALLOCATE(dia_inc(nphases),vfrac_inc(nphases), xc_gener(nbody,ndim), rad_gener(nbody))
       
       partstart = 1
       sum_volfrac = zero
       min_dia = LARGE_NUMBER
       Do iphs = 1, nphases
          if(phase_array(iphs)%dia.lt.min_dia)then
             min_dia = phase_array(iphs)%dia
          end if
       End Do
       if(my.eq.undefined_I)then
          if(TRIM(psd_type).eq.'bidisp')then
             my = (DOML(2)/DIA_PHYS)*DBYDX
          else
             my = (DOML(2)/MIN_DIA)*DBYDX
          end if
          PRINT*,'SETTING MY BASED ON THE DBYDX FOR THE SMALLEST SPHERE, WHICH IS:', MY
          if(MOD(my,2).ne.0) then
             my = my+1
             PRINT*,'MY IS NOT EVEN:', MY, 'ADDING 1 : ', MY
          endif
          mx = my+1
          mz = my
          
          
          mxf = mx 
          my2 =  my/2+1
          mx1 = mx - 1 
       end if
       dx = doml(2)/real(my,prcn)
       dy = dx
       dz = dx 

       Do iphs = 1,nphases
          percent_buf(iphs) = (min_part_sep*dx/phase_array(iphs)%dia)*100.d0
          percent_buf(iphs) = percent_buf(iphs) + one
          
          dia_inc(iphs) = phase_array(iphs)%dia*(one + percent_buf(iphs)/100.d0)
          partend = partstart + phase_array(iphs)%npart-1
          do m = partstart,partend
             dbdy(m) = dia_inc(iphs)
          end do
          partstart = partend + 1
          sum_volfrac = sum_volfrac + phase_array(iphs)%volfrac
          vfrac_inc(iphs) = phase_array(iphs)%volfrac*(one+percent_buf(iphs)/100.d0)**3.d0
          !PRINT*,'percent_buf IN = ', percent_buf(iphs)
       End Do
       
       Write(*,'(3(2x,A25,g17.8,/))')&
            'VOL ORIG = ', sum_volfrac, &
            'VOL FRAC INC =', SUM(vfrac_inc(:))
       
       IF(GOF_AVG) then 
          WRITE(*,*)'GOF AVG IS TRUE, THEREFORE GENERATING', mis_gof, 'initializations, averaging and then stopping' 
          WRITE(ounit,*)'GOF AVG IS TRUE, THEREFORE GENERATING', mis_gof, 'initializations, averaging and then stopping' 
          nsim = mis_hs
       else 
          nsim = 1
       ENDIF

		if (.not.allocated(gofr_avg)) then
			allocate(gofr_avg(nbins), gofr_mis(nsim, nbins))
			allocate(rad_bin(nbins))
			allocate(int_dist(nsim))

			gofr_avg = zero
			gofr_mis = zero
			rad_bin  = zero
			int_dist = zero
		endif

       do imis_gof = 1, nsim
          WRITE(*,*) 'GENERATING RANDOM CONFIG FOR IMIS =', imis_gof		 
          CALL gener_config_random(dbdy(1:nbody),dia_inc(1:nphases),vfrac_inc(1:nphases))
          gofr_mis(imis_gof,:) = gofr(:)

!		call interstitial_dist(xc_gener, rad_gener, int_dist(imis_gof), nbody)
       end do

       gofr_avg = zero 
       do j=1,nbins
          gofr_avg(j) = gofr_avg(j) + one/dble(nsim)*sum(gofr_mis(1:nsim,j))
       end do
     
       do j=1,nbins
          IF(nsim.ge.2) THEN 
             conf=1.96/sqrt(dble(nsim))*sqrt(1./dble(nsim-1)* sum((gofr_mis(1:nsim,j)- gofr_avg(j))**2))
          end IF

          write(gofavgunit,'(10(E20.10,1x))') rad_bin(j), rad_bin(j)*lybyd, gofr_avg(j), conf
       end do
       close(gofavgunit, status = "keep") 

		int_dist_avg = sum(int_dist)/nsim

		if (nsim>2) then
			call get_confin(nsim, confint)
			int_dist_var = sum( (int_dist(1:nsim)- int_dist_avg) **2) / nsim
			int_dist_sd  = sqrt(int_dist_var)

			conf = int_dist_sd / sqrt(float(nsim)) * confint
		else
			conf = zero
		endif

!		call int_dist_gofr(nbins, rad_bin(1:nbins), gofr_avg(1:nbins), int_dist2)

		filename = "NMIS_interparticle.dat"
		open (unit=1,file=trim(filename),status="replace",action="write")
		write (1,"(4D15.7)") maxvolfrac, int_dist_avg, conf, int_dist2
		close (1)

       IF(GOF_AVG) then 
          if(I_AM_NODE_ZERO)WRITE(*,*) 'GOF AVG IS TRUE, SO STOPPING THE SIMULATION'
          PARALLEL_FINISH()
          STOP 
       end IF
       
       if(TRIM(psd_type).eq.'csd')then
          nphases = nbody
          CALL DEALLOC_PHASE_RELATED
          CALL ALLOC_PHASE_RELATED
          do iphs = 1, nphases
             phase_array(iphs)%npart = 1
          end do
       end if
    END if
    
    BROADCAST_INT(nphases,1,NODE_ZERO,decomp_group)
    
    if(.not.I_AM_NODE_ZERO)CALL ALLOC_PHASE_RELATED
    
    sum_part = 0
    do iphs = 1, nphases
       BROADCAST_INT(phase_array(iphs)%npart,1,NODE_ZERO,decomp_group)
       sum_part = sum_part + phase_array(iphs)%npart
    end do

    nbody = sum_part

    CALL ALLOC_BND_RELATED

    if(I_AM_NODE_ZERO)then
       do m = 1, NBODY
          XC(m,1:ndim) = xc_gener(m,1:ndim)
          RADBDY(m) = RAD_GENER(m)
       end do
       DEALLOCATE(XC_GENER, RAD_GENER)
    end if

    BROADCAST_INT(my,1,NODE_ZERO,decomp_group)
    BROADCAST_DOUBLE(DOML(1),ndim,NODE_ZERO,decomp_group)
    if(.not.I_AM_NODE_ZERO)then
       dx = doml(2)/real(my,prcn)
       dy = dx
       dz = dx 
       mx = my+1
       mz = my
       mxf = mx 
       my2 =  my/2+1
       mx1 = mx - 1 
    end if

    if(I_AM_NODE_ZERO)then
       WRITE(*,*) 'BOX SIZE   = ', mx, my, mz
       WRITE(*,*) 'mxf, mx1, my2   = ', mxf, mx1, my2
       if(TRIM(psd_type).eq.'bidisp')then
          sigd3 = zero
          sigd2 = zero
          do m = 1, nbody
             sigd3 = sigd3 + (two*radbdy(m)*dx)**3.d0
             sigd2 = sigd2 + (two*radbdy(m)*dx)**2.d0
          end do
          char_length = sigd3/sigd2
          Write(*,'(A30, 2x,g17.8)')'CHAR LENGTH FOR BIDISP CASE = ', char_length

          if(I_AM_NODE_ZERO)then
             OPEN(unit = 1001, FILE = TRIM(RUN_NAME)//'_yis.dat', status='unknown')
             Write(1001,*)phase_array(1)%dia/char_length, phase_array(2)%dia/char_length
             CLOSE(1001, status='keep')
          end if
       end if
    end if 
    BROADCAST_DOUBLE(char_length, 1, NODE_ZERO, decomp_group)
  END SUBROUTINE generate_psd_config



  SUBROUTINE gener_config_random(DBDY,DIA_IN, PHI_IN)
    USE dem_mod
    USE collision_mod
    USE dependent_functions
    implicit none 
    
    REAL(prcn), INTENT(IN),DIMENSION(:) :: DIA_IN, PHI_IN,DBDY
    REAL(prcn) :: vbox, xtemp, ly_lat, dmax, davg,&
         & DAVG_CALC, SIGMAD_CALC, TSTOP1, TSTOP2, YLEN,&
         & DOML_TMP(3), A, B, ymin, ymax , DIST&
         &, davgmax, sigmabydavgmax,  DAVG_BAR, RHOP_BAR, TMP&
         &,dia_sm , DIA_AIM, RHO_AIM, DRATIO, mfp, tmfp
    REAL(prcn) :: total_phi, min_dia, max_dia, drms_calc, ndens
    REAL(prcn), DIMENSION(:), ALLOCATABLE :: frandn, radbdy_tmp
    Integer ::  nratio, MCOUNT, MINDEX, NSPEC_AIM,&
         & i,j,k, NP, iphs, partstart, partend
    
    INTEGER SIZE_FRANDN, IDIM, NSPEC_INIT(nphases), trialcount, npart(nphases)
    CHARACTER*80 ::  collision_type_orig    
  
    Integer, Dimension(1):: mloc

	logical, allocatable :: contact(:,:)
	real(prcn) :: max_overlap


    CALL screen_separator(80,'R')
    
    vbox = doml(1)*doml(2)*doml(3)       
    collision_type_orig = collision_type
    do iphs = 1, nphases
       npart(iphs) = phase_array(iphs)%npart
    end do

    nspec_init(1:nphases) = npart(1:nphases)
    
    
    SIZE_FRANDN = 4*NBODY
    
    
    ALLOCATE(radbdy_TMP(nbody), frandn(SIZE_FRANDN))
    
    DMAX = MAXVAL(DIA_IN(1:nphases))

    trialcount = 1
!!$    
!!$    IF(NSPEC1.LT.NSPEC2) THEN
!!$       dbdy(:) = dia2_in
!!$       NSPEC_AIM = NSPEC1
!!$       DIA_AIM = dia1_in
!!$    ELSE
!!$       dbdy(:) = dia1_in
!!$       NSPEC_AIM = NSPEC2
!!$       DIA_AIM = dia2_in
!!$    end IF
!!$    
!!$    IF(.not.ibidisperse) goto 2000
!!$    
!!$1000 CONTINUE
!!$    !       IF(DRATIO.GT.ONE) THEN 
!!$    MCOUNT = 0
!!$    CALL uni_dist(FRANDN(1:SIZE(FRANDN)))
!!$    FRANDN(:) = FRANDN(:)*NBODY
!!$
!!$    DO I = 1, SIZE(FRANDN)
!!$       
!!$       MINDEX = NINT(FRANDN(I))!NBODY
!!$       IF(MINDEX.EQ.0) CYCLE
!!$!       PRINT*,'MINDEX = ', MINDEX, DIA_AIM, NSPEC_AIM
!!$       IF(DRATIO.GT.ONE) THEN 
!!$          IF(dbdy(MINDEX).NE.DIA_AIM) THEN 
!!$             DBDY(MINDEX) = DIA_AIM
!!$             MCOUNT = MCOUNT + 1
!!$             IF(MCOUNT.EQ.NSPEC_AIM) GOTO 2000
!!$          end IF
!!$       ELSE
!!$          MCOUNT = MCOUNT + 1
!!$          IF(MCOUNT.EQ.NSPEC_AIM) GOTO 2000
!!$       end IF
!!$    end DO
!!$    
!!$    IF(MCOUNT.LT.NSPEC_AIM) THEN 
!!$       trialcount = trialcount+1
!!$       PRINT*,'DOING UNIDIST AGAIN'
!!$       PRINT*,'NSPEC_AIM = ', NSPEC_AIM, ' AND MCOUNT = ', MCOUNT
!!$       if(trialcount.gt.50)then 
!!$          WRITE(*,'(A,/,A,/,A)')'EXCEEDED 50 TRIALS OF TRYING TO GENERATE MINDEX',' GIVING UP AND STOPPING HERE', ' TRY A DIFFERENT INITIAL SEED.'
!!$          STOP
!!$       else
!!$          GOTO 1000
!!$       end IF
!!$    end IF
!!$          
!!$2000 CONTINUE
!!$    !       end IF
!!$    DAVG = ZERO
!!$    DO I = 1, NBODY
!!$       DAVG  = DAVG + DBDY(I)/NBODY
!!$    end DO
!!$

!!$
!!$    DO I = 1, NBODY
!!$       TOTAL_PHI = TOTAL_PHI + PI*(DBDY(I)**THREE)/6.d0
!!$       DRMS_CALC  = DRMS_CALC + DBDY(I)**2.d0
!!$       DAVG_CALC  = DAVG_CALC + DBDY(I) 
!!$       
!!$    end DO
!!$
!!$    
!!$    TOTAL_PHI = TOTAL_PHI/(DOML(1)*DOML(2)*DOML(3))
!!$
!!$    PRINT*,'TOTAL VOL FRACTION', TOTAL_PHI
!!$    WRITE(*,*) 'MAX DIA RATIO = ', MAXVAL(DBDY)/MINVAL((DBDY))
    
    DOML_TMP(:) = DOML(:)
    MAX_DIA = MAXVAL(DBDY)
    MIN_DIA = MINVAL(DBDY)
    PRINT*,'DMAX = ', MAX_DIA, MIN_DIA
    CALL gener_lattice_mod(nbody,doml(1:3),xc_gener(1:nbody,1:ndim),dmax,&
         & dbdy(1:nbody))
    PRINT*,'OUT OF GENER_LATTICE_MOD'
    IF(MAXVAL(DOML(1)-XC_gener(1:nbody,1))-MAX_DIA*HALF.LT.ZERO) PRINT*,'PARTICLE COMING OUT OF +X'
    IF(MINVAL(XC_gener(1:nbody,1))-MAX_DIA*HALF.LT.ZERO) THEN 
       PRINT*,'PARTICLE COMING OUT OF -X', MINVAL(XC_gener(1:nbody,1)), MAX_DIA, MINVAL(XC_gener(:,1))- MAX_DIA
    end IF
    

    
    IF((DOML(2)-MAXVAL(XC_GENER(1:nbody,2)))-MAX_DIA*HALF.LT.ZERO) PRINT*,'PARTICLE COMING OUT OF +Y'
    IF(MINVAL(XC_GENER(1:nbody,2))-MAX_DIA*HALF.LT.ZERO) PRINT*,'PARTICLE COMING OUT OF -Y'
    IF(MAXVAL(DOML(3)-XC_GENER(1:nbody,3))-MAX_DIA*HALF.LT.ZERO) PRINT*,'PARTICLE COMING OUT OF +Z'
    IF(MINVAL(XC_GENER(1:nbody,3))-MAX_DIA*HALF.LT.ZERO) PRINT*,'PARTICLE COMING OUT OF -Z'
    
    
    XLENGTH = DOML(1)
    ZLENGTH = DOML(3)
    YLENGTH  =  MAXVAL(XC_GENER(1:nbody,2))+DMAX*HALF
    
    RAD_GENER(1:NbodY) = HALF*DBDY(1:NBODY)
    !goto 20000
    GENER_CONFIG_CASE = .TRUE.    
    DT = LARGE_NUMBER
    
    IF(YLENGTH.LT.DOML(2)) THEN 
       WRITE(*,'(A7,2x,g17.8,A7,2x,g17.8)')'YLEN  = ', YLENGTH, 'LT DOML(2) = ', DOML(2)
       
       SHRINK = .FALSE.
       
       ymin = MINVAL(XC_GENER(1:nbody,2)) - MINVAL(dbdy(1:nbody))*half
       ymax = MAXVAL(XC_GENER(1:nbody,2)) + MAXVAL(dbdy(1:nbody))*half
       
       A = doml(2)/(ymax-ymin)
       B = -A*ymin
       PRINT*,'A, B = ', A,b
       XC_GENER(1:nbody,2) = XC_GENER(1:nbody,2)*A + B
       !CALL des_time_march(.true.)    
       
    ELSE 
       WRITE(*,'(A7,2x,g17.8,A7,2x,g17.8)')'YLEN  = ', YLENGTH, 'GT DOML(2) = ', DOML(2)
       collision_type="softsphere"
       TSTOP = 20.d0*SQRT((two*YLENGTH)/980.d0)
       SHRINK = .TRUE.
       TEST_YMAXVAL = .TRUE.
       YMAXVAL = DOML(2)
       DES_EN_INPUT = 0.3
       DES_EN_WALL_INPUT = 1.0

       CALL  des_time_march(.true.)
       CALL  des_time_march(.false.)
       
       SHRINK = .FALSE.
       TEST_YMAXVAL = .FALSE.
       
       XC_GENER(1:NBODY,1:3) = DES_POS_NEW(1:NBODY, 1:3)
    end IF
    
    collision_type="softsphere"
    
    XLENGTH = DOML(1)
    ZLENGTH = DOML(3)
    YLENGTH = DOML(2)
    DES_EN_INPUT(:) = 1.0
    DES_EN_WALL_INPUT(:) = 1.0
    !PRINT*,'SI =', SIZE(DES_EN_INPUT), DES_EN_INPUT
!!$    ndens = (nspec1+nspec2)/(vbox)
    ndens = SUM(npart(1:nphases))/(vbox)
    !PRINT*,'NUMBER DENSITY = ', ndens
    mfp = (one/ndens)**(one/three)
    !mfp = dchar/(phi1+phi2)
    tmfp = mfp/dsqrt(pvel_var)
    tstop = 10.d0*tmfp
    
    WRITE(*,*)'mfp, tmfp, tstop = ', mfp, tmfp, tstop
    !    20000 continue
    
    CALL  des_time_march(.true.)
    
    
    CALL  des_time_march(.false.)
    
    DO IDIM = 1, 3
       DES_POS_NEW(1:NBODY,IDIM) = DES_POS_NEW(1:NBODY, IDIM)/DOML(IDIM)
    end DO
    
    radbdy_tmp(1:NBODY)=rad_gener(1:NBODY)/DOML(2)
    
    DO I = 1, NBODY
       mloc = MINLOC(radbdy_tmp)
       rad_gener(I) = radbdy_tmp(MLOC(1))
       
       xc_gener(I,:) = des_pos_new(MLOC(1),:)
       
       radbdy_tmp(mloc(1)) = LARGE_NUMBER
    ENDDO
    
    CALL scale_to_grid_units(nbody,npart(1:nphases),nphases,my,mxf,xperiodic&
         &,percent_buf(1:nphases),xc_gener(1:nbody,1:3),rad_gener(1:nbody),&
         & min_part_sep, toscale=.TRUE.) 
    
    nbody = SUM(npart(1:nphases))
    
    do iphs = 1, nphases
       phase_array(iphs)%npart = npart(iphs)
       if(TRIM(psd_type).ne.'csd')then
          write(*,'(A80,g17.8,2x,i4)')'FINAL VOL FRAC AND NSPECS AFTER LATTICE, SHRINKAGE, AND RESCALING = ', npart(iphs)*fourthirdpi*(phase_array(iphs)%dia/(two*dx*real(my,prcn)))**3, npart(iphs)
!!$       write(*,'(A80,g17.8,2x,i4)')'RADBDY = ', radbdy(1)
          WRITE(*,'(A80,2(2x,i2))') 'NUMBER OF PARTICLES LOST FOR PAHSE = ', iphs, nspec_init(iphs)-npart(iphs)
          !read(*,*)
!!$    IF(nspec_init(iphs)-npart(iphs)1.gt.0.OR.nspec2_init-nspec2.gt.0) then 
!!$	OPEN(1001, file=TRIM(RUN_NAME)//"_DELETION_TRUE.dat", form="formatted")
!!$    WRITE(*,'(A80,2(2x,i2))') 'NUMBER OF PARTICLES LOST IN 1st AND 2nd SPECIES = ', nspec1_init-nspec1, nspec2_init-nspec2
!!$    CLOSE(1001, status="keep")
!!$    END IF
       end if
    END do
    
!    nbins = 200
    
    rescaling = .true.
    
    
    
    OPEN(1001, file=TRIM(RUN_NAME)//"_xc_post_grid_scal.dat", form="formatted")
    
    write(1001,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" ',' "MARK" '
    partstart = 1
    DO iphs = 1, nphases
       partend = partstart + npart(iphs) -1 
       do NP = PARTSTART,PARTEND
          WRITE(1001,'(10(2x,g15.8))')( XC_GENER(NP, i), i = 1, ndim), rad_gener(NP),iphs
       ENDDO
       partstart = partend + 1
    END DO

    CLOSE(1001,status="keep")
    
    if(xperiodic)then
		if (.not.allocated(contact)) allocate(contact(nbody,nbody))
       CALL calculate_gofr_homog(nbody,xc_gener(1:nbody,1:3), contact, my, mxf, nbins, rescaling, gofr(1:nbins), rad_bin(1:nbins), max_overlap)
    else
       CALL calc_gofr(nbody, xc_gener(1:nbody,1:3), rad_gener(1:nbody), my, mxf,xperiodic,nbins, 3, rescaling, gofr(1:nbins),rho_est(1:nbins), rad_bin(1:nbins)) 
    end if
    
!!$    do j=1,nbins
!!$       write(gofunit,'(10(E20.10,1x))')rad_bin(j),rad_bin(j)/(one/(two*lybyd))&
!!$            &,rho_est(j),gofr(j)
!!$    end do
!!$    close(gofunit, status = "keep")
    
    GENER_CONFIG_CASE = .FALSE.
    
    
    collision_type = collision_type_orig
    CALL screen_separator(80,'R')
    
  end SUBROUTINE gener_config_random



  SUBROUTINE jacinv(nx,ny,nz,jj,cp,sp,st,ct,r1,r2)

    IMPLICIT NONE

    REAL(prcn) ::  nx,ny,nz
    INTEGER ::  i,j
    REAL(prcn) ::  jj(3,3)
    REAL(prcn) :: st,ct,sp,cp

    REAL(prcn) :: r1,r2


    !	r1=dsqrt(ny*ny+nz*nz)
    !	r2=dsqrt(nx*nx+ny*ny+nz*nz)

    !	ct=ny/r1
    !	st=nz/r1

    IF(r1.NE.(0.d0)) THEN

       !       ct=ny/r1
       !       st=nz/r1

       !       cp=nx/r2
       !       sp=r1/r2
       jj(1,1)=cp
       jj(1,2)=sp*ct
       jj(1,3)=sp*st

       jj(2,1)=-sp
       jj(2,2)=cp*ct
       jj(2,3)=cp*st

       jj(3,1)=0.d0
       jj(3,2)=-st/sp
       jj(3,3)=ct/sp

    ELSE

       !special case for phi=zero
       !otherwise matrix is singular

       DO i=1,3
          DO j=1,3
             jj(i,j)=0.d0
          ENDDO
       ENDDO

       jj(1,1)=nx/ABS(nx)
       !write(*,*) nx,jj(1,1)
       jj(2,2)=1.d0
       jj(3,3)=1.d0

    ENDIF

    RETURN

  END SUBROUTINE jacinv




  
  SUBROUTINE ALLOC_BND_RELATED
    IMPLICIT NONE 
    
    INTEGER :: NBODY_TMP

    ALLOCATE(force(nbody,ndim),pres(nbody,ndim),visc(nbody,ndim),torq(nbody,ndim), force_loc(nbody,ndim), force_chem(nbody,ndim), contact_force(nbody,ndim))
    ALLOCATE(presloc(nbody,ndim), viscloc(nbody,ndim),torqloc(nbody,ndim))

    ALLOCATE(ap(nbody,ndim,2),up(nbody,ndim,2),angv(nbody,ndim,2),anga(nbody,ndim,2))

    ALLOCATE(mp(nbody),mpart(nbody),mompart(nbody), xc(nbody, ndim), tmpfor(nbody,ndim), velbdy(nbody,ndim))

	allocate(color(nbody))

	angv = zero
	anga = zero

    ALLOCATE(radbdy(nbody),radibdy(nbody),radobdy(nbody),rado2bdy(nbody))
    ALLOCATE(part_array(nbody))
    IF(iscalon.eq.1)THEN
       ALLOCATE(phisurfall(nbody, nspmx))
       ALLOCATE(phirmean(nspmx), fphirmean(nspmx),sum_flux_nm1(nbody)&
            &, sum_flux_nm2(nbody), sum_flux(nbody), flux_body(nbody&
            &,nspmx), flux_body2(nbody,nspmx))
       phisurfall = zero
       phisurfall_nm1 = zero
       phisurfall_nm2 = zero
       sum_flux = zero
       sum_flux_nm1 = zero
       sum_flux_nm2 = zero
    END IF
    
  end SUBROUTINE ALLOC_BND_RELATED

  SUBROUTINE DEALLOC_BND_RELATED
    IMPLICIT NONE 
    INTEGER :: m

    DEALLOCATE(force,pres,visc,torq,force_loc,force_chem,contact_force)
    DEALLOCATE(presloc, viscloc,torqloc)
    DEALLOCATE(ap,up,angv,anga)
    DEALLOCATE(mp,mpart,mompart, xc, tmpfor, velbdy)
        
    DEALLOCATE(radbdy,radibdy,radobdy,rado2bdy)
    do m = 1, nbody
       IF (ASSOCIATED(part_array(m)%if_rev)) DEALLOCATE(part_array(m)%if_rev)
    end do
    DEALLOCATE(part_array)

    IF(iscalon.eq.1)THEN
       DEALLOCATE(phisurfall)
       DEALLOCATE(phirmean, fphirmean,sum_flux_nm1, sum_flux_nm2, sum_flux, flux_body, flux_body2)
    END IF
  end SUBROUTINE DEALLOC_BND_RELATED


  
  SUBROUTINE ALLOC_PHASE_RELATED
    IMPLICIT NONE
    ALLOCATE(phase_array(nphases), percent_buf(nphases),norm_drag_spec(nphases),norm_drag_chem_spec(nphases))
  END SUBROUTINE ALLOC_PHASE_RELATED

  SUBROUTINE DEALLOC_PHASE_RELATED
    IMPLICIT NONE
    INTEGER :: iphs

    do iphs = 1, nphases
       IF (ASSOCIATED(phase_array(iphs)%bndpts)) DEALLOCATE(phase_array(iphs)%bndpts)
       IF (ASSOCIATED(phase_array(iphs)%ferror_array)) DEALLOCATE(phase_array(iphs)%ferror_array)
    end do

    DEALLOCATE(phase_array, percent_buf, norm_drag_spec, norm_drag_chem_spec)
  END SUBROUTINE DEALLOC_PHASE_RELATED


  SUBROUTINE DOMAIN_1D_DECOMP( n, numprocs, myid, s, e )
    Integer,Intent(in):: n, numprocs, myid
    Integer, Intent(out) :: s, e
    Integer :: nlocal, deficit
    
    nlocal  = n / numprocs
    s	      = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s	      = s + min(myid,deficit)
    
    if (myid .lt. deficit) then
       nlocal = nlocal + 1
    endif
    
    e = s + nlocal - 1
    
    if (e .gt. n .or. myid .eq. numprocs-1) e = n
    
  end SUBROUTINE DOMAIN_1D_DECOMP


END MODULE initialize_flo
