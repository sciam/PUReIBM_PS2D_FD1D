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


PROGRAM dtibm_main
  !-----------------------------------------------------------------------
  !	3D code. 
  !	Finite difference in the x-direction
  !	Fourier pseudospectral in the y-direction
  !	Fourier pseudospectral in the z-direction
  ! AUTHOR: RAHUL GARG 
  !	mx	number of interior (u,v,w)-gridpoints in the x-direction
  !	mx1	number of interior (p)-gridpoints in the x-direction!
  !	mxf	number of yz- planes where forcing exists
  !	my2 	number of fourier modes in the y-direction
  !	mz 	number of fourier modes in the z-direction
  !	my	number of (y-) gridpoints on the grid 
  !	mz	number of (z-) gridpoints on the grid 

  !-----------------------------------------------------------------------

#include "ibm.h"

	USE precision  
	USE restart_funcs
	USE constants  
	USE global_data
	USE initialize_flo 
	USE nlcalc
	USE steptheflow
	USE initialize  
	USE dem_mod
	USE collision_mod
	USE writeoutput 
	USE machine 
	USE general_funcs
	USE stepthescalar
	USE outputscalar

	use mypost_process

#if PARALLEL
	USE nlarrays, ONLY : ur11, ur22,uatminus1, uatnxp2
	USE nlmainarrays, ONLY : ubcp
#endif
	IMPLICIT NONE

	real(prcn) :: cpu0, cpu1, cpu2, cpu3, cpu4, cpu5, cpu6, cputime_used, cputime_hrs, cputime_min
	INTEGER :: run_name_len, orig_rank

#if PARALLEL
	communicator_done = .false.
#endif

	PARALLEL_START()

	GET_NPROCS(comm_group,nproc)
	GET_PROCESSOR_RANK(comm_group,myid)

	! MYID = 0 FOR SERIAL VERSION
	if (I_AM_NODE_ZERO)then
		CALL GETARG(1 , RUN_NAME)
		run_name_len = LEN_TRIM(RUN_NAME)
	end if

	BROADCAST_INT(run_name_len,1,NODE_ZERO,comm_group)
	BROADCAST_CHARARR(RUN_NAME,run_name_len,NODE_ZERO,comm_group)

	! BROADCAST_STRING(RUN_NAME, run_name_len, NODE_ZERO,comm_group)
  
	CALL calculate_constants
	CALL init_params !set the parameters

	! Virtual processor topology for parallel decomposition
	CREATE_CART_TOPOLOGY(comm_group,1,nproc,xperiodic,.true.,decomp_group)  
	orig_rank = myid

#if PARALLEL
	communicator_done = .true.
#endif

	!decomp_group = 1
	GET_PROCESSOR_RANK(decomp_group,myid)
	GET_SHIFT_PROCS(decomp_group,0,1,fromproc,toproc)

	if(I_AM_NODE_ZERO) then
		!--------
		! open error file
		!--------
		eunit = getnewunit(minunitno,maxunitno)
		IF (eunit.LT.0) THEN   ! have to handle this locally
			WRITE(*,*)'Cannot find unused unit for error file'
			WRITE(*,*)'Quitting...'
			PARALLEL_FINISH()
			STOP
		ENDIF
		errfile = TRIM(RUN_NAME)//'_'//TRIM(errfile)
		IF(IRESTART.EQ.0) THEN 
			OPEN(unit=eunit,file=errfile,form="formatted",status="replace")
		ELSE
			OPEN(unit=eunit,file=errfile,form="formatted",POSITION="append")
		endIF
       
		!--------
		! open file for standard out
		!--------
		ounit = getnewunit(minunitno,maxunitno)
		!ounit = -1
		IF (ounit.LT.0) CALL printerror("newunit","ounit")
		!Print*,'UNit for outfile =', ounit, eunit
		outputfile = TRIM(RUN_NAME)//'_'//TRIM(outputfile)
		IF(IRESTART.EQ.0) THEN 
			OPEN(unit=ounit,file=outputfile,form="formatted",status="replace")
		ELSE
			OPEN(unit=ounit,file=outputfile,form="formatted",POSITION="append")
		endIF
     
		CALL GET_RUN_ID 
		CALL screen_separator(80,'*')

		call separator(ounit,62,'-')
		IF(IRESTART.EQ.0) THEN
#if PARALLEL
			WRITE(*,'(A,i2,A)') 'STARTING A NEW PARALLEL RUN ON ',nproc,' NODES : RUN_NAME IS ', RUN_NAME
			WRITE(ounit,'(A,A)') 'STARTING A NEW PARALLEL RUN ON ',nproc,' NODES : RUN_NAME IS ', RUN_NAME
#else
			WRITE(*,'(A,A)') 'STARTING A NEW RUN: RUN_NAME IS ', RUN_NAME
			WRITE(ounit,'(A,A)') 'STARTING A NEW RUN: RUN_NAME IS', RUN_NAME
#endif
		ELSE
#if PARALLEL
			WRITE(*,'(A,A,A,i2,A)') 'RESTARTING AN OLD RUN HAVING RUN_NAME ', RUN_NAME,' ON ',nproc,'  NODES'
!        WRITE(ounit,'(A,i2,A)') 'RESTARTING AN OLD RUN HAVING RUN_NAME = ', RUN_NAME,' ON ', nproc, ' NODES'

#else
			WRITE(*,'(A,A)') 'RESTARTING A OLD RUN HAVING RUN_NAME = ', RUN_NAME
			WRITE(ounit,'(A,A)') 'RESTARTING A OLD RUN HAVING RUN_NAME = ', RUN_NAME
#endif
		end IF
     
		WRITE(*, '(A40,i2,A,i2,A,i4)') 'SIMULATION STARTED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
		WRITE(*, '(A40,i2,A,i2,A,i2)') 'SIMULATION STARTING TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND
		WRITE(*,'(A40,A,A,I1)') 'RUNNING ON MACHINE ', ID_NODE,' WITH RANK ', MYID

		call separator(ounit,62,'-')

		WRITE(ounit, '(A40,i2,A,i2,A,i4)') 'SIMULATION STARTED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
		WRITE(ounit, '(A40,i2,A,i2,A,i2)') 'SIMULATION STARTING TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND
		WRITE(ounit,'(A40,A)') 'RUNNING ON MACHINE ', ID_NODE

		CALL screen_separator(80,'*')

		call separator(ounit,62,'-')
		CALL CPU_TIME (CPU0) 
	end if
  
	BROADCAST_DOUBLE(cpu0,1,NODE_ZERO,decomp_group)
#if PARALLEL
	!WRITE(*,'(3(2x,A25,i4))')'RANK = ', myid, 'SORUCE = ', fromproc, 'DEST = ', toproc
#endif
	BARRIER(decomp_group)
	CALL main
	!WRITE(*, *) 'Simulation started on date:', date,' at time:', time
  
CONTAINS 
	SUBROUTINE main

!	use mypost_process

		IMPLICIT NONE  

		INTEGER ::  sflag,intstep,n,i,j,k,ifilstep,idim,lcstart, res_temp, ierr
		INTEGER idum,ifirstcheck,ic,ii, m,snapunit, count_resfiles, iphs
#if PARALLEL
		INTEGER :: nprocold,fnamelen
		LOGICAL :: rstok
#endif
		Integer, SAVE :: unit_t_his

		REAL(prcn),DIMENSION(:,:,:,:,:),ALLOCATABLE ::  velr,velr_mxf

		REAL*8 :: etime, elapsed(2), total

		REAL(prcn) :: dtold,s_temp, re_umax, abtparam, time_since_last_res, CPUTIME_USED_LAST_HRS
		CHARACTER*80 :: conv_file, res_file, filename1
		COMPLEX(prcn) :: usumloc, usum
		LOGICAL:: filexist, isopen, CALL_FLOW, killjob, stop_criterion
		REAL(prcn) ::  rhoofs(ndim),mean_vel(ndim), lag_str_func(ndim), delta_meshpos(ndim)
		!Allocate(velr(mxf,my,mz,ndim,2))
		REAL(prcn) :: iterstart_time, iterend_time, avg_iter_time
#if PARALLEL
		REAL(prcn) :: global_avg_iter_time, global_max_iter_time
#endif
		real(prcn) :: taup, tauf, froude
 
		ifilstep = 50
		if(I_AM_NODE_ZERO)then
			unit_t_his  = getnewunit(minunitno,maxunitno)

			conv_file = TRIM(RUN_NAME)//"_CONVERGED"
			res_file = TRIM(RUN_NAME)//"_RESTART"
			!First check if Convergence indicator file exists or not.
			!if it does, then delete it. 
			INQUIRE(FILE=conv_file,EXIST=filexist,OPENED=isopen)
			IF (filexist) THEN
				OPEN(unit = 1000, file=conv_file, status="old")
				close(1000, status="delete")
			end IF
		end if

		IF(irestart.eq.0) then 
			if(I_AM_NODE_ZERO)then
				OPEN(unit = 1000, file=res_file, status="unknown")
				WRITE(FILENAME1, '(I1)')0
				write(1000, *) FILENAME1
#if PARALLEL
				write(1000, *) nproc
#endif
				close(1000,status="keep")
			end if

		!BROADCAST_STRING(FILENAME1,fnamelen,NODE_ZERO,decomp_group)
		else if(irestart.eq.1)then
			if(I_AM_NODE_ZERO)then
				!check for this case if something screwed up during the
				! original run and as a result restart files were not properly
				!  written out. Can't happen with the logic in our code. But
				! once hpc4 guys deleted some files during run time  and all this mess
				! hapenned. MIS*_RESTART had 0 written in it and all the
				! restart files were also not written out. RG 11/14/08

				res_temp = 0 
				OPEN(unit = 1000, file=res_file, status="unknown", IOSTAT=IERR)
				READ(1000, *, IOSTAT=IERR) res_temp
#if PARALLEL
				read(1000, *, IOSTAT=IERR) nprocold
				if(nprocold.eq.nproc)rstok = .TRUE.
#endif
				WRITE(*,*)'FIELD IN RESTART INDICATOR FILE =', res_temp, irestart  
				CLOSE(1000, status='keep')
				!WRITE(*,*)'IERR =  ', IERR
				if(res_temp.eq.0) then 
					CALL screen_separator(80,'E')
					WRITE(*,*) 'WARNING'
					WRITE(*,'(A,/,A,/,A,/,A)') 'EVENTHOUGH IRESTART = 1, BUT FILE NAMED', res_file,'has&
							& 0 or nothing written in it. Something screwed up',' So starting &
							&the run from scratch'
					WRITE(ounit,*) 'WARNING'
					WRITE(ounit,'(A,/,A,/,A,/,A)') 'EVENTHOUGH IRESTART = 1, BUT FILE NAMED', res_file,'has&
							& 0 written in it. Something screwed up',' So starting &
							&the run from scratch'
							IRESTART = 0
							ISCAL_RESTART = 0
#if PARALLEL
					RSTOK = .TRUE.
#endif
					!Now delete this file and write 0 in it. This is to make sure
					! that in case nothing was writtin in this file, then at the
					! restart stage we do not get any read error. RG 11/16/08

					OPEN(unit = 1000, file=res_file, status="unknown", IOSTAT=IERR)
					CLOSE(1000, status='delete')
					OPEN(unit = 1000, file=res_file, status="unknown")
					WRITE(FILENAME1, '(I1)')0
					write(1000, *) FILENAME1
					close(1000,status="keep")

					CALL screen_separator(80,'E')
				end if
			end IF
#if PARALLEL
			BROADCAST_LOGICAL(rstok,1,NODE_ZERO,decomp_group)
			BROADCAST_INT(irestart, 1, NODE_ZERO, decomp_group)
			BROADCAST_INT(iscal_restart, 1, NODE_ZERO, decomp_group)
			if(.not.rstok)then
				if(I_AM_NODE_ZERO)then
					WRITE(*,'(A,i2)')'RESTARTING RUN WITH DIFFERENT NUMBER OF NODES: ', nproc
					WRITE(*,'(A,i2,A)')'PREVIOUSLY RUN WITH ', nprocold, ' NODES. RESTART AGAIN.'
				end if
				PARALLEL_FINISH()
				STOP
			end if
#endif
		end IF
		if(I_AM_NODE_ZERO)then
			IF(irestart.eq.0) then 
				OPEN(unit = unit_t_his,file=TRIM(RUN_NAME)//'_dt_history.dat',form='formatted', status='unknown')
				WRITE(unit_t_his,'(10(2x,A17))') 'ISTEP', 'DT', '#STEPS REQD', 'CFL'
			ELSE
				OPEN(unit = unit_t_his,file=TRIM(RUN_NAME)//'_dt_history.dat',form='formatted', POSITION='APPEND')
			end IF
		END if

		FROM_POST = .FALSE.
		POST_NO_FLOW_MEM_ALLOC = .FALSE.

		CALL initflo

		DO intstep=1,itrmax
		END DO

		final_restart = .FALSE.

		!-------------------------------------------------------------------
		!	timestep through simulation
		!------------------------------------------
		t=0.d0
		abtparam  = 0 
		idumstep = 0
		iglobstep = 0

		t = tstart  !in case of restart, tstart will be equal to the last
		! time step. Restart reading is done from initflo.f90

		if(I_AM_NODE_ZERO) PRINT*,' TOLERANCE FOR FERROR_HIST =',  TOL_FERROR
		CALL_FLOW = .TRUE.

		!IF(Re.EQ.ZERO) THEN 
		!   ferror_hist = 0.1d0*tol_ferror
		!end IF


		flow_converged = ferror_hist.lt.tol_ferror
		scal_converged = nu_error_hist.lt.tol_ferror

		flow_converged = .TRUE.
!!$    grant_converged = .FALSE.
!!$    if(imove.eq.1)grant_converged = .TRUE.

		do iphs = 1, nphases
			flow_converged = flow_converged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
		end do
!!$    if(imove.eq.1)then
!!$       do iphs = 1, nphases
!!$          grant_converged = grant_converged.and.(phase_array(iphs)%gran_error_hist.lt.tol_gran_error)
!!$       end do
!!$    end if


		CALL_FLOW = .not.flow_converged 
		if(iscalon.eq.0) scal_converged = .TRUE.
		if(I_AM_NODE_ZERO)then 
			WRITE(*,*)'FLOW_CONVERGED  = ', FLOW_CONVERGED
			WRITE(*,*)'SCAL_CONVERGED  = ', SCAL_CONVERGED
			!WRITE(*,*)'GRAN_CONVERGED  = ', GRANT_CONVERGED
			PRINT*,'MAX WALL TIME ALLOWED =', wtime_maxhrs
		end if
		!dO WHILE (t.LE.tend)
		cputime_hrs = zero
		time_since_last_res = zero
		count_resfiles = 0

		BARRIER(decomp_group)
	!!	stop_criterion = flow_converged.and.scal_converged
              
		if(imove.eq.1) stop_criterion = .FALSE. ! For a moving particle case run until wall time exceeds the limit.
		if(imove.eq.1) CALL_FLOW = .TRUE. ! For a moving particle case always call the flow.
               !    stop_criterion = .true.

		avg_iter_time = zero
		do While(.not.stop_criterion)
			if(I_AM_NODE_ZERO) CALL CPU_TIME (CPU2)

			iter_u = 0
			iter_p = 0
			iter_phi = 0
			iter_scal = 0

			if(I_AM_NODE_ZERO) then
				WRITE(*,'(A,2(2x,g12.5))')'FERROR_HIST AND TOL = ', ferror_hist, TOL_FERROR
				do iphs = 1, nphases
					WRITE(*,'(A,I8,2(2x,g12.5))')'FERROR_HIST AND TOL = ', iphs,phase_array(iphs)%ferror_hist, TOL_FERROR
					if(move_particles) WRITE(*,'(A,I8,2(2x,g12.5))')'GRAN_ERROR_HIST AND TOL = ', &
					&	iphs,phase_array(iphs)%gran_error_hist, TOL_GRAN_ERROR
				end do
			end if
			IF(ISCALON.EQ.1) then
				if(I_AM_NODE_ZERO) WRITE(*,'(A25,g12.5)')'NU ERROR HIST = ', nu_error_hist 
			END IF
       
			if(I_AM_NODE_ZERO) then
				CALL screen_separator(80,'-')
				CALL screen_separator(80,'-')
			end if

			idumstep = idumstep + 1
			iglobstep = iglobstep + 1
       
#if !PARALLEL
			CALL CPU_TIME (iterstart_time)
#else
			iterstart_time = MPI_WTIME()
#endif
			DO intstep=1,itrmax
				if(imove.eq.1.and.fixed_moving)then
					if((.not.move_particles))  move_particles = flow_converged
				else if(imove.eq.1.and.(.not.fixed_moving))then
!             if((.not.move_particles))  move_particles = .true.
					move_particles = .true.
				end if
          
				if(I_AM_NODE_ZERO) WRITE(*,'(A,2(2x,i5))') 'RKS AND GLOBAL STEP # = ', intstep, iglobstep
				if(I_AM_NODE_ZERO) WRITE(*,'(A)') 'COMPUTING THE NEW TIME STEP'
				CALL compute_new_timestep(intstep)

				if((mean_vel_to_particles).and.(imove.ne.1).and.(.not.movingcv))then
					delta_meshpos(1:ndim) = mesh_vel(1:ndim)*dt/dx
					Write(*,*)'Delta_meshpos : ', delta_meshpos(1:ndim)
					CALL interpolate_fields_to_new_mesh(delta_meshpos(1:ndim))
				end if
!!$
!!$          IF(MOD(iglobstep,saveitns).EQ.0) THEN 
!!$             
!!$             CALL output
!!$             STOP
!!$          END IF

				IF(CALL_FLOW) THEN
					CALL nonlinear(intstep)
					!In non linear now the dt is reset to newly calculated dt based on 
					!cfl criteria
					!Now t is incremented in nl to ensure correct value goes in bcset

					if(I_AM_NODE_ZERO) CALL CPU_TIME (CPU3)
					if(.not.only_dem) CALL velstep(sflag,intstep)
					if(I_AM_NODE_ZERO) then
						CALL CPU_TIME (CPU4)
						WRITE(unit_t_his, '(20(2x,g17.8))') REAL(iglobstep,prcn), DT, (tend-t)/dt + REAL(iglobstep, PRCN), uchar(1)*dt/dx
					endif

					if(debug_check)then
						usumloc = SUM(u(1:nx,1,1,1))
						GLOBAL_COMPLEX_SUM(usumloc,usum,1,decomp_group)
						if(I_AM_NODE_ZERO)then
							Write(*,'(A25,2(2x,g12.5))') 'UAVG = ', usum/mx1 !SUM(u(1:mx1,1,1,1))/mx1 
							Write(*,'(A25,2(2x,g12.5))') 'MAX DIVERGENCE = ', divmax
						endif
					endif
					!if(t.gt.tendused)
					flow_converged = ferror_hist.lt.tol_ferror
					flow_converged = .TRUE.
					do iphs = 1, nphases
						flow_converged = flow_converged.and.(phase_array(iphs)%ferror_hist.lt.tol_ferror)
					end do
		          
					CALL_FLOW = .not.flow_converged 
					if(imove.eq.1)CALL_FLOW = .TRUE.
				ELSE
					if(I_AM_NODE_ZERO) WRITE(*,*)'NOT CALLING FLOW ROUTINES BECAUSE EITHER CALL_FLOW = ', &
					&	CALL_FLOW, 'OR FLOW_CONVERGED = ', flow_converged
				endif

				IF(iscalon.EQ.1)THEN
					if(I_AM_NODE_ZERO) CALL CPU_TIME (CPU5)
#if PARALLEL
					if(flow_converged)then
						if(.not.vel_converted)then
							CALL ff2cr(u(2,1:my2,1:mz,1),ur11(1:my,1:mz))
							CALL ff2cr(u(nx-1,1:my2,1:mz,1),ur22(1:my,1:mz))
							RSENDRECV(ur11(1,1),my*mz,fromproc,1,uatnxp2(1,1),my*mz,toproc,1,decomp_group,status)
							RSENDRECV(ur22(1,1),my*mz,toproc,1,uatminus1(1,1),my*mz,fromproc,1,decomp_group,status)
							uatnxp2(:,:) = uatnxp2(:,:)+umean(1)
							uatminus1(:,:) = uatminus1(:,:)+umean(1)
							vel_converted = .TRUE.
						endif
					else
						uatminus1(1:my,1:mz) = ubcp(-1,1:my,1:mz,1)
						uatnxp2(1:my,1:mz) = ubcp(nx+2,1:my,1:mz,1)
					endif
#endif
					!nlphi now called in scalstep 
					CALL scalstep(sflag,intstep)

					!if(t.gt.tendused)
					scal_converged = nu_error_hist.lt.tol_ferror
					!scal_converged  = .true.
					if(I_AM_NODE_ZERO) CALL CPU_TIME (CPU6)
				ENDIF

				stop_criterion = flow_converged.and.scal_converged
				if(imove.eq.1) stop_criterion = .FALSE.
				IF (ISCALON.EQ.1) then
					usumloc = SUM(phif(1:nx,1,1,1))
					GLOBAL_COMPLEX_SUM(usumloc,usum,1,decomp_group)
					if (I_AM_NODE_ZERO) Write(*,'(A25,2(2x,g12.5))') 'PHI_AVG = ', usum/mx1
				end IF

!          if(move_particles)then
!             !CALL part_snapshot
!          end if
          
				if(I_AM_NODE_ZERO)CALL screen_separator(80,'*')
!				FIRST_PASS = .FALSE.   
				!cf = cforig
			ENDDO

			if(I_AM_NODE_ZERO)then
				CALL CPU_TIME (CPU1) 
				!BROADCAST_DOUBLE(cpu1,1,NODE_ZERO,decomp_group)

				CPUTIME_USED_LAST_HRS = CPUTIME_HRS

				CPUTIME_USED  = CPU1 - CPU0

				CPUTIME_HRS = CPUTIME_USED/(3600.d0)
				CPUTIME_MIN = CPUTIME_USED/(60.d0)
				time_since_last_res = time_since_last_res + (-cputime_used_last_hrs+cputime_hrs)

				PRINT*,'CPU TIMEHRS = ', CPUTIME_HRS

				hydro_time_current = (cpu1-cpu2)-dem_time_current
				hydro_time = hydro_time + hydro_time_current

				vp_time_current = (cpu4-cpu3)
				vp_time = vp_time + vp_time_current

				scalar_time_current = (cpu6-cpu5)
				scalar_time = scalar_time + scalar_time_current

				WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED IN HYDRO (S), CURRENT, TOTAL, AVG : ', hydro_time_current , hydro_time, hydro_time/iglobstep
				WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED IN NL (S)   , CURRENT, TOTAL, AVG : ', nl_time_current , nl_time, nl_time/iglobstep
				WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED IN BC (S)   , CURRENT, TOTAL, AVG : ', bc_time_current , bc_time, bc_time/iglobstep
				WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED IN VP (S)   , CURRENT, TOTAL, AVG : ', vp_time_current , vp_time, vp_time/iglobstep

				if (move_particles) then
					write (*,*)
					WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED IN DEM (S),   CURRENT, TOTAL, AVG : ', dem_time_current, dem_time, dem_time/iglobstep
					WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED IN H+DEM(S),  CURRENT, TOTAL, AVG : ', dem_time_current+hydro_time_current, dem_time+hydro_time, (dem_time+hydro_time)/iglobstep
				endif

				if (iscalon==1) then
					write (*,*)
					WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED IN SC (S),    CURRENT, TOTAL, AVG : ', scalar_time_current, scalar_time, scalar_time/iglobstep
					WRITE(*,'(A,3(2x,g15.8))') 'CPUTIME USED H+D+S (S),    CURRENT, TOTAL, AVG : ', dem_time_current+hydro_time_current+scalar_time_current, dem_time+hydro_time, (dem_time+hydro_time+scalar_time)/iglobstep
				endif

				write (*,*)
				WRITE(*,'(A,4(2x,g15.8))') 'CPUTIME USED (H:M:S)', CPUTIME_HRS,CPUTIME_MIN, CPUTIME_USED
				WRITE(*,'(A,4(2x,g15.8))') 'TIME SINCE LAST RES, RES_TIME', time_since_last_res, saveforres_time

				IF(CPUTIME_HRS.GT.WTIME_MAXHRS) THEN 
					PRINT*,'KIILING THE JOB BECAUSE WTIME_MAXHRS EXCEEDED'
					PRINT*, 'CONVERGENCE NOT YET REACHED'
					killjob = .TRUE.
				else
					killjob = .FALSE.
				endif
			end if

			BROADCAST_LOGICAL(killjob,1,NODE_ZERO,decomp_group)
			BROADCAST_DOUBLE(time_since_last_res,1,NODE_ZERO,decomp_group)
			!       PRINT*,'KILL JOB =', killjob       
			if(killjob)goto 999

#if !PARALLEL
			CALL CPU_TIME (iterend_time)
			avg_iter_time = avg_iter_time + (iterend_time - iterstart_time)
#else
			iterend_time = MPI_WTIME()
			avg_iter_time = (iterend_time - iterstart_time)
			GLOBAL_DOUBLE_MAX(avg_iter_time,global_max_iter_time,1,decomp_group)
			global_avg_iter_time = global_avg_iter_time + global_max_iter_time
#endif

!write (*,"(1a,3d15.7)") "MPG = ", mpg
!grav(:) = - mpg(:) / (rhof*(one-maxvolfrac) + rhos*maxvolfrac)
!archno = sqrt( dot_product(grav,grav)  ) * dia_phys**3 * (rhos-rhof)/rhof / vis**2
!write (*,"(1a,3d15.7)") "gravty = ", grav
!write (*,"(1a,3d15.7)") "ARCHNO = ", archno





!			taup = rhos/rhof * dia_phys**2/18/vis
!			tauf = dia_phys / sqrt(dot_product(ufmean(:)-usmean(:), ufmean(:)-usmean(:)))
!			froude = mpg(1) * taup**2 / doml(1)
!			if (I_AM_NODE_ZERO) then
!				write (*,"(1a,3d15.7)") "UFMEAN = ", ufmean
!				write (*,"(1a,3d15.7)") "USMEAN = ", usmean
!				write (*,"(1a,3d15.7)") "MPG    = ", mpg
!				write (*,"(1a,1d15.7)") "RE     = ", RE
!				write (*,"(1a,3d15.7)") "DOML   = ", doml
!				write (*,"(1a,1d15.7)") "FROUDE = ", froude
!				write (*,"(1a,3d15.7)") "TAUP, TAUF, ST = ", taup, tauf, taup/tauf
!			endif
!			PARALLEL_FINISH()
!			STOP
		!!	call reynolds_stress_tensor

			if (mod(iglobstep,saveitns)==1) then
			!!	call compute_sijsij
				if (write_snapshot) call flow_snapshot2
!			!!	call compute_AiVi
			endif

			!output intermediate restart file for every 100 steps
!!$

			IF ((saveforrestart.eq.1).AND.(time_since_last_res.gt.saveforres_time)) then 
				count_resfiles = count_resfiles+1
				if(I_AM_NODE_ZERO)then
					WRITE(*,'(A,g17.8,A,g17.8,A)')'WRITING RESTART FILES BECAUSE time_since_last_res (= ', &
						&	 time_since_last_res,') is > saveforres_time (= ', saveforres_time, ')'
					WRITE(ounit,'(A,g17.8,A,g17.8,A)')'WRITING RESTART FILES BECAUSE time_since_last_res (= ', &
						&	time_since_last_res,') is > saveforres_time (= ', saveforres_time, ')'
					WRITE(*,'(A,i4)')'NUMBER OF TIMES RESTART FILES WRITTEN IN THIS RUN = ', count_resfiles
					WRITE(ounit,'(A,i4)')'NUMBER OF TIMES RESTART FILES WRITTEN IN THIS RUN = ', count_resfiles
				endif
				CALL save_part_restart           
				time_since_last_res = zero
			endIF
			FIRST_PASS = .FALSE.   

                   if(iscalon.eq.1.and.mod(iglobstep,10).eq.0) call output_scal

                 !!  call output
 !!   go to 5555		
                ENDDO

		if(I_AM_NODE_ZERO)then 
			WRITE(*,*) 'HERE BECAUSE CONVERGENCE REACHED'
			WRITE(*,'(A,2(2x,g12.5))')'FERROR_HIST AND TOL = ', ferror_hist, TOL_FERROR
			IF(ISCALON.EQ.1) WRITE(*,'(A25,g12.5)')'NU ERROR HIST = ', nu_error_hist 
			OPEN(unit = 1000, file=conv_file, status="replace")
			close(1000, status="keep")
		end if
           
               !!  call output
            if(iscalon.eq.1.and.mod(iglobstep,1).eq.0) call output_scal


999	continue 

		IF(saveforrestart.EQ.1) THEN 
			count_resfiles = count_resfiles+1
			if(I_AM_NODE_ZERO)then
				WRITE(*,'(A,g17.8)')'WRITING RESTART FILES AT THE END: TIME_SINCE_LAST_RESTART FILE WRITTEN = ', time_since_last_res
				WRITE(ounit,'(A,g17.8)')'WRITING RESTART FILES AT THE END: TIME_SINCE_LAST_RESTART FILE WRITTEN = ', time_since_last_res

				WRITE(*,'(A,i4)')'NUMBER OF TIMES RESTART FILES WRITTEN IN THIS RUN = ', count_resfiles
				WRITE(ounit,'(A,i4)')'NUMBER OF TIMES RESTART FILES WRITTEN IN THIS RUN = ', count_resfiles
			end if
			final_restart = .FALSE.
			CALL save_part_restart
		ENDIF

	!!	call output 
!5555     continue 
		if(iscalon.eq.1) call output_scal

111	FORMAT(i4)
		if(I_AM_NODE_ZERO)then
			PRINT*,'Done with writing the output'
			CLOSE(unit_t_his,status='keep')
			CLOSE(unitnormdrag,status='keep')
			CLOSE(unitnormdragchem,status='keep')
			CLOSE(unitdragtavg,status='keep')
			CLOSE(unitdrag_comps,status='keep')
			CLOSE(sphrunit,status='keep')
			if(move_particles)then
				CLOSE(unitpartinfo,status='keep')
			end if
			!if(base_state_saved) CLOSE(partunit, status='keep')
			!t = 20.0
			!Print*,'time before entering test flux = ', t
			!call calc_anal

			CALL GET_RUN_ID 
			CALL screen_separator(80,'*')
			WRITE(*, '(A40,i2,A,i2,A,i4)') 'SIMULATION ENDED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
			WRITE(*, '(A40,i2,A,i2,A,i2)') 'AT TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND

			call separator(ounit,62,'-')

			WRITE(ounit, '(A40,i2,A,i2,A,i4)') 'SIMULATION ENDED ON DATE (MM/DD/YYYY):', ID_MONTH,'/', ID_DAY,'/', ID_YEAR
			WRITE(ounit, '(A40,i2,A,i2,A,i2)') 'AT TIME (HH:MM:SS)  ', ID_HOUR,':', ID_MINUTE,':', ID_SECOND
#if !PARALLEL       
			WRITE(*, '(A40,2x,g17.5, 2x,A10 )') 'COST/GRIDCELL/ITER:', avg_iter_time/real(iglobstep, prcn), 'SECONDS'
#else
			WRITE(*, '(A40,2x,g17.5, 2x,A10 )') 'COST/GRIDCELL/ITER:', global_avg_iter_time/real(iglobstep, prcn), 'SECONDS'
#endif
			CALL screen_separator(80,'*')
			call separator(ounit,62,'-')
		endif
		PARALLEL_FINISH()
		STOP
	END SUBROUTINE main
END PROGRAM dtibm_main

SUBROUTINE part_snapshot
	USE precision
	USE global_data
	USE general_funcs
	USE dependent_functions

	IMPLICIT NONE
	Integer  :: m
	LOGICAL, SAVE :: first_time_here=.TRUE.
	REAL(prcn) :: ucg
	CHARaCTER*80 :: FILENAME
	CHARACTER(LEN=80) :: formfile

	formfile='formatted'

	!if(irestart.eq.1)first_time = .FALSE.
	IF(first_time_here)THEN
		FILENAME = TRIM(RUN_NAME)//'_part_snapshot.dat'
		CALL  RUN_TIME_FILE_OPENER(partunit,FILENAME, formfile)
		first_time_here = .FALSE.
	END IF

	WRITE(partunit,*)'ZONE'
	WRITE(partunit,*)t
	DO m=1,nbody
		WRITE(partunit,'(10(2x,f12.8))')  xc(m,1), xc(m,2), xc(m,3), radbdy(m),velbdy(m,1:ndim),force(m,1:ndim)
	enddo
END SUBROUTINE part_snapshot

SUBROUTINE u_periodic_bc
	USE precision  
	USE constants  
	USE global_data
	USE scalar_data
	IMPLICIT NONE  
	INTEGER :: j, k 
	REAL(prcn) :: a(4), b(4), c(4), r1(4), sol(4)
	DO k = 1, mz
		DO j = 1, my2 
			uin(j,k,:)  = u(mx1, j,k,:)
			uout(j,k,:)  = u(2, j,k,:)
			IF(iscalon.eq.1) THEN 
				phiin(j,k,:) = phif(mx1,j,k,:)
				phiout(j,k,:) = phif(2,j,k,:)
			end IF
		END DO
	END DO
	a(1:4) = one
	c(1:4) = one
	if(I_AM_NODE_ZERO)then
		b(1) = 3.d0
		b(2) = 3.d0
		r1(1) = 1.d0
		r1(2) = 1.d0
	else
		b(1) = 2.d0
		b(2) = 6.d0 
		r1(1) = 3.d0
		r1(2) = 0.5d0 
	endif
!    CALL tridag(a,b,c,r,sol,4)
!    CALL mpi_tridag(a(1:2),b(1:2),c(1:2),r1(1:2),sol,2)
END SUBROUTINE u_periodic_bc
