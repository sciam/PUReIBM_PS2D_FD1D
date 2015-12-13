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

!-------
! Dependent module
!-------
! Routines related to reading/writing restart files
! COMMENTS: 
!-------
! Code:   DTIBM
! Author: Rahul Garg
!         ISU
! Date:   October 2005
!-------
MODULE restart_funcs
#include "ibm.h"
  USE precision          ! Independent modules
  USE general_funcs
  
  USE global_data        ! Dependent modules
  USE scalar_data 
  USE errormesgs

  Implicit none
  Private
  Public:: save_part_restart, read_restart, delete_file, filename_gen
  
  !-------
CONTAINS
   !-------
  
  !-------
  !-------
  ! Save restart file
  !-------
  !-------
   SUBROUTINE save_part_restart
     IMPLICIT NONE

     INTEGER  :: count_routine
     CHARACTER*10 :: FILENAME1
     CHARACTER*80 :: FILENAME,stat
     LOGICAL :: filexist, isopen
     INTEGER :: runit,strlen,node
!     CALL hello
 !    PARALLEL_FINISH()
  !   STOP

     If (rstsave.eq."formatted") then
        !call save_formatted
     Elseif (rstsave.eq."unformatted") then
     
        runit = getnewunit(minunitno,maxunitno)
        !BROADCAST_INT(runit,1,NODE_ZERO,decomp_group)
        if (runit.lt.0) call printerror("newunit","runit")
        
        if(I_AM_NODE_ZERO)then
           FILENAME = TRIM(RUN_NAME)//'_RESTART'
           open(unit=runit,file=FILENAME,form="formatted"      &
                ,status="unknown")    
           read(runit, *) count_routine
           close(runit, status='keep')
        end if
        BROADCAST_INT(count_routine,1,NODE_ZERO,decomp_group)
        count_routine = count_routine+1
        IF(count_routine.eq.3) count_routine = 1

        if(I_AM_NODE_ZERO)then
           WRITE(*,*) 'COUNT FOR DUMPING RESTART FILES ', COUNT_ROUTINE
           WRITE(ounit,*) 'COUNT FOR DUMPING RESTART FILES ', COUNT_ROUTINE
        end if

	if (igeometry==1) return

        CALL FILENAME_GEN(FILENAME,'u',count_routine)
        
        
        INQUIRE(FILE=TRIM(FILENAME),EXIST=filexist,OPENED=isopen)

        IF (.NOT.filexist) THEN
           stat="new"
           call save_unformatted(count_routine,stat)
        ELSEIF(filexist.AND..NOT.isopen) THEN
           stat="replace"
           call save_unformatted(count_routine,stat)
        ENDIF

     Endif


     count_restart = count_routine
     !-------
     ! write to output unit
     !-------

   END SUBROUTINE save_part_restart
  
   SUBROUTINE FILENAME_GEN(FILENAME,VAR_NAME,RES_COUNT)
     
     IMPLICIT NONE

     Character(LEN=*),Intent(out) :: filename
     Character(LEN=*),Intent(in) :: var_name
     INTEGER, INTENT(in) :: RES_COUNT

     CHARACTER*10 :: FILENAME1
     Character*80 :: FILENAMELOC
     Integer :: node, strlen

     if(I_AM_NODE_ZERO)then
#if PARALLEL

		!^^^^ MODIFIED FOR MORE PROCS
		if (nproc<=100) then
			write (filename1,fmt="('NODE',i2.2,'_',i1)") myid, res_count
		else
			write (filename1,fmt="('NODE',i3.3,'_',i1)") myid, res_count
		endif

#else
        WRITE(FILENAME1, '(I1)') res_count
#endif
        FILENAME = TRIM(RUN_NAME)//'_'//TRIM(var_name)//'_'//TRIM(FILENAME1)//'.rst'
        FILENAMELOC = ""
     else
        FILENAME = ""
     end if
     if(I_AM_NODE_ZERO)then
        do node=1,nproc-1

				!^^^^ MODIFIED FOR MORE PROCS
				if (nproc<=100) then
					write (filename1,fmt="('NODE',i2.2,'_',i1)") node, res_count
				else
					write (filename1,fmt="('NODE',i3.3,'_',i1)") node, res_count
				endif


           FILENAMELOC = TRIM(RUN_NAME)//'_'//TRIM(var_name)//'_'//TRIM(FILENAME1)//'.rst'
           !              PRINT*,'NODE ZERO', node, filenameloc
           SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
        end do
     else  
        RECV_STRING(filename,strlen,node_zero,0,1,decomp_group,status)
     end if

   END SUBROUTINE FILENAME_GEN


  !-------
  !-------
  ! If restart file is unformatted
  !-------
  !-------

   SUBROUTINE save_unformatted(count_routine,stat)
     Implicit none

     INTEGER, Intent(in) :: count_routine

     CHARACTER*80,Intent(in) :: stat
     INTEGER :: count_tmp, strlen, node, iphs, nerr_steps_tmp
     CHARACTER*80 :: FILENAME
     CHARACTER*10 :: FILENAME1
     LOGICAL:: filexist, isopen


     runit = getnewunit(minunitno,maxunitno)
     if (runit.lt.0) call printerror("newunit","runit")
     if(I_AM_NODE_ZERO)then
        WRITE(FILENAME1, '(I1)') count_routine 	
     end if

     !-------
     ! Basic grid data
     !-------
     !    PRINT*,'sat = ', myid, stat
     if(I_AM_NODE_ZERO)then
        FILENAME = TRIM(RUN_NAME)//'_sphere_config.rst'
        open(unit=runit,file=FILENAME,form=rstsave      &
             ,status="unknown")

        write(runit) nbody
        if(input_type.eq."random".or.input_type.eq."default")then
           write(runit)nphases
           do iphs = 1, nphases
              write(runit) phase_array(iphs)%npart
           end do
           do iphs = 1, nphases
              write(runit) phase_array(iphs)%dia
           end do
           write(runit) DOML(1:ndim)
        endif
        write(runit) char_length
        write(runit) xc(1:nbody, 1:3)
        write(runit) radbdy(1:nbody)
        write(runit) velbdy(1:nbody, 1:3)
        close(runit, status='keep')

        if(imove.eq.1)then
           FILENAME = TRIM(RUN_NAME)//'_sphere_config_'//TRIM(FILENAME1)//'.rst'
           open(unit=runit,file=FILENAME,form=rstsave      &
                ,status="unknown")

           write(runit) nbody
           if(input_type.eq."random".or.input_type.eq."default")then
              write(runit)nphases
              do iphs = 1, nphases
                 write(runit) phase_array(iphs)%npart
              end do
              do iphs = 1, nphases
                 write(runit) phase_array(iphs)%dia
              end do
              write(runit) DOML(1:ndim)
              write(runit) move_particles
           endif
           write(runit) char_length
           write(runit) xc(1:nbody, 1:3)
           write(runit) radbdy(1:nbody)
           write(runit) velbdy(1:nbody, 1:3)
           write(runit) frame_vel(1:ndim)
           write(runit) frame_pos(1:ndim)
           close(runit, status='keep')
        end if

     end if

     CALL FILENAME_GEN(FILENAME,'u',count_routine)
     
     open(unit=runit,file=FILENAME,form=rstsave      &
          ,status=stat)
     
     write(runit) t, dt
     write(runit)umean(1), umean(2), umean(3)
     write(runit)usmean_des(1), usmean_des(2), usmean_des(3)
     write(runit)ufmean_des(1), ufmean_des(2), ufmean_des(3)
     write(runit) cf 
     write(runit) nerr_steps

     write(runit) ferror, fold

     write(runit) ferror_array(1:nerr_steps)

     do iphs = 1, nphases
        write(runit) phase_array(iphs)%ferror, phase_array(iphs)%fold
        if(imove.eq.1)write(runit) phase_array(iphs)%grant_error, phase_array(iphs)%grant_old
     end do
     
     do iphs = 1, nphases
        write(runit) phase_array(iphs)%ferror_array(1:nerr_steps)
        if(imove.eq.1)write(runit) phase_array(iphs)%grant_array(1:nerr_steps)
     end do

#if PARALLEL
     write(runit)u(0:nx+1,1:my2,1:mz,1:ndim)
#else
     write(runit)u
#endif
     !write(runit)fluid_atijk(1:mxf,1:my,1:mz)
     close(runit, status='keep')
     
     CALL FILENAME_GEN(FILENAME,'nl',count_routine)
     
     open(unit=runit,file=FILENAME,form=rstsave      &
          ,status=stat)
     write(runit)nl
     
     close(runit, status='keep')
     
     CALL FILENAME_GEN(FILENAME,'p',count_routine)

     open(unit=runit,file=TRIM(FILENAME),form=rstsave      &
          ,status=stat)    
     write(runit)mpg(1),mpg(2),mpg(3)
#if PARALLEL
     write(runit)p(0:nx+1,1:my2,1:mz)
#else
     write(runit)p
#endif
     close(runit, status='keep')


     if(iscalon.eq.1) then

        CALL FILENAME_GEN(FILENAME,'scal',count_routine)
        
        open(unit=runit,file=TRIM(FILENAME),form=rstsave      &
             ,status=stat)    
        write(runit) nu_error, nu_old
        write(runit) nu_error_array(1:nerr_steps)
#if PARALLEL
        write(runit)phif(0:nx+1,1:my2,1:mz,1:nspmx)
        write(runit)nlphif(1:nx,1:my2,1:mz,1:nspmx)
#else
        write(runit)phif
        write(runit)nlphif(1:nx+1,1:my2,1:mz,1:nspmx)
#endif
        write(runit) phisurfall(1:nbody,1:nspmx)
        write(runit) phirmean(1:nspmx), heat_ratio(1:nspmx), heat_ratio_old(1:nspmx)!fphirmean(1:nspmx)
        close(runit, status='keep')
     endif


     if(I_AM_NODE_ZERO)then
        FILENAME = TRIM(RUN_NAME)//'_RESTART'
        WRITE(FILENAME1, '(I1)')count_routine 	    
        open(unit=runit,file=TRIM(FILENAME),form="formatted"      &
             ,status="replace")    
        write(runit, *) FILENAME1
#if PARALLEL
        write(runit, *) nproc
#endif
        close(runit, status='keep')
     end if
     IF(count_routine.eq.1) count_tmp = 2
     IF(count_routine.eq.2) count_tmp = 1
!!$    PRINT*,'count_tmp = ', myid, count_tmp
     CALL delete_restart_files(count_tmp)


     if(I_AM_NODE_ZERO)then 
        call separator(ounit,62,'-')
        write(ounit,12,ADVANCE='NO')t
        write(ounit,*)TRIM(FILENAME)
        call separator(ounit,62,'-')
     endif
12   format('Saving restart data at time = ',E12.5,' in &
          &file with extension ')

   END SUBROUTINE save_unformatted

  SUBROUTINE delete_restart_files(count_routine)
    INTEGER, intent(in):: count_routine
    CHARACTER*80 :: FILENAME, form, FILENAMELOC
    CHARACTER*8 :: FILENAME1, FILENAME2
    INTEGER :: runitno, node, strlen
    LOGICAL:: filexist, isopen
    
    CALL FILENAME_GEN(FILENAME,'u',count_routine)
    
    INQUIRE(FILE=FILENAME,EXIST=filexist,OPENED=isopen)
    IF(filexist) THEN 
       runitno = getnewunit(minunitno, maxunitno)
       CALL delete_file(FILENAME,rstsave,runitno)
       
       CALL FILENAME_GEN(FILENAME,'nl',count_routine)
       CALL delete_file(FILENAME,rstsave,runitno)
       
       CALL FILENAME_GEN(FILENAME,'p',count_routine)
       CALL delete_file(FILENAME,rstsave,runitno)
       
       IF(ISCALON.EQ.1) then 
          CALL FILENAME_GEN(FILENAME,'scal',count_routine)             
          CALL delete_file(FILENAME,rstsave,runitno)
       end IF
       IF(I_AM_NODE_ZERO)then
          IF(IMOVE.EQ.1)then
             CALL FILENAME_GEN(FILENAME,'sphere_config',count_routine)             
             CALL delete_file(FILENAME,rstsave,runitno)
          END IF
       END IF
    ENDIF
  end SUBROUTINE delete_restart_files
  
  SUBROUTINE delete_file(filename,formtype,runitno)
    CHARACTER*80, INTENT(in) :: FILENAME
    CHARACTER(LEN=11), INTENT(in):: FORMTYPE
    INTEGER, Intent(in) :: runitno
!    PRINT*,'*******', myid, FILENAME
    OPEN(runitno, FILE=FILENAME, form=formtype, status="replace")
    close(runitno, status='delete')
  end SUBROUTINE delete_file
  

  !------
  ! Read from restart file
  !-------
  
  SUBROUTINE read_restart
    Implicit none
    Logical:: fileexists

    !INQUIRE(file=restartfile, exist=fileexists)
    !If (fileexists) then
    If (rstread.eq."formatted") then
       ! call read_formatted
    Elseif (rstread.eq."unformatted") then
       call read_unformatted
    Endif
    !Else
    if(I_AM_NODE_ZERO)then
       call separator(ounit,52,'-')
       write(ounit,12)t
12     format('Read field from restartfile at time = ', F12.5)
       call separator(ounit,52,'-')
    end if
  END SUBROUTINE read_restart
  
  SUBROUTINE read_unformatted
    Implicit none
    INTEGER :: count_tmp, i, nerr_steps_tmp, node,strlen, iphs
    CHARACTER*80 :: FILENAME, FILENAMELOC
    CHARACTER*10 :: FILENAME1,FILENAME2
    LOGICAL :: filex,isopen
    if(I_AM_NODE_ZERO)then
       CALL screen_separator(80,'*')
       WRITE(*,*) 'IN READ UNFORMATTED RESTART'
    end if
    
    runit = getnewunit(minunitno,maxunitno)
    
    if (runit.lt.0) call printerror("newunit","runit")
    
    count_tmp = count_restart
    
    CALL FILENAME_GEN(FILENAME,'u',count_tmp)

    open(unit=runit,file=FILENAME,form=rstsave)
    
    read(runit)tstart, dt
    read(runit)umean(1), umean(2), umean(3)
    read(runit)usmean_des(1), usmean_des(2), usmean_des(3)
    read(runit)ufmean_des(1), ufmean_des(2), ufmean_des(3)

    read(runit)cf
    
    cforig = cf
    Write(*,*)'cf from restart = ', cf
    read(runit) nerr_steps_tmp
    
    read(runit)ferror, fold
    
    if(SIZE(ferror_array,1).ne.nerr_steps)then
       DEALLOCATE(ferror_array)
       ALLOCATE(ferror_array(nerr_steps))
    end if
    
    if(nerr_steps_tmp.gt.nerr_steps) then
       read(runit) ferror_array(1:nerr_steps)
    ELSE
       ferror_array = 1.d0
       read(runit) ferror_array(1:nerr_steps_tmp)
    end if
    
    do iphs = 1, nphases
       read(runit)phase_array(iphs)%ferror, phase_array(iphs)%fold
		!^^^^^^^ 06-24-2010 Mohammad: RESTARTING A FIXED BED FOR MOVING CASE ^^^^
		if (imove==1) then
			if (moving_from_fixed_bed) then
				phase_array(iphs)%grant_error = one
				phase_array(iphs)%grant_old = gran_temp
			else
				read(runit)phase_array(iphs)%grant_error, phase_array(iphs)%grant_old
			endif
		endif
		!------------------------------------------------------------------------
    enddo

    do iphs = 1, nphases
       if(SIZE(phase_array(iphs)%ferror_array,1).ne.nerr_steps)then
          if(ASSOCIATED(phase_array(iphs)%ferror_array)) DEALLOCATE(phase_array(iphs)%ferror_array)
          ALLOCATE(phase_array(iphs)%ferror_array(nerr_steps))
       end if
       if(imove.eq.1)then
          if(SIZE(phase_array(iphs)%grant_array,1).ne.nerr_steps)then
             DEALLOCATE(phase_array(iphs)%grant_array)
             ALLOCATE(phase_array(iphs)%grant_array(nerr_steps))
          end if
       end if
    end do

    do iphs = 1, nphases
       if(nerr_steps_tmp.gt.nerr_steps) then
          read(runit) phase_array(iphs)%ferror_array(1:nerr_steps)
       ELSE
          phase_array(iphs)%ferror_array = 1.d0
          read(runit) phase_array(iphs)%ferror_array(1:nerr_steps_tmp)
       end if
       if(imove.eq.1)then
		!^^^^^^^ 06-24-2010 Mohammad: RESTARTING A FIXED BED FOR MOVING CASE ^^^^          
		if (moving_from_fixed_bed) then
          	phase_array(iphs)%grant_array(1:nerr_steps) = one
		else
		     if(nerr_steps_tmp.gt.nerr_steps) then
		        read(runit) phase_array(iphs)%grant_array(1:nerr_steps)
		     ELSE
		        phase_array(iphs)%grant_array = 1.d0
		        read(runit) phase_array(iphs)%grant_array(1:nerr_steps_tmp)
		     end if
		endif
		!--------------------------------------------------------------------------
       end if
    end do

#if PARALLEL
    read(runit)u(0:nx+1,:,:,:)
#else
    read(runit)u
#endif
    close(runit, status='keep')
    
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A40,g12.5)')'READING THE RESTART FILE WRITTEN AT ', tstart
       WRITE(*,'(A40,g12.5)')'DT READ FROM RESTART FILES = ', dt
    end if
    ferror_hist = SUM(ferror_array(1:nerr_steps))/nerr_steps
    if(I_AM_NODE_ZERO)then
       WRITE(*,'(A25,2(2x,g12.5))')'FOLD, FERROR = ', FOLD, FERROR
       
       WRITE(*,'(A25,g12.5)') 'FERROR_HIST = ', ferror_hist
    end if
    
    do iphs = 1, nphases
       phase_array(iphs)%ferror_hist = SUM(phase_array(iphs)%ferror_array(1:nerr_steps))/nerr_steps
       if(I_AM_NODE_ZERO)then
          WRITE(*,'(A10,I4,A25,2(2x,g12.5))')'PHASE = ', iphs, 'FOLD, FERROR = ', phase_array(iphs)%fold, phase_array(iphs)%ferror
       
          WRITE(*,'(A10,I4,A25,g12.5)') 'PHASE = ', iphs, 'FERROR_HIST = ', phase_array(iphs)%ferror_hist
       end if
    end do
    
    CALL FILENAME_GEN(FILENAME,'nl',count_tmp)
    open(unit=runit,file=FILENAME,form=rstsave)
    
    read(runit)nl
    
    close(runit, status='keep')
    
    CALL FILENAME_GEN(FILENAME,'p',count_tmp)
    open(unit=runit,file=FILENAME,form=rstsave)
    
    read(runit)mpg(1),mpg(2),mpg(3)
#if PARALLEL
    read(runit)p(0:nx+1,:,:)
#else
    read(runit)p
#endif
    close(runit, status='keep')
    
    if(iscalon.eq.1.and.iscal_restart.eq.1) then
      write(*,*)'read_scal_file'
       CALL FILENAME_GEN(FILENAME,'scal',count_tmp)
       open(unit=runit,file=FILENAME,form=rstsave)
       
       read(runit) nu_error, nu_old
       
       if(nerr_steps_tmp.gt.nerr_steps) then
          read(runit) nu_error_array(1:nerr_steps)
       ELSE
          nu_error_array = 1.d0
          read(runit) nu_error_array(1:nerr_steps_tmp)
       end if
#if PARALLEL
       read(runit)phif(0:nx+1,:,:,1:nspmx)
       read(runit)nlphif(1:nx,1:my2,1:mz,1:nspmx)
#else
       read(runit)phif
       read(runit)nlphif(1:nx+1,1:my2,1:mz,1:nspmx)
#endif
    !!   read(runit) nlphif
write(*,*)'read_scal_filE2'
       read(runit) phisurfall(1:nbody,1:nspmx)
       read(runit)phirmean(1:nspmx),heat_ratio(1:nspmx), heat_ratio_old(1:nspmx) !fphirmean(1:nspmx)
       close(runit, status = "keep")
       
       if(I_AM_NODE_ZERO)then
          WRITE(*,*) 'IN SCALAR READ UNFORMATTED RESTART'
          
          DO i=1,nbody
             !WRITE(*,'(A,g12.5,2x,i2)')'phisurf of bodies =', phisurfall(i,1),i
             WRITE(ounit,'(A,g12.5,2x,i2)')'phisurf of bodies =', phisurfall(i,1),i
          END DO
          WRITE(*,'(A30,2(2x,g12.5))')'phirmean and fphirmean = ', phirmean, fphirmean 
          
          WRITE(*,'(A25,2(2x,g12.5))')'nu_old, nu_error = ', nu_old, nu_error
       end if
       nu_error_hist = SUM(nu_error_array(1:nerr_steps))/nerr_steps
       if(I_AM_NODE_ZERO) WRITE(*,'(A25,g12.5)')'NU_ERROR_HIST = ', NU_ERROR_HIST
    endif
    if(I_AM_NODE_ZERO)CALL screen_separator(80,'*')
    
  END SUBROUTINE read_unformatted

end MODULE restart_funcs
