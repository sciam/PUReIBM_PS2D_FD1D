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

module writeoutput
#include "ibm.h"
  Use precision  
  Use constants 
  Use global_data
  USE errormesgs
  USE general_funcs
  Use nlmainarrays 
  USE nlarrays, ONLY : ff1=>uf1,ff2=>uf2,ff3=>uf3
  USE nlmainarrays, velr=>ubcp, pr=>pbcp
  !Use bcsetarrays :  verl_mxf=>omega, pr=>dudy
  Use dependent_functions
  USE postproc_funcs
  !------------------------------------------------------------------------
  !       Edited, commented and partially re-written by Pravin & Dr.Shankar
  !       Formatted output to suit Tecplot Specifications.
  !------------------------------------------------------------------------
Contains 
  subroutine output

    implicit none 


    !-----------------------------------------------------------------------
    !	local variables
    !-----------------------------------------------------------------------

    !complex(prcn) ::  ff1(my2,mz),ff2(my2,mz),ff3(my2,mz) 
    INTEGER, SAVE  :: count_routine=0
    !real(prcn), dimension(:,:,:),allocatable, Target :: pr
    real(prcn) ::  vort(ndim)
    !real(prcn),dimension(:,:,:,:),allocatable, Target ::  velr_mxf
    !        real(prcn) ::  fr(mx,my,mz,ndim)
    CHARACTER*80 :: FILENAME1, FILENAME2, FILENAME3,  FILENAME4,FILENAMELOC
    CHARACTER*8:: FILENAME,FlTEMP
    LOGICAL :: filexist, isopen
    INTEGER :: unit1, count_tmp ,unit2, idim,jdim,node,strlen
    integer::i,j,k,m,n, unitno, meshunit, ii, nhbins, nrbins, ibin, dimi,dimj, iphs

    real(prcn) ::  umag, wt(nbody),fmin,fmax,ftemp(nbody) ,divur(mx,my,mz), aivj(ndim,ndim), vchar, fchar
    real(prcn), DIMENSION(:), allocatable :: hist, rad_bin
    real(prcn), DIMENSION(:,:,:), allocatable :: fij
    LOGICAL :: rescaling

    vchar = dsqrt(gran_temp)
    fchar = 6*pi*vis*dia_phys*ucharmod
    count_routine  = count_routine + 1

    IF(count_routine.eq.3) count_routine = 1
    if(I_AM_NODE_ZERO) then
       WRITE(*,*) 'COUNT FOR OUTPUT FILES= ', COUNT_ROUTINE
       WRITE(ounit,*) 'COUNT FOR OUTPUT FILES= ', COUNT_ROUTINE
    end if
    if(I_AM_NODE_ZERO)then
#if PARALLEL    
       write (filename,fmt="('NODE',i2.2,'_',i1)") myid,count_routine
       FILENAMELOC=""
#else
       WRITE(FILENAME,'(I1)')count_routine
#endif
    end if

    if(I_AM_NODE_ZERO)then
       FILENAME1 = TRIM(RUN_NAME)//'_U_'//TRIM(FILENAME)//'.dat'
    else
       FILENAME1 = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (fltemp,fmt="('NODE',i2.2,'_',i1)") node,count_routine
          FILENAMELOC = TRIM(RUN_NAME)//'_U_'//TRIM(FlTEMP)//'.dat'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(filename1,strlen,node_zero,0,1,decomp_group,status)
    end if

    if(I_AM_NODE_ZERO)then
       WRITE(FILENAME2, '(I1)')count_routine	
       WRITE(FILENAME3, '(I1)')count_routine	
       WRITE(FILENAME4, '(I1)')count_routine	
       FILENAME2 = TRIM(RUN_NAME)//'_FORCE_'//TRIM(FILENAME2)//'.dat'
       FILENAME3 = TRIM(RUN_NAME)//'_FPDF_'//TRIM(FILENAME3)//'.dat'
       FILENAME4 = TRIM(RUN_NAME)//'_FCORR_'//TRIM(FILENAME4)//'.dat'
    end if

    INQUIRE(FILE=FILENAME1,EXIST=filexist,OPENED=isopen)

    unit1 = getnewunit(minunitno, maxunitno)
    IF (.NOT.filexist) THEN
       OPEN(unit1,FILE=FILENAME1, status='new')
    ELSEIF(filexist.AND..NOT.isopen) THEN
       OPEN(unit1,FILE=FILENAME1, status="replace")
    ENDIF

    WRITE(unit1, '("# generated at time = ",g12.5)') t
    IF(write_output)then
       call calc_velreal(velr(1:nx,1:my,1:mz,1:ndim))
       call calc_pressure
    END IF
    !       Transfer pressure to subroutine that computes pressure 
    !       on the periphery of the sphere
    !       Note: here the entire pressure matrix is not required..
    !       this can be improved later.
    
!!$    call deterpressph!(pr,velr_mxf)
    

!!$    
    if(write_output.and.iscalon.ne.1)then
       mesh_vel(:) = velbdy(1,:)
       write(unit1,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" ', ' "UY" &
            &',' "UZ" ',' "P" ' !, ' "nl2" ', ' "nl3" '
       write(unit1,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
       
       write(*,*) 'Mesh velocity in output : ', mesh_vel(1:ndim)
       do k=1,1 !mz
          do j=1,my
             do i=1,nx !mx
                
!!$                umag = DSQRT(velr(i,j,k,1)**2. + velr(i,j,k,2)**2. + velr(i,j,k,3)**2.)
!!$                write(unit1,*)(GLOBAL_INDEX(i)),(j),(k),(velr(i,j,k,1))&
!!$                     &,(velr(i,j,k,2)),(velr(i,j,k,3)),pr(i,j,k)!/(half*upi(1)*upi(1))!,divur(i,j,k)!,nlbcp(i,j,k,1),nlbcp(i,j,k,2),nlbcp(i,j,k,3)!!velr_mxf(i,j,k,3)
!!$                !write(14,*)(i),(j),(k),divur(i,j,k)
!              !  if(fluid_atijk(i,j,k))then
                   write(unit1,*)(GLOBAL_INDEX(i)),(j),(k),(velr(i,j,k,1)-mesh_vel(1))&
                        &,(velr(i,j,k,2)-mesh_vel(2)),(velr(i,j,k,3)-mesh_vel(3)) !,pr(i,j,k)
                  if(j.eq.1) write(*,*)u(i,j,k,1)
!                else
!                   write(unit1,*)(GLOBAL_INDEX(i)),(j),(k),zero,zero,zero,  zero
!                end if
             enddo
          enddo
       enddo
    end if

    close(unit1, status='keep')
    !close(14, status='keep')
    
  call velocity_pdf

    if(I_AM_NODE_ZERO)then
		INQUIRE(FILE=FILENAME2,EXIST=filexist,OPENED=isopen)
		IF (.NOT.filexist) THEN
			OPEN(unit1,FILE=FILENAME2, status='new')
			write(unit1, *)nbody
			do iphs = 1, nphases
				write(unit1,*)phase_array(iphs)%npart
			end do
		ELSEIF(filexist.AND..NOT.isopen) THEN
			OPEN(unit1,FILE=FILENAME2, status="old",position="append")
		ENDIF

!!$       WRITE(unit1, '("# generated at time = ",g12.5)') t
!!$       write(unit1,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "FX" ',' "FY" ',' "FZ" '
		write(unit1,*) t/t_conv
		do m = 1, nbody
			write(unit1,'(9(2x,g17.8))') xc(m,1:ndim), velbdy(m,1:ndim), force(m,1:ndim)
		end do
		close(unit1, status='keep')

stop


!!$----------------------------------------------------------------------------------------------------!
!!$                                     Compute the pdf of forces
!!$----------------------------------------------------------------------------------------------------!
       
       if(nbody.gt.1)then
		INQUIRE(FILE=FILENAME3,EXIST=filexist,OPENED=isopen)
          IF (.NOT.filexist) THEN
             OPEN(unit1,FILE=FILENAME3, status='new')
          ELSEIF(filexist.AND..NOT.isopen) THEN
             OPEN(unit1,FILE=FILENAME3, status="replace")
          ENDIF
!!$
!!$
!!$       nhbins = 20

!!$       ALLOCATE(hist(nhbins))

!!$       do m = 1, nbody
!!$          wt(m) = 1.d0/real(nbody, prcn)
!!$       end do
!!$       
!!$       do idim=1,ndim
!!$          ftemp(1:nbody) = force(1:nbody,idim)/fchar
!!$          CALL histogram(ftemp,wt,nbody,nhbins,fmin,fmax,hist)
!!$          PRINT*,"i = ", idim, "fmin = ", fmin, "fmax =", fmax
!!$          CALL plothist(hist(1:nhbins),fmin,fmax,nhbins,unit1,t,t_conv)
!!$       end do
!!$       
!!$       DEALLOCATE(hist)

          close(unit1, status='keep')

!!$----------------------------------------------------------------------------------------------------!
!!$                                 Compute the force-force correlation of particles
!!$----------------------------------------------------------------------------------------------------!
		INQUIRE(FILE=FILENAME4,EXIST=filexist,OPENED=isopen)
          IF (.NOT.filexist) THEN
             OPEN(unit1,FILE=FILENAME4, status='new')
          ELSEIF(filexist.AND..NOT.isopen) THEN
             OPEN(unit1,FILE=FILENAME4, status="replace")
          ENDIF

          nrbins = 200
          ALLOCATE(fij(ndim,ndim,nrbins), rad_bin(nrbins))
          rescaling = .TRUE.
          CALL calc_vector_correlation(nbody, xc(1:nbody,1:ndim), radbdy(1:nbody), my, mxf, xperiodic, nrbins,&
               & ndim, rescaling, fij, rad_bin, force(1:nbody,1:ndim))
          do dimj =1,ndim
             write(unit1,*) 'Zone'
             do ibin=1,nrbins
                write(unit1,'(4(2x,f12.8))') rad_bin(ibin), fij(1:ndim,dimj,ibin)
             end do
          end do
          DEALLOCATE(fij,rad_bin)
          close(unit1,status='keep')
       end if
    
       !-----------------------------------------------------------------------
!!$    CLOSE(unit1, status="keep")    
       
       
    endif
    
    IF(count_routine.eq.1) count_tmp = 2
    IF(count_routine.eq.2) count_tmp = 1
    
    if(I_AM_NODE_ZERO)then
#if PARALLEL    
       write (filename,fmt="('NODE',i2.2,'_',i1)") myid,count_tmp
       FILENAMELOC=""
#else
       WRITE(FILENAME,'(I1)')count_tmp
#endif
    end if
    
    if(I_AM_NODE_ZERO)then
       FILENAME1 = TRIM(RUN_NAME)//'_U_'//TRIM(FILENAME)//'.dat'
    else
       FILENAME1 = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (fltemp,fmt="('NODE',i2.2,'_',i1)") node,count_tmp
          FILENAMELOC = TRIM(RUN_NAME)//'_U_'//TRIM(FlTEMP)//'.dat'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(filename1,strlen,node_zero,0,1,decomp_group,status)
    end if


    if(I_AM_NODE_ZERO)then
       WRITE(FILENAME2, '(I1)')count_tmp
       WRITE(FILENAME3, '(I1)')count_tmp	
       WRITE(FILENAME4, '(I1)')count_tmp	

       FILENAME2 = TRIM(RUN_NAME)//'_FORCE_'//TRIM(FILENAME2)//'.dat'
       FILENAME3 = TRIM(RUN_NAME)//'_FPDF_'//TRIM(FILENAME3)//'.dat'
       FILENAME4 = TRIM(RUN_NAME)//'_FCORR_'//TRIM(FILENAME4)//'.dat'
    end if
    
    INQUIRE(FILE=FILENAME1,EXIST=filexist,OPENED=isopen)
    IF(filexist) THEN 
       OPEN(unit1, FILE=FILENAME1, form="formatted", status="replace")
       close(unit1, status='delete')
       if(I_AM_NODE_ZERO)then
          OPEN(unit1, FILE=FILENAME2, form="formatted", status="replace")
          close(unit1, status='delete')
          if(nbody.gt.1)then
             OPEN(unit1, FILE=FILENAME3, form="formatted", status="replace")
             close(unit1, status='delete')
             OPEN(unit1, FILE=FILENAME4, form="formatted", status="replace")
             close(unit1, status='delete')
          end if
       end if
    ENDIF

    !WRITE(*,*)'UMAX = ', MAXVAL(velr(:,:,:,1))
111 format(e12.4)
112 format(e12.4,1x,e12.4)
113 format(3e12.4)
114 format(2e12.4,i4)
221 format(i4)
222 format(2i4)
223 format(3i4)
225 format(5i4)
211 format(i4,2e12.4)

  end subroutine output

  subroutine write_velrmid(velr)
    implicit none 
    real(prcn), DImension(:,:,:,:), Intent(out) :: velr
    integer :: i,j,k, n 
    open(unit=155,file='umid.dat',form='formatted',status='unknown')
    write(155,*)'VARIABLES= ',' "X" ',' "Y" ',' "UX" ',   &
         &    ' "UY" ',' "UZ" ' 
    write(155,*)'ZONE F=POINT, I=', mx,  ', J=', my
    k = mz/2
    do j=1,my
       do i=1,mx
          write(155,*)(i-1)*dx,(j-1)*dy,velr(i,j,k,1)&
               &,velr(i,j,k,2),velr(i,j,k,3) 
       enddo
    enddo
    close(155, status = 'keep')
  end subroutine write_velrmid
        subroutine calc_avr_var(nvar, var, avr, variance)
                implicit none
                integer, intent(in) :: nvar
                real(prcn), intent(in) :: var(nvar)
                real(prcn), intent(out) :: avr, variance
                integer :: ivar
                avr = sum(var(:))/nvar
                variance = zero
                do ivar=1, nvar
                        variance = variance + (var(ivar)-avr)**2
                enddo
        end subroutine calc_avr_var







        subroutine velocity_pdf
!               use boundary_condition
                use bcsetarrays, only : ppr, diffn
                implicit none
                real(prcn), allocatable :: f_acc(:,:), p_acc(:,:), wt(:),hist(:), radbin(:), tmp_array(:)
                real(prcn) :: u_min, u_max
                integer :: ibody, idim, i, j, k, count, unit1, nvar, nvar1,nvar2, ip, im,nbins
                character*50 filename1, filename2
                call screen_separator(30,'^')
                write (*,*) "IN FLUID_PARTICLE_ACCELERATION"

               nbins =500
                ! if (.not.post_no_flow_mem_alloc) then
                        allocate(f_acc(ndim, count_fluid))
                        f_acc = zero
                        p_acc = zero
!                       call calc_pgrad
!                       call calc_visc
                        count = 0
                        do k=1, mz
                                do j=1, my
                                        do i=1, nx
                                                if (fluid_atijk(i,j,k)) then
                                                        count = count+1
                                                        do idim=1, ndim
                                                                f_acc(idim,count) = velr(i,j,k,idim) !- ufmean(idim)
                                                        enddo
                                                endif
                                        enddo
                                enddo
                        enddo
                        unit1 = 1
                        filename1 = trim(run_name)//"_pdf_u.dat"
                        open  (unit=unit1, file=trim(filename1),status="replace")
                        allocate(wt(count_fluid))
                        wt = 1d0/count_fluid
                        allocate(hist(nbins), radbin(nbins))
                        allocate(tmp_array(count_fluid))
                        do idim=1, ndim
                                tmp_array(:) = f_acc(idim,:)
                                call make_histogram(tmp_array(:), count_fluid,nbins, radbin, hist)
                                write (unit1,*) "zone"
                                do i=1, nbins
                                        if (hist(i)>small_number) write(1,"(2d15.7)") radbin(i), hist(i)
                                enddo
                        enddo
                        close (unit1)
                        deallocate(tmp_array)
              !  else
              !          nvar  = 2
              !          nvar1 = 1
              !          nvar2 = 1
              !          line = nbins
              !          filename1 = "_pdf_u.dat"
              !          call mis_average(nvar, nvar1, nvar2, filename1, nbins,.false.)
              !  endif
                call screen_separator(30,'-')
        end subroutine velocity_pdf

        subroutine make_histogram(array, n, nbins, radbin, hist)
                implicit none
                integer, intent(in) :: n, nbins
                real(prcn), intent(inout) :: array(n)
                real(prcn), intent(out) :: radbin(nbins), hist(nbins)
!               real(prcn) :: tmp_array(n)
                real(prcn) :: mean, var, sd, left, right, dr
                integer :: i, ibin
      !          call calc_avr_var(n, array, mean, var)
     !           sd = sqrt(var/n)
    !            array(:) = (array(:) - mean) / sd
  !!write (*,*) mean, sd
                call calc_avr_var(n, array, mean, var)
                sd = sqrt(var/n)
write (*,*) mean, sd
                left = minval(array(:))
                right = maxval(array(:))
                write (*,*) "LEFT, RIGHT = ", left, right
                dr = (right-left) / nbins
                do i=1, nbins
                        radbin(i) = left + (i-.5)*dr
                enddo
                hist = zero
                do i=1, n
                        ibin = (array(i)-left) / dr + 1
                        if (ibin>nbins) ibin = nbins
                        hist(ibin) = hist(ibin) + 1
                enddo
                hist(:) = hist(:) / n / dr
                write (*,*) "SUM OF HIST = ", sum(hist(:)) * dr
                radbin(:) = radbin(:)
        end subroutine make_histogram
end module writeoutput





    
!!$    meshunit  = getnewunit(minunitno,maxunitno)
!!$ 
!!$    open(unit=meshunit,file='mesh.dat',form='formatted',status='unknown')
!!$    write(meshunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "VAR" '
!!$    write(meshunit,*)'ZONE T = "UX"'
!!$        
!!$    Write(meshunit, *) 'ZONE DATAPACKING=POINT, I=',mxf,', J=',my,', K=', mz
!!$    do k=1,mz
!!$       do j=1,my
!!$          do i=1,mxf
!!$             write(meshunit, '(4(2x,e17.8))')DREAL((i-1)),DREAL((j-1)),DREAL((k-1)), velr_mxf(i,j,k,1) 
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    write(meshunit,*)'ZONE T = "UY"'
!!$
!!$    Write(meshunit, *) 'ZONE DATAPACKING=POINT, I=',mxf,', J=',my,', K=', mz
!!$    write(meshunit,*)', VARSHARELIST= ([1-3]=1)'
!!$    do k=1,mz
!!$       do j=1,my
!!$          do i=1,mxf
!!$             
!!$             write(meshunit,('3(2x,e17.8)'))velr_mxf(i,j,k,2)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    
!!$    close(meshunit, status="keep")
