module outputscalar
#include "../FLO/ibm.h"
  use scalar_data 
  use fftw_interface
  Use dependent_functions
  USE restart_funcs
  USE general_funcs
  USE nlmainarrays, velr=>ubcp, phir=>nlbcp
  Use nlarrays , only : ur1,uf1
  USE postproc_funcs
  USE bcsetarrays, ONLY : gradphi=>ppr,fr,diffn
  USE bc_scalar

  implicit none 
  Integer :: ix
  !real(prcn), Private,dimension(:,:,:,:), allocatable :: phir, velr

  PUBLIC :: compute_local_Nuss_x
contains
  subroutine output_scal
    implicit none 
    Integer :: i, j, k, count_tmp, runitno, nrbins, ibin, node, strlen, isp, idim,count_phi
    INTEGER, SAVE :: count_routine=0
    CHARACTER*80 :: SCALANDU_FOR, THETAVSNU3,SCALANDU,SCALONLY,&
         & junk_char, stat, filenameloc 
    LOGICAL :: filexist , isopen
    CHARACTER*10 :: filename2,filename1
    real(prcn), dimension(:), allocatable :: phi_corr, rad_bin
    real(prcn):: mean_phi(mx),mean_phi_v(nx),theta_u_fluc(mx),theta_u(mx),phi_phi(mx),theta_theta(mx),theta_mean(mx),umean_vel,axial,axial_con(nx)

    !CALL calc_velreal(velr(:,:,:,:))
    !Write(*,*)'SIZE : ', SIZE(velr,1)

    mean_phi = zero
    mean_phi_v =zero

 
    

 OPEN(200,FILE=TRIM(RUN_NAME)//'_meanphi.dat',form="formatted",status='unknown')    
    do idim = 1, ndim
#if PARALLEL
       do i=0,nx+1 !mx,1
#else
       do i = 1, nx+1
#endif   
          do j=1,my2,1
             do k=1,mz,1
                uf1(j,k) = u(i,j,k,idim)
             enddo
          enddo
          
          call ff2cr(uf1,ur1)	 
          do j=1,my,1
             do k=1,mz,1
                velr(i,j,k,idim) = ur1(j,k) + umean(idim)
             enddo
          enddo
       enddo
    enddo
umean_vel =zero
         do j =1,my
          do k=1,mz
          
           umean_vel= umean_vel + velr(1,j,k,1)
          enddo
         enddo
          umean_vel =umean_vel/my/mz 
mean_phi =zero
    do isp=1,nspmx
       do i=1,nx+1 !mx,1
          do j=1,my2,1
             do k=1,mz,1
                uf1(j,k) = phif(i,j,k,isp)
             enddo
          enddo
          
          call ff2cr(uf1,ur1)	 
          do j=1,my,1
             do k=1,mz,1
                phir(i,j,k,isp) = ur1(j,k) 
             enddo
          enddo
    
       !call compute_phim_heat_ratio_1(.true.)

        
        count_phi =zero
          do j=1,my,1
             do k=1,mz,1
              if(fluid_atijk(i,j,k)) then

                mean_phi(i)=mean_phi(i)+phir(i,j,k,isp)
                mean_phi_v(i)=mean_phi_v(i)+phir(i,j,k,isp)    
              endif
              count_phi = count_phi +1 
             enddo
          enddo

         !  write(200,*)(i-1)*dx,mean_phi(i)/count_phi
            mean_phi(i)= mean_phi(i)/count_phi !!mean phif
             mean_phi_v(i)= mean_phi_v(i)/my/mz !!mean phif 
            write(*,*)"m",i,mean_phi(i)

         enddo
       enddo
      !!   write(*,*)"m",i,mean_phi(i)
         call compute_phim_heat_ratio_1(.true.)      
        
        
        do isp=1,nspmx
theta_u_fluc =zero
theta_u =zero
phi_phi =zero
theta_theta =zero
theta_mean =zero
axial_con=zero
          do i=1,nx
              count_phi=0
          do j=1,my,1
          do k=1,mz,1
            if(fluid_atijk(i,j,k)) then
              count_phi = count_phi +1
              theta_u_fluc(i)  =theta_u_fluc(i)+(velr(i,j,k,1)-ufmean_des(1)) *phir(i,j,k,isp)/phim(i,isp) 
              theta_u(i)=theta_u(i)+(velr(i,j,k,1)-ufmean_des(1))*(phir(i,j,k,isp)-mean_phi(i))
              phi_phi(i) = phi_phi(i) +(phir(i,j,k,isp)-mean_phi(i))*(phir(i,j,k,isp)-mean_phi(i))
             theta_theta(i) = theta_theta(i)+phir(i,j,k,isp)/phim(i,isp)*phir(i,j,k,isp)/phim(i,isp)
              theta_mean(i) = theta_mean(i)+phir(i,j,k,isp)/phim(i,isp)
            !  theta_mean(i) = theta_mean(i)+(velr(i,j,k,1)-ufmean_des(1))*(phir(i,j,k,isp))
             axial_con(i) = axial_con(i)+(phir(i+1,j,k,isp)-phir(i-1,j,k,isp))/(dx*two)
            else
              theta_u_fluc(i)  =theta_u_fluc(i)+(velr(i,j,k,1)-umean(1)) *zero/phim(i,isp)
              theta_u(i)=theta_u(i)+(velr(i,j,k,1)-ufmean_des(1))*zero !! (zero-mean_phi(i))
              phi_phi(i) = phi_phi(i)+(zero-mean_phi(i))*(zero-mean_phi(i))
            endif
          enddo
          enddo
                  
                theta_u_fluc(i)= theta_u_fluc(i)/count_phi
                theta_u(i)=theta_u(i)/count_phi
                phi_phi(i)=phi_phi(i)/my/mz
                theta_theta(i) =theta_theta(i)/my/mz
                theta_mean(i) =theta_mean(i)/count_phi
                axial_con(i) = axial_con(i)/my/mz
              write(200,21)(i-1)*dx, mean_phi(i),phim(i,1)  !! theta_u_fluc(i)/ufmean_des(1),theta_u(i)/ufmean_des(1), axial_con(i) !-axial_con(i-1))/dx !!,theta_mean(i) !, &
                                !  & phi_phi(i),theta_theta(i),theta_mean(i), & 
                                !  & phim(i,isp),mean_phi(i)
              
          enddo
         enddo
          ! write(*,*)'theta_u_flucdone',sum(theta_u_fluc(1:nx))/nx/ufmean_des(1),"mean_pthf",(theta_u(1)-theta_u(nx))/doml(2),doml(2),"mean_phiu",sum(theta_u(1:nx)),(mean_phi_v(1)-mean_phi_v(nx)),"alpha",&
          !      &sum(theta_u(1:nx))/(mean_phi_v(1)-mean_phi_v(nx)),"umean",umean(1),umean_vel,umean(1)/(1-maxvolfrac),ufmean_des(1), &
          !     & "mean_theta",sum(theta_mean(1:nx))/nx
      
       !   axial =zero
       !   do i =1,nx-1
       !    axial=axial+(axial_con(i+1)-axial_con(i))/dx
       !    write(*,*)(axial_con(i+1)-axial_con(i))/dx

       !   enddo 
           
        !   write(*,*)'fluc_u_theta', t, sum(theta_u_fluc(1:nx))/nx

21  FORMAT(10(1xe17.4))

   close(200)
      count_theta = 80
       CALL calc_diffn_2

        call compute_local_Nuss_x

        OPEN(100,FILE=TRIM(RUN_NAME)//'_meanphi_time.dat')
      if(mod(iglobstep,20).eq.0) then
       write(100,*)'ZONE T= "', t, ' " '
       do i =1,nx,3
       
        write(100,*) (i-1)*dx, mean_phi(i),theta_u_fluc(i)/umean(1) 
       enddo
      endif 
      
    count_theta = 80
    nrbins = NINT(real(my,prcn)/1.4d0)

  ! 	CALL calc_diffn_2

 ! 	call compute_local_Nuss_x
!stop    
#if 0
    !!!!!! pdf
      call calc_generate
 !       call calc_array
        call calc_pdf
   !   
          
           call calc_generate_vel
     !!!     call calc_array_vel
           call calc_pdf_vel
#endif
 ! stop 
    !!!!!!
    count_routine  = count_routine + 1
    IF(count_routine.eq.3) count_routine = 1
    
    if(I_AM_NODE_ZERO) WRITE(*,*) 'COUNT IN SCALAR = ', COUNT_ROUTINE
    !READ(*,*)
    if(I_AM_NODE_ZERO)then
#if PARALLEL
       write (filename1,fmt="('NODE',i2.2,'_',i1)") myid,count_routine
       FILENAMELOC = ""
#else
       WRITE(FILENAME1, '(I1)')count_routine 	
#endif
    end if
    if(I_AM_NODE_ZERO)then
       SCALANDU = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME1)//'.dat&
            &'
       SCALONLY = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME1)//'.dat&
            &'
    else
       SCALANDU = ""
       SCALONLY = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (filename2,fmt="('NODE',i2.2,'_',i1)") node&
               &,count_routine
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME2)//'.dat&
               &'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME2)//'.dat&
               &'

          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(SCALANDU,strlen,node_zero,0,1,decomp_group,status)
       RECV_STRING(SCALONLY,strlen,node_zero,0,1,decomp_group,status)
    end if

        
    INQUIRE(FILE=SCALANDU,EXIST=filexist,OPENED=isopen)
    IF (.NOT.filexist) THEN
       stat = "new"
    ELSEIF(filexist.AND..NOT.isopen) THEN
       stat ="replace"
    ENDIF
    call write_mid_plane (SCALANDU,SCALONLY,stat)
 !!   call write_inlet_plane (SCALANDU,SCALONLY,stat)
    
    if(I_AM_NODE_ZERO)then
       SCALANDU_FOR = TRIM(RUN_NAME)//'_SCALANDU_FOR_'&
            &//TRIM(FILENAME1)//'.dat'    
    else
       SCALANDU_FOR = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (filename2,fmt="('NODE',i2.2,'_',i1)") node&
               &,count_routine
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_FOR_'//TRIM(FILENAME2)//'.dat&
               &'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(SCALANDU_FOR,strlen,node_zero,0,1,decomp_group,status)
    end if
    
    INQUIRE(FILE=SCALANDU_FOR,EXIST=filexist,OPENED=isopen)
    IF (.NOT.filexist) THEN
       stat = "new"
    ELSEIF(filexist.AND..NOT.isopen) THEN
       stat = "replace"
    ENDIF
    call write_forcing_wake(SCALANDU_FOR,stat)
    

!!$    THETAVSNU3 = TRIM(RUN_NAME)//'_THETAVSNU3_'//TRIM(FILENAME1)//'.dat'
!!$    INQUIRE(FILE=THETAVSNU3,EXIST=filexist,OPENED=isopen)
!!$    IF (.NOT.filexist) THEN
!!$       stat = "new"
!!$    ELSEIF(filexist.AND..NOT.isopen) THEN
!!$       stat ="replace"
!!$    ENDIF
    
    !call write_nu_no(THETAVSNU3,stat)

    IF(count_routine.eq.1) count_tmp = 2
    IF(count_routine.eq.2) count_tmp = 1


    runitno = getnewunit(minunitno, maxunitno)
    junk_char = "formatted"
    if(I_AM_NODE_ZERO)then
#if PARALLEL
       write (filename1,fmt="('NODE',i2.2,'_',i1)") myid,count_tmp
       FILENAMELOC = ""
#else
       WRITE(FILENAME1, '(I1)')count_tmp 	
#endif
    end if
    if(I_AM_NODE_ZERO)then
       SCALANDU = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME1)//'.dat&
            &'
       SCALONLY = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME1)//'.dat&
            &'
       SCALANDU_FOR = TRIM(RUN_NAME)//'_SCALANDU_FOR_'//TRIM(FILENAME1)//'.dat&
            &'
    else
       SCALANDU = ""
       SCALONLY = ""
       SCALANDU_FOR = ""
    end if
    if(I_AM_NODE_ZERO)then
       do node=1,nproc-1
          write (filename2,fmt="('NODE',i2.2,'_',i1)") node&
               &,count_tmp
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_'//TRIM(FILENAME2)//'.dat&
               &'
          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALONLY_'//TRIM(FILENAME2)//'.dat&
               &'

          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
          FILENAMELOC = TRIM(RUN_NAME)//'_SCALANDU_FOR__'//TRIM(FILENAME2)//'.dat&
               &'

          SEND_STRING(filenameloc,strlen,node,0,1,decomp_group)
       end do
    else  
       RECV_STRING(SCALANDU,strlen,node_zero,0,1,decomp_group,status)
       RECV_STRING(SCALONLY,strlen,node_zero,0,1,decomp_group,status)
       RECV_STRING(SCALANDU_FOR,strlen,node_zero,0,1,decomp_group,status)
    end if
    INQUIRE(FILE=SCALANDU,EXIST=filexist,OPENED=isopen)
    IF(filexist) THEN 
       CALL delete_file(SCALANDU_FOR,junk_char, runitno)
       CALL delete_file(SCALANDU,junk_char, runitno)
       CALL delete_file(SCALONLY,junk_char, runitno)
!!       CALL delete_file(THETAVSNU3,junk_char, runitno)
    end IF

  end subroutine output_scal
  
  subroutine write_mid_plane(scalandu,scalonly,stat)

    implicit none 
    Integer :: i,j,k, ii, isp
    Real(prcn) :: dist
    INTEGER, SAVE  :: count_routine=0
    CHARACTER*80 :: scalandu,scalonly,stat
    LOGICAL :: filexist, isopen
    INTEGER :: count_tmp , scalunit, scalnuunit
    
    
    scalnuunit = getnewunit(minunitno,maxunitno)
    open(scalnuunit,FILE=scalandu,form = "formatted",status=stat)
    
    WRITE(scalnuunit, '("# generated at time = ",g12.5)') t    

    write(scalnuunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "UX" ',   &
         &    ' "UY" ',' "UZ" ' ,' "phi"', '"temp" '
    write(scalnuunit,*)'ZONE F=POINT, I=', nx+1,  ', J=', my
    
    k = mz/2
    do j=1,my
       do i=1,nx+1
          write(scalnuunit,'(10(3x,g12.5))')real(GLOBAL_INDEX(i)-1)*dx,real(j-1)*dx,velr(i,j,k,1), velr(i,j,k,2),&
               & velr(i,j,k,3), phir(i,j,k,1),phir(i,j,k,1)*(-100.d0)+400.d0 !(phir(i,j,k,1)-phistream)/(phisurf&
          ! &-phistream) 
          
       enddo
    enddo

    close(scalnuunit,status='keep')
  end subroutine write_mid_plane
  
  subroutine write_inlet_plane(scalandu,scalonly,stat)

    implicit none
    Integer :: i,j,k, ii, isp
    Real(prcn) :: dist
    INTEGER, SAVE  :: count_routine=0
    CHARACTER*80 :: scalandu,scalonly,stat
    LOGICAL :: filexist, isopen
    INTEGER :: count_tmp , scalunit, scalnuunit


    scalnuunit = getnewunit(minunitno,maxunitno)
    open(scalnuunit,FILE=scalandu,form = "formatted",status=stat)

    WRITE(scalnuunit, '("# generated at time = ",g12.5)') t

    write(scalnuunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "UX" ',   &
         &    ' "UY" ',' "UZ" ' ,' "phi"', '"temp" '
    write(scalnuunit,*)'ZONE F=POINT, I=', mz,  ', J=', my

    i = 1
    do j=1,my
       do k=1,mz
          write(scalnuunit,'(10(3x,g12.5))')k,j,velr(i,j,k,1), velr(i,j,k,2),&
               & velr(i,j,k,3), phir(i,j,k,1),phir(i,j,k,1)*(-100.d0)+400.d0
!(phir(i,j,k,1)-phistream)/(phisurf&
          ! &-phistream)

       enddo
    enddo

    close(scalnuunit,status='keep')
  end subroutine write_inlet_plane


  subroutine write_forcing_wake(scalandu_for,stat)
    USE errormesgs
    USE general_funcs
    USE nlmainarrays
    USE nl_allphi
    implicit none 
    INTEGER :: i,j,k, unitno,unitno1,unitno2,unitno3, isp, ii, jj, kk, refi,iii
    INTEGER, SAVE  :: count_routine=0
    CHARACTER*80 :: scalandu_for,stat
    REAL(prcn) :: perfac_scal, phimi, phimo, umi, perfacscal,phim_total(mx),phim_total_z(mx),Um_total_z(mx),Um_total(mx)
    character*20 ch
    integer ich    

    unitno = getnewunit(minunitno, maxunitno)
   
unitno2 = getnewunit(minunitno, maxunitno)
unitno3 = getnewunit(minunitno, maxunitno)

    
    CALL compute_phim_heat_ratio(.TRUE.)!(velr,phir)

    !!!Tm

    phim_total = zero
    Um_total =zero
    phim_total_z = zero
    Um_total_z =zero


     if(nproc-1.eq.myid) then
      do i=1,nx+1
       phim_total(i+myid*nx)=phim(i,1)
       Um_total(i+myid*nx)=um(i)
      end do

     else
       do i=1,nx
         phim_total(i+myid*nx)=phim(i,1)
         Um_total(i+myid*nx)=um(i)
       end do

     end if

   do i=1,mx
    GLOBAL_DOUBLE_SUM(phim_total(i),phim_total_z(i),1,decomp_group)
    GLOBAL_DOUBLE_SUM(Um_total(i),Um_total_z(i),1,decomp_group)
  end do

    if(I_AM_NODE_ZERO)then
      
       OPEN(unitno,FILE=TRIM(RUN_NAME)//'_Tm_and_Um.dat', form="formatted",status='unknown')
     
       
        do i = 1,mx
          ii = i !+ mx1/2
          refi = 1
          perfac_scal = one
#if 0
          if(TRIM(RUN_NAME).eq.'SIMPLE')then
             ii = i
          else
             refi = 1+mx1/2
          end if

          if(ii.gt.mx) then
             ii = ii-mx+1
             perfac_scal = perfacscal
          end if
#endif
          Write(unitno,21)(GLOBAL_INDEX(i)-1)*dx, phim_total_z(i),Um_total_z(i)
            !perfac_scal!/phim(refi,1)
       end do
       close(unitno,status='keep')
      
    end if
   

    Write(*,*)'phim reference: ', phim(refi,1)
    OPEN(unitno,FILE=scalandu_for, form="formatted",status=stat)
    
    !WRITE(unitno, '("# generated at time = ",g12.5)') t
    
    write(unitno,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "phi_n" ' !,' &
        ! &"UX" ',' "UY" ',' "UZ" ' 
    write(unitno,*)'ZONE F=POINT, I=', nx+1,  ', J=', my, ', K=', mz
    
    if(write_output) then
       do k=1,mz
          do j=1,my
             do i=1,nx+1
                ii = i !+ mx1/2
                perfac_scal = one
                refi = i
#if 0
                if(TRIM(RUN_NAME).eq.'SIMPLE') then
                   ii = i
                else
                   refi = 1+mx1/2
                end if

                if(ii.gt.mx)then
                   ii = ii-mx+1
                   perfac_scal = perfacscal
                end if
#endif
              if(fluid_atijk(i,j,k)) then
                write(unitno,21) REAL(GLOBAL_INDEX(i)),REAL(j),REAL(k), velr(ii,j,k,1)
                   !  &,velr(ii,j,k,2),velr(ii,j,k,3) 
              else 
                write(unitno,21) REAL(GLOBAL_INDEX(i)),REAL(j),REAL(k), 0 !phir(ii,j&
                    ! &,k,1),0 &
                    ! &,0,0
 
              endif
!!$                write(unitno,21) REAL(GLOBAL_INDEX(i)),REAL(j)
!!$                ,REAL(k),velr(ii,j,k,1) ,velr(ii,j,k,2),velr(ii,j,k
!!$                ,3) 
             enddo
          enddo
       enddo
    end if

21  FORMAT(10(1xe17.4))
    CLOSE(unitno,status= 'keep')

  end subroutine write_forcing_wake
  
  subroutine write_nu_no(THETAVSNU3,stat)
    USE errormesgs
    USE general_funcs
    IMPLICIT NONE 
    INTEGER :: i,j,k,m,l,isp, unitno, iphs
    INTEGER, SAVE  :: count_routine=0
    CHARACTER*80 :: thetavsnu3,stat
    LOGICAL :: filexist, isopen

    unitno = getnewunit(minunitno, maxunitno)
    OPEN(unitno,FILE=thetavsnu3, form="formatted",status=stat)
    
!!$    
!!$    IF (unitno.LT.0) CALL printerror("newunit","ounit")
!!$
!!$    OPEN(unit=unitno,file='thetavsNU3.dat',form='formatted',status='re&
!!$         &place')

    DO m=1,nbody
       iphs = part_array(m)%iphs
       nbnd = phase_array(iphs)%nbnd
       nrpr = phase_array(iphs)%nrpr
       
       bndarray => phase_array(iphs)%bndpts
       
       WRITE(unitno,*)'Zone'
       DO isp = 1,nspmx
          DO l=1,nrpr
             IF (bndarray(2,l).GE.zero.AND.bndarray(3,l).EQ.zero) THEN 
                IF(bndarray(1,l).GE.zero) THEN 
                   WRITE(unitno,21)180.-((180.*ATAN(bndarray(2,l)/bndarray(1,l)))/pi)&
                        &,-Nu3(m,l,isp)*two*radbdy(m)*dx
                ELSE
                   WRITE(unitno,21)((180. *ATAN(-bndarray(2,l)/bndarray(1,l)))/pi)&
                        &,-Nu3(m,l,isp)*two*radbdy(m)*dx
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    END DO
21  FORMAT(10(1xe17.4))
    CLOSE(unitno,status= 'keep')

  end subroutine write_nu_no

  SUBROUTINE flow_snapshot
    Use nlarrays , only : ur1,uf1
    USE dem_mod, only : is_mobile, des_pos_new, des_radius
    IMPLICIT NONE
    Integer  :: sunit,i,j,k,m,isp, mark, idim
    INTEGER, SAVE :: zone_count = 0
    LOGICAL, SAVE :: first_time=.TRUE.
    REAL(prcn) :: ucg, fluct_vel(ndim), fluct_force(ndim),&
         & mean_vel(ndim), mean_force(ndim), position(ndim)
    CHARaCTER*80 :: FILENAME 
    CHARACTER(LEN=80) :: formfile

    formfile='formatted' 
    
    sunit  = getnewunit(minunitno,maxunitno)
    !if(irestart.eq.1)first_time = .FALSE.
    !CALL calc_velreal(velr)
    if(iscalon.eq.1)then
       do isp=1,nspmx
          do i=1,nx !mx,1
             do j=1,my2,1
                do k=1,mz,1
                   uf1(j,k) = phif(i,j,k,isp)
                enddo
             enddo
             
             call ff2cr(uf1,ur1)	 
             do j=1,my,1
                do k=1,mz,1
                   phir(i,j,k,isp) = ur1(j,k)+phirmean(isp)
                enddo
             enddo
          enddo
       enddo
    else
       do isp=1,nspmx
          do i=1,nx !mx,1
             do j=1,my,1
                do k=1,mz,1
                   phir(i,j,k,isp) = zero
                end do
             end do
          end do
       end do
    end if

    IF(first_time)THEN
       OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted',      &
            &       status='unknown')
       !write(sunit,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "UX" ',' "UY" ',' "UZ" '!,' "PHI" '
       write(sunit,*)'ZONE T = "', t, '",'
       write(sunit,*)'DATAPACKING=POINT, I =', nx,  ', J=', my, ', K=', 3    
       first_time = .FALSE.
       k = 1
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz/2
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       FILENAME = TRIM(RUN_NAME)//'_sphr_motion.dat'
       CALL  RUN_TIME_FILE_OPENER(sphrunit,FILENAME, formfile)
    ELSE
       OPEN(unit = sunit,file=TRIM(RUN_NAME)//'_MOVIE.dat',form='formatted',      &
            &       POSITION='append')
       write(sunit,*)'ZONE T = "', t, '",'
       write(sunit,*)'DATAPACKING=POINT, I =', nx,  ', J=', my, ', K=', 3
       !write(sunit,*)', VARSHARELIST= ([1-3]=1)'
       k = 1
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz/2
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
       k = mz
       do j=1,my
          do i=1,nx
             write(sunit,*)(i),(j),(k),(velr(i,j,k,1)+frame_vel(1)),(velr(i,j,k,2)+frame_vel(2)),(velr(i,j,k,3)+frame_vel(3))!, phir(i,j,k,1)
          enddo
       enddo
       
    END IF

    close(sunit,status='keep')
!!$    DEALLOCATE(velr)
    
    WRITE(sphrunit,*)'ZONE T= "', t, ' " '
!!$    ucg = DSQRT((rhos/rhof - one)*9.8*dia_phys)
    do idim = 1, ndim
       mean_force(idim) = SUM(force(1:nbody,idim))/real(nbody,prcn)
       mean_vel(idim) = SUM(velbdy(1:nbody,idim))/real(nbody,prcn)
    end do

    DO m=1,nbody
!!$       WRITE(sphrunit,'(4(2x,f12.8))')  t*ucg/dia_phys,velbdy(m,1:ndim)/ucg
!!$       WRITE(*,'(4(2x,f12.8))')  t*ucg/dia_phys,velbdy(m,1:ndim)/ucg

       fluct_vel(1:ndim) = velbdy(m,1:ndim)-mean_vel(1:ndim)
       fluct_force(1:ndim) = force(m,1:ndim)-mean_force(1:ndim)
       mark = -1
       if(DOT_PRODUCT(fluct_force(1:ndim),fluct_vel(1:ndim)).gt.zero) mark = 1

       do idim = 1, ndim
          position(idim) = XC(m,idim)+frame_pos(idim)
          if(position(idim).lt.one) position(idim) = position(idim)+&
               & real(my,prcn)
          if(position(idim).ge.real(my+1,prcn))position(idim) =&
               & position(idim)- real(my,prcn)
       end do

       WRITE(sphrunit,'(10(2x,f12.8))')  position(1), position(2), position(3), radbdy(m),radbdy(m),real(mark)
    enddo
!!$    close(sphrunit,status='keep')
    !Write(*,*) 'WRITITNG SNAPSHOT FILE, ZONE COUNT = ', zone_count
  END SUBROUTINE flow_snapshot

  subroutine calc_pdf
                
    implicit none
       
       real(prcn), DIMENSION(:),allocatable ::  hist, radbin, tmp_array,tmp_u, array,arraycopy
       real(prcn) :: u_min, u_max, phirmax,phi_model,kmean
       integer ::  idim, i, j, k, count, unit1,imax,jmax,kmax,ii,iii
       integer :: nbins   , countfl_plane ,count_ne
       character*50 filename1, filename2
    
       call screen_separator(30,'^')
                
    write (*,*) "IN pdf"
    
    nbins = 500       
     
     !   filename2 = trim(run_name)//"_pdf.dat"
     !   open (2, file=trim(filename2), status="replace")

                allocate(hist(nbins), radbin(nbins))
                allocate(tmp_array(count_fluid),tmp_u(count_fluid),arraycopy(count_fluid))
      
       open(2,file='MIS1_array1_scal.dat')
          read(2,*) count_ne
             allocate(array(count_ne))
    
         do i=1,count_ne


            read(2,*)array(i)

         enddo

       close(2,status='keep')

          call unequal_histogram(array(1:count_ne),count_ne,nbins,radbin(1:nbins),hist(1:nbins))

            !!  call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))

         filename2 = trim(run_name)//"_pdf_scal_1.dat"
         open (2, file=trim(filename2), status="replace")
           do iii=1, nbins
                  ! if (hist(iii)>=small_number) then
                      write(2,"(2d15.7)")radbin(iii),hist(iii)  
                  ! else
                  !    hist(iii) = zero
                  ! endif 
           enddo
          deallocate(array)
         close(2) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               open(2,file='MIS1_array2_scal.dat')
          read(2,*) count_ne
             allocate(array(count_ne))
    
         do i=1,count_ne


            read(2,*)array(i)

         enddo

       close(2,status='keep')

          call unequal_histogram(array(1:count_ne),count_ne,nbins,radbin(1:nbins),hist(1:nbins))

            !!  call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))

         filename2 = trim(run_name)//"_pdf_scal_2.dat"
         open (2, file=trim(filename2), status="replace")
           do iii=1, nbins
                  ! if (hist(iii)>=small_number) then
                      write(2,"(2d15.7)")radbin(iii),hist(iii)  
                  ! else
                  !    hist(iii) = zero
                  ! endif 
           enddo

          deallocate(array)
         close(2) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       open(2,file='MIS1_array3_scal.dat')
          read(2,*) count_ne
             allocate(array(count_ne))
    
         do i=1,count_ne


            read(2,*)array(i)

         enddo

       close(2,status='keep')

          call unequal_histogram(array(1:count_ne),count_ne,nbins,radbin(1:nbins),hist(1:nbins))

            !!  call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))

         filename2 = trim(run_name)//"_pdf_scal_3.dat"
         open (2, file=trim(filename2), status="replace")
           do iii=1, nbins
                  ! if (hist(iii)>=small_number) then
                      write(2,"(2d15.7)")radbin(iii),hist(iii)  
                  ! else
                  !    hist(iii) = zero
                  ! endif 
           enddo
         deallocate(array)
         close(2) 


#if 0
!! ---------------------------------------
            phi_fluid_mean(1) = 0.32
            write(*,*)'Re,PHIAVG,phi_fluid_mean(1)=',Re, PHIAVG,phi_fluid_mean(1)
              phi_model = 4.50d0*PHIAVG/(Re+2.718/2.d0)/0.70d0*10.4 !!6.84 !!change

           
                count = 0
        um = zero           
       phim = zero               
       countfl_plane = 0
       
          PHIrmax = zero
           kmean = zero
                u_min = minval(radbin(:))
                u_max = maxval(radbin(:))
          
            do i=1,nx
             do j = 2,nbins-1
                    phirmax = exp(-phi_model*(i-1)*dx)
                 kmean = kmean + (2.8E+19)*exp(-242400.0d0/8.314/(-10.0*phirmax*radbin(j)+800.0))* & 
                            phi_model*hist(j)*dx*(radbin(j+1)-radbin(j-1))/2.0d0
       !  if(i==1)write(*,*)phirmax,radbin(j),hist(j),exp(-242400.0d0/8.314/(-10.0*phirmax*radbin(j)+800.0)), & 
!!(2.8E+19)*exp(-242400.0d0/8.314/(-10.0*phirmax*radbin(j)+800.0))*phi_model*hist(j)*dx*(radbin(j+1)-radbin(j-1))/2.0d0
 
             enddo
            enddo 
             write(*,*)'kpdf=',kmean,'kmean=',(2.8E+19)*exp(-242400/8.314/(-10.0*0.32+800.0))
#endif
#if 0       
                        do i=1, nx,nx/3
                           if(i.eq.1)filename2 = trim(run_name)//"_pdf_1.dat"
                           if(i.eq.1+nx/3)filename2 = trim(run_name)//"_pdf_2.dat"
                           if(i.eq.1+2.0d0*nx/3)filename2 = trim(run_name)//"_pdf_3.dat"
                            open (2, file=trim(filename2), status="replace")
                              countfl_plane = 0
                              count_ne = 0
                              phirmax = zero
                              count = 0
                           write(*,*)'i=',i,1+2.0d0*nx/3
                        do j=1, my
                        do k=1, mz
                            if (fluid_atijk(i,j,k)) then
                                     count = count+1
                                 if(phir(i,j,k,1).gt.phirmax) then 
                                   phirmax= phir(i,j,k,1)
                                     imax= i 
                                     jmax= j
                                    kmax =k
                                 endif
                                 tmp_u(count)= velr(i,j,k,1)   
                                 tmp_array(count) =phir(i,j,k,1)
                         if(phir(i,j,k,1).ge.0.0d0) then 
                           count_ne = count_ne+1
                           arraycopy(count_ne) =phir(i,j,k,1)
                         endif
                         phim(i,1) = phim(i,1) + velr(i,j,k,1)*phir(i,j,k,1)
                      
                         um(i) = um(i) + velr(i,j,k,1)
             !!           countfl_plane(i) = countfl_plane(i) + 1
                 
         
           !!  phim(i,isp) = phim(i,isp)/(um(i)*real(countfl_plane(i),prcn))
         
                !                 write (2,"(1e15.7)") tmp_array(count)   
                         if(count.eq.1) write(*,*)'phir,tmp',phir(i,j,k,1),tmp_array(count) 
                         if(mod(count,50000).eq.0)     write(*,*) count, count_fluid, tmp_array(count)
                 
                          endif
                        enddo
                        enddo
                
                         ! phim(i,1) = phim(i,1)/um(i)
                           allocate(array(count_ne))
                           do ii= 1,count_ne
                             array(ii)= arraycopy(ii)
                           enddo
                    call unequal_histogram(array(1:count_ne), count_ne,nbins,radbin(1:nbins),hist(1:nbins))
   
                      !    call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))
                          
                            do iii=1, nbins
                            if (hist(iii)>small_number) write(2,"(2d15.7)") radbin(iii),hist(iii)   
                            enddo
       

                          deallocate(array)
                           close(2) 
                        enddo
#endif               
                  
             write(*,*)'i,j,k max',imax,jmax,kmax,phir(imax,jmax,kmax,1)
           !  write(*,*)'u_max and min=',minval(tmp_u(:)), maxval(tmp_u(:))        
             write(*,*)'count',count,'nbins',nbins, 'count_ne',count_ne
           !  write(*,*)'phim',phim(:,1)
 ! stop

#if 0                                        
               u_min = minval(tmp_array(:))
                u_max = maxval(tmp_array(:))
                 do i = 1,count
                    tmp_array(i)= (tmp_array(i)-u_min)/ (u_max- u_min)
                 enddo                             
                  do i= 1,count_ne
                      array(i)= arraycopy(i)
                  enddo     
                  u_min = minval(array(:))
                  u_max = maxval(array(:))
                 !  do i= 1,count_ne
                     ! array(i)= (array(i)-u_min)/ (u_max- u_min)
                 ! enddo 
      !    call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))
   
       !   call make_histogram(tmp_array(1:count), count, nbins,radbin(1:nbins), hist(1:nbins))
                              
    
                !  do i=1, nbins
                !       if (hist(i)>small_number) write(2,"(2d15.7)") radbin(i), hist(i)
               !   enddo
                        
#endif                     
                 deallocate(tmp_array,tmp_u)
                  
   end subroutine calc_pdf

subroutine calc_generate

    implicit none

       real(prcn), DIMENSION(:),allocatable ::  hist, radbin, tmp_array,tmp_u, array,arraycopy
       real(prcn) :: u_min, u_max, phirmax
       integer ::  idim, i, j, k, count, unit1,imax,jmax,kmax,ii,iii
       integer :: nbins   , countfl_plane ,count_ne
       character*50 filename1, filename2

       call screen_separator(30,'^')

    write (*,*) "IN pdf"

    nbins = 200

     !   filename2 = trim(run_name)//"_pdf.dat"
     !   open (2, file=trim(filename2), status="replace")

                allocate(hist(nbins), radbin(nbins))
                allocate(tmp_array(count_fluid),tmp_u(count_fluid),arraycopy(count_fluid))

             radbin = 0.0d0
                count = 0
        um = zero
       phim = zero
       countfl_plane = 0
       count_ne = 0
          phirmax = zero
                         filename2 = trim(run_name)//"_array1_scal.dat"
            
                        
                        !  filename2 = trim(run_name)//"_array2_scal.dat"

                         ! filename2 =trim(run_name)//"_array3_scal.dat"

                            open (2, file=trim(filename2), status="replace")
                              countfl_plane = 0
                              count_ne = 0
                              phirmax = zero
                              count = 0
                           
                        do i=1,nx
                        do j=1, my
                        do k=1, mz
                           if (fluid_atijk(i,j,k)) then
                              phim(i,1) = phim(i,1) + velr(i,j,k,1)*phir(i,j,k,1)

                              um(i) = um(i) + velr(i,j,k,1)
             !!               countfl_plane(i) = countfl_plane(i) + 1
                           endif
                        enddo
                        enddo

                           phim(i,1) = phim(i,1)/(um(i))  !! *real(countfl_plane(i),prcn))
                        enddo 
                         write(*,*) phim(:,1)

                        i=1
                        do j=1, my
                        do k=1, mz 
                            if (fluid_atijk(i,j,k)) then
                                     count = count+1
                                 if(phir(i,j,k,1).gt.phirmax) then
                                   phirmax= phir(i,j,k,1)
                                     imax= i
                                     jmax= j
                                    kmax =k
                                 endif
                                 tmp_u(count)= velr(i,j,k,1)
                                 tmp_array(count) =phir(i,j,k,1)
                        ! if(phir(i,j,k,1).ge.0.0d0) then
                        !   count_ne = count_ne+1
                         !  arraycopy(count_ne) =phir(i,j,k,1)/phim(i,1)
                        ! endif
                         

                !                 write (2,"(1e15.7)") tmp_array(count)
                           

                          endif
                       ! enddo
                        enddo
                        enddo
                       
                          ! allocate(array(count_ne))
                          ! do ii= 1,count_ne
                          !   array(ii)= arraycopy(ii)
                          ! enddo
                            
                            write(2,*)count
                            
                            do iii=1, count
                               write(2,*) tmp_array(iii)
                            enddo
			
				close(2)
                         
			  filename2 =trim(run_name)//"_array2_scal.dat"

                            open (2, file=trim(filename2), status="replace")
							count = 0
                        i=75
                        do j=1, my
                        do k=1, mz 
                            if (fluid_atijk(i,j,k)) then
                                     count = count+1
                                 
                                 tmp_u(count)= velr(i,j,k,1)
                                 tmp_array(count) =phir(i,j,k,1) 
                             endif  
                        enddo
						enddo
						 write(2,*)count
						  do iii=1, count
                               write(2,*) tmp_array(iii)
                          enddo
                        close(2)

			  filename2 =trim(run_name)//"_array3_scal.dat"

                            open (2, file=trim(filename2), status="replace")
			count = 0
                        i=150
                        do j=1, my
                        do k=1, mz 
                            if (fluid_atijk(i,j,k)) then
                                     count = count+1
                                 
                                 tmp_u(count)= velr(i,j,k,1)
                                 tmp_array(count) =phir(i,j,k,1) 
                             endif  
                        enddo
						enddo
			 write(2,*)count
			  do iii=1, count
                               write(2,*) tmp_array(iii)
                          enddo
                        close(2)



                         !!! deallocate(array)
                          
                       

             write(*,*)'i,j,k max',imax
    
         deallocate(tmp_array,tmp_u,arraycopy) 

   end subroutine calc_generate


   subroutine make_histogram(array, n, nbins, radbins, hist)
          implicit none
                integer, intent(in) :: n, nbins
                real(prcn), intent(inout) :: array(n)
                real(prcn), intent(inout) :: radbins(nbins), hist(nbins)
                real(prcn) :: mean, var, sd, left, right, dr
                integer :: i, ibin
                
              call calc_avr_var(n, array, mean, var)
                 sd = sqrt(var/n)
               ! array(:) = (array(:) - mean) / sd
                write (*,*)"mean,sd", mean, sd
           
                left = minval(array(:))
                right = maxval(array(:))
                write (*,*) "LEFT, RIGHT = ", left, right,nbins,hist(1:2),n
                dr = (right-left) / nbins
       
                do i=1, nbins
                        radbins(i) = left + (i-0.5d0)*dr
                enddo
   
                hist = zero

                do i=1, n
                        ibin = (array(i)-left) / dr + 1
                        if (ibin>nbins) ibin = nbins
                        hist(ibin) = hist(ibin) + 1
                enddO
                hist(:) = hist(:) / n / dr
                write (*,*) "SUM OF HIST = ", sum(hist(:)) * dr
                

   end subroutine make_histogram
   
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

  subroutine calc_array ! (nvar, var, avr, variance)
                implicit none
                integer :: nvar
                real(prcn), allocatable  :: var(:)
               !! real(prcn) :: avr, variance
                integer :: ivar,i,avr
     
      !!!!
         open(2,file='MIS1_array1.dat')   
          read(2,*) nvar
           avr = nvar
          allocate(var(150*150*150*5))
        do i=1,nvar

                   
          read(2,*)var(i)
   
        enddo
         
        
    close(2)

          open(2,file='MIS2_array1.dat')   
          read(2,*) nvar
           avr= avr+nvar
         
          do i=1,nvar
             

          read(2,*)var(i+avr-nvar)

         enddo
  
         close(2)
                open(2,file='MIS3_array1.dat')
          read(2,*) nvar
           avr= avr+nvar
          
          do i=1,nvar


          read(2,*)var(i+avr-nvar)

         enddo
         
         close(2)

        open(2,file='MIS4_array1.dat')

          read(2,*) nvar
           avr= avr+nvar
          
          do i=1,nvar


          read(2,*)var(i+avr-nvar)

         enddo
         
         close(2)

         open(2,file='MIS5_array1.dat')
          read(2,*) nvar
           avr= avr+nvar
          
          do i=1,nvar


          read(2,*)var(i+avr-nvar)

         enddo
         
         close(2)

               write(*,*) 'avr',avr,150*150*5

        open(2,file='MIS1_pdfarray_1.dat')
        
          write(2,*) avr
          
          do i=1,avr


          write(2,*)var(i)

         enddo
         
         close(2)

        end subroutine calc_array
   
   subroutine unequal_histogram(array, n, nbins, radbins, hist)
          
         implicit none
                integer, intent(in) :: n, nbins
                real(prcn), intent(inout) :: array(n)
                real(prcn), intent(inout) :: radbins(nbins), hist(nbins)
                real(prcn) :: mean, var, sd, left, right, dr
                integer :: i, ibin
                integer :: j,k
                real(prcn) :: b,fac,p,x(nbins+1),y(nbins+1),tot,delx,uneq_xdiff(nbins),pist(nbins),ratio(nbins),bal_l,bal_r
            
            !!   call calc_avr_var(n, array, mean, var)
            !!     sd = sqrt(var/n)
             !!   array(:) = (array(:) - mean) / sd
             !!   write (*,*)"mean,sd", mean, sd
       
         
                left = minval(array(:))
                right = maxval(array(:))
                write (*,*) "LEFT, RIGHT = ", left, right,nbins,hist(1:2),n

            ! bldd = 10^25
            ! brdd = -10^25
           !! bldd = min(min(free_norm(:)),bldd)
           !! brdd = max(max(free_norm(:)),brdd)
!!% bldd = 0
!!% brdd = max(free_norm(:))
!!%brdd = 10^-1;
  !!%%%%%%%%%%%%%%%unequal intervals
            p = 1.0d0
          !!  b = 1.0000001
!!%optimumb = 1.0001;
            b = 1.0001
          fac = (b+1)/(b-1)
          x = 0
          y = 0
          tot = right-left
          delx = p/nbins
      
        do j = 1,nbins
          x(j+1) = x(j) + delx
        enddo
        do j = 2,nbins+1
          y(j) = (p*(((b+1)-(b-1)*(fac**(1-x(j))))/(fac**(1-x(j))+1)))*tot   
        enddo
          y(:) = y(:) + left
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
     pist = 0
     uneq_xdiff = 0

    do j = 1,nbins
      uneq_xdiff(j) = y(j+1)- y(j)
    enddo

!!% xdiff = (brdd-bldd)/(M-1);
    do i = 1,n
     do j = 1,nbins
      if (j==1) then
        bal_l = array(i) - y(1)
        bal_r = array(i) - y(2)
        if((bal_l>=0).and.(bal_r<0)) then
!!% ibin = round((free_norm(i)-bldd)/uneq_xdiff(i)) + 1;
         pist(j) = pist(j) + 1.0d0/n !!;%(N*uneq_xdiff(j));
          goto 1000
        endif
      elseif((j.ne.1).and.(j.ne.nbins)) then
        bal_l = array(i) - y(j)
        bal_r = array(i) - y(j+1)
        if((bal_l>=0).and.(bal_r<0)) then
!% ibin = round((free_norm(i)-bldd)/uneq_xdiff(i)) + 1;
         pist(j) = pist(j) + 1.0d0/n  !!;%/(N*uneq_xdiff(j));
         goto 1000
        endif
      elseif(j==nbins) then
         bal_l = array(i) - y(nbins)
         bal_r = array(i) - y(nbins+1)
         if((bal_l>=0).and.(bal_r<=0)) then
!!% ibin = round((free_norm(i)-bldd)/uneq_xdiff(i)) + 1;
           pist(j) = pist(j) + 1.0d0/n  !;%/(N*uneq_xdiff(j));
           goto 1000
         endif
      endif
     enddo

1000 continue

    enddo
       ratio = 0
     do k = 1,nbins
       ratio(k) = pist(k)/uneq_xdiff(k)
     enddo

#if 0
                dr = (right-left) / nbins

                do i=1, nbins
                        radbins(i) = left + (i-0.5d0)*dr
                enddo

                hist = zero

                do i=1, n
                        ibin = (array(i)-left) / dr + 1
                        if (ibin>nbins) ibin = nbins
                        hist(ibin) = hist(ibin) + 1
                enddo

                hist(:) = hist(:) / n / dr
#endif       
              do i=1, nbins
                        radbins(i) = y(i)
              enddo

               hist = ratio
                write (*,*) "SUM OF HIST = ",uneq_xdiff !!sum(hist(:)) * dr


   end subroutine unequal_histogram

 SUBROUTINE compute_phim_heat_ratio_1(output)!(velr,phi)
       Use nlmainarrays, Only : velr=>ubcp, phi=>nlbcp
       IMPLICIT NONE
       !real(prcn), DImension(:,:,:,:), Intent(in) :: velr,phi
       LOGICAL, Intent(in) :: output
       INTEGER :: i,j,k, isp


       um = zero
       phim = zero
       countfl_plane = 0

       do i = 1, nx+1 !mx
          do k = 1, mz
             do j = 1, my
                if(output)then
                   if(fluid_atijk(i,j,k))then
                      do isp = 1, nspmx
                         phim(i,isp) = phim(i,isp) +velr(i,j,k,1)*phi(i,j,k,isp)
                      end do
                      um(i) = um(i) + velr(i,j,k,1)
                      countfl_plane(i) = countfl_plane(i) + 1
                   end if
                else
                   do isp = 1, nspmx
                      phim(i,isp) = phim(i,isp) + velr(i,j,k,1)*phi(i,j,k,isp)
                   end do
                   um(i) = um(i) + velr(i,j,k,1)
                   countfl_plane(i) = countfl_plane(i) + 1
                end if
             end do
          end do
   um(i) = um(i)/real(countfl_plane(i),prcn)

          do isp = 1, nspmx
             phim(i,isp) = phim(i,isp)/(um(i)*real(countfl_plane(i),prcn))
          end do
       end do
       if(I_AM_LAST_PROC)then
          do isp = 1, nspmx
             heat_ratio_new(isp) =one/phim(nx+1,isp)!phim(1,isp)/phim(nx+1,isp)!
          end do
       end if
       BROADCAST_DOUBLE(heat_ratio_new(1),nspmx,END_PROC,decomp_group)
       if(I_AM_NODE_ZERO)then
          Write(*,*)'HEAT RATIO@N-1 : ', heat_ratio(:), 'HEAT_RATIO@N: ',&
                     & heat_ratio_new(:)
          !Write(*,*)'umo: ', um(mx), 'umi: ', um(1)
          Write(*,*) 'phimo: ', phim(1:mx,1)!, 'phimi: ', phim(1)
          !Write(*,*)'NLMEAN REQ:', (phim(mx)-phim(1))*um(1)
          !READ(*,*)
       end if
     END SUBROUTINE compute_phim_heat_ratio_1

	 
SUBROUTINE compute_local_Nuss_x
  USE bc_scalar
  USE dependent_functions
  use scalar_data
 
   implicit none 
    !!INTEGER, Intent(in) ::imis
        Integer :: m
!    LOGICAL, SAVE :: routine_called = .false. 
    REAL(prcn) ::  dfll(nspmx), dtheta, rad_proj,xltemp
!#if 0
    REAL(prcn) ::  nlphil(nspmx),onlphil(nspmx), nx1, theta, ny1,nz1
    REAL(prcn) ::  xl(ndim), phil(nspmx),philo(nspmx),philo2(nspmx)&
         &,phili(nspmx), cphi, phiang, thetaang
    INTEGER :: sp, i, j, k, n, l 
    REAL(prcn) :: rad, tempor(nspmx), gradphibnd(3,nspmx), normal(3)

    INTEGER :: ib, ie, jb, je, kb, ke, onew, ii, jj, kk, unitno, unitno1,unitno2  
!#if 0
    Integer :: pcell(3), funit, iphi
    integer ::  is(ndim),iii(ndim),io(ndim), io2(ndim), index_phi,&
         & index_theta, itheta, count_model,model_total_count
    REAL(prcn) ::  tempt, areabdy(nbody),  &
         &  area_spec1, flux_tmp, xcor, dx_x, &
         &  nusselt_loc(nspmx),nusselt_loc_1(nspmx),  nusselt_bdy,&
         & total_nuss_theta(nbody), total_nuss_phi(nbody),  &
         & nuss_phi_avg(count_phi),&
         & nuss_theta_avg(count_theta), x_globa(count_phi),&
         & rad_proj_area(mx),Um_total_x(mx),Um_total(mx),phim_total(mx),phim_total_x(mx), &
         & right_first(nx+1),right_third(nx+1),grad_right_first(mx),grad_right_third(mx),&
         & meanvelo(0:nx+1,3),meantemp(0:nx+1), right_first_total(mx),right_third_total(mx),&
          & right_first_x(mx),right_third_x(mx), &
          & left_meantemp(mx),left_meanvelo(mx,3),grad_left_meantemp(mx),left_meantemp_x(mx),left_meanvelo_x(mx,3), &  
          &  grad_right(nx+1), grad_fluc(0:nx+1,3), &
          & nl_total(nx+1,my,mz),right_third_yz(nx+1), &
          & grad_left_1(nx+1),grad_left_2(nx+1), grad_meanu(nx+1),grad_meanphi(nx+1),grad_meanmean(nx+1), &
          & grad_meanu_x(mx),grad_meanphi_x(mx),grad_fluc_x(mx,3),grad_meanu_total(mx),grad_meanphi_total(mx),grad_fluc_total(mx,3), &
          & grad_meanmean_total(mx),grad_meanmean_x(mx), &
          & nlphir_one(0:mx,0:my,0:mz,1),Nu_num_avg,term_one,term_two,term_three, axial_con, out_plane, &
          & model_grad(mx), model_flux(mx),model_mean,model_total,model_mean_grad(mx),model_grad_uphi(mx), &
          & model_mean_grad_uphi(mx),model_uphi(mx),model_meanmean(mx),phim_mean
       
     complex(prcn) :: nl_zf(my2,mz),nl_yf(my2,mz),nl_yzf(my2,mz)
!#if 0
    REAL(prcn) :: nuss_vsxavg(mx), nuss_vsx(mx),nuss_vsxavg2(mx), nuss_vsx2(mx), xcl, xcr, nussfac, rad_proj_line(mx)
    INTEGER :: norm_phi(count_phi), norm_theta(count_theta),nbdymax,phicelltemp, pl, & 
               & count_term,count_m_in,count_total,count_node,count_point
    CHARACTER*80 :: NUVSX, budget_temp_eqn, nusselt_number_total, budget_Tm 
    LOGICAL :: phiterm
!#endif
   
    write(*,*)'MYID : ', 'In SUBROUTINE calc_Nuss_loc_x'

!!#if 0
!!!##### Tm and Um ######

    phim_total = zero
    Um_total =zero
    phim_total_x = zero
    Um_total_x =zero

     if(nproc-1.eq.myid) then
      do i=1,nx+1
       ii=i
       jj=GLOBAL_INDEX(ii)
       phim_total(jj)=phim(i,1)
       Um_total(jj)=um(i)
      
      end do

     else
       do i=1,nx
          ii=i
         jj=GLOBAL_INDEX(ii)
         phim_total(jj)=phim(i,1)
         Um_total(jj)=um(i)
       end do

     end if

!#if 0
   do i=1,mx
    GLOBAL_DOUBLE_SUM(phim_total(i),phim_total_x(i),1,decomp_group)
    GLOBAL_DOUBLE_SUM(Um_total(i),Um_total_x(i),1,decomp_group)
  end do

  !!! write(*,*)"done diff", gradphi(1:nx,1:my/2,1:10,1)

  !!!!!#### Nusselt number VS x direction #######
       
       if(I_AM_NODE_ZERO)then
       
       NUVSX = TRIM(RUN_NAME)//'_NUVSX.dat'

       unitno = getnewunit(minunitno, maxunitno)

       OPEN(unitno,FILE=NUVSX, form="formatted",status="unknown") 

       end if

    write(*,*)"theta",count_theta
    dtheta  = twopi/real(count_theta-1,prcn)
    nusselt_loc= zero
    nusselt_loc_1= zero

    nuss_vsxavg = zero
    nuss_vsx = zero
    nuss_vsxavg2 = zero
    nuss_vsx2 = zero
    rad_proj_area =zero
    rad_proj_line =zero

   
    do pl = 1, mx
     count_m_in = 0
       do m = 1, nbody
          xcl = xc(m,1)-radbdy(m)
          xcr = xc(m,1)+radbdy(m)
       !  write(*,*)"xcl",xcl,xcr,pl
        
            if((xcl.le.zero).and.(pl.ge.(xcl+mx-1))) then
		      xl(1)= pl
                     count_m_in=count_m_in + 1
	             nx1= real(xl(1)-(xc(m,1)+mx-1),prcn)/radbdy(m)
	             rad_proj = radbdy(m)**2.d0 - (xc(m,1)+mx-1-pl)**2.d0
                     rad_proj = DSQRT(rad_proj)
		 
	    else if((xcr.ge.mx).and.(pl.le.(xcr-(mx-1)))) then
		 
		      xl(1)= pl
                      count_m_in=count_m_in + 1
		      nx1= real(xl(1)-(xc(m,1)-mx+1),prcn)/radbdy(m)
	     	      rad_proj = radbdy(m)**2.d0 - (xc(m,1)-mx+1-pl)**2.d0
                      rad_proj = DSQRT(rad_proj)
          
            else if((xcl.le.pl).and.(xcr.ge.pl))then
                      xl(1)= pl !global coordinate
                      count_m_in=count_m_in + 1 

                      nx1= real(xl(1)-xc(m,1),prcn)/radbdy(m)
             
             !   cphi=(nx1)
             !   phiang = ACOS(cphi)
            !    rad_proj = sin(phiang) !*radbdy(m)
                
              rad_proj = radbdy(m)**2.d0 - (xc(m,1)-pl)**2.d0
              rad_proj = DSQRT(rad_proj)
            else 
             goto 2222
            end if 
               
               if(I_AM_NODE_ZERO)then
                    
                rad_proj_area(pl) =rad_proj_area(pl) + pi*(rad_proj*dx)**2.0d0 
                            !!pi*(rad_proj*radbdy(m)*dx)**2.d0
                rad_proj_line(pl) =rad_proj_line(pl) + 2.0d0*pi*(rad_proj*dx)
                            

               end if         
             !! Write(*,*)'PL: ', pl,' rad_proj: ', rad_proj, radbdy(1)
             do itheta = 1, count_theta 
                theta = dtheta*real(itheta-one,prcn)
                thetaang = 180.0d0*theta/pi
                ny1 = rad_proj*cos(theta)
                nz1 = rad_proj*sin(theta) 
                xl(2)=xc(m,2)+ny1 !*radbdy(m) !global coordinate
                xl(3)=xc(m,3)+nz1 !*radbdy(m) !global coordinate

                phil(:)=zero
                dfll(:) = zero
                DO n=1,ndim
                   
                   if(xl(n).lt.zero) then 
                      is(n) = int(xl(n)-1)
                   else 
                      is(n) = int(xl(n))
                   end if
                   
                ENDDO

                rad = radbdy(m)

                
                normal(1) = nx1/(rad_proj/radbdy(m))
                normal(2) = cos(theta) !ny1
                normal(3) = sin(theta) !nz1
                
                xltemp = xl(1)
         
                phicelltemp = is(1)
                if(.not.CELL_IN_PROC(phicelltemp))then
                   WEST_PERIODIC_IMAGE(is(1),phicelltemp,xl(1),xltemp)
#if PARALLEL
                   EAST_PERIODIC_IMAGE_MOD(is(1),phicelltemp,xl(1),xltemp)
#else
                   EAST_PERIODIC_IMAGE(is(1),phicelltemp,xl(1),xltemp)
#endif
                   if(.not.CELL_IN_PROC(phicelltemp))goto 2222
                end if
                xl(1) = xltemp
                is(1) = phicelltemp
#if PARALLEL
                phiterm = CELL_IN_VEL_GRID(phicelltemp)
#else
                phiterm = .TRUE.
#endif
                
               call interpolate_phidata(is,xl,ib,&
                              & ie,jb,je,kb,ke,phil,nlphil,onlphil,dfll, 1,onew,1) 
                                
                 do k = 1, onew
                   kk = kb+k-1
                if(kk.lt.1) kk = mz+kk
                   if(kk.gt.mz) kk = kk-mz 
                   
                   do j = 1, onew
                      jj = jb+j-1
                      if(jj.lt.1) jj = my+jj
                      if(jj.gt.my) jj = jj-my
                      
                      do i = 1, onew
                         DO sp = 1, nspmx
                            ii = ib+i-1
#if !PARALLEL
!!$                         if(ii.lt.1) ii = mxf+ii-1
!!$                         if(ii.gt.mxf-1) ii = ii-mxf +1
                            if((ii.lt.1).or.(ii.gt.mxf))Write(*,*)'SOME PROBLEM WI&
                                 &TH INNER REV POINTS. CHECK'
#endif
                            

                            LOCAL_INDEX(ii)
                            
                            gradphisten(i,j,k,:,sp) = gradphi(ii,jj,kk,:)
                             !!write(*,*)gradphisten(i,j,k,:,sp),weightp(i,j,k)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
                
                do sp = 1,nspmx
                   
                   do n = 1, ndim 
                      gradphibnd(n,sp) = &
                           & array_dot_product(gradphisten(1:onew,1:onew&
                           &,1:onew,n,sp),weightp(1:onew,1:onew,1:onew))
                   end do

                   nusselt_loc(sp) = array_dot_product(gradphibnd(2:3,sp),&
                         & normal(2:3)) !!!!!!!!!change 
                   nusselt_loc_1(sp)=array_dot_product(gradphibnd(1:1,sp),&
                         & normal(1:1))
                   nusselt_loc(sp) = nusselt_loc(sp) *(rad_proj*dx)*dtheta

                        !&/(phim(pl)) !!!!!!!!!!!!!
                    nusselt_loc_1(sp) = nusselt_loc_1(sp) *(rad_proj*dx)*dtheta
                   nuss_vsxavg(pl) = nuss_vsxavg(pl) +&
                        &(nusselt_loc(sp)) !/(two*pi*rad_proj*dx)
                   nuss_vsx(pl) = nuss_vsx(pl) +&
                        &(nusselt_loc_1(sp))!/(two*pi*rad_proj)

                end do
             end do
                  ! nuss_vsxavg(pl)=nuss_vsxavg(pl)/real(count_theta,prcn)
                    
          
          !!   goto 2222
          
2222      continue
       end do
           ! nuss_vsxavg(pl)= nuss_vsxavg(pl)/real(count_m_in,prcn)
    end do

#if PARALLEL

    do sp=1,nspmx
       do i = 1, mx
          GLOBAL_DOUBLE_SUM(nuss_vsx(i),nuss_vsx2(i),1,decomp_group)
          GLOBAL_DOUBLE_SUM(nuss_vsxavg(i),nuss_vsxavg2(i),1,decomp_group)
       enddo
    end do
    
    if(I_AM_NODE_ZERO)then
       nuss_vsx=nuss_vsx2
       nuss_vsxavg=nuss_vsxavg2
    endif

#endif

    if(I_AM_NODE_ZERO)then
       nussfac = (one-maxvolfrac)/(6.d0*maxvolfrac)*char_length
       
       do i = 1, mx
        
          write(unitno,21) (i-1)*dx,rad_proj_line(i), nuss_vsx(i),nuss_vsxavg(i),rad_proj_area(i) !,nuss_vsx(i)*nussfac/phim_total_x(i)/rad_proj_area(i),&
 !                          & rad_proj_area(i),gamma(1)*nuss_vsx(i)/phim_total_x(i),nuss_vsx(i)/phim_total_x(i) !,rad_proj_line(i)
       end do

  21     FORMAT(10(2xe17.4))
       CLOSE(unitno,status= 'keep')
       
    end if
  !!  nussfac = (one-maxvolfrac)/(6.d0*maxvolfrac)*char_length
    Write(*,*)'AXIALLY AVG NUSS NO: ', SUM(nuss_vsxavg(1:mx))/real(mx,prcn), SUM(nuss_vsx(1:mx))/real(mx,prcn)
    Write(*,*)'AXIALLY AVG NUSS NO: ', SUM(nuss_vsxavg(1:mx))/(real(mx,prcn))*nussfac, SUM(nuss_vsx(1:mx))/(real(mx,prcn))*nussfac
  
  !!!####   budget of the mean temperature equation in NSF ##########
   
 !    write(*,*) "nl",nlphif(10,:,:,1)


#if 0
!!#if PARALLEL
 !!   do i=1,nx
!!#else
    do i=1,mx
!!#endif
      !! write(*,*)nlphif(i,:,:,1) !! why it doesn't work
          do j=1,my2,1
             do k=1,mz,1
                uf1(j,k) = nlphif(i,j,k,1)
             enddo
          enddo

          call ff2cr(uf1,ur1)  
        ! do j=1,my,1
        !     do k=1,mz,1
              !  nlphir_one(i,j,k,1) = ur1(j,k)
        !     enddo
        !  enddo

     !! call ff2cr(nlphif(i,:,:,1),nlphir(i,:,:,1))
    end do 
#endif   
    !! if(myid==0)write(*,*) "vemphi_trybudget=",velr(1,my/2,mz/2,1),phir(1,my/2,mz/2,1)
!!#if 0 
   right_first = zero
   right_third =zero
   right_first_x = zero
   right_third_x = zero
   right_first_total = zero
   right_third_total =zero
   meantemp =zero
   meanvelo =zero
   grad_right_first= zero
   grad_right_third= zero
   left_meantemp =zero
   left_meanvelo =zero
   right_third_yz =zero
   nl_yf =czero
   nl_zf= czero
   nl_yzf= czero
   model_mean_grad=zero 
   model_grad_uphi=zero
   model_uphi =zero
   count_total= 0
 
      !!####### mean value ######   
   do i=1,nx+1
      count_term = 0 
      do j=1,my
       do k=1,mz
        if(fluid_atijk(i,j,k)) then
           count_term = count_term +1 
            model_uphi(i)= model_uphi(i)+ velr(i,j,k,1)*phir(i,j,k,1)*dx*dx
            model_mean_grad(i)= model_mean_grad(i) + gradphi(i,j,k,1)
            
          !  model_grad_uphi(i)=model_grad_uphi(i)+ nlphir_one(i,j,k,1)*dx*dx
            
             right_first(i) = right_first(i) +  gradphi(i,j,k,1)*dx*dx !! q
           
           ! right_third(i) = right_third(i) +  (nlphir(i,j,k,1)+diffn(i,j,k,2)+diffn(i,j,k,3))*dx*dx  !! change  !!!dq/d x
              right_third(i) = right_third(i) + diffn(i,j,k,1)*dx*dx  !! change  !!!dq/d x
         
        end if
       end do 
         
      end do
       ! model_uphi(i)= model_uphi(i)/count_term
         model_mean_grad(i)= model_mean_grad(i)/count_term
       ! model_grad_uphi(i)=model_grad_uphi(i)/count_term
   end do
   ! if(myid==0) write(*,*)"grad_check_in******",fluid_atijk(16,my/2,mz/2),fluid_atijk(18,my/2,mz/2), & 
    !                                         &  diffn(16,my/2,mz/2,1)/gamma(1),diffn(17,my/2,mz/2,1)/gamma(1), & 
 
!!#if 0
    
     do i=0,nx+1
       count_term = 0
     
      do j=1,my
       do k=1,mz
        if(fluid_atijk(i,j,k)) then
            do ii=1,3
            meanvelo(i,ii) = meanvelo(i,ii)+velr(i,j,k,ii)
            end do
            meantemp(i) = meantemp(i)+phir(i,j,k,1)

            count_term= count_term+1

        end if
       end do

      end do
       meanvelo(i,1:3)= meanvelo(i,1:3)/real(count_term,prcn)
       meantemp(i)= meantemp(i)/real(count_term,prcn)

     end do
    

 !!!##### fluctuate value #####    
!#if 0
   
!! #### nonlinear term split ######
 grad_right=zero
 grad_left_1 =zero
 grad_left_2 =zero
 grad_fluc=zero
 grad_meanu =zero
 grad_meanphi=zero
 
 model_flux=zero
 model_grad = zero
 model_total_count =zero
 model_total = zero
  
    do i=1,nx
      count_model =zero
      do j=1,my
       do k=1,mz
         if(fluid_atijk(i,j,k)) then
          count_model= count_model +1
          grad_left_1(i)=grad_left_1(i)+(velr(i,j,k,1)- meanvelo(i,1))*(phir(i,j,k,1)- meantemp(i))*dx*dx
         model_flux(i)= model_flux(i)+(velr(i,j,k,1)- meanvelo(i,1))*(phir(i,j,k,1)- meantemp(i))
       
         grad_fluc(i,1)= grad_fluc(i,1) +( (velr(i+1,j,k,1)- meanvelo(i+1,1))*(phir(i+1,j,k,1)- meantemp(i+1))- &   
                                        & (velr(i-1,j,k,1)- meanvelo(i-1,1))*(phir(i-1,j,k,1)- meantemp(i-1)))/(2.0d0*dx)*dx*dx  !! phi''*u'' in one point
                                     

           grad_meanu(i)= grad_meanu(i)+ (meanvelo(i+1,1)*(phir(i+1,j,k,1)- meantemp(i+1))- &
                             & meanvelo(i-1,1)*(phir(i-1,j,k,1)- meantemp(i-1)))/(2.0*dx)*dx*dx
           grad_meanphi(i)= grad_meanphi(i) + ((velr(i+1,j,k,1)- meanvelo(i+1,1))*meantemp(i+1)- &
                             & (velr(i-1,j,k,1)- meanvelo(i-1,1))*meantemp(i-1))/(2.0*dx)*dx*dx
      
        ! if(ii.eq.1)  grad_left_1(i) =grad_left_1(i) + nlphir(i,j,k,1)*dx*dx     

         end if

       end do
      end do

         grad_left_2(i) =grad_left_2(i) + ( meanvelo(i+1,1)* meantemp(i+1)- &
                                & meanvelo(i-1,1)* meantemp(i-1))/(2.0d0*dx)*(DOML(2)*DOML(2)-rad_proj_area(i)) !! mean u * mean phi
         model_flux(i)= model_flux(i)/count_model
         model_total =model_total + model_flux(i)
         model_total_count= model_total_count + count_model
     end do !! i          
         model_total = model_total/model_total_count
!!#if 0      
grad_meanu_total= zero
grad_meanphi_total= zero
grad_fluc_total= zero
grad_meanu_x= zero
grad_meanphi_x= zero
grad_fluc_x= zero
grad_meanmean_x =zero
grad_meanmean_total =zero 
!!#if PARALLEL

    if(nproc-1.eq.myid) then
      do i=1,nx+1
         ii=GLOBAL_INDEX(i)
         right_third_total(ii)=right_third(i)  !!! dq/dx
         right_first_total(ii)=right_first(i)  !! q
       !  left_meanvelo(ii,1:3)=meanvelo(i,1:3)
       !  left_meantemp(ii)=meantemp(i)        

      end do

    else
       do i=1,nx
         ii=GLOBAL_INDEX(i)
         right_third_total(ii)=right_third(i)
         right_first_total(ii)=right_first(i)

      end do

    end if
     
       do i=1,nx
         ii=GLOBAL_INDEX(i)
         grad_fluc_total(ii,1)=grad_fluc(i,1)
         grad_meanmean_total(ii)= grad_left_2(i)
         grad_meanphi_total(ii)= grad_meanphi(i)
         grad_meanu_total(ii)= grad_meanu(i)           
       end do

   do sp=1,nspmx
       do i = 1, mx
          GLOBAL_DOUBLE_SUM(right_first_total(i),right_first_x(i),1,decomp_group)
          GLOBAL_DOUBLE_SUM(right_third_total(i),right_third_x(i),1,decomp_group)
      !!    GLOBAL_DOUBLE_SUM(left_meantemp(i),left_meantemp_x(i),1,decomp_group)
              
          GLOBAL_DOUBLE_SUM(grad_meanu_total(i),grad_meanu_x(i),1,decomp_group)
          GLOBAL_DOUBLE_SUM(grad_meanphi_total(i),grad_meanphi_x(i),1,decomp_group)
          GLOBAL_DOUBLE_SUM(grad_meanmean_total(i),grad_meanmean_x(i),1,decomp_group)
          GLOBAL_DOUBLE_SUM(grad_fluc_total(i,1),grad_fluc_x(i,1),1,decomp_group)
     
       enddo
    end do

!!#endif
 
!!#if 0
!    if(I_AM_NODE_ZERO)then
    do i=2,mx-1
     
     grad_right_first(i)=(right_first_x(i+1)-right_first_x(i-1))/(2.0d0*dx)
     
    !! grad_right_third(i)=(right_third_x(i+1)-right_third_x(i-1))/(2.0d0*dx)
     
    !! grad_right_third(i)=(grad_fluc_x(i+1)-grad_fluc_x(i-1))/(2.0d0*dx)
   

    end do
       left_meantemp_x=zero 
    do i=1,nx
    do k=1,mz
    do j=1,my
        if(fluid_atijk(i,j,k)) left_meantemp_x(i) =left_meantemp_x(i) + sum(diffn(i,j,k,1:3))*dx*dx   
    end do 
    end do
    end do
     write(*,*) 'ufmean',ufmean(1),ufmean_des(1)

#if 0
      budget_temp_eqn = TRIM(RUN_NAME)//'_budget_temp_eqn.dat'

       unitno1 = getnewunit(minunitno, maxunitno)
  
       OPEN(unitno1,FILE=budget_temp_eqn, form="formatted",status="unknown") 

         do i=2,mx-1
        
      write(unitno1,41)(GLOBAL_INDEX(i)-1)*dx, & !! 1
                
        & right_first_x(i), &  !! 2
        & nuss_vsx(i)/((DOML(2)**2.0d0 - rad_proj_area(i))*Um_total_x(i)*phim_total_x(i)-gamma(1)*right_first_x(i)),  & !! 3
        & right_third_x(i)/phim_total_x(i), & !! 4
        & gamma(1)*nuss_vsxavg(i)/phim_total_x(i), &  !! 5
        & nuss_vsxavg(i)/phim_total_x(i)/rad_proj_line(i),  &  !! 6
        & right_third_x(i)/phim_total_x(i)/gamma(1), & !! diffn  7
        & grad_right_first(i)/phim_total_x(i), &      !!d/dx intege q 8 
        & (grad_right_first(i)-right_third_x(i)/gamma(1))/phim_total_x(i), &  !!9
        & (nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/phim_total_x(i), &  !! 10
          
        & (nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/phim_total_x(i)/rad_proj_line(i), & !! 11
        & abs(right_third_x(i)/gamma(1)+nuss_vsx(i))/(nuss_vsxavg(i)+nuss_vsx(i)), & !! 12
        & nuss_vsx(i)/(nuss_vsxavg(i)+nuss_vsx(i)), & !! 13  
        & abs(grad_right_first(i))/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1)), & !! 14  axial
        & (grad_right_first(i)-right_third_x(i)/gamma(1))/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1)), & !! 15
                      
        &  nuss_vsx(i)/phim_total_x(i), & !! 16
        & (nuss_vsxavg(i)+nuss_vsx(i))/phim_total_x(i), & !! 17
        & (nuss_vsxavg(i)+nuss_vsx(i))/phim_total_x(i)/rad_proj_line(i), & !! 18
        
        & -grad_fluc_x(i,1)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/gamma(1), & !19
        & -grad_meanmean_x(i)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/gamma(1), & ! 20 
        & -grad_meanphi_x(i)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/gamma(1), & ! 21
        & -grad_meanu_x(i)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/gamma(1), & ! 22
                  !! & left_meantemp_x(i)*left_meanvelo_x(i,1)/uchar(1),left_meantemp_x(i)*left_meanvelo_x(i,2)/uchar(1)  ! 
                  !!  & grad_fluc_x(i,1)/(left_meantemp_x(i)*left_meanvelo_x(i,1))
        & -grad_fluc_x(i,1)/ufmean_des(1)/(DOML(2)*DOML(2)-rad_proj_area(i)), -grad_meanmean_x(i)/ufmean_des(1), -grad_meanphi_x(i)/ufmean_des(1),-grad_meanu_x(i)/ufmean_des(1), & !23-26
        & -(grad_meanphi_x(i)+grad_meanmean_x(i)+ grad_meanu_x(i)+grad_fluc_x(i,1))/((nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))*gamma(1)), & !27 
        & grad_left_1(i)/(DOML(2)**2.0d0 - rad_proj_area(i)), meanvelo(i,1)*meantemp(i),grad_left_1(i)/(meanvelo(i,1)*meantemp(i))/(DOML(2)**2.0d0 - rad_proj_area(i)), & !!28-30
        & (nuss_vsxavg(i)*gamma(1) -right_third_x(i)+(grad_meanphi_x(i)+grad_meanmean_x(i)+ grad_meanu_x(i)+grad_fluc_x(i,1)))/( nuss_vsxavg(i)*gamma(1) - right_third_x(i))
         end do

     41  FORMAT(40(2xe17.4))

       close(unitno1,status='keep')
#endif
!#if 0
       nusselt_number_total = TRIM(RUN_NAME)//'_nusselt_number_total.dat'

       unitno2 = getnewunit(minunitno, maxunitno)

     !!  OPEN(unitno2,FILE=nusselt_number_total, form="formatted",status="unknown")
           OPEN(500,FILE=TRIM(RUN_NAME)//'_nusselt_number_time.dat', status="old", position="append")

       Nu_num_avg=zero
       phim_mean = zero
       axial_con =zero
       out_plane =zero
       count_point=zero
       model_mean_grad_uphi=zero

           Write(*,*)'AVG NUSS NO: begin'
       Nu_num_avg = zero
       out_plane = zero
       do i = 1, mx-1
        !  if(rad_proj_line(i).ne.zero) then
             count_point=count_point+1
            ! phim_mean = phim_mean + meantemp(i)*(DOML(2)*DOML(2)-rad_proj_area(i))
           !!  Nu_num_avg= Nu_num_avg+ (grad_right_first(i)-right_third_x(i)/gamma(1)+nuss_vsxavg(i))/phim(i,1)/rad_proj_line(i)
            Nu_num_avg= Nu_num_avg+ (nuss_vsx(i)+nuss_vsxavg(i))/phim(i,1)/rad_proj_line(i) 
             !! out_plane = out_plane +(grad_right_first(i)-right_third_x(i)/gamma(1)+nuss_vsxavg(i))

            ! write(unitno2,71) (GLOBAL_INDEX(i)-1)*dx, (nuss_vsxavg(i)+nuss_vsx(i))/doml(2)**2 , (grad_right_first(i)-right_third_x(i)/gamma(1))

            axial_con= axial_con + grad_right_first(i)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))
            !! out_plane = out_plane -grad_fluc_x(i,1)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/gamma(1)
        !  end if
       end do
     !!   Write(*,*)'MFIX_estimate_NU', Nu_num_avg/phim_mean/6/maxvolfrac,'volfac',maxvolfrac,'mean',meantemp(1),meantemp(2),meantemp(3),meanvelo(1,1)
!!#if 0     
    model_mean =zero
    model_meanmean =zero

       do i=2,mx-2
          model_grad(i)= (meantemp(i+1)-meantemp(i-1))/(2.0d0*dx)
          model_mean_grad_uphi(i)= (model_uphi(i+1)-model_uphi(i-1))/(2.0d0*dx)
          model_meanmean(i)= ( meanvelo(i+1,1)* meantemp(i+1)*(DOML(2)*DOML(2)-rad_proj_area(i+1))- &
                                & meanvelo(i-1,1)* meantemp(i-1)*(DOML(2)*DOML(2)-rad_proj_area(i-1)))/(2.0d0*dx)

         !! Write(unitno2,71)model_flux(i)/ufmean_des(1),meantemp(i),(GLOBAL_INDEX(i)-1)*dx,model_flux(i),model_mean_grad(i), &
                         !!  & model_grad(i),model_grad_uphi(i)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/gamma(1),-model_mean_grad_uphi(i)/(nuss_vsxavg(i)+grad_right_first(i)-right_third_x(i)/gamma(1))/gamma(1),-model_meanmean(i),-grad_left_2(i)
          model_mean= model_mean + model_flux(i)/model_grad(i) 
       end do
          model_mean= model_mean/(mx-3)
        
        !   model_mean = model_total*(DOML(2)-2*dx)/(model_grad(mx-2)-model_grad(2))
          Nu_num_avg = Nu_num_avg /real(count_point,prcn)

          Write(500,*)t,Nu_num_avg
           
            Write(*,*)'Ave Nusselt number:and axial con ',t ,Nu_num_avg, axial_con/real(count_point,prcn), out_plane/real(count_point,prcn) 
         !   Write(*,*)'model_T',model_mean,ufmean_des(1)
    

        71  FORMAT(20(2xe17.4))

     !!  close(unitno2,status='keep')

        close(500) 
      !###########
#if 0   
       budget_Tm = TRIM(RUN_NAME)//'_Tm_budget.dat'

       unitno2 = getnewunit(minunitno, maxunitno)

       OPEN(unitno2,FILE=budget_Tm, form="formatted",status="unknown")
       
      
           Write(*,*)'AVG budg NO: begin'

       do i = 2, mx-1
            
            term_one= (phim_total_x(i+1)*Um_total_x(i+1)-phim_total_x(i-1)*Um_total_x(i-1))/2.0d0/dx
            term_two= 1.0/(DOML(2)*DOML(2)-rad_proj_area(i))*(phim_total_x(i)*Um_total_x(i)-gamma(1)*(phim_total_x(i+1)-phim_total_x(i-1))/2.0d0/dx)*(-rad_proj_area(i+1)+rad_proj_area(i-1))/dx/2.0d0
            term_three = gamma(1)*(phim_total_x(i-1)-2.0d0*phim_total_x(i)+phim_total_x(i+1))/(dx2)
            Write(unitno2,71) (i-1)*dx,Um_total_x(i)*(DOML(2)*DOML(2)-rad_proj_area(i)),term_one,term_two,term_three,gamma(1)*nuss_vsxavg(i)/(DOML(2)*DOML(2)-rad_proj_area(i)),& 
                             & (term_one+term_two+term_three+gamma(1)*nuss_vsxavg(i)/(DOML(2)*DOML(2)-rad_proj_area(i)))/(term_one+term_two)
       end do
           

        91  FORMAT(20(2xe17.4))

       close(unitno2,status='keep')
#endif

   !! endif

!#endif
  
  END SUBROUTINE compute_local_Nuss_x

    SUBROUTINE calc_diffn_2 
    USE nlarrays, ONLY : ff1=>uf1, ff2=>uf2, ff3=>uf3, ff4=>uf11,fr1&
         &=>ur1,fr2=>ur2,fr3=>ur3,ff5 => uf22, ff6 => uf33
    IMPLICIT NONE 
    INTEGER :: sp, i, j, k, ioffset,lend,lstart
    Real(prcn) :: sumdiff1, sumdiff2, sumdiff3
    sumdiff1 = zero
    sumdiff2 = zero
    sumdiff3 = zero
    lstart = 1
#if PARALLEL
    lend = nx
#else
    lend = nx+1
#endif
    DO sp =1, nspmx ! Initialize the loop for all the species (nspmx)
#if PARALLEL
       diffn(0,:,:,sp) = zero
       gradphi(0,:,:,:)= zero
#endif
       DO i=lstart,lend 
          ioffset=i
          DO k=1,mz
             DO j=1,my2
                if(xperiodic) then
#if PARALLEL
                   ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                        &,k,sp)+phif(ioffset+1,j,k,sp)) 
#else                   
                   if(ioffset.eq.1) THEN
                      ff1(j,k)=(1./dx2)*(phif(nx,j,k,sp)*heat_ratio(sp)-2.*phif(ioffset,j&
                           &,k,sp)+ phif(ioffset+1,j,k,sp)) 
                   else if(ioffset.eq.nx) THEN
                      ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+phif(1,j,k,sp)/heat_ratio(sp))
                   else if(ioffset.eq.nx+1) THEN
                      ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+phif(2,j,k,sp)/heat_ratio(sp))
                   else
                      
                      ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.*phif(ioffset,j&
                           &,k,sp)+phif(ioffset+1,j,k,sp)) 
                   endif
#endif
	        else
                   ff1(j,k)=(1./dx2)*(phif(ioffset-1,j,k,sp)-2.&
                        &*phif(ioffset,j,k,sp)+phif(ioffset+1,j,k,sp)) 
                endif
                !ff1(j,k)=ff1(j,k)-w2(j,k)*phif(ioffset,j,k,sp)
                
                ff1(j,k)=ff1(j,k)*gamma(sp)
                
                ff5(j,k)= gamma(sp)*dreal(wy(j)*wy(j))*phif(ioffset,j,k,sp)
                ff6(j,k) = gamma(sp)*dreal(wz(k)*wz(k))*phif(ioffset,j,k,sp)
                !end of diffusion terms calculation 
                !beginnig of grad phi tems calculation
#if PARALLEL
                ff2(j,k) = (phif(i+1,j,k,sp)-phif(i-1,j,k,sp))/(two&
                     &*dx)
#else
		
                   if(i.eq.1) then
                      ff2(j,k) = (phif(2,j,k,sp)-phif(mx1,j,k,sp)*heat_ratio(sp))/(two*dx) 
                   else if(i.eq.(nx))then
                      ff2(j,k) = (phif(1,j,k,sp)/heat_ratio(sp)&
                           &-phif(i-1,j,k,sp))/(two*dx) 
                   else if(i.eq.(nx+1))then
                      ff2(j,k) = (phif(2,j,k,sp)/heat_ratio(sp)&
                           &-phif(i-1,j,k,sp))/(two*dx) 
                   ELSE
                      ff2(j,k) = (phif(i+1,j,k,sp)-phif(i-1,j,k,sp))&
                           &/(two*dx)
                   end if
#endif
                   ff3(j,k)=phif(i,j,k,sp)*wy(j)  !starts at foffset+1
                   ff4(j,k)=phif(i,j,k,sp)*wz(k)
                   
                ENDDO
             ENDDO
             
           
          ! CALL ff2cr(ff1(:,:),diffn(i,:,:,sp)
	 
            CALL ff2cr(ff1(:,:),diffn(i,:,:,1))
            CALL ff2cr(ff5(:,:),diffn(i,:,:,2))
            CALL ff2cr(ff6(:,:),diffn(i,:,:,3))
 
             
             CALL ff2cr(ff2(:,:),gradphi(i,:,:,1))
             CALL ff2cr(ff3(:,:),gradphi(i,:,:,2))
             CALL ff2cr(ff4(:,:),gradphi(i,:,:,3))
             
          END DO
#if PARALLEL
             diffn(nx+1,:,:,sp) = zero
             gradphi(nx+1,:,:,:)= zero
#endif
       ENDDO
    !write(90,'(12(2x,e20.10))') t, sumdiff1/(mxf*my*mz), sumdiff2/(mxf&
    !     &*my*mz), sumdiff3/(mxf*my*mz)

       write(*,*)"done diff" !!, gradphi(1:nx,1:my/2,1:10,1)
     END SUBROUTINE calc_diffn_2

 subroutine calc_pdf_vel
                
    implicit none
       
       real(prcn), DIMENSION(:),allocatable ::  hist, radbin, tmp_array,tmp_u, array,arraycopy
       real(prcn) :: u_min, u_max, phirmax,phi_model,kmean
       integer ::  idim, i, j, k, count, unit1,imax,jmax,kmax,ii,iii
       integer :: nbins   , countfl_plane ,count_ne
       character*50 filename1, filename2
    
       call screen_separator(30,'^')
                
    write (*,*) "IN pdf"
    
    nbins = 500      
     
     !   filename2 = trim(run_name)//"_pdf.dat"
     !   open (2, file=trim(filename2), status="replace")

                allocate(hist(nbins), radbin(nbins))
                allocate(tmp_array(count_fluid),tmp_u(count_fluid),arraycopy(count_fluid))
      
     !  open(2,file='MIS1_pdfarray_vel.dat')
        open(2,file='MIS1_array1_vel.dat')
          read(2,*) count_ne
           
          allocate(array(count_ne))
    
         do i=1,count_ne

            read(2,*)array(i)

         enddo

       close(2,status='keep')

        write(*,*)"s",array(count_ne),count_ne

           call unequal_histogram(array(1:count_ne),count_ne,nbins,radbin(1:nbins),hist(1:nbins))

         !!    call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))

         filename2 = trim(run_name)//"_pdf_vel.dat"
         open (2, file=trim(filename2), status="replace")
           do iii=1, nbins
                  ! if (hist(iii)>=small_number) then
                      write(2,"(2d15.7)")radbin(iii),hist(iii)  
                  ! else
                  !    hist(iii) = zero
                  ! endif 
           enddo
                  deallocate(array)
         close(2) 
!!stop

#if 0
!! ---------------------------------------
            phi_fluid_mean(1) = 0.32
            write(*,*)'Re,PHIAVG,phi_fluid_mean(1)=',Re, PHIAVG,phi_fluid_mean(1)
              phi_model = 4.50d0*PHIAVG/(Re+2.718/2.d0)/0.70d0*10.4 !!6.84 !!change

           
                count = 0
        um = zero           
       phim = zero               
       countfl_plane = 0
       
          PHIrmax = zero
           kmean = zero
                u_min = minval(radbin(:))
                u_max = maxval(radbin(:))
          
            do i=1,nx
             do j = 2,nbins-1
                    phirmax = exp(-phi_model*(i-1)*dx)
                 kmean = kmean + (2.8E+19)*exp(-242400.0d0/8.314/(-10.0*phirmax*radbin(j)+800.0))* & 
                            phi_model*hist(j)*dx*(radbin(j+1)-radbin(j-1))/2.0d0
       !  if(i==1)write(*,*)phirmax,radbin(j),hist(j),exp(-242400.0d0/8.314/(-10.0*phirmax*radbin(j)+800.0)), & 
!!(2.8E+19)*exp(-242400.0d0/8.314/(-10.0*phirmax*radbin(j)+800.0))*phi_model*hist(j)*dx*(radbin(j+1)-radbin(j-1))/2.0d0
 
             enddo
            enddo 
             write(*,*)'kpdf=',kmean,'kmean=',(2.8E+19)*exp(-242400/8.314/(-10.0*0.32+800.0))
#endif
#if 0       
                        do i=1, nx,nx/3
                           if(i.eq.1)filename2 = trim(run_name)//"_pdf_1.dat"
                           if(i.eq.1+nx/3)filename2 = trim(run_name)//"_pdf_2.dat"
                           if(i.eq.1+2.0d0*nx/3)filename2 = trim(run_name)//"_pdf_3.dat"
                            open (2, file=trim(filename2), status="replace")
                              countfl_plane = 0
                              count_ne = 0
                              phirmax = zero
                              count = 0
                           write(*,*)'i=',i,1+2.0d0*nx/3
                        do j=1, my
                        do k=1, mz
                            if (fluid_atijk(i,j,k)) then
                                     count = count+1
                                 if(phir(i,j,k,1).gt.phirmax) then 
                                   phirmax= phir(i,j,k,1)
                                     imax= i 
                                     jmax= j
                                    kmax =k
                                 endif
                                 tmp_u(count)= velr(i,j,k,1)   
                                 tmp_array(count) =phir(i,j,k,1)
                         if(phir(i,j,k,1).ge.0.0d0) then 
                           count_ne = count_ne+1
                           arraycopy(count_ne) =phir(i,j,k,1)
                         endif
                         phim(i,1) = phim(i,1) + velr(i,j,k,1)*phir(i,j,k,1)
                      
                         um(i) = um(i) + velr(i,j,k,1)
             !!           countfl_plane(i) = countfl_plane(i) + 1
                 
         
           !!  phim(i,isp) = phim(i,isp)/(um(i)*real(countfl_plane(i),prcn))
         
                !                 write (2,"(1e15.7)") tmp_array(count)   
                         if(count.eq.1) write(*,*)'phir,tmp',phir(i,j,k,1),tmp_array(count) 
                         if(mod(count,50000).eq.0)     write(*,*) count, count_fluid, tmp_array(count)
                 
                          endif
                        enddo
                        enddo
                
                         ! phim(i,1) = phim(i,1)/um(i)
                           allocate(array(count_ne))
                           do ii= 1,count_ne
                             array(ii)= arraycopy(ii)
                           enddo
                    call unequal_histogram(array(1:count_ne), count_ne,nbins,radbin(1:nbins),hist(1:nbins))
   
                      !    call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))
                          
                            do iii=1, nbins
                            if (hist(iii)>small_number) write(2,"(2d15.7)") radbin(iii),hist(iii)   
                            enddo
       

                          deallocate(array)
                           close(2) 
                        enddo
#endif               
#if 0                  
             write(*,*)'i,j,k max',imax,jmax,kmax,phir(imax,jmax,kmax,1)
             write(*,*)'u_max and min=',minval(tmp_u(:)), maxval(tmp_u(:))        
             write(*,*)'count',count,'nbins',nbins, 'count_ne',count_ne
           !  write(*,*)'phim',phim(:,1)
#endif 

#if 0                                        
               u_min = minval(tmp_array(:))
                u_max = maxval(tmp_array(:))
                 do i = 1,count
                    tmp_array(i)= (tmp_array(i)-u_min)/ (u_max- u_min)
                 enddo                             
                  do i= 1,count_ne
                      array(i)= arraycopy(i)
                  enddo     
                  u_min = minval(array(:))
                  u_max = maxval(array(:))
                 !  do i= 1,count_ne
                     ! array(i)= (array(i)-u_min)/ (u_max- u_min)
                 ! enddo 
      !    call make_histogram(array(1:count_ne), count_ne, nbins,radbin(1:nbins),hist(1:nbins))
   
       !   call make_histogram(tmp_array(1:count), count, nbins,radbin(1:nbins), hist(1:nbins))
                              
    
                !  do i=1, nbins
                !       if (hist(i)>small_number) write(2,"(2d15.7)") radbin(i), hist(i)
               !   enddo
                        
#endif                     
                 deallocate(tmp_array,tmp_u)
                  
   end subroutine calc_pdf_vel

subroutine calc_generate_vel

    implicit none

       real(prcn), DIMENSION(:),allocatable ::  hist, radbin, tmp_array,tmp_u, array,arraycopy
       real(prcn) :: u_min, u_max, phirmax
       integer ::  idim, i, j, k, count, unit1,imax,jmax,kmax,ii,iii
       integer :: nbins   , countfl_plane ,count_ne
       character*50 filename1, filename2

       call screen_separator(30,'^')

    write (*,*) "IN pdf"

    nbins = 150

     !   filename2 = trim(run_name)//"_pdf.dat"
     !   open (2, file=trim(filename2), status="replace")
       
                allocate(tmp_array(count_fluid),tmp_u(count_fluid),arraycopy(count_fluid))

             radbin = 0.0d0
                count = 0
        um = zero
       phim = zero
       countfl_plane = 0
       count_ne = 0
          phirmax = zero
                 filename2 = trim(run_name)//"_array1_vel.dat"
                      
                        
                          !! if(i.eq.1+nx/3)filename2 = trim(run_name)//"_array2.dat"

                 
                            open (2, file=trim(filename2), status="replace")
                              countfl_plane = 0
                              count_ne = 0
                              phirmax = zero
                              count = 0
                           
                            write(*,*)"umean(1)",umean(1)

                       do i=1,nx
                        do j=1, my
                        do k=1, mz 
                            if (fluid_atijk(i,j,k)) then
                                     count = count+1      
                           
                              !! arraycopy(count) = sqrt(velr(i,j,k,1)**2+velr(i,j,k,2)**2+velr(i,j,k,3)**2)/umean(1)
                               arraycopy(count) = velr(i,j,k,1) ! /umean(1)

                            endif
                        enddo
                        enddo
                        enddo
                           count_ne = count

                           allocate(array(count_ne))
                           do ii= 1,count_ne
                             array(ii)= arraycopy(ii)
                           enddo
                            
                            write(2,*)count_ne
                            
                            do iii=1, count_ne
                               write(2,*) array(iii)
                            enddo
                       

                          deallocate(array)
                           close(2)
                        

           !  write(*,*)'i,j,k max',imax
    
         deallocate(tmp_array,tmp_u,arraycopy) 

   end subroutine calc_generate_vel

end module outputscalar
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
