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


MODULE boundary_condition
#include "ibm.h"
  USE precision 
  use constants 
  USE global_data
  USE fftw_interface   
  USE dependent_functions
  !USE  initialize_flo
  USE collision_mod
  !uncomment the following line omega needs to be calculated 
  Use nlmainarrays, Only : ur=>ubcp,pr=> pbcp!, nlr=>nlbcp, onlr=>onlbcp 
  USE bcsetarrays, ONLY :  fr, ppr, diffn 
  USE nlarrays, ONLY : ff1=>uf1,ff2=>uf2,ff3=>uf3, fr1=>ur1, fr2=>ur2&
       &, fr3=>ur3, dudy=>ur11, dudz=>ur22, dvdz=>ur33, dwdy=>ur12
  implicit none 
  PRIVATE 
  REAL(prcn) ::  vort(ndim), da(2)
  REAL(prcn) ::  usph(ndim),ucar(ndim) 
  REAL(prcn) ::   furf

  REAL(prcn) ::  pl,nll(ndim),onll(ndim),ppll(ndim),dfll(ndim)
  REAL(prcn) ::  snorm(ndim),rad
  REAL(prcn) ::  unorm(ndim),utang(ndim)
  REAL(prcn) ::  sumforcepoint(3),sumforcenode(3)
  REAL(prcn) ::  pnmag,dfmag,unmago,unmagb,unmagi, normal(ndim)
  INTEGER ioffset
  REAL(prcn) ::tempfor1(ndim),  tempt, tempfor2(ndim),unmagl,ucarl,unorml(ndim),&
       & utangl(ndim),ucarf(ndim), tempdist
  REAL(prcn) ::  drm,utangm(ndim)&
       &,unormm(ndim),unmag, tmp_rnp, upimod, force_tmp(ndim), ubnorm(ndim), &
       & ubtang(ndim), unormi(ndim), utangi(ndim), ubnmag, utemp
  INTEGER :: isp, i,j,k,l,m,n, is(ndim),iii(ndim)&
       &,io(ndim), c,d,rcount,bcount,ib, ie, jb, je, kb, ke, onew,&
       & ii, jj, kk , pcell(3), mxftemp
  CHARaCTER*80 :: FILENAME1, FILENAME2, FILENAME3, FILENAME4, FILENAME5, FILENAME6, FILENAME7
  LOGICAL:: filexist, isopen, REVERSE
  REAL(prcn) :: ulb(ndim),plb,plo,pli,deriv_p,deriv_u,deriv_un,&
       & frmeanfluidloc(ndim), frmeanloc(ndim), norm_factor, mpg_without_unst(ndim)

  !REAL(prcn), ALLOCATABLE,DIMENSION(:,:,:) :: fromleft,fromright
  REAL(prcn) :: unsteady_term(3)
#if PARALLEL
  REAL(prcn) :: visc_totalloc(ndim), pres_totalloc(ndim)
#endif
  Public :: bcset, calc_pgrad, calc_visc, compute_omega
CONTAINS

#if 0
	SUBROUTINE compute_mpg_at_n2(rks)
		IMPLICIT NONE 
		integer, Intent(in) ::  rks
		integer ::  n , idim, m, l, iphs
		REAL(prcn) :: force_total(ndim), x

		if(I_AM_NODE_ZERO)WRITE(*,*)'CALCULATING MPG @ NTH TIME-STEP'
!		CALL compute_omega
		visc_total_old = visc_total
		pres_total = zero 
		visc_total = zero 
#if PARALLEL
		pres_totalloc = zero 
		visc_totalloc = zero 
#endif

		DO m=1,nbody
			if(rks.eq.1) then 
				DO n=1,ndim
					pres(m,n)=zero
					visc(m,n)=zero
					torq(m,n)=zero
#if PARALLEL
					presloc(m,n)=zero
					viscloc(m,n)=zero
					torqloc(m,n)=zero
					force_loc(m,n) = zero
#endif
				ENDDO
			end if
       
			CALL calc_pres_visc_drag2(m,rks)
		end DO

		do idim=1,ndim
			pres_avg(idim) = SUM(pres(1:nbody,idim))!/real(nspec1,prcn)
			visc_avg(idim) = SUM(visc(1:nbody,idim))!/real(nspec1,prcn)
		end do
    
		pres_drag= DSQRT(dot_product(pres_avg(1:ndim), pres_avg(1:ndim)))
		visc_drag = DSQRT(dot_product(visc_avg(1:ndim), visc_avg(1:ndim)))

		pres_drag = pres_drag/norm_factor/(one-maxvolfrac)
		visc_drag = visc_drag/norm_factor/(one-maxvolfrac)
		total_drag = (pres_drag+visc_drag)

		IF(FROM_POST) goto 1000

		IF(.not.set_mpg) then 
			if(move_particles)then
				mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + cf*(usmean(:)-usmean_des(:)) 
				mpg(:) = mpg(:) + (pres_total(:) -visc_total(:)) * (one/(rhof*(one-mean_volfrac)) + one/(mean_volfrac*rhos))/voldom
				mpg(:) = mpg(:)/(one/rhof - one/rhos)
				frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf*(usmean_des(:)-usmean(:))
			else
				mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + (pres_total(:)/voldom -visc_total(:)/voldom)/(one-maxvolfrac)  
			end if

			mpg_without_unst(:) = pres_total(:)/voldom -(visc_total(:))/voldom  

			if(include_frmeanfluid)  then 
				WRITE(*,*)'INCLUDING FRMEAN FLUID'
				mpg(:) = mpg(:) + frmeanfluid(:)
			end if
			mpg_without_unst(:) = mpg_without_unst(:)/(one - maxvolfrac)
			mpg_without_unst(:) = mpg_without_unst(:)/(coef(rks,1) + coef(rks,2))
		else
			if(move_particles)then
				frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf *(usmean_des(:)-usmean(:))
			endif
		endif
1000	continue 
    
		unsteady_term(:) = cf*(ufmean_des(:)-ufmean(:))*(one-maxvolfrac)
		force_total = zero 

		Do m = 1, nbody
			force(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:)) -(coef(rks,1)+coef(rks,2))*mpg(:)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0 
			force_total(:) = force_total(:) + force(m,:)

			force_chem(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:))
		End Do

		if(I_AM_NODE_ZERO)then
			if(move_particles)then
				WRITE(*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:)*voldom/norm_factor
				WRITE(*,'(A30,3(2x,g17.8))') 'PRES TERM = ', (coef(rks,1)+coef(rks,2))*pres_total(:)*(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)*norm_factor)
		       
				WRITE(*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:)*(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)*norm_factor) 

				WRITE(*,'(A30,3(2x,g17.8))') 'mpg TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom/norm_factor
			else
				WRITE(*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:)*voldom/(one-maxvolfrac)/norm_factor

				WRITE(*,'(A30,3(2x,g17.8))') 'PRES TERM = ', (coef(rks,1)+coef(rks,2))*pres_total(:)/((one-maxvolfrac))/norm_factor

				WRITE(*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:)/((one-maxvolfrac))/norm_factor 

				WRITE(*,'(A30,3(2x,g17.8))') 'mpg TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom/norm_factor
				WRITE(*,'(A30,3(2x,g17.8))') 'mpg WITHOUT UNST = ', (coef(rks,1)+coef(rks,2))*mpg_without_unst(:)*voldom/norm_factor

				WRITE(*,'(A30,3(2x,g17.8))') 'FORCE_TOTAL = ', (coef(rks,1)+coef(rks,2))*force_total(:)/norm_factor

			end if
		end if
		RETURN 
	end SUBROUTINE compute_mpg_at_n2


	SUBROUTINE calc_pres_visc_drag2(m,rks)
		use mypost_process
		IMPLICIT NONE
		integer, intent(in) :: m,rks

		integer :: i, j, k, ii, jj, kk
		integer :: istart, iend, jstart, jend, kstart, kend
		integer :: local_solid, local_fluid, ibody
		real(prcn) :: local_length, local_ufmean(ndim), local_usmean(ndim), local_volfrac, local_Rem, part_force, local_meanslip(ndim), local_meanslip_mag


		!-----------------------------------------------------------------------


		if (lybyd_small<lybyd-small_number) then
			local_length = lybyd_small*(2*radbdy(1))

			istart = xc(m,1) - local_length/2
			iend   = xc(m,1) + local_length/2
			jstart = xc(m,2) - local_length/2
			jend   = xc(m,2) + local_length/2
			kstart = xc(m,3) - local_length/2
			kend   = xc(m,3) + local_length/2

			local_solid = 0
			local_fluid = 0
			local_ufmean = zero
			do kk=kstart, kend
				do jj=jstart, jend
					do ii=istart, iend
						if (ii<1) then
							i=ii+nx
						elseif (ii>nx) then
							i=ii-nx
						else
							i=ii
						endif

						if (jj<1) then
							j=jj+my
						elseif (jj>my) then
							j=jj-my
						else
							j=jj
						endif

						if (kk<1) then
							k=kk+mz
						elseif (kk>mz) then
							k=kk-mz
						else
							k=kk
						endif

						if (fluid_atijk(i,j,k)) then
							local_fluid = local_fluid+1
							local_ufmean(:) = local_ufmean(:) + ur(i,j,k,:)
						else
							local_solid = local_solid+1
							call find_body(i,j,k,ibody, doml(2)/2)
							local_usmean(:) = local_usmean(:) + velbdy(ibody,:)
						endif
					enddo
				enddo
			enddo

			local_volfrac = dble(local_solid) / dble(local_fluid+local_solid)
			local_ufmean(:) = local_ufmean(:) / local_fluid
		else
			local_volfrac = maxvolfrac
			local_ufmean(:) = ufmean(:)
		endif

		local_usmean(:) = velbdy(m,:)

		local_meanslip(:) = local_ufmean(:) - local_usmean(:)
		local_meanslip_mag = sqrt(dot_product(local_meanslip, local_meanslip))

		local_Rem = (1-local_volfrac) * local_meanslip_mag * (radbdy(m)*2*dx) / vis

		call compute_ibm_drag(local_volfrac, local_Rem, part_force)

!write (*,"(1i,4d15.7)") m, local_meanslip(:), local_meanslip_mag
!write (*,"(3d15.7)") local_volfrac, local_Rem, part_force
!write (*,*)

		part_force = part_force * 3 * pi * (radbdy(m)*2*dx) * vis * (1-local_volfrac)
		visc(m,:) = part_force * local_meanslip(:) + (coef(rks,1)+coef(rks,2)) * mpg(:) * pi * (two*radbdy(m)*dx)**3.d0 /6.d0 
		pres(m,:) = zero

		visc_total(:) = visc_total(:) + visc(m,:)
		pres_total(:) = pres_total(:) + pres(m,:)

	END SUBROUTINE calc_pres_visc_drag2

#endif












  SUBROUTINE compute_mpg_at_n(rks)

    IMPLICIT NONE 

    integer, Intent(in) ::  rks
    integer ::  n , idim, m, l, iphs
    
    REAL(prcn) :: force_total(ndim), x
    
    if(I_AM_NODE_ZERO)WRITE(*,*)'CALCULATING MPG @ NTH TIME-STEP'

    
    CALL compute_omega

    visc_total_old = visc_total

    pres_total = zero 
    visc_total = zero 
#if PARALLEL
    pres_totalloc = zero 
    visc_totalloc = zero 
#endif

    DO m=1,nbody            ! loop over bodies

!if (I_AM_NODE_ZERO) write (*,*) "particle = ", m

       if(rks.eq.1) then 
          DO n=1,ndim
             pres(m,n)=zero
             visc(m,n)=zero
             torq(m,n)=zero
#if PARALLEL
             presloc(m,n)=zero
             viscloc(m,n)=zero
             torqloc(m,n)=zero
             force_loc(m,n) = zero
#endif
          ENDDO
       end if
       
       iphs = 1!part_array(m)%iphs
       nbnd = phase_array(iphs)%nbnd
       NULLIFY(bndarray)
       bndarray => phase_array(iphs)%bndpts
       
       da(1)=4.*pi*(radbdy(m)*dx)**2./real(nbnd,prcn)
!!$       PRINT*,'  m = ', m, 'phase = ', iphs, 'nbnd = ', nbnd
!!$       PRINT*,'  l = ', 1, 'bndpt = ', bndarray(1:ndim,100)
!!$       READ(*,*)
       CALL calc_pres_visc_drag(m,rks)
#if PARALLEL

!if (m==9426) write (*,"(1i,3d15.7)") myid, presloc(m,:)
!if (m==9426) write (*,"(1i,3d15.7)") myid, viscloc(m,:)
!if (m==9426) write (*,"(1i,3d15.7)") myid, torqloc(m,:)

       do n = 1, ndim
          GLOBAL_DOUBLE_SUM(presloc(m,n),pres(m,n),1,decomp_group)
          GLOBAL_DOUBLE_SUM(viscloc(m,n),visc(m,n),1,decomp_group)
          if(move_particles) GLOBAL_DOUBLE_SUM(torqloc(m,n),torq(m,n),1,decomp_group)
       end do
#endif
    end DO!CLOSE LOOP OVER ALL BODIES


#if PARALLEL    
    GLOBAL_DOUBLE_SUM(pres_totalloc(1),pres_total(1),3,decomp_group)
    GLOBAL_DOUBLE_SUM(visc_totalloc(1),visc_total(1),3,decomp_group)
#endif
    
    do idim=1,ndim
       pres_avg(idim) = SUM(pres(1:nbody,idim))!/real(nspec1,prcn)
       visc_avg(idim) = SUM(visc(1:nbody,idim))!/real(nspec1,prcn)
    end do
    
    pres_drag= DSQRT(dot_product(pres_avg(1:ndim), pres_avg(1:ndim)))
    
    visc_drag = DSQRT(dot_product(visc_avg(1:ndim), visc_avg(1:ndim)))
    
    pres_drag = pres_drag/norm_factor/(one-maxvolfrac)
    visc_drag = visc_drag/norm_factor/(one-maxvolfrac)
    total_drag = (pres_drag+visc_drag)
    
!    IF(FROM_POST) goto 1000
    
    
    IF(.not.set_mpg) then 
       if(move_particles)then
          mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + cf*(usmean(:)-usmean_des(:)) 
!!$          mpg(:) = mpg(:) + (pres_total(:) -(visc_total(:)))*(one&
!!$               &/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))&
!!$               &/voldom
          mpg(:) = mpg(:) + (pres_total(:) -(visc_total(:)))*(one/(rhof*(one-mean_volfrac)) + one/(mean_volfrac*rhos))/voldom
          mpg(:) = mpg(:)/(one/rhof - one/rhos)
          frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf *(usmean_des(:)-usmean(:))
       else
          mpg(:) = cf*(ufmean_des(:)-ufmean(:)) + (pres_total(:)/voldom -(visc_total(:))/voldom)/(one-maxvolfrac)  
       end if
       
       mpg_without_unst(:) = pres_total(:)/voldom -(visc_total(:))/voldom  
       
       if(include_frmeanfluid)  then 
          WRITE(*,*)'INCLUDING FRMEAN FLUID'
          mpg(:) = mpg(:) + frmeanfluid(:)
       end if
       !if(.not.move_particles) 
       !mpg(:) = mpg(:)/(one - maxvolfrac)
       mpg_without_unst(:) = mpg_without_unst(:)/(one - maxvolfrac)

       ! mpg(:) = mpg(:)/(coef(rks,1) + coef(rks,2))
       mpg_without_unst(:) = mpg_without_unst(:)/(coef(rks,1) + coef(rks,2))
	else
		if(move_particles)then
	        frame_accln(:) = -mpg(:)/rhos - (pres_total(:) - visc_total(:))/(mean_volfrac*rhos*voldom) + cf *(usmean_des(:)-usmean(:))

!			write (*,"(3d15.7)") usmean_des(:)
!			write (*,"(3d15.7)") usmean(:)
!			write (*,"(3d15.7)") frame_accln(:)
!			read (*,*)
		endif
    endif

1000 continue 
    
    unsteady_term(:) = cf*(ufmean_des(:)-ufmean(:))*(one-maxvolfrac)
    force_total = zero 
    
    Do m = 1, nbody
       force(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:)) -(coef(rks,1)+coef(rks,2))*mpg(:)*pi*((two*radbdy(m)*dx)**3.d0)/6.d0 
       force_total(:) = force_total(:) + force(m,:)
       
       force_chem(m,:) =  (coef(rks,1)+coef(rks,2))*(pres(m,:)+visc(m,:))
    End Do
    
    if(I_AM_NODE_ZERO)then
       if(move_particles)then
          WRITE(*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:)*voldom/norm_factor
         
          WRITE(*,'(A30,3(2x,g17.8))') 'PRES TERM = ', (coef(rks,1)+coef(rks,2))*pres_total(:)*(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)*norm_factor)
          
          WRITE(*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:)*(one/(rhof*(one-maxvolfrac)) + one/(maxvolfrac*rhos))/((one/rhof - one/rhos)*norm_factor) 
          
          WRITE(*,'(A30,3(2x,g17.8))') 'mpg TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom/norm_factor
       else
          
          WRITE(*,'(A30,3(2x,g17.8))') 'UNSTEADY TERM = ', unsteady_term(:)*voldom/(one-maxvolfrac)/norm_factor
          !WRITE(*,'(A30,3(2x,g17.8))') 'SOURCE HYDRO = ', -source_hydro(:)*voldom/(one)/norm_factor
          
          WRITE(*,'(A30,3(2x,g17.8))') 'PRES TERM = ', (coef(rks,1)+coef(rks,2))*pres_total(:)/((one-maxvolfrac))/norm_factor
          
          WRITE(*,'(A30,3(2x,g17.8))') 'VISC TERM = ', -(coef(rks,1)+coef(rks,2))*visc_total(:)/((one-maxvolfrac))/norm_factor 
          
          WRITE(*,'(A30,3(2x,g17.8))') 'mpg TERM = ', (coef(rks,1)+coef(rks,2))*mpg(:)*voldom/norm_factor
          WRITE(*,'(A30,3(2x,g17.8))') 'mpg WITHOUT UNST = ', (coef(rks,1)+coef(rks,2))*mpg_without_unst(:)*voldom/norm_factor
          
          WRITE(*,'(A30,3(2x,g17.8))') 'FORCE_TOTAL = ', (coef(rks,1)+coef(rks,2))*force_total(:)/norm_factor
          
       end if
    end if
    RETURN 

  end SUBROUTINE compute_mpg_at_n
  

  SUBROUTINE bcset(rks)!(ur,pr,nlr,onlr)

    IMPLICIT NONE 

    integer, Intent(in) ::  rks
    REAL(prcn) ::  c_drag_st(nbody), c_drag_tr(nbody), norm_drag_chem, norm_drag_poly_spec(nphases), norm_drag_poly
    REAL(prcn) ::  avg_force(ndim),avg_force_chem(ndim),avg_force_spec(nphases,ndim),avg_force_chem_spec(nphases,ndim), tmp_ferror_array(nerr_steps),LHS_UF(ndim), MOD_LHS_UF, cpu0, cpu1, cpu2, cpu3

    INTEGER :: idim, iphs, pstart, pend
    REAL(prcn) :: slip_diff
#if PARALLEL
    REAL(prcn) :: frsumloc(ndim), frsum(ndim)
    COMPLEX(prcn) :: fsumloc,fsum, frtemp
#endif
    !ag=7.8d-4
    CHARACTER(LEN=80) :: formfile

    formfile='formatted'

	if (I_AM_NODE_ZERO) CALL CPU_TIME (CPU2) 

    IF(I_AM_NODE_ZERO.and.first_pass) THEN 
       FILENAME1 = TRIM(RUN_NAME)//'_norm_drag_chem'//'.dat'
       FILENAME2 = TRIM(RUN_NAME)//'_dragcoeffy'//'.dat'
       FILENAME3 = TRIM(RUN_NAME)//'_dragcoeffz'//'.dat'
!!$       FILENAME2 = TRIM(RUN_NAME)//'_dragcoeffsum'//'.dat'
       FILENAME4 = TRIM(RUN_NAME)//'_norm_drag'//'.dat'
       FILENAME5 = TRIM(RUN_NAME)//'_force_part'//'.dat'
       FILENAME6 = TRIM(RUN_NAME)//'_drag_components'//'.dat'
       FILENAME7 = TRIM(RUN_NAME)//'_normdrag_poly'//'.dat'
       
       CALL  RUN_TIME_FILE_OPENER(unitnormdragchem,FILENAME1, formfile)
       
       CALL  RUN_TIME_FILE_OPENER(unitdragtavg,FILENAME2, formfile)
       
       CALL  RUN_TIME_FILE_OPENER(unitnormdrag,FILENAME4,formfile)

       CALL  RUN_TIME_FILE_OPENER(unitnormdragpoly,FILENAME7,formfile)

       CALL  RUN_TIME_FILE_OPENER(unitforce,FILENAME5,formfile)

       CALL  RUN_TIME_FILE_OPENER(unitdrag_comps,FILENAME6,formfile)

       !first_pass=.false.
    ENDIF

    norm_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*char_length)*real(nbody,prcn)
!!$    if(TRIM(input_type).eq.'lubtest')norm_factor = (3.d0*pi*vis*(ucharmod)*char_length)*real(nbody,prcn)
!    norm_factor = (3.d0*pi*vis*(one-maxvolfrac)*char_length)*real(nbody,prcn)
    if(TRIM(input_type).eq.'lubtest')norm_factor = (3.d0*pi*vis*(one-maxvolfrac)*char_length)*real(nbody,prcn)

    if(.not.only_dem)  CALL calc_visc
    
    IF(rks.eq.1)then
       Do m = 1, nbody
          DO n=1,ndim
             force(m,n)=zero
             !zero the drag on the body --> the total, viscous and pressure drags.
          ENDDO
       ENDDO
    ENDIF

    !CALL compute_mpg_at_n(rks)
    
		if (.not.only_dem) then
			if (use_drag_law) then
!				CALL compute_mpg_at_n2(rks)
			else
				CALL compute_mpg_at_n(rks)
			endif
		endif

    if(move_particles)then
       if(I_AM_NODE_ZERO) CALL calc_part_statistics(rks)

       if (I_AM_NODE_ZERO) CALL CPU_TIME (CPU0) 
       CALL des_time_march(.FALSE.)
       if (I_AM_NODE_ZERO) then
			CALL CPU_TIME (CPU1)
			dem_time_current = cpu1-cpu0
			dem_time = dem_time + dem_time_current
       endif

       CALL grid_nodes_insphere
       do m = 1, nbody
          CALL update_nrpr_array(m)
       end do
    end if
    
    if(only_dem)RETURN
    
    drm = dr !1.0
    radm = 1.0
    
    furf = 1.0
    CALL calc_pgrad
    
    !ALLOCATE(fromleft(my,mz,ndim),fromright(my,mz,ndim))
    DO i=1,3
       sumforcepoint(i)=zero
       sumforcenode(i)=zero
    ENDDO
    
    apmax = 0.d0
    norm_drag_spec = 0.d0
    norm_drag = 0.d0
    norm_drag_chem_spec = zero
    norm_drag_chem = zero
    norm_drag_poly = zero
    norm_drag_poly_spec = zero

    frmean(1:ndim) = zero
    frmeanloc(1:ndim) = zero
    frmeanfluid(1:ndim) = zero
    frmeanfluidloc(1:ndim) = zero
    DO l=1,ndim
       DO k=1,mz
          DO j=1,my
#if PARALLEL          
             fr(0,j,k,l)=zero
#endif
             !fromleft(j,k,l) = zero
             !fromright(j,k,l) = zero
             
             DO i=1,nx+1
                fr(i,j,k,l)=zero
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    !-----------------------------------------------------------------------
    !     loop over all bodies

   
    DO m=1,nbody            ! loop over bodies
       !---------------------------------------------------------------------------
       !     zero the drag on the body --> the total, viscous and pressure drags.
       
       iphs = 1!part_array(m)%iphs
       nbnd = phase_array(iphs)%nbnd
       nrpr = phase_array(iphs)%nrpr
       NULLIFY(bndarray)
       bndarray => phase_array(iphs)%bndpts
       
       da(1)=4.*pi*(radbdy(m)*dx)**2./real(nbnd,prcn)
       da(2)=4.*pi*(radibdy(m)*dx)**2./real(part_array(m)%nrpr_active, prcn)

       if(dobnd) CALL calc_bnd_forcing(m,rks)
       
       IF(dorpr) CALL calc_inner_reversal_forcing(m,rks)
!#if PARALLEL
!!$       if(set_umean)then
!!$          do n = 1, ndim
!!$             GLOBAL_DOUBLE_SUM(force_loc(m,n),force(m,n),1,decomp_group)
!!$          end do
!!$       end if
!#endif
       
    end DO!CLOSE LOOP OVER ALL BODIES
    
#if PARALLEL
    GLOBAL_DOUBLE_SUM(frmeanloc(1),frmean(1),ndim,decomp_group)
    if(include_frmeanfluid) GLOBAL_DOUBLE_SUM(frmeanfluidloc(1),frmeanfluid(1),ndim,decomp_group)
#endif
    frmean(1:ndim) = frmean(1:ndim)/real(mx1*my*mz,prcn)
    frmeanfluid(1:ndim) = frmeanfluid(1:ndim)/real(count_fluid,prcn)
    if(debug_check)then
       if(I_AM_NODE_ZERO)WRITE(*,'(A25,3(2x,g17.8))')'frmean(1:ndim)=', frmean(1:ndim)
    end if
    !if(I_AM_NODE_ZERO)WRITE(*,*)'frmean(1:ndim)=', frmean(1:ndim)
    
#if PARALLEL
    do n = 1,ndim
       VECSENDRECV(fr(nx+1,1,1,n),1,twodrslice,toproc,1,diffn(0,1,1,n),1,fromproc,1,decomp_group,status)
       VECSENDRECV(fr(0,1,1,n),1,twodrslice,fromproc,1,diffn(nx+1,1,1,n),1,toproc,1,decomp_group,status)
       !fromleft(j,k,n) = diffn(0,j,k,n)
       !fromright(j,k,n) = fr(1,j,k,n)
       !frtemp = fr(1,j,k,n) + diffn(0,j,k,n)
       !fr(1,j,k,n) = fromleft(j,k,n)+fromright(j,k,n)
       
       do k = 1, mz
          do j = 1, my 
             fr(1,j,k,n) = fr(1,j,k,n) + diffn(0,j,k,n)
             fr(nx,j,k,n) = fr(nx,j,k,n) + diffn(nx+1,j,k,n)
          end do
       end do
       if(debug_check)then
          frsumloc(n) = SUM(fr(1:nx,:,:,n))
          GLOBAL_DOUBLE_SUM(frsumloc(n),frsum(n),1,decomp_group)
       end if
    end DO
    if(debug_check)then
       if(I_AM_NODE_ZERO)WRITE(*,*)'frmean(1:ndim) FROM  DOMAIN = ', frsum(1:ndim)/(mx1&
            &*my*mz)
    end if
#endif
    
    
#if 0 
    if(iglobstep.eq.1)then
#if PARALLEL
       if(myid.eq.1) CALL write_real_forcing
#else
       CALL write_real_forcing
#endif
       PARALLEL_FINISH()
       STOP
    endif
#endif

    DO n = 1, ndim 
       DO i = 1, nx
          do j = 1, my 
             do k = 1, mz
                fr(i,j,k,n) = fr(i,j,k,n) - frmean(n)
             end do
          end do
          
          CALL ff2rc(fr(i,1:my,1:mz,n), ff(i,1:my2,1:mz,n))
          
          if(i.eq.1)then
             VECSENDRECV(ff(i,1,1,n),1,twodcslice,fromproc,1,ff(nx+1,1,1,n),1,toproc,1,decomp_group,status)
          end if
          !WARNING: FF IS ALSO USED BY SCALAR ARRAYS. THEREFORE, ISCALON=1, THEN DO NOT USE UNDER RELAXATION HERE
          !ff(i,:,:,n) = ff(i,:,:,n) + furf*(ff1(:,:)-ff(i,:,:,n))
!          if(n.eq.1) Write(*,'(I,3(2x,g17.8))')GLOBAL_INDEX(i), SUM(fr(i,1:my,1:mz,1))/(real(my*mz,prcn)),dreal(ff(i,1,1,1))
       END DO
       
    END DO
    !This is very important or else the gauss siedel will fail as the
    ! div of ff at i = mx1 is calculated from ff at mxf and mxf-2
#if !PARALLEL
    ff(nx+1, :, :,1:ndim) = ff(1,:,:,1:ndim)
#endif
   
    
    
#if PARALLEL
#if 0
    PARALLEL_FINISH()
    STOP
    
!!$             Write(*,'(I,3(2x,g17.8))')GLOBAL_INDEX(i), fr(i,1,1,1),dreal(ff(i,1,1,1)), dreal(ff2(1,1))
!!$    fsumloc = SUM(ff(1:nx,1,1,1))
!!$    GLOBAL_COMPLEX_SUM(fsumloc,fsum,1,decomp_group)
!!$    if(I_AM_NODE_ZERO)PRINT*,'FLO: AVG FLUC FORCE  = ',fsum/(mx1)
#endif
#else
    if(I_AM_NODE_ZERO)PRINT*,'FLO: AVG FLUC FORCE  = ',SUM(ff(1:mx1,1,1,1))/(mx1)
#endif
    

    !-----------------------------------------------------------------------
    !      POST PROCESSING AND DIAGNOSTICS
    !-----------------------------------------------------------------------
    pstart = 1
    do iphs = 1, nphases
       pend = pstart + phase_array(iphs)%npart - 1
       do idim=1,ndim
          avg_force_spec(iphs,idim) = SUM(force(pstart:pend,idim))/phase_array(iphs)%npart
          avg_force_chem_spec(iphs,idim) = SUM(force_chem(pstart:pend,idim))/phase_array(iphs)%npart
       end do
       pstart = pend + 1
    end do

    avg_force = zero
    avg_force_chem = zero

    do iphs = 1, nphases
       avg_force(1:ndim) = avg_force(1:ndim) + phase_array(iphs)%volfrac*avg_force_spec(iphs,1:ndim)
       
       avg_force_chem(1:ndim) = avg_force_chem(1:ndim) + phase_array(iphs)%volfrac*avg_force_chem_spec(iphs,1:ndim)
    end do
    
    avg_force(1:ndim) = avg_force(1:ndim)/mean_volfrac
    avg_force_chem(1:ndim) = avg_force_chem(1:ndim)/mean_volfrac
!!$    do iphs = 1, nphases
!!$       if(I_AM_NODE_ZERO)WRITE(*,'(A,4(2x,g17.8))') 'FORCES: avg_force_spec : ', phase_array(iphs)%volfrac,avg_force_spec(iphs,1:ndim)
!!$    end do
!!$    if(I_AM_NODE_ZERO)WRITE(*,'(A,3(2x,g17.8))') 'FORCES: avg_force : ', avg_force(1:ndim)
    do iphs = 1, nphases
       norm_drag_spec(iphs) = DSQRT(avg_force_spec(iphs,1)**2.d0 + avg_force_spec(iphs,2)**2.d0 + avg_force_spec(iphs,3)**2.d0)
       norm_drag_chem_spec(iphs) = DSQRT(avg_force_chem_spec(iphs,1)**2.d0 + avg_force_chem_spec(iphs,2)**2.d0 + avg_force_chem_spec(iphs,3)**2.d0)

    end do
    
    norm_drag = DSQRT(avg_force(1)**2.d0 + avg_force(2)**2.d0 + avg_force(3)**2.d0)
    norm_drag_chem = DSQRT(avg_force_chem(1)**2.d0 + avg_force_chem(2)**2.d0 + avg_force_chem(3)**2.d0)
    

    mpg_avg = DSQRT(dot_product(mpg(1:ndim), mpg(1:ndim)))
    
    total_drag_mpg = mpg_avg*voldom/norm_factor

    LHS_UF(:) = dufmeandt(:)*voldom -pres_total(:)/((one-maxvolfrac)) + visc_total(:)/((one-maxvolfrac))

    MOD_LHS_UF = DSQRT(dot_product(LHS_UF(1:3), LHS_UF(1:3)))

    MOD_LHS_UF = MOD_LHS_UF/norm_factor 
    
    if(I_AM_NODE_ZERO)WRITE(*,'(A,5(2x,g17.8))') 'DRAG: PRES, VISC, TOTAL, MOD_LHS_UF/N, TOTAL_MPG', pres_drag,visc_drag, total_drag, MOD_LHS_UF, total_drag_mpg 
    
    if(TRIM(input_type).eq.'lubtest')then
       do iphs = 1, nphases
          norm_drag_spec(iphs) = norm_drag_spec(iphs)/(3.d0*pi*vis*(ucharmod)*phase_array(iphs)%dia)
          norm_drag_poly_spec(iphs) = norm_drag_chem_spec(iphs)
          norm_drag_chem_spec(iphs) = norm_drag_chem_spec(iphs)/(3.d0*pi*vis*(ucharmod)*phase_array(iphs)%dia)
          norm_drag_poly_spec(iphs) = norm_drag_poly_spec(iphs)!/(vis**2.d0)
       end do
       norm_drag  = norm_drag/(3.d0*pi*vis*(ucharmod)*char_length)
       norm_drag_poly  = norm_drag_chem
       norm_drag_chem  = norm_drag_chem/(3.d0*pi*vis*(ucharmod)*char_length)
       norm_drag_poly  = norm_drag_poly!/(vis**2.d0)
    else
		do iphs = 1, nphases
			if (impose_grav.or.set_mpg) then
				norm_factor = (3.d0*pi*vis*(mixmeanslipmod+SMALL_NUMBER)*phase_array(iphs)%dia)
			else
				norm_factor = (3.d0*pi*vis*ucharmod*phase_array(iphs)%dia)
			endif

			norm_drag_spec(iphs) = norm_drag_spec(iphs)/norm_factor
			norm_drag_poly_spec(iphs) = norm_drag_chem_spec(iphs)
			norm_drag_chem_spec(iphs) = norm_drag_chem_spec(iphs)/norm_factor
			norm_drag_poly_spec(iphs) = norm_drag_poly_spec(iphs)!/(vis**2.d0)
		enddo
		norm_drag  = norm_drag/norm_factor
		norm_drag_poly  = norm_drag_chem
		norm_drag_chem  = norm_drag_chem/norm_factor
		norm_drag_poly  = norm_drag_poly!/(vis**2.d0)
    end if

    do iphs = 1, nphases
       IF(norm_drag_spec(iphs).gt.ZERO)then
          phase_array(iphs)%ferror = ABS(norm_drag_spec(iphs) - phase_array(iphs)%fold)/norm_drag_spec(iphs)
       ELSE
          phase_array(iphs)%ferror = ONE
       END IF
    end do
    IF(norm_drag.gt.ZERO)THEN
       ferror = ABS(norm_drag - fold)/norm_drag
    ELSE
       ferror = one
    END IF
    fold = norm_drag
    
    do iphs = 1, nphases
       phase_array(iphs)%fold = norm_drag_spec(iphs)
    end do

!!$       IF(Re.EQ.Zero.and.ReT.gt.zero)Then
!!$          norm_drag1 = (norm_drag1)/(6.d0*pi*vis*SQRT(gran_temp)*radbdy(1)*dx)
!!$       ELSE
!!$          norm_drag1 = (norm_drag1)/(6.d0*pi*vis*ucharmod*radbdy(1)*dx)
!!$          norm_drag2 = norm_drag2/(6.d0*pi*vis*ucharmod*radbdy(nbody)*dx)
!!$          norm_drag  = norm_drag/(3.d0*pi*vis*ucharmod*dia_phys)
!!$       END IF

    IF(FROM_POST) RETURN
    if(rks.eq.itrmax) then 
       !Rearrange the ferror_array array so that the last entry is flushed out
       tmp_ferror_array(1:nerr_steps) = ferror_array(1:nerr_steps)
       ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
       ferror_array(1) = ferror
!!$       PRINT*,'FERROR_A =', FERROR_ARRAY
       ferror_hist = SUM(ferror_array(1:nerr_steps))/nerr_steps
       do iphs = 1, nphases
          tmp_ferror_array(1:nerr_steps) = phase_array(iphs)%ferror_array(1:nerr_steps)
          phase_array(iphs)%ferror_array(2:nerr_steps) = tmp_ferror_array(1:nerr_steps-1)
          phase_array(iphs)%ferror_array(1) = phase_array(iphs)%ferror
!!$       PRINT*,'FERROR_A =', FERROR_ARRAY
          phase_array(iphs)%ferror_hist = SUM(phase_array(iphs)%ferror_array(1:nerr_steps))/nerr_steps
       end do
    end if
    
    
    if(I_AM_NODE_ZERO) WRITE(*,'(A25,4(2x,g17.8))') 'NORM DRAGS, FERRO&
         &R', norm_drag, ferror
    if(I_AM_NODE_ZERO)then
       if(.not.((TRIM(input_type).eq.'random').and.(TRIM(psd_type).eq.'psd')))WRITE(*,'(A25,4(2x,g17.8))') 'NORM DRAGS PHASES:', norm_drag_spec(1:nphases)
    end if
    if(I_AM_NODE_ZERO) then 
       IF(rks.eq.itrmax) then 
          !c_drag_st(1:nbody) = force(1:nbody,1)/(3.d0*pi*vis*dia_phys*ucharmod)
          WRITE(unitnormdragchem,'(500(2x, e20.12))') t/t_conv/(1-maxvolfrac), t/t_vis, t&
               &/t_diff, norm_drag_chem, norm_drag_chem_spec(1:nphases)
          
          !c_drag_st(1:nbody) = force(1:nbody,2)/(3.d0*pi*vis*dia_phys*ucharmod)
          
          !WRITE(unitdragtavg,'(15(2x,e20.12))')  t/t_conv/(1-maxvolfrac), t/t_vis, t/t_diff,(mean_drag_tavg(idim)/t,idim=1,ndim),(drag_tavg(1,idim)/t,idim=1,ndim), (drag_tavg(2,idim)/t,idim=1,ndim)
          
          !c_drag_st(1:nbody) = force(1:nbody,3)/(3.d0*pi*vis*dia_phys*ucharmod)
          
          WRITE(unitnormdrag,'(500(2x, e20.12))') t/t_conv/(1-maxvolfrac), t/t_vis, t&
               &/t_diff, norm_drag, norm_drag_spec(1:nphases), usmean_act(1:ndim), (ufmean(idim)/ucharmod,idim=1,ndim)

          
          WRITE(unitdrag_comps,'(20(2x, e20.12))') pres_drag, visc_drag, total_drag, total_drag_mpg, norm_drag, ABS(total_drag - total_drag_mpg)/(total_drag + SMALL_NUMBER)
          
          WRITE(unitnormdragpoly,'(500(2x, e20.12))') t/t_conv/(1-maxvolfrac), t/t_vis, t&
               &/t_diff, norm_drag_poly, norm_drag_poly_spec(1:nphases)
          
       end IF
    end if
    
    ! Checking for blown-up simulations
    if(norm_drag.gt.1E+06) THEN 
       if(I_AM_NODE_ZERO)then 
          WRITE(*,'(A,2x,g17.8,2x,A)') 'NORM DRAG', norm_drag,' is greater than 1E+06'
          WRITE(*,'(A)') 'STOPPING THIS CASE AFTER WRITING THE BLOW UP INDICATOR FILE'
          OPEN(2000, file=TRIM(RUN_NAME)//'_BLOWUP.dat', form='formatted')
          
          close(2000, status="keep")
       end if
       PARALLEL_FINISH()
       STOP
    end if

1000 FORMAT(18(E14.6,1x))

2000 continue


	if (I_AM_NODE_ZERO) then
		CALL CPU_TIME (CPU3) 
		bc_time_current = cpu3-cpu2 - dem_time_current
		bc_time = bc_time + bc_time_current
	endif


    RETURN

  END SUBROUTINE bcset

  SUBROUTINE calc_bnd_forcing(m,rks)
    IMPLICIT NONE
    Integer, Intent(in) :: m,rks
    Integer :: l

    REAL(prcn) ::  xl(ndim),xlo(ndim),xli(ndim)
    REAL(prcn) ::  xpb(ndim),xpo(ndim),xpi(ndim), force_tmp_vis(ndim)
    REAL(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, df(nbnd,ndim), CROSSP(ndim),linvel(ndim)
    LOGICAL :: velterm, pressterm
    INTEGER :: count, vcell(ndim)
#if PARALLEL
    REAL(prcn) :: xltemp
    INTEGER :: bndcount,bndtot, vcelltemp, pcelltemp,focus_point,focus_particle
    
    focus_point = -1
    focus_particle = -1
#endif
    !-----------------------------------------------------------------------
    !WRITE(*,*) 'IN FLOW BND'
    bcount = 0
    frombnd = .TRUE.
    fromrpr = .FALSE.
    count = 0
    
    BNDLOOP: DO l=1,nbnd

       rad = zero
       DO n=1,ndim

          !     Note: x co-ordinate of xc() is absolute center (center with
          !     respect to global origin) minus foffset
          !     Hence, is(1) is the currently the cell location w.r.t foffset

          xl(n)=xc(m,n)+bndarray(n,l)*radbdy(m)

          rad=rad+(bndarray(n,l)*radbdy(m))**2.0


          is(n)=INT(xl(n))
          
          ul(n)=zero
          nll(n)=zero
          onll(n)=zero

          ppll(n)=zero
          dfll(n)=zero
       ENDDO
       rad = DSQRT(rad)
       xpb(1) = xl(1)-0.5
       xpb(2:3)=xl(2:3)
       do n = 1, ndim
          if(xpb(n).lt.zero) then 
             pcell(n) = int(xpb(n)-1)
             !because of int of example -1.1 is -1, but we want -2. So,
             ! the adjustment 
          else 
             pcell(n) = int(xpb(n))
          end if
          if(xl(n).lt.zero) then
             vcell(n) = int(xl(n)-1)
          else 
             vcell(n) = int(xl(n))
          end if

       end do



#if PARALLEL
       vcelltemp = vcell(1)
       xltemp  = xl(1)
       if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
          PRINT*,' INNER REV PT = ', xli(1),vcelltemp, myid, m
          PRINT*,' EXTERNAL REV PT = ', xlo(1), myid,m
       end if
       
       if(.not.CELL_IN_PROC(vcelltemp))then
          WEST_PERIODIC_IMAGE(vcell(1),vcelltemp,xl(1),xltemp)
          EAST_PERIODIC_IMAGE_MOD(vcell(1),vcelltemp, xl(1), xltemp)
          if(.not.CELL_IN_PROC(vcelltemp))then
             if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
                PRINT*,' INNER REVERSAL PT = ', xl(1),vcelltemp, myid, m
                
             end if

             goto 2600
          end if
       end if

       velterm = CELL_IN_VEL_GRID(vcelltemp)
       
       vcell(1) = vcelltemp
       xl(1) = xltemp

       pressterm = .TRUE.
       pcelltemp = pcell(1)
       xltemp = xpb(1)
       if(.not.CELL_IN_PRESS_GRID(pcelltemp))then
          WEST_PERIODIC_IMAGE_PRES(pcell(1),pcelltemp, xpb(1), xltemp)
          EAST_PERIODIC_IMAGE_PRES(pcell(1),pcelltemp, xpb(1), xltemp)
          if(l.eq.FOCUS_POINT)then
             PRINT*,' INNER PRES PT = ', xpi(1),pcelltemp, myid, m
          end if
          pressterm = CELL_IN_PRESS_GRID(pcelltemp)
       end if
       pcell(1) = pcelltemp
       xpb(1) = xltemp


#else
       velterm = .TRUE.
       pressterm = .TRUE.
#endif
       if(pressterm)then

          pl=zero
          call interpolate_pdata(pcell,xpb,ppll,pl,l)

       end if

       call interpolate_udata(vcell,xl,ib&
            &,ie,jb,je,kb,ke,ul,nll,onll,dfll, 1,m, l, onew) 
       
       CALL CROSS_PRODUCT3(CROSSP(1:3),ANGV(M,1:3,1),SNORM(1:3))
       linvel(1:ndim) = velbdy(m,1:ndim) + CROSSP(1:ndim)*RADBDY(m)*dx 

       
       DO n=1,3
          force_tmp(n) = zero
          if(velterm)force_tmp(n) = cf*(-linvel(n)+ ul(n))
          
          force_tmp(n) = force_tmp(n) - coef(rks,3)*nll(n)-coef(rks,4)&
               &*onll(n)
          force_tmp(n) = force_tmp(n) + (coef(rks,1)+coef(rks,2))*(ppll(n)-dfll(n))
          
          force_tmp(n) = force_tmp(n)*da(1)*drm*dx

          sumforcepoint(n)=sumforcepoint(n)+force_tmp(n)
!!$#if PARALLEL          
!!$          IF(set_umean) force_loc(m,n)= force_loc(m,n) - force_tmp(n)
!!$#else
!!$          IF(set_umean) force(m,n)= force(m,n) - force_tmp(n)
!!$#endif
       ENDDO

              
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

                DO n=1,ndim
                   tmppa = weightp(i,j,k)*force_tmp(n)
                   fr(ii,jj,kk,n) = fr(ii,jj,kk,n) + tmppa/(dx**3.d0)
                   if(fluid_atijk(ii,jj,kk)) then
                      if(include_frmeanfluid)  frmeanfluid(n) = frmeanfluid(n) + tmppa/(dx**3.d0)
                   end if
#if PARALLEL
                   frmeanloc(n) = frmeanloc(n) + tmppa/(dx**3.d0)
#else
                   frmean(n) = frmean(n) + tmppa/(dx**3.d0)
#endif
                                      
#if 0
                   if((GLOBAL_INDEX(ii).eq.1))then
                      if((pcell(1).eq.my/2))then
                         fromleft(jj,kk,n) = fromleft(jj,kk,n) + tmppa/(dx**3.d0)
                      else if(pcell(1).eq.my/2+1)then
                         fromright(jj,kk,n) = fromright(jj,kk,n) + tmppa/(dx**3.d0)
                      end if
                   endif
#endif

                      
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          

#if PARALLEL
2600      continue
#endif          

       END DO BNDLOOP

    frombnd = .FALSE.
  END SUBROUTINE calc_bnd_forcing

  SUBROUTINE calc_inner_reversal_forcing(m,rks)
    IMPLICIT NONE
    Integer, Intent(in) :: m,rks
    Integer :: l, count_fl, count_so
    INTEGER :: pcelli(ndim), vcelli(ndim), vcello(ndim)
    REAL(prcn) ::  xl(ndim),xlo(ndim),xli(ndim), force_fl(ndim), force_dist(ndim)
    REAL(prcn) ::  xpb(ndim),xpo(ndim),xpi(ndim), force_tmp_vis(ndim)
    REAL(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, xltemp,CROSSP(ndim),linvel(ndim)
    LOGICAL :: velterm, pressterm
#if PARALLEL
    INTEGER :: rprcountloc, rprcount,rprcom,rprcomloc, FOCUS_POINT, VELGL, vcelltemp, pcelltemp, FOCUS_PARTICLE
#endif

    !WRITE(*,*) 'IN FLOW REVERSAL, nrpr :', nrpr

    frombnd = .FALSE.
    fromrpr = .TRUE.
#if PARALLEL
    FOCUS_POINT = -1
    FOCUS_PARTICLE = -1
#endif

!open (unit=1, file="vel.dat", action="write",position="append", status="old")
!write (1,*) "zone"
!
!open (unit=2, file="vel2.dat", action="write",position="append", status="old")
!write (2,*) "zone"
!
!open (unit=3, file="force.dat", action="write",position="append", status="old")
!write (3,*) "zone"

!write (*,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
!write (*,*) "rbidy  = ", radibdy(m)
!write (*,*) "rbdy   = ", radbdy(m)
!write (*,*) "rbody  = ", radobdy(m)
!write (*,*) "rbo2dy = ", rado2bdy(m)
!write (*,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

    DO l=1,nrpr
       if(.NOT.PART_ARRAY(M)%if_rev(L)) GOTO 666
       
       rad = zero
       DO n=1,ndim
          
          !     location of internal points
          
          xli(n)=xc(m,n)+ bndarray(n,l)*radibdy(m)
          
          iii(n)=INT(xli(n))
          uli(n)=zero

          !     location of external points

          xlo(n)=xc(m,n)+ bndarray(n,l)*radobdy(m)


          io(n)=INT(xlo(n))
          ulo(n)=zero

          rad=rad+(bndarray(n,l)*radibdy(m))**2.0

          onll(n)=zero
          nll(n)=zero
          ppll(n)=zero
          dfll(n)=zero

       ENDDO

       plb = zero
       plo = zero
       pli = zero
       rad=dsqrt(rad)
       DO n=1,ndim                       
          snorm(n)=(bndarray(n,l)*radibdy(m))/rad
       ENDDO
       xpi(1) = (xli(1)-0.5)
       xpi(2:3) = xli(2:3)

       do n = 1, ndim
          if(xlo(n).lt.zero) then 
             vcello(n) = int(xlo(n)-1)
          else 
             vcello(n) = int(xlo(n))
          end if
          if(xpi(n).lt.zero) then 
             pcelli(n) = int(xpi(n)-1)
          else 
             pcelli(n) = int(xpi(n))
          end if
          if(xli(n).lt.zero) then 
             vcelli(n) = int(xli(n)-1)
          else 
             vcelli(n) = int(xli(n))
          end if
          
       end do
#if PARALLEL
       vcelltemp = vcelli(1)
       xltemp  = xli(1)
       if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
          PRINT*,' INNER REV PT = ', xli(1),vcelltemp, myid, m
          PRINT*,' EXTERNAL REV PT = ', xlo(1), myid,m
       end if
       
       if(.not.CELL_IN_PROC(vcelltemp))then
          WEST_PERIODIC_IMAGE(vcelli(1),vcelltemp,xli(1),xltemp)
          EAST_PERIODIC_IMAGE_MOD(vcelli(1),vcelltemp, xli(1), xltemp)
          if(.not.CELL_IN_PROC(vcelltemp))then
             if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
                PRINT*,' INNER REVERSAL PT = ', xli(1),vcelltemp, myid, m
                
                !PARALLEL_FINISH()
                !STOP
             end if
             !             PRINT*,' XLTEMP = ', myid, xltemp, xli(1), l
             goto 666
          end if
       end if
       if(EAST_NO_MANS_LAND(vcelli(1)).or.EAST_NO_MANS_LAND(vcelltemp)) then 
          velterm = .not.CONCAVE(xli,1,m)
       else if(WEST_NO_MANS_LAND(vcelli(1)).or.WEST_NO_MANS_LAND(vcelltemp))then
          velterm = CONCAVE(xli,1,m)
       else
          velterm = .TRUE.
       end if
       vcelli(1) = vcelltemp
       xli(1) = xltemp

       pressterm = .TRUE.
       pcelltemp = pcelli(1)
       xltemp = xpi(1)
       if(.not.CELL_IN_PRESS_GRID(pcelltemp))then
          WEST_PERIODIC_IMAGE_PRES(pcelli(1),pcelltemp, xpi(1), xltemp)
          EAST_PERIODIC_IMAGE_PRES(pcelli(1),pcelltemp, xpi(1), xltemp)
          if(l.eq.FOCUS_POINT)then
             PRINT*,' INNER PRES PT = ', xpi(1),pcelltemp, myid, m
          end if
          pressterm = CELL_IN_PRESS_GRID(pcelltemp)
       end if
       pcelli(1) = pcelltemp
       xpi(1) = xltemp

       if(velterm)then
          vcelltemp = vcello(1)
          xltemp = xlo(1)
          if(.not.RPR_CELL_IN_PROC(vcelltemp))then
             WEST_PERIODIC_IMAGE(vcello(1),vcelltemp, xlo(1),xltemp)
             EAST_PERIODIC_IMAGE_MOD(vcello(1),vcelltemp,xlo(1),xltemp)
             if(.not.RPR_CELL_IN_PROC(vcelltemp))then
                if(I_AM_NODE_ZERO)then
                   if(vcelltemp.eq.mxf-3)then
                      vcelltemp = vcelltemp-(mxf-1)+1
                      xltemp = xltemp-(mxf-1)
                   endif
                else
                   PRINT*,' ERROR WITH EXTERNAL POINT IN THIS PROCESSOR : ', myid, m, l, xlo(1), vcelltemp,vcello(1),xli(1)
                end if
             end if
          end if
          vcello(1) = vcelltemp
          xlo(1) = xltemp
       end if

       if(l.eq.FOCUS_POINT.and.m.eq.FOCUS_PARTICLE)then
          PRINT*,' PRESTERM = ', l,pressterm, myid,m
          PRINT*,' VELTERM = ', l,velterm, myid, m
          PARALLEL_FINISH()
          STOP
       endif
#else
       velterm = .TRUE.
       pressterm = .TRUE.
#endif

       if(velterm)then
          call interpolate_udata(vcello,xlo,ib&
               &,ie,jb,je,kb,ke,ulo,nll,onll,dfll, 0, m, l, onew) 
          
          unmago=zero
          DO n=1,ndim
             unmago = unmago + snorm(n)*ulo(n)
          ENDDO
       end if
       !-----------------------------------------------------------------------
       !     calculate internal local velocities
       !-------------------------------------------------------


       if(pressterm)then
          call interpolate_pdata(pcelli,xpi,ppll,pl,l)
          
       end if

       call interpolate_udata(vcelli,xli,ib&
            &,ie,jb,je,kb,ke,uli,nll,onll,dfll, 1, m, l, onew) 
       

!write (1,"(1I8,3D15.7)") l,(uli(:)+ulo(:))/umeanslip
!write (*,"(A,3D15.7)") "ro", ulo(:)/umeanslip
!write (*,*)

       if(velterm)then
          unmagi=zero
          DO n=1,ndim
             unmagi = unmagi + snorm(n)*ulo(n)
          ENDDO

          unmag=zero
          ubnmag = zero ! Velocity of the body (dot) normal vector
          CALL CROSS_PRODUCT3(CROSSP(1:3),ANGV(M,1:3,1),SNORM(1:3))
          linvel(1:ndim) = velbdy(m,1:ndim) + CROSSP(1:ndim)*RADBDY(m)*dx 

          do d=1,ndim
             ucar(d)=zero
             unmag=unmag+snorm(d)*ulo(d)
             ubnmag = ubnmag + snorm(d)*linvel(d)
          enddo

          !-----------------------------------------------------------------------
          !     set internal forcing to reverse external velocity
          !     reverse both tangential velocity components
          !     zero radial component
          !     scale tangential velocities by ratio of radii
          !-----------------------------------------------------------------------

          do d=1,ndim
             unorm(d)=snorm(d)*unmag
             utang(d)=ulo(d)-unorm(d)
             ubnorm(d) = snorm(d)*ubnmag  ! velocity of the body in the normal direction
!!$          ubtang(d) = velbdy(m,d)-ubnorm(d) ! velocity of the body in the tangential direction
             ubtang(d) = linvel(d)-ubnorm(d) ! velocity of the body in the tangential direction
             utangi(d) = ubtang(d)*(radobdy(m)-radibdy(m))-utang(d)*(radbdy(m)-radibdy(m))
             utangi(d) = utangi(d)/(radobdy(m)-radbdy(m))

             ! Find tangential velocity at the internal reversal point so that the no slip condition is satisfied

             IF(revernorm.EQ.1)THEN

                unormi(d) = ubnorm(d)*(radobdy(m)-radibdy(m))-unorm(d)*(radbdy(m)-radibdy(m))
                unormi(d) = unormi(d)/(radobdy(m)-radbdy(m))
             ELSE
                unormi(d) = zero
             ENDIF
             ucar(d) = utangi(d) + unormi(d)
          ENDDO
       end if

       DO n=1,ndim
          force_tmp(n) = zero

          if(velterm)force_tmp(n)=cf*(-ucar(n)+uli(n))

!write (2,"(1I8,3D15.7)") l,(ucar(:)+ulo(:))/umeanslip
!write (*,"(A,3D15.7)") "ri", ucar(:)/umeanslip
!write (*,"(A,3D15.7)") "ro", ulo(:)/umeanslip
!read(*,*)
          
          force_tmp(n) = force_tmp(n) - coef(rks,3)*nll(n)-coef(rks,4)&
               &*onll(n)

          force_tmp(n) = force_tmp(n) + (coef(rks,1)+coef(rks,2))*(ppll(n)-dfll(n))

          force_tmp(n) = force_tmp(n)*da(2)*drm*dx

          force_tmp_vis(n) = force_tmp_vis(n)*da(1)*drm*dx
          
          sumforcepoint(n)=sumforcepoint(n)+force_tmp(n)
!!$#if PARALLEL
!!$          if(set_umean) force_loc(m,n)= force_loc(m,n) - force_tmp(n)
!!$#else
!!$          if(set_umean) force(m,n)= force(m,n) - force_tmp(n)
!!$#endif
       ENDDO

!!$       if(pcell(1).eq.mx1.or.pcell(1).eq.0)PRINT*,' DING = ', myid, l, force_tmp(1)
       count_fl = 0
       force_fl = zero
       
       do k = 1, onew
          kk = kb+k-1
          if(kk.lt.1) kk = mz+kk
          if(kk.gt.mz) kk = kk-mz 
          
          do j = 1, onew
             jj = jb+j-1
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             
             do i = 1, onew
                ii = ib+i-1
                
#if !PARALLEL
                if(ii.lt.1) ii = mxf+ii-1
                if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
                LOCAL_INDEX(ii)
                if(fluid_atijk(ii,jj,kk)) then 
                   count_fl = count_fl+1
                   DO n=1,ndim
                      tmppa = weightp(i,j,k)*force_tmp(n)
                      force_fl(n) = force_fl(n) + tmppa/(dx**3.d0)
                   end DO
                end if
             end do
          end do
       end do
       count_so = onew*onew*onew - count_fl
       if(count_fl.gt.0)then
          !Write(*,*) ' count_fl gt 0....node = ', myid, count_fl, l, ib,force_fl(1)
       end if
       force_dist(:) = force_fl(:)/real(count_so)


      
       do k = 1, onew
          kk = kb+k-1
          if(kk.lt.1) kk = mz+kk
          if(kk.gt.mz) kk = kk-mz 
          
          do j = 1, onew
             jj = jb+j-1
             if(jj.lt.1) jj = my+jj
             if(jj.gt.my) jj = jj-my
             
             do i = 1, onew
                
                ii = ib+i-1
#if !PARALLEL
                if(ii.lt.1) ii = mxf+ii-1
                if(ii.gt.mxf-1) ii = ii-mxf +1
#endif
                LOCAL_INDEX(ii)
                DO n=1,ndim                   
                   tmppa = weightp(i,j,k)*force_tmp(n)
                   if(include_frmeanfluid)  then 
                      if(fluid_atijk(ii,jj,kk)) then
#if PARALLEL
                         frmeanfluidloc(n) = frmeanfluidloc(n) + tmppa/(dx**3.d0)
#else
                         frmeanfluid(n) = frmeanfluid(n) + tmppa/(dx**3.d0)
#endif
                      end if
                      fr(ii,jj,kk,n) = fr(ii,jj,kk,n) + tmppa/(dx**3.d0)
                   else
                      if(.not.fluid_atijk(ii,jj,kk)) then 
                         fr(ii,jj,kk,n) = fr(ii,jj,kk,n) + tmppa/(dx**3.d0)&
                              & + force_dist(n)
                      end if
                   end if
                   
#if PARALLEL
                   frmeanloc(n) = frmeanloc(n) + tmppa/(dx**3.d0)
#else
                   frmean(n) = frmean(n) + tmppa/(dx**3.d0)
#endif
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       !------------------------------------------
       !     Close loop over all reversal points
666    CONTINUE             

!write (3,"(1I8,3D15.7)") l,force_tmp(:)


88  ENDDO                ! loop over reversal points
    fromrpr = .FALSE.

!	close (1)
!	close (2)
!	close (3)
!	read(*,*)

  END SUBROUTINE calc_inner_reversal_forcing

  SUBROUTINE calc_pres_visc_drag(m,rks)
    USE bcsetarrays, ONLY :  omega => fr
    IMPLICIT NONE
    Integer, Intent(in) :: m,rks
    Integer :: l, pcelltemp,vcelltemp, pcellb(ndim), vcellb(ndim)

    REAL(prcn) ::  xl(ndim),xlo(ndim),xli(ndim)
    REAL(prcn) ::  xpb(ndim),xpo(ndim),xpi(ndim), force_tmp_vis(ndim)
    REAL(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, df(nbnd,ndim)
    REAL(prcn) :: tempforce(ndim), crossp(ndim), xltemp, xptemp
#if PARALLEL
    INTEGER :: FOCUS_POINT
    FOCUS_POINT = -1
#endif

    !-----------------------------------------------------------------------
    !WRITE(*,*) 'IN FLOW BND'
    bcount = 0
    
    DO l=1,nbnd

       rad = zero
       DO n=1,ndim
          
          xl(n)=xc(m,n)+ bndarray(n,l)*radbdy(m)
          
          rad=rad+(bndarray(n,l)*radbdy(m))**2.0
          
          
          is(n)=INT(xl(n))
          
          ul(n)=zero
          
          ppll(n)=zero
          
       ENDDO

       rad = DSQRT(rad)

       xpb(1) = xl(1)-0.5
       xpb(2:3)=xl(2:3)
       do n = 1, ndim
          if(xpb(n).lt.zero) then 
             pcellb(n) = int(xpb(n)-1)
             !because of int of example -1.1 is -1, but we want -2. So,
             ! the adjustment 
          else 
             pcellb(n) = int(xpb(n))
          end if
          if(xl(n).lt.zero) then 
             vcellb(n) = int(xl(n)-1)
          else 
             vcellb(n) = int(xl(n))
          end if
       end do

#if PARALLEL

       xltemp = xl(1)
       xptemp = xpb(1)
       pcelltemp = pcellb(1)
       vcelltemp = vcellb(1)
       if(l.eq.FOCUS_POINT)then
          PRINT*,' xl = ', myid, xltemp, xptemp, pcelltemp
       end if
       if(.not.CELL_IN_VEL_GRID(vcelltemp))then
          WEST_PERIODIC_IMAGE(vcellb(1),vcelltemp,xl(1),xltemp)
          WEST_PERIODIC_IMAGE(pcellb(1),pcelltemp,xpb(1),xptemp)
          EAST_PERIODIC_IMAGE(vcellb(1),vcelltemp,xl(1),xltemp)
          EAST_PERIODIC_IMAGE_PRES(pcellb(1),pcelltemp,xpb(1),xptemp)
          if(l.eq.FOCUS_POINT)then
             PRINT*,' xl IMAGES = ', myid,xltemp, xptemp, pcelltemp
          end if

          if(.not.CELL_IN_VEL_GRID(vcelltemp)) goto 777
       end if
       vcellb(1) = vcelltemp
       pcellb(1) = pcelltemp
       xl(1) = xltemp
       xpb(1) = xptemp
#endif

       pl=zero
       
       call interpolate_pdata(pcellb,xpb,ppll,pl,l)

    
       call interpolate_udata(vcellb,xl,ib&
            &,ie,jb,je,kb,ke,ul,nll,onll,dfll, 0,m, l, onew) 
       
       vort(:) = zero 
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
                DO n=1,ndim
                   vort(n)=vort(n)+ weightp(i,j,k)*omega(ii,jj,kk,n) 
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       
       !VISC_CALC: IF(FROM_POST) then
!!$       
       df(l,1)=vort(3)*cd(2,l)
       df(l,1)=df(l,1)-vort(2)*cd(3,l)
       
       df(l,2)=vort(1)*cd(3,l)
       df(l,2)=df(l,2)-vort(3)*cd(1,l)

       df(l,3)=vort(2)*cd(1,l)
       df(l,3)=df(l,3)-vort(1)*cd(2,l)

       df(l,1)=-df(l,1)
       df(l,2)=-df(l,2)
       df(l,3)=-df(l,3)

       !---------------------------------------------------------------
       !     calculate the viscous and pressure components separately
       !---------------------------------------------------------------
       
       DO d=1,ndim,1
#if PARALLEL
          presloc(m,d)= presloc(m,d)-pl*cd(d,l)*da(1)
          pres_totalloc(d) = pres_totalloc(d) + pl*cd(d,l)*da(1)

          viscloc(m,d)=viscloc(m,d)+df(l,d)*vis*da(1)
          visc_totalloc(d) = visc_totalloc(d) + df(l,d)*vis*da(1)
#else
          pres(m,d)=pres(m,d)-pl*cd(d,l)*da(1)
          pres_total(d) = pres_total(d) + pl*cd(d,l)*da(1)
          visc(m,d)=visc(m,d)+df(l,d)*vis*da(1)
          visc_total(d) = visc_total(d) + df(l,d)*vis*da(1)
#endif
       
          tempforce(d) = df(l,d)*vis*da(1)!+ df(l,d)*vis*da(1)
       ENDDO
       CALL CROSS_PRODUCT3(CROSSP(1:3),bndarray(1:3,l),TEMPFORCE(1:3))
#if PARALLEL
       TORQLOC(m,1:3) = TORQLOC(m,1:3) + CROSSP(1:3)*RADBDY(m)*dx
#else
       TORQ(m,1:3) = TORQ(m,1:3) + CROSSP(1:3)*RADBDY(m)*dx
#endif
       !---------------------------------------------
       !     close loop over all boundary points
       !---------------------------------------------
       !end IF VISC_CALC
#if PARALLEL
777    continue      
#endif
    ENDDO


  END SUBROUTINE calc_pres_visc_drag


  SUBROUTINE compute_omega
    USE bcsetarrays, ONLY :  omega => fr
    IMPLICIT NONE
    Integer i,j,k
    
#if PARALLEL
    Do i = 0, nx+1
#else
    DO i=1,nx+1 !mxf
#endif
       ioffset=i+foffset
       DO k=1,mz
          DO j=1,my2
             ff1(j,k)=wy(j)*u(ioffset,j,k,1) ! dudy
          ENDDO
       ENDDO
       CALL ff2cr(ff1,dudy(:,:))
       
       DO k=1,mz
          DO j=1,my
#if PARALLEL
             omega(i,j,k,3)=(ur(i+1,j,k,2)-ur(i-1,j,k,2))/(two*dx)-dudy(j,k)
#else
             if(i.eq.1) THEN 
               
                IF(xperiodic)  THEN 
                   omega(i,j,k,3)=(ur(i+1,j,k,2)-ur(mxf-1,j,k,2))/(two*dx)-dudy(j,k)
                ELSE
                   omega(i,j,k,3)=(ur(i+1,j,k,2)-ur(i,j,k,2))/(dx)-dudy(j,k)
                end IF

             ELSE IF(i.eq.mxf) THEN
                IF(xperiodic)  THEN 
                   omega(i,j,k,3)=(ur(2,j,k,2)-ur(mxf-1,j,k,2))/(two*dx)-dudy(j,k)
                ELSE
                   omega(i,j,k,3)=(ur(i,j,k,2)-ur(i-1,j,k,2))/(dx)-dudy(j,k)
                end IF

             ELSE 
                IF(xperiodic)THEN
                   omega(i,j,k,3)=(ur(i+1,j,k,2)-ur(i-1,j,k,2))/(two*dx)-dudy(j,k)
                ELSE
                   omega(i,j,k,3)=(ur(i+1,j,k,2)-ur(i,j,k,2))/(dx)-dudy(j,k)
                ENDIF
             ENDIF
#endif
             
          ENDDO
       ENDDO
    ENDDO

#if PARALLEL
    DO i = 0, nx+1
#else
    DO i=1,nx+1 !mxf
#endif
       ioffset=i+foffset
       DO k=1,mz
          DO j=1,my2
             ff2(j,k)=wz(k)*u(ioffset,j,k,1) ! dudz
          ENDDO
       ENDDO
       CALL ff2cr(ff2,dudz(:,:))
       
       DO k=1,mz
          DO j=1,my
#if PARALLEL
             omega(i,j,k,2)=dudz(j,k)- (ur(i+1,j,k,3)-ur(i-1,j,k,3))/(two*dx)
#else
             if(i.eq.1) THEN 
                IF(xperiodic)  THEN 
                   omega(i,j,k,2)=dudz(j,k)- (ur(i+1,j,k,3)-ur(mxf-1,j,k,3))/(two*dx)
                ELSE
                   omega(i,j,k,2)=dudz(j,k)- (ur(i+1,j,k,3)-ur(i,j,k,3))/dx
                end IF

             ELSE IF(i.eq.mxf) THEN 
                IF(xperiodic)  THEN 
                   omega(i,j,k,2)=dudz(j,k)- (ur(2,j,k,3)-ur(mxf-1,j,k,3))/(two*dx)
                ELSE
                   omega(i,j,k,2)=dudz(j,k)- (ur(i,j,k,3)-ur(i-1,j,k,3))/dx
                end IF

             ELSE
                IF(xperiodic)THEN
                   omega(i,j,k,2)=dudz(j,k)- (ur(i+1,j,k,3)-ur(i-1,j,k,3))/(two*dx)
                ELSE
                   omega(i,j,k,2)=dudz(j,k)- (ur(i+1,j,k,3)-ur(i,j,k,3))/dx
                ENDIF
             ENDIF
#endif
          ENDDO
       ENDDO
    ENDDO!i = 1, mxf

#if PARALLEL
    Do i = 0, nx+1
#else
    DO i=1,nx+1 !mxf
#endif
       ioffset=i+foffset
       DO k=1,mz
          DO j=1,my2
             ff1(j,k)=wz(k)*u(ioffset,j,k,2) ! dvdz
             ff2(j,k)=wy(j)*u(ioffset,j,k,3) ! dwdy
          ENDDO
       ENDDO

       CALL ff2cr(ff1,dvdz(:,:))
       CALL ff2cr(ff2,dwdy(:,:))
       
       DO k=1,mz 
          DO j=1,my
             omega(i,j,k,1)=dwdy(j,k)-dvdz(j,k) 
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE compute_omega

  SUBROUTINE calc_visc
    IMPLICIT NONE
    INTEGER :: n, i, j, k, ioffset
    DO n=1,ndim
#if PARALLEL
       diffn(0,1:my,1:mz,n) = zero
       DO i=1,nx
#else
          DO i = 1, nx+1
#endif
             ioffset=foffset+i
             DO k=1,mz
                DO j=1,my2

                   if(xperiodic) then
#if PARALLEL
                      ff1(j,k)=(1./dx2)*(u(ioffset-1,j,k,n)-2.*u(ioffset,j&
                           &,k,n)+u(ioffset+1,j,k,n)) 
#else                   
                      if(ioffset.eq.1) THEN
                         ff1(j,k)=(1./dx2)*(u(nx,j,k,n)-2.*u(ioffset,j&
                              &,k,n)+u(ioffset+1,j,k,n)) 
                      else if(ioffset.eq.nx+1) THEN
                         ff1(j,k)=(1./dx2)*(u(ioffset-1,j,k,n)-2.*u(ioffset,j&
                              &,k,n)+u(2,j,k,n))
                         
                      else

                         ff1(j,k)=(1./dx2)*(u(ioffset-1,j,k,n)-2.*u(ioffset,j&
                              &,k,n)+u(ioffset+1,j,k,n)) 
                      endif
#endif
                   else
                      ff1(j,k)=(1./dx2)*(u(ioffset-1,j,k,n)-2.*u(ioffset,j&
                           &,k,n)+u(ioffset+1,j,k,n)) 
                   endif
                   ff1(j,k)=ff1(j,k)-w2(j,k)*u(ioffset,j,k,n)

                   ff1(j,k)=ff1(j,k)*vis
                ENDDO
             ENDDO
             
             CALL ff2cr(ff1,diffn(i,:,:,n))
          ENDDO
#if PARALLEL
          diffn(nx+1,1:my,1:mz,n) = zero
#endif
       ENDDO
     END SUBROUTINE calc_visc
     
     SUBROUTINE calc_pgrad
       IMPLICIT NONE
       INTEGER :: i,j,k,n
       
       !     To Check the contribution of pressure terms in the forcing (at the grid locations)
       !     Compute pressure gradient and store at velocity grid points
       
       !c     To Check the contribution of pressure terms in the forcing
       !c     Compute pressure gradient and store at --pressure-- grid points
       !c     earlier was being stored at velocity grid points, and implicit
       !c     smoothing of the pressure field was being done in the x-direction
       
       !    goto 9901
       
       if(xperiodic.and.nlcons)then

          DO i=1,nx
             DO k=1,mz
                DO j = 1, my2 
#if PARALLEL
                   ff1(j,k) = (p(i+1,j,k)-p(i-1,j,k))/(two*dx)
#else
                   if(i.eq.1) then
                      ff1(j,k) = (p(i+1,j,k)-p(nx,j,k))/(2.*dx) 
                      
                   else if(i.eq.nx)then
                      ff1(j,k) = (p(1,j,k)-p(i-1,j,k))/(2.*dx) 
                   ELSE

                      ff1(j,k) = (p(i+1,j,k)-p(i-1,j,k))/(two*dx)
                   end if
#endif

                   ff2(j,k)=p(i,j,k)*wy(j)  !starts at foffset+1
                   ff3(j,k)=p(i,j,k)*wz(k)
                end DO
             end DO

             CALL ff2cr(ff1(:,:),ppr(i,:,:,1))
             CALL ff2cr(ff2(:,:),ppr(i,:,:,2))
             CALL ff2cr(ff3(:,:),ppr(i,:,:,3))
             
             do k = 1, mz
                do j = 1, my
                   ppr(i,j,k,1) = ppr(i,j,k,1) + mpg(1)!total_mean_forcing(1)/real((mxf-1), prcn)
                   ppr(i,j,k,2) = ppr(i,j,k,2) + mpg(2)!total_mean_forcing(2)/real((mxf-1), prcn)
                   ppr(i,j,k,3) = ppr(i,j,k,3) + mpg(3)!total_mean_forcing(3)/real((mxf-1), prcn)
                   
                   !   pr(i,j,k) = ppr(i,j,k,1)*(i-1)*dx + ppr(i,j,k,2)*(j-1)*dy + ppr(i,j,k,3)*(k-1)*dz
                end DO
             end do
             if(i.eq.nx) then 
                do n = 1, ndim
                   VECSENDRECV(ppr(i,1,1,n),1,twodrslice,toproc,1,ppr(0,1,1,n),1,fromproc,1,decomp_group,status)
                end do
             end if
          end DO
#if PARALLEL
          do n = 1,ndim
             ppr(nx+1,:,:,n) = zero
          end do
#endif

       elseif(.not.xperiodic.and.nlcons) then
          DO i=1,mxf
             ioffset=foffset+i-1
             DO k=1,mz
                DO j=1,my2
                   ff1(j,k) = (p(ioffset+1,j,k)-p(ioffset-1,j,k))/(2.*dx) 
                   ff2(j,k)=p(ioffset,j,k)*wy(j)  !starts at foffset+1
                   ff3(j,k)=p(ioffset,j,k)*wz(k)
                ENDDO
             ENDDO
             CALL ff2cr(ff1,ppr(i,:,:,1))
             CALL ff2cr(ff2,ppr(i,:,:,2))
             CALL ff2cr(ff3,ppr(i,:,:,3))
          ENDDO
       endif
       
     END SUBROUTINE calc_pgrad

     
     SUBROUTINE calc_local_pres_visc_drag_plane
       USE bcsetarrays, ONLY :  omega => fr
       IMPLICIT NONE
       LOGICAL, SAVE :: routine_called = .false. 
       Integer :: m
       Integer :: l, iphi, itheta, unitno, unitno1, unitno2, unitno3, unitno4, unitno5
       REAL(prcn) ::  xl(ndim),xlo(ndim),xli(ndim)
       REAL(prcn) ::  xpb(ndim),xpo(ndim),xpi(ndim), force_tmp_vis(ndim)
       REAL(prcn) ::  ul(ndim),ulo(ndim),uli(ndim), tmppa, df(nbnd,ndim)
       REAL(prcn) :: tempforce(ndim), crossp(ndim), pres_loc, visc_loc,&
            & pres_bnd(3), visc_bnd(3), force_bdy_tmp(3), force_bdy_mag,&
            & axial_direction(3), thetaang, theta
       REAL(prcn) :: cphi, phiang, xcor, xsmin, xsmax, xmin, xmax, dx_x,&
            & nx, ny,nz, rad_proj, dtheta, normal(3), total_Fvis(nbody&
            &,2), total_Pres(nbody,2), areabdy(nbody), norm_factor2,&
            & avg_Fvis(2), avg_Pres(2), Fvis_theta_avg(count_theta),&
            & Fvis_phi_avg(count_phi),  Pres_phi_avg(count_phi),&
            & Pres_theta_avg(count_theta), conf1, conf2
       INTEGER :: norm_phi(count_phi), norm_theta(count_theta)

       norm_factor2 = (3.d0*pi*vis*(meanslipmod+SMALL_NUMBER)*dia_phys) 
       total_Fvis = zero 
       total_pres = zero 

       WRITE(*,*) 'count_phi = ', count_phi, count_theta
       axial_direction(1:3) = uchar(1:3)/ucharmod
       unitno = getnewunit(minunitno, maxunitno)

       OPEN(unitno,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_PHI_THETAZERO.dat', form="form&
            &atted",status="unknown") 
       unitno1 = getnewunit(minunitno, maxunitno)
       OPEN(unitno1,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_THETA_PHIZERO.dat', form="form&
            &atted",status="unknown") 

       unitno2 = getnewunit(minunitno, maxunitno)
       OPEN(unitno2,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_PHI_NBODY.dat', form="form&
            &atted",status="unknown") 

       unitno3 = getnewunit(minunitno, maxunitno)
       OPEN(unitno3,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_THETA_NBODY.dat', form="fo&
            &rmatted",status="unknown")  


       unitno4 = getnewunit(minunitno, maxunitno)
       OPEN(unitno4,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_PHI_AVG.dat', form="fo&
            &rmatted",status="unknown") 

       unitno5 = getnewunit(minunitno, maxunitno)
       OPEN(unitno5,FILE=TRIM(RUN_NAME)//'_FVIS_PRES_THETA_AVG.dat', form="fo&
            &rmatted",status="unknown")  

       dtheta  = twopi/real(count_theta-1,prcn)

       NULLIFY(bndarray)
       bndarray => phase_array(1)%bndpts

       xsmin = MINVAL(bndarray(1,1:nbnd))
       xsmax = MAXVAL(bndarray(1,1:nbnd))

       IF(.not.routine_called) then 
          ALLOCATE(Fvis_theta(nbody, count_theta), Fvis_phi(nbody,&
               & count_phi), Pres_theta(nbody, count_theta), Pres_phi(nbody, count_phi))
          routine_called = .true.
       end IF
       L = 1
       force_bdy_tmp = zero 
       !-----------------------------------------------------------------------
       WRITE(*,*) 'CALCULATING FVIS AND PRES DRAG COMPONENTS ALONG THETA &
            &AND PHI' 

       BODYLOOP: DO m = 1, nbody !Loop over bodies  
          areabdy(m) = 4.*pi*(radbdy(m)*dx)**2
          norm_phi = zero 
          norm_theta = zero 
          Fvis_phi(m,:) = zero 
          Fvis_theta(m,:) = zero 
          Pres_phi(m,:) = zero 
          Pres_theta(m,:) = zero 

          WRITE(unitno,*)'Zone'
          WRITE(unitno1,*)'Zone'
          WRITE(unitno2,*)'Zone'
          WRITE(unitno3,*)'Zone'

          xmin = xsmin!*radbdy(m)
          xmax = xsmax!*radbdy(m)

          dx_x = (xmax-xmin)/real(count_phi-one,prcn)

          PHILOOP: do iphi = 1, count_phi
             nx = xmin + dx_x*real((iphi-1),prcn)

             xl(1)=xc(m,1)+nx*radbdy(m) !global coordinate

             cphi=(nx)

             phiang = ACOS(cphi)
             rad_proj = sin(phiang)
             phiang = oneeighty*phiang/pi
             !WRITE(*,'(A,5(2x,g17.8))')'xcor= ', iphi, nx, phiang, rad_proj
             phi_array(iphi) = phiang
             THETALOOP: do itheta = 1, count_theta 
                theta = dtheta*real(itheta-one,prcn)
                thetaang = oneeighty*theta/pi
                !IF(iphi.eq.1) WRITE(*,*) 'theta =', thetaang
                ny = rad_proj*cos(theta)
                nz = rad_proj*sin(theta) 
                
                xl(2)=xc(m,2)+ny*radbdy(m) !global coordinate
                xl(3)=xc(m,3)+nz*radbdy(m) !global coordinate
                
                theta_array(itheta) = theta*oneeighty/pi
                norm_phi(iphi) = norm_phi(iphi) + 1
                norm_theta(itheta) = norm_theta(itheta) + 1
                
                !print*,iphi, itheta, norm_phi(iphi), norm_theta(itheta)
                DO n=1,ndim
                   
                   is(n)=INT(xl(n))
                   
                   ul(n)=zero

                   ppll(n)=zero

                ENDDO
                rad = radbdy(m)

                pl=zero
                isp=INT(xl(1)-0.5)
                xpb(1) = xl(1)-0.5
                xpb(2:3)=xl(2:3)
                do n = 1, ndim
                   if(xpb(n).lt.zero) then 
                      pcell(n) = int(xpb(n)-1)
                      !because of int of example -1.1 is -1, but we want -2. So,
                      ! the adjustment 
                   else 
                      pcell(n) = int(xpb(n))
                   end if
                end do


                call interpolate_pdata(pcell,xpb, ppll,pl,l)

                do n = 1, ndim
                   if(xl(n).lt.zero) then 
                      pcell(n) = int(xl(n)-1)
                   else 
                      pcell(n) = int(xl(n))
                   end if
                end do

                call interpolate_udata(pcell,xl,ib&
                     &,ie,jb,je,kb,ke,ul,nll,onll,dfll, 0,m, l, onew) 

                vort(:) = zero 
                do k = 1, onew 
                   do j = 1, onew
                      do i = 1, onew
                         ii = ib+i-1
                         jj = jb+j-1
                         kk = kb+k-1
                         if(ii.lt.1) ii = mxf+ii-1
                         if(ii.gt.mxf-1) ii = ii-(mxf-1)
                         if(jj.lt.1) jj = my+jj
                         if(jj.gt.my) jj = jj-my
                         if(kk.lt.1) kk = mz+kk
                         if(kk.gt.mz) kk = kk-mz 

                         DO n=1,ndim
                            vort(n)=vort(n)+ weightp(i,j,k)*omega(ii,jj,kk,n) 
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO

                normal(1) = nx
                normal(2) = ny
                normal(3) = nz

                !WRITE(*,'(6(2x,g17.8))') 'norm = ', normal(1:3), vort(1:3)
                df(L,1)=vort(3)*normaL(2)
                df(L,1)=df(L,1)-vort(2)*normal(3)

                df(L,2)=vort(1)*normal(3)
                df(L,2)=df(L,2)-vort(3)*normal(1)

                df(L,3)=vort(2)*normal(1)
                df(L,3)=df(L,3)-vort(1)*normal(2)

                df(L,1)=-df(L,1)
                df(L,2)=-df(L,2)
                df(L,3)=-df(L,3)

                !IF(THETAANG.GT.oneeighty) WRITE(*,*)'p1 = ', pl
                !---------------------------------------------------------------
                !     calculate the viscous and pressure components separately
                !---------------------------------------------------------------

                DO d=1,ndim,1
                   pres_bnd(d) = -pl
                   visc_bnd(d) = df(l,d)
                ENDDO

                visc_loc = array_dot_product(visc_bnd(1:3),&
                     & axial_direction(1:3))
                pres_loc = -pl*array_dot_product(normal(1:3),&
                     & axial_direction(1:3))

                visc_loc = vis*visc_loc*areabdy(m)/norm_factor2
                pres_loc = pres_loc*areabdy(m)/norm_factor2

                !WRITE(*,'(6(2x,g17.8))') 'norm = ', pres_loc, visc_loc
                !pres_loc = pres_loc/(0.5d0*meanslipmod*meanslipmod)
                !visc_loc = visc_loc*dchar/meanslipmod
                Fvis_phi(m,iphi) = Fvis_phi(m,iphi) + visc_loc

                Fvis_theta(m,itheta) = Fvis_theta(m,itheta) + visc_loc

                pres_phi(m,iphi) = pres_phi(m,iphi) + pres_loc
                pres_theta(m,itheta) = pres_theta(m,itheta) + pres_loc
                
                IF (thetaang.eq.zero) THEN 
                   WRITE(unitno,21) oneeighty-phiang,pres_loc,visc_loc

                ENDIF
                

                IF (phiang.eq.half*oneeighty) THEN 
                   WRITE(unitno1,21) thetaang,pres_loc,visc_loc
                ENDIF
             ENDDO THETALOOP
             

          end do PHILOOP
          !READ(*,*)
          DO j = 1, count_phi 
             Fvis_phi(m,j) = Fvis_phi(m,j)/real((norm_phi(j)),prcn)
             Pres_phi(m,j) = Pres_phi(m,j)/real((norm_phi(j)),prcn)

             total_Fvis(m,1) =  total_Fvis(m,1) + Fvis_phi(m,j)
             total_pres(m,1) =  total_Pres(m,1) + Pres_phi(m,j)

             WRITE(unitno2,21) oneeighty-phi_array(j), Fvis_phi(m,j),&
                  & pres_phi(m,j) 

          end DO
          total_Fvis(m,1) = total_Fvis(m,1)/count_phi
          total_Pres(m,1) = total_Pres(m,1)/count_phi

          DO j = 1, count_theta
             Fvis_theta(m,j) = Fvis_theta(m,j)/real((norm_theta(j)),prcn)
             Pres_theta(m,j) = Pres_theta(m,j)/real((norm_theta(j)),prcn)

             total_Fvis(m,2) = total_Fvis(m,2) + Fvis_theta(m,j)
             total_pres(m,2) = total_Pres(m,2) + Pres_theta(m,j)

             WRITE(unitno3,21) theta_array(j), Fvis_theta(m,j),&
                  & pres_theta(m,j)
          end DO

          total_Fvis(m,2) = total_Fvis(m,2)/count_theta
          total_Pres(m,2) = total_Pres(m,2)/count_theta

          WRITE(*,'(A40,/,2x, i4, 2(2x,g17.8))') 'AVG VISC DRAG  ALONG PHI AND THET&
               &A FOR M = ',M, total_Fvis(m,1), & 
               & total_Fvis(m,2)


          WRITE(*,'(A40,/,2x, i4, 2(2x,g17.8))') 'AVG PRES DRAG  ALONG PHI AND THET&
               &A FOR M = ',M, total_Pres(m,1), & 
               & total_Pres(m,2)

       end DO BODYLOOP

       avg_Fvis(1) =  SUM(total_Fvis(:,1))/real(nbody,prcn)
       avg_Fvis(2) =  SUM(total_Fvis(:,2))/real(nbody,prcn)
       avg_Pres(1) =  SUM(total_Pres(:,1))/real(nbody,prcn)
       avg_Pres(2) =  SUM(total_Pres(:,2))/real(nbody,prcn)

       do j = 1, count_phi 
          Fvis_phi_avg(j) = SUM(Fvis_phi(1:nbody, j))/real(nbody,prcn)
          Pres_phi_avg(j) = SUM(Pres_phi(1:nbody, j))/real(nbody,prcn)
       end do

       do j=1, count_phi
          conf1 = zero; conf2 = zero
          IF(nbody.ge.2) THEN  
             do m = 1, nbody
                conf1 = conf1 + (Fvis_phi(m,j) - Fvis_phi_avg(j))&
                     &**two
                conf2 = conf2 + (Pres_phi(m,j) - Pres_phi_avg(j))&
                     &**two
             end do
             conf1 = dsqrt(conf1/real(nbody-1,prcn))
             conf1=1.96/sqrt(dble(nbody))*conf1 
             conf2 = dsqrt(conf2/real(nbody-1,prcn))
             conf2=1.96/sqrt(dble(nbody))*conf1 
          end IF
          WRITE(unitno4,21) oneeighty-phi_array(j), Fvis_phi_avg(j),&
               & pres_phi_avg(j) , conf1, conf2 
       end do

       do j = 1, count_theta
          Fvis_theta_avg(j) = SUM(Fvis_theta(1:nbody, j))/real(nbody,prcn)
          Pres_theta_avg(j) = SUM(Pres_theta(1:nbody, j))/real(nbody,prcn)
       end do

       do j=1, count_theta
          conf1 = zero; conf2 = zero
          IF(nbody.ge.2) THEN 
             do m = 1, nbody
                conf1 = conf1 + (Fvis_theta(m,j) - Fvis_theta_avg(j))&
                     &**two
                conf2 = conf2 + (Pres_theta(m,j) - Pres_theta_avg(j))&
                     &**two

             end do
             conf1 = dsqrt(conf1/real(nbody-1,prcn))
             conf1=1.96/sqrt(dble(nbody))*conf1 
             conf2 = dsqrt(conf2/real(nbody-1,prcn))
             conf2=1.96/sqrt(dble(nbody))*conf1 
          end IF

          WRITE(unitno5,21) theta_array(j), Fvis_theta_avg(j),&
               & pres_theta_avg(j) , conf1, conf2 
       end do

       WRITE(*,'(A40,/, 4(2x,g17.8))') 'AVG VISC DRAG  ALONG PHI AND THET&
            &A  = ', avg_Fvis(1:2), SUM(Fvis_phi_avg(1:count_phi))&
            &/real(count_phi),  SUM(Fvis_theta_avg(1:count_theta))&
            &/real(count_theta)


       WRITE(*,'(A40,/, 4(2x,g17.8))') 'AVG PRES DRAG  ALONG PHI AND THET&
            &A', avg_Pres(1:2), SUM(Pres_phi_avg(1:count_phi))&
            &/real(count_phi),  SUM(Pres_theta_avg(1:count_theta))&
            &/real(count_theta) 
21     FORMAT(10(2xe17.4))

       CLOSE(unitno,status= 'keep')
       CLOSE(unitno1,status= 'keep')
       CLOSE(unitno2,status= 'keep')
       CLOSE(unitno3,status= 'keep')
       CLOSE(unitno4,status= 'keep')
       CLOSE(unitno5,status= 'keep')
     END SUBROUTINE calc_local_pres_visc_drag_plane


  SUBROUTINE write_complex_forcing
    IMPLICIT NONE
    INTEGER :: i, j, k

    if(I_AM_NODE_ZERO)then
       
       OPEN(1000,FILE=TRIM(RUN_NAME)//'_complex_force.dat',status='unknown')
       write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "FX" ', ' "FY" &
            &',' "FZ" '!,' "P" ' !, ' "nl2" ', ' "nl3" '
#if PARALLEL
       write(1000,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
#else
       write(1000,*)'ZONE F=POINT, I=', nx/2,  ', J=', my, ', K=', mz
#endif
       do k=1,mz
          do j=1,my2
#if PARALLEL
             do i=1,nx !mx
#else
             do  i = 1,nx/2
#endif
                write(1000,*)(GLOBAL_INDEX(i)),(j),(k),dreal(ff(i,j,k,1))&
                           &,dreal(ff(i,j,k,2)),dreal(ff(i,j,k,3))
!!$
!!$                      write(1000,*)(GLOBAL_INDEX(i)),(j),(k),fr(i,j,k,1)&
!!$                           &,fr(i,j,k,2),fr(i,j,k,3)!,pr(i,j,k)/(half*upi(1)*upi(1))!,divur(i,j,k)!,nlbcp(i,j,k,1),nlbcp(i,j,k,2),nlbcp(i,j,k,3)!!velr_mxf(i,j,k,3)
!!$                      !write(14,*)(i),(j),(k),divur(i,j,k)
             enddo
          enddo
       enddo
       close(1000,status='keep')
    end if

    end subroutine write_complex_forcing

  SUBROUTINE write_real_forcing
    USE nlmainarrays , ONLY : nlbc, ubc
    IMPLICIT NONE
    INTEGER :: i, j, k
        

#if PARALLEL       
    OPEN(1000,FILE='_real_force.dat',status='unknown')
#else
    OPEN(1000,FILE=TRIM(RUN_NAME)//'_real_force.dat',status='unknown')
#endif
    write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "FX" ', ' "FY" &
               &',' "FZ" '
#if PARALLEL
    write(1000,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
#else
    write(1000,*)'ZONE F=POINT, I=', nx,  ', J=', my, ', K=', mz
#endif
       do k=1,mz
          do j=1,my
#if PARALLEL
             do i=1,nx !mx
#else             
             do  i =  1, nx
#endif
!             write(1000,*)(GLOBAL_INDEX(i)),(j),(k),(fromleft(j,k,1))&
!                  &,(fromright(j,k,1)),(fromleft(j,k,1)+fromright(j,k,1)), nlbc(i,j,k,1)
                
!             write(1000,*)(GLOBAL_INDEX(i)),(j),(k),(fromleft(j,k,1))&
!                  &,(fromright(j,k,1)),(fromleft(j,k,1)+fromright(j,k,1)), nlbc(i,j,k,1)

             write(1000,*)(GLOBAL_INDEX(i)),(j),(k),(fr(i,j,k,1))&
                  &,(ppr(i,j,k,1)),(nlbc(i,j,k,1))


!!$
!!$                      write(1000,*)(GLOBAL_INDEX(i)),(j),(k),fr(i,j,k,1)&
!!$                           &,fr(i,j,k,2),fr(i,j,k,3)!,pr(i,j,k)/(half*upi(1)*upi(1))!,divur(i,j,k)!,nlbcp(i,j,k,1),nlbcp(i,j,k,2),nlbcp(i,j,k,3)!!velr_mxf(i,j,k,3)
!!$                      !write(14,*)(i),(j),(k),divur(i,j,k)
             enddo
          enddo
       enddo
       close(1000,status='keep')
!       end if
     end subroutine write_real_forcing

END MODULE boundary_condition
    
