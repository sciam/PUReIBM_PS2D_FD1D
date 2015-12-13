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

MODULE machine
  USE global_data 

  !              record length used in open statement for unformatted,
  !              direct access file, with 512 bytes per record
  INTEGER  OPEN_N1
  !
  !              number of DOUBLE PRECISION words in 512 bytes
  INTEGER  NWORDS_DP
  !
  !              number of REAL words in 512 bytes
  INTEGER  NWORDS_R
  !
  !              number of INTEGER words in 512 bytes
  INTEGER  NWORDS_I
  !
  LOGICAL :: JUST_FLUSH = .TRUE.

CONTAINS

  !
  SUBROUTINE MACHINE_CONS
    !
    !
    IMPLICIT NONE
    !
    OPEN_N1   = 512
    NWORDS_DP =  64
    NWORDS_R  = 128
    NWORDS_I  = 128
    JUST_FLUSH = .TRUE.
    !
    RETURN
  END SUBROUTINE MACHINE_CONS
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !                                                                      C
  !  Module name: GET_RUN_ID                                             C
  !  Purpose: get the run id for this run                                C
  !                                                                      C
  !  Author: P. Nicoletti                               Date: 16-DEC-91  C
  !  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
  !                                                                      C
  !  Revision Number: 1                                                  C
  !  Purpose: add ndoe name                                              C
  !  Author: P.Nicoletti                                Date: 07-FEB-92  C
  !  Reviewer:                                          Date: dd-mmm-yy  C
  !                                                                      C
  !  Literature/Document References:                                     C
  !                                                                      C
  !  Variables referenced: None                                          C
  !  Variables modified: ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, ID_MINUTE   C
  !                      ID_SECOND, ID_NODE                              C
  !                                                                      C
  !  Local variables: TIME_ARRAY                                         C
  !                                                                      C
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  !
  SUBROUTINE GET_RUN_ID
    !
    IMPLICIT NONE
    !
    !             temporary array to hold time data
    INTEGER DAT(8)
    CHARACTER*10 DATE, TIM, ZONE


    ! Intel Linux compiler supports this function thru it's portability library
    CALL DATE_AND_TIME(DATE, TIM, ZONE, DAT)
    ID_YEAR   = DAT(1)
    ID_MONTH  = DAT(2)
    ID_DAY    = DAT(3)
    ID_HOUR   = DAT(5)
    ID_MINUTE = DAT(6)
    ID_SECOND = DAT(7)



    ! Intel Linux compiler supports this function thru it's portability library
    call hostnm(ID_NODE)      
    !
    RETURN
  END SUBROUTINE GET_RUN_ID
  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !                                                                      C
  !  Module name: CPU_TIME (CPU)                                         C
  !  Purpose: get the CPU time for the run                               C
  !                                                                      C
  !  Variables referenced: None                                          C
  !  Variables modified: None                                            C
  !                                                                      C
  !  Local variables: TA, XT                                             C
  !                                                                      C
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  !
  SUBROUTINE CPU_TIME(CPU)
    !
    IMPLICIT NONE
    !
    ! passed arguments
    !
    !                      cpu time since start of run
    DOUBLE PRECISION CPU

    INTEGER, SAVE :: COUNT_OLD=0, WRAP=0
    !
    ! local variables
    !

    !                       clock cycle
    INTEGER           COUNT

    !                       number of cycles per second
    INTEGER           COUNT_RATE

    !                       max number of cycles, after which count is reset to 0
    INTEGER           COUNT_MAX

    !
    ! Intel Linux compiler supports this function thru it's portability library
    CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
    IF(COUNT_OLD .GT. COUNT) THEN
       WRAP = WRAP + 1
    ENDIF
    COUNT_OLD = COUNT

    CPU           = DBLE(COUNT)/DBLE(COUNT_RATE) &
         + DBLE(WRAP) * DBLE(COUNT_MAX)/DBLE(COUNT_RATE)


    RETURN
  END SUBROUTINE CPU_TIME
end MODULE machine
