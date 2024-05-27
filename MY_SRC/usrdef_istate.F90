MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                  ===  ISOMIP configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2016-11 (S. Flavoni)             Original code
   !!                 ! 2017-02 (P. Mathiot, S. Flavoni) Adapt code to ISOMIP case
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce
   USE phycst         ! physical constants
   USE oce,    ONLY: rhd
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called by istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by domqco.F90
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10074 2018-08-28 16:15:49Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
!   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here ISOMIP configuration 
      !!
      !! ** Method  : - set temperature field
      !!              - set salinity    field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
!      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  ::   ji, jj, jk     ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : ISOMIP configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with a constant salinity and temperature.      '
      pu  (:,:,:) = 0._wp * umask(:,:,:)       ! ocean at rest
      pv  (:,:,:) = 0._wp * vmask(:,:,:) 
      pu  (:,:,3:jpkm1) = 0._wp * umask(:,:,3:jpkm1)        ! ocean at rest
      pv  (:,:,3:jpkm1) = 0._wp * vmask(:,:,3:jpkm1)
!      pssh(:,:)   = 0._wp
      !
      !                          ! T & S profiles
!      pts(:,:,:,jp_tem) = 2  * ptmask(:,:,:)          ! ISOMIP configuration : start from constant T+S fields
      pts(:,:,:,jp_tem) = 4._wp * ptmask(:,:,:)          ! ISOMIP configuration : start from constant T+S fields
      pts(:,:,:,jp_sal) = 34._wp * ptmask(:,:,:)


      DO ji=1,jpi
      DO jj=1,jpj
      DO jk=1,jpk
!         pts(ji,jj,jk,jp_tem) = ABS(glamt(ji,jj))/10000._wp
!         pts(ji,jj,jk,jp_sal) = ABS(gphit(ji,jj))/10000._wp
         pts(ji,jj,jk,jp_sal) = (34.25_wp + rhd(ji,jj,jk)/0.79_wp)*ptmask(ji,jj,jk)
         pv  (ji,jj,jk) = -((0.1120_wp)* vmask(ji,jj,jk) &
      &                 - 3.601_wp*((jk-1._wp)/(jpk-2._wp))* &
      &                   ((- grav*(+0.04_wp/2.E4_wp))/(-14E-5_wp))* vmask(ji,jj,1)* vmask(ji,jj,jk))
      ENDDO
      ENDDO
      ENDDO


!      pts(:,:,3:jpkm1,jp_tem) = 1  * ptmask(:,:,3:jpkm1)          ! ISOMIP configuration : start from constant T+S fields
      !   
   END SUBROUTINE usr_def_istate


   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !!
      !! ** Purpose :   Initialization of ssh
      !!
      !! ** Method  :   Set ssh as null, ptmask is required for test cases
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask   [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    !sea-surface height   [m]
      INTEGER  ::   ji, jj, jk     ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : GYRE configuration,analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~   Ocean at rest, ssh is zero'
      !
      ! Sea level:
      pssh(:,:) = 0._wp
      !

      DO ji=1,jpi
      DO jj=1,jpj
      DO jk=1,jpk
!         pts(ji,jj,jk,jp_tem) = ABS(glamt(ji,jj))/10000._wp
!         pts(ji,jj,jk,jp_sal) = ABS(gphit(ji,jj))/10000._wp
        pssh(ji,jj) = (-0.02_wp+0.04_wp*ABS(glamt(ji,jj))/2.E4_wp)*ptmask(ji,jj,1)
      ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE usr_def_istate_ssh

   !!======================================================================
END MODULE usrdef_istate
