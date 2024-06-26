MODULE usrdef_zgr
   !!======================================================================
   !!                   ***  MODULE  usrdef_zgr  ***
   !!
   !!                   ===      ICE_AGRIF case     ===
   !!
   !! Ocean domain : user defined vertical coordinate system 
   !!======================================================================
   !! History :  4.0  ! 2016-08  (G. Madec, S. Flavoni)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system (required)
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain
   USE usrdef_nam     ! User defined : namelist variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr   ! called by domzgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 12597 2020-03-25 08:57:21Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v , pe3f ,               &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw,                      &   !     -      -      -
      &                    k_top  , k_bot    )                             ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors 
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   jk, k_dz  ! dummy indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : ICE_AGRIF configuration '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (slab ocean - advection of an ice patch in a biperiodic square box domain)'
      !
      !
      ! type of vertical coordinate  ==>>>   here ICE_AGRIF : slab ocean always
      ! ---------------------------
      ld_zco    = .TRUE.       ! z-full-step coordinate
      ld_zps    = .FALSE.      ! z-partial-step coordinate
      ld_sco    = .FALSE.      ! s-coordinate
      ld_isfcav = .FALSE.      ! ISF Ice Shelves Flag
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !                       !==  UNmasked meter bathymetry  ==!
      !
      ! 
      k_dz = 10
      pdepw_1d(1) = 0.0
      pdept_1d(1) = k_dz / 2.0
      pe3w_1d = k_dz
      pe3t_1d = k_dz
      !
      DO jk = 2, jpk
         pdepw_1d(jk) = pdepw_1d(jk-1) +   k_dz
         pdept_1d(jk) = pdept_1d(jk-1) +   k_dz
      END DO
      !                       !==  top masked level bathymetry  ==!  (all coordinates)
      !
      ! no ocean cavities : top ocean level is ONE, except over land
      k_top(:,:) = 1
      !
      !                       !==  z-coordinate  ==!   (step-like topography)
      !                                !* bottom ocean compute from the depth of grid-points
      jpkm1 = jpk-1
      k_bot(:,:) = 10         ! here use k_top as a land mask
!      WHERE (( glamt < 10000 .OR. glamt > 30000 ))
!         k_bot(:,:) = 4
!      END WHERE
      !                                !* horizontally uniform coordinate (reference z-co everywhere)

! Obstacle 2
!      WHERE (( gphit < 5000 .AND. glamt > 20000))
!         k_bot(:,:) = 9
!      END WHERE

! Obstacle 3
!      WHERE (( glamt < 13000)) !.AND. glamt > 10000))
!         k_bot(:,:) = 9
!      END WHERE
!      WHERE (( glamt < 12000)) ! .AND. glamt > 11000))
!         k_bot(:,:) = 8
!      END WHERE
!      WHERE (( glamt < 11000)) ! .AND. glamt > 12000))
!         k_bot(:,:) = 7
!      END WHERE

      DO jk = 1, jpk
         pdept(:,:,jk) = pdept_1d(jk)
         pdepw(:,:,jk) = pdepw_1d(jk)
         pe3t (:,:,jk) = pe3t_1d (jk)
         pe3u (:,:,jk) = pe3t_1d (jk)
         pe3v (:,:,jk) = pe3t_1d (jk)
         pe3f (:,:,jk) = pe3t_1d (jk)
         pe3w (:,:,jk) = pe3w_1d (jk)
         pe3uw(:,:,jk) = pe3w_1d (jk)
         pe3vw(:,:,jk) = pe3w_1d (jk)
      END DO
      !
   END SUBROUTINE usr_def_zgr

   !!======================================================================
END MODULE usrdef_zgr
