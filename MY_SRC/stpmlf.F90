MODULE stpmlf
   !!======================================================================
   !!                       ***  MODULE stpMLF  ***
   !! Time-stepping   : manager of the ocean, tracer and ice time stepping
   !!                   using Modified Leap Frog for OCE
   !!======================================================================
   !! History :  OPA  !  1991-03  (G. Madec)  Original code
   !!             -   !  1991-11  (G. Madec)
   !!             -   !  1992-06  (M. Imbard)  add a first output record
   !!             -   !  1996-04  (G. Madec)  introduction of dynspg
   !!             -   !  1996-04  (M.A. Foujols)  introduction of passive tracer
   !!            8.0  !  1997-06  (G. Madec)  new architecture of call
   !!            8.2  !  1997-06  (G. Madec, M. Imbard, G. Roullet)  free surface
   !!             -   !  1999-02  (G. Madec, N. Grima)  hpg implicit
   !!             -   !  2000-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
   !!   NEMO     1.0  !  2002-06  (G. Madec)  free form, suppress macro-tasking
   !!             -   !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-01  (C. Ethe) Add the KPP closure scheme
   !!             -   !  2005-11  (G. Madec)  Reorganisation of tra and dyn calls
   !!             -   !  2006-01  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   !  2006-07  (S. Masson)  restart using iom
   !!            3.2  !  2009-02  (G. Madec, R. Benshila)  reintroduicing z*-coordinate
   !!             -   !  2009-06  (S. Masson, G. Madec)  TKE restart compatible with key_cpl
   !!            3.3  !  2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal
   !!            3.6  !  2012-07  (J. Simeon, G. Madec. C. Ethe)  Online coarsening of outputs
   !!            3.6  !  2014-04  (F. Roquet, G. Madec) New equations of state
   !!            3.6  !  2014-10  (E. Clementi, P. Oddo) Add Qiao vertical mixing in case of waves
   !!            3.7  !  2014-10  (G. Madec)  LDF simplication
   !!             -   !  2014-12  (G. Madec) remove KPP scheme
   !!             -   !  2015-11  (J. Chanut) free surface simplification (remove filtered free surface)
   !!            4.0  !  2017-05  (G. Madec)  introduction of the vertical physics manager (zdfphy)
   !!            4.1  !  2019-08  (A. Coward, D. Storkey) rewrite in preparation for new timestepping scheme
   !!            4.x  !  2020-08  (S. Techene, G. Madec)  quasi eulerian coordinate time stepping
   !!----------------------------------------------------------------------
#if defined key_qco   ||   defined key_linssh
   !!----------------------------------------------------------------------
   !!   'key_qco'                        Quasi-Eulerian vertical coordinate
   !!                          OR
   !!   'key_linssh                       Fixed in time vertical coordinate
   !!----------------------------------------------------------------------
   !!
   !!----------------------------------------------------------------------
   !!   stp_MLF       : NEMO modified Leap Frog time-stepping with qco or linssh
   !!----------------------------------------------------------------------
   USE step_oce       ! time stepping definition modules
   !
   USE domqco         ! quasi-eulerian coordinate
   USE traatf_qco     ! time filtering                 (tra_atf_qco routine)
   USE dynatf_qco     ! time filtering                 (dyn_atf_qco routine)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp_MLF   ! called by nemogcm.F90

   !                                          !**  time level indices  **!
   INTEGER, PUBLIC ::   Nbb, Nnn, Naa, Nrhs   !: used by nemo_init

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   RECURSIVE SUBROUTINE stp_MLF( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp_MLF( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp_MLF  ***
      !!
      !! ** Purpose : - Time stepping of OCE  (momentum and active tracer eqs.)
      !!              - Time stepping of SI3 (dynamic and thermodynamic eqs.)
      !!              - Time stepping of TRC  (passive tracer eqs.)
      !!
      !! ** Method  : -1- Update forcings and data
      !!              -2- Update ocean physics
      !!              -3- Compute the t and s trends
      !!              -4- Update t and s
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, hdiv,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk, jtile   ! dummy loop indice
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zgdept
      !! ---------------------------------------------------------------------
         WRITE(numout,*) 'The step MLF subroutine stp_MLF is called'

#if defined key_agrif
      IF( nstop > 0 ) RETURN   ! avoid to go further if an error was detected during previous time step (child grid)
      kstp = nit000 + Agrif_Nb_Step()
      Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs   ! agrif_oce module copies of time level indices
      IF( lk_agrif_debug ) THEN
         IF( Agrif_Root() .and. lwp)   WRITE(*,*) '---'
         IF(lwp)   WRITE(*,*) 'Grid Number', Agrif_Fixed(),' time step ', kstp, 'int tstep', Agrif_NbStepint()
      ENDIF
      IF( kstp == nit000 + 1 )   lk_agrif_fstep = .FALSE.
# if defined key_xios
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap( cxios_context )
# endif
#endif
      !
      IF( ln_timing )   CALL timing_start('stp_MLF')
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! model timestep
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF( l_1st_euler ) THEN     ! start or restart with Euler 1st time-step
         rDt   = rn_Dt   
         r1_Dt = 1._wp / rDt
      ENDIF
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! update I/O and calendar
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF( kstp == nit000 ) THEN                       ! initialize IOM context (must be done after nemo_init for AGRIF+XIOS+OASIS)
                             CALL iom_init( cxios_context, ld_closedef=.FALSE. )   ! for model grid (including possible AGRIF zoom)
         IF( lk_diamlr   )   CALL dia_mlr_iom_init    ! with additional setup for multiple-linear-regression analysis
                             CALL iom_init_closedef
         IF( ln_crs      )   CALL iom_init( TRIM(cxios_context)//"_crs" )  ! for coarse grid
      ENDIF
      IF( kstp == nitrst .AND. lwxios ) THEN
                             CALL iom_swap(                     cw_ocerst_cxt )
                             CALL iom_init_closedef(            cw_ocerst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_ocerst_cxt )
#if defined key_top
                             CALL iom_swap(                     cw_toprst_cxt )
                             CALL iom_init_closedef(            cw_toprst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_toprst_cxt )
#endif
      ENDIF
      IF( kstp + nn_fsbc - 1 == nitrst .AND. lwxios ) THEN
#if defined key_si3
                             CALL iom_swap(                     cw_icerst_cxt )
                             CALL iom_init_closedef(            cw_icerst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_icerst_cxt )
#endif
         IF( ln_abl      ) THEN
                             CALL iom_swap(                     cw_ablrst_cxt )
                             CALL iom_init_closedef(            cw_ablrst_cxt )
                             CALL iom_setkt( kstp - nit000 + 1, cw_ablrst_cxt )
         ENDIF
      ENDIF
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1,      cxios_context          )   ! tell IOM we are at time step kstp
      IF( ln_crs         )   CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" )   ! tell IOM we are at time step kstp

                         CALL sbc     ( kstp, Nbb, Nnn )              ! Sea Boundary Condition (including sea-ice)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_floats  )   CALL flo_stp   ( kstp, Nbb, Nnn )      ! drifting Floats
      IF( ln_diacfl  )   CALL dia_cfl   ( kstp,      Nnn )      ! Courant number diagnostics
                         CALL dia_hth   ( kstp,      Nnn )      ! Thermocline depth (20 degres isotherm depth)
      IF( ln_diadct  )   CALL dia_dct   ( kstp,      Nnn )      ! Transports 
                         CALL dia_ar5   ( kstp,      Nnn )      ! ar5 diag
                         CALL dia_ptr   ( kstp,      Nnn )      ! Poleward adv/ldf TRansports diagnostics
                         CALL dia_wri   ( kstp,      Nnn )      ! ocean model: outputs
      IF( ln_crs     )   CALL crs_fld   ( kstp,      Nnn )      ! ocean model: online field coarsening & output
      IF( lk_diadetide ) CALL dia_detide( kstp )                ! Weights computation for daily detiding of model diagnostics
      IF( lk_diamlr  )   CALL dia_mlr                           ! Update time used in multiple-linear-regression analysis

      IF( lrst_oce   )   CALL rst_write    ( kstp, Nbb, Nnn )   ! write output ocean restart file
      IF( ln_sto_eos )   CALL sto_rst_write( kstp )   ! write restart file for stochastic parameters

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! File manipulation at the end of the first time step
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( kstp == nit000 ) THEN                          ! 1st time step only
                                        CALL iom_close( numror )   ! close input  ocean restart file
         IF( lrxios )                   CALL iom_context_finalize( cr_ocerst_cxt )
         IF(lwm)                        CALL FLUSH    ( numond )   ! flush output namelist oce
         IF(lwm .AND. numoni /= -1 )    CALL FLUSH    ( numoni )   ! flush output namelist ice (if exist)
      ENDIF

      CALL dom_qco_r3c( ssh(:,:,Nnn), r3t(:,:,Nnn), r3u(:,:,Nnn), r3v(:,:,Nnn), r3f(:,:) )
      !  CALL dom_qco_r3c( ssh(:,:,Naa), r3t(:,:,Naa), r3u(:,:,Naa), r3v(:,:,Naa) )
         WRITE(numout,*) 'Sampled_r3t=', r3t(10,10,Nnn)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Coupled mode
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!      IF( lk_oasis .AND. nstop == 0 )   CALL sbc_cpl_snd( kstp, Nbb, Nnn )     ! coupled mode : field exchanges
      !
#if defined key_xios
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Finalize contextes if end of simulation or error detected
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( kstp == nitend .OR. nstop > 0 ) THEN
                      CALL iom_context_finalize(      cxios_context          ) ! needed for XIOS+AGRIF
         IF( ln_crs ) CALL iom_context_finalize( trim(cxios_context)//"_crs" ) !
      ENDIF
#endif
      !
      IF( l_1st_euler ) THEN         ! recover Leap-frog timestep
         rDt   = 2._wp * rn_Dt
         r1_Dt = 1._wp / rDt
         l_1st_euler = .FALSE.
      ENDIF
      !
!            ts  (:,:,:,:,Nnn) = ts (:,:,:,:,Nbb)       ! set present values to previous ones
!            uu    (:,:,:,Nnn) = uu   (:,:,:,Nbb)
!            vv    (:,:,:,Nnn) = vv   (:,:,:,Nbb)
!            ts  (:,:,:,:,Naa) = ts (:,:,:,:,Nbb)       ! set future values to previous ones
!            uu    (:,:,:,Naa) = uu   (:,:,:,Nbb)
!            vv    (:,:,:,Naa) = vv   (:,:,:,Nbb)
!                         uu(:,:,:,Nrhs) = 0._wp        ! set trends to zero
!                         vv(:,:,:,Nrhs) = 0._wp
!                         ts(:,:,:,:,Nrhs) = 0._wp

      IF( ln_timing )   CALL timing_stop('stp_MLF')
      !
   END SUBROUTINE stp_MLF


   SUBROUTINE mlf_baro_corr( Kmm, Kaa, puu, pvv )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mlf_baro_corr  ***
      !!
      !! ** Purpose :   Finalize after horizontal velocity.
      !!
      !! ** Method  : * Ensure after velocities transport matches time splitting
      !!             estimate (ln_dynspg_ts=T)
      !!
      !! ** Action :   puu(Kmm),pvv(Kmm)   updated now horizontal velocity (ln_bt_fw=F)
      !!               puu(Kaa),pvv(Kaa)   after horizontal velocity
      !!----------------------------------------------------------------------
      USE dynspg_ts, ONLY : un_adv, vn_adv   ! updated Kmm barotropic transport 
      !!
      INTEGER                             , INTENT(in   ) ::   Kmm, Kaa   ! before and after time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv   ! velocities
      !
      INTEGER  ::   ji,jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::   zue, zve
      !!----------------------------------------------------------------------

      ! Ensure below that barotropic velocities match time splitting estimate
      ! Compute actual transport and replace it with ts estimate at "after" time step
      DO_2D( 0, 0, 0, 0 )
         zue(ji,jj) = e3u(ji,jj,1,Kaa) * puu(ji,jj,1,Kaa) * umask(ji,jj,1)
         zve(ji,jj) = e3v(ji,jj,1,Kaa) * pvv(ji,jj,1,Kaa) * vmask(ji,jj,1)
      END_2D
      DO jk = 2, jpkm1
         DO_2D( 0, 0, 0, 0 )
            zue(ji,jj) = zue(ji,jj) + e3u(ji,jj,jk,Kaa) * puu(ji,jj,jk,Kaa) * umask(ji,jj,jk)
            zve(ji,jj) = zve(ji,jj) + e3v(ji,jj,jk,Kaa) * pvv(ji,jj,jk,Kaa) * vmask(ji,jj,jk)
         END_2D
      END DO
      DO jk = 1, jpkm1
         DO_2D( 0, 0, 0, 0 )
            puu(ji,jj,jk,Kaa) = ( puu(ji,jj,jk,Kaa) - zue(ji,jj) * r1_hu(ji,jj,Kaa) + uu_b(ji,jj,Kaa) ) * umask(ji,jj,jk)
            pvv(ji,jj,jk,Kaa) = ( pvv(ji,jj,jk,Kaa) - zve(ji,jj) * r1_hv(ji,jj,Kaa) + vv_b(ji,jj,Kaa) ) * vmask(ji,jj,jk)
         END_2D
      END DO
      !
      IF( .NOT.ln_bt_fw ) THEN
         ! Remove advective velocity from "now velocities"
         ! prior to asselin filtering
         ! In the forward case, this is done below after asselin filtering
         ! so that asselin contribution is removed at the same time
         DO jk = 1, jpkm1
            puu(:,:,jk,Kmm) = ( puu(:,:,jk,Kmm) - un_adv(:,:)*r1_hu(:,:,Kmm) + uu_b(:,:,Kmm) )*umask(:,:,jk)
            pvv(:,:,jk,Kmm) = ( pvv(:,:,jk,Kmm) - vn_adv(:,:)*r1_hv(:,:,Kmm) + vv_b(:,:,Kmm) )*vmask(:,:,jk)
         END DO
      ENDIF
      !
   END SUBROUTINE mlf_baro_corr


   SUBROUTINE finalize_lbc( kt, Kbb, Kaa, puu, pvv, pts )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE finalize_lbc  ***
      !!
      !! ** Purpose :   Apply the boundary condition on the after velocity
      !!
      !! ** Method  : * Apply lateral boundary conditions on after velocity
      !!             at the local domain boundaries through lbc_lnk call,
      !!             at the one-way open boundaries (ln_bdy=T),
      !!             at the AGRIF zoom   boundaries (lk_agrif=T)
      !!
      !! ** Action :   puu(Kaa),pvv(Kaa)   after horizontal velocity and tracers
      !!----------------------------------------------------------------------
#if defined key_agrif
      USE agrif_oce_interp
#endif
      USE bdydyn         ! ocean open boundary conditions (define bdy_dyn)
      !!
      INTEGER                                  , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                                  , INTENT(in   ) ::   Kbb, Kaa   ! before and after time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt)     , INTENT(inout) ::   puu, pvv   ! velocities to be time filtered
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts        ! active tracers
      !!----------------------------------------------------------------------
      !
      ! Update after tracer and velocity on domain lateral boundaries
      !
# if defined key_agrif
            CALL Agrif_tra                     !* AGRIF zoom boundaries
            CALL Agrif_dyn( kt )
# endif
      !                                        ! local domain boundaries  (T-point, unchanged sign)
      CALL lbc_lnk( 'finalize_lbc', puu(:,:,:,       Kaa), 'U', -1., pvv(:,:,:       ,Kaa), 'V', -1.   &
                       &          , pts(:,:,:,jp_tem,Kaa), 'T',  1., pts(:,:,:,jp_sal,Kaa), 'T',  1. )
      !
      ! lbc_lnk needed for zdf_sh2 when using nn_hls = 2, moved here to allow tiling in zdf_phy
      IF( nn_hls == 2 .AND. l_zdfsh2 ) CALL lbc_lnk( 'stp', avm_k, 'W', 1.0_wp )

      ! dom_qco_r3c defines over [nn_hls, nn_hls-1, nn_hls, nn_hls-1]
      IF( nn_hls == 2 .AND. .NOT. lk_linssh ) THEN
         CALL lbc_lnk( 'finalize_lbc', r3u(:,:,Kaa), 'U', 1._wp, r3v(:,:,Kaa), 'V', 1._wp, &
            &                          r3u_f(:,:),   'U', 1._wp, r3v_f(:,:),   'V', 1._wp )
      ENDIF
      !                                        !* BDY open boundaries
      IF( ln_bdy )   THEN
                               CALL bdy_tra( kt, Kbb, pts,      Kaa )
         IF( ln_dynspg_exp )   CALL bdy_dyn( kt, Kbb, puu, pvv, Kaa )
         IF( ln_dynspg_ts  )   CALL bdy_dyn( kt, Kbb, puu, pvv, Kaa, dyn3d_only=.true. )
      ENDIF
      !
   END SUBROUTINE finalize_lbc

#else
   !!----------------------------------------------------------------------
   !!   default option             EMPTY MODULE           qco not activated
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================
END MODULE stpmlf
