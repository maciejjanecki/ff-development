      module modal_aero_data

!--------------------------------------------------------------
! ... Basic aerosol mode parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod,  only: r8 => shr_kind_r8
      use constituents,  only: pcnst
      use radconstants,  only: nswbands, nlwbands

      implicit none
      save
#if (defined MODAL_AERO_7MODE)
      integer, parameter :: maxd_amode = 7,  &
#elif (defined MODAL_AERO_3MODE)
      integer, parameter :: maxd_amode = 3,  &
#endif
                            maxd_aspectype = 14

      integer                                               &   !
          ntot_amode,                                       &   !
          ntot_aspectype,                                   &   !
          msectional,                                       &   !
          nspec_amode( maxd_amode ),                        &   !
          lspectype_amode( maxd_aspectype, maxd_amode ),    &   !
          lmassptr_amode( maxd_aspectype, maxd_amode ),     &   !
          lmassptrcw_amode( maxd_aspectype, maxd_amode ),   &   !
!         lwaterptr_amode( maxd_amode ),                    &   !
          numptr_amode( maxd_amode ),                       &   !
          numptrcw_amode( maxd_amode ),                     &   !
          mcalcwater_amode( maxd_amode ),                   &   !
          mprognum_amode( maxd_amode ),                     &   !
          mdiagnum_amode( maxd_amode ),                     &   !
          mprogsfc_amode( maxd_amode )
! following not needed
!         nspec_amode_nontracer( maxd_amode ),              &   !
!         lkohlercptr_amode( maxd_amode ),                  &   !
!         lsfcptr_amode( maxd_amode ),                      &   !
!         lsfcptrcw_amode( maxd_amode ),                    &   !
!         lsigptr_amode( maxd_amode ),                      &   !
!         lsigptrcw_amode( maxd_amode ),                    &   !
!         lsigptrac_amode( maxd_amode ),                    &   !

      real(r8) ::                                 &   !
          dgnum_amode( maxd_amode ),              &   !
          dgnumlo_amode( maxd_amode ),            &   !
          dgnumhi_amode( maxd_amode ),            &   !
          sigmag_amode( maxd_amode ),             &   !
          alnsg_amode( maxd_amode ),              &   !
          voltonumb_amode( maxd_amode ),          &   !
          voltonumblo_amode( maxd_amode ),        &   !
          voltonumbhi_amode( maxd_amode ),        &   !
          alnv2n_amode( maxd_amode ),             &   !
          alnv2nlo_amode( maxd_amode ),           &   !
          alnv2nhi_amode( maxd_amode ),           &   !
          specdens_amode( maxd_aspectype ),       &   !
          specmw_amode( maxd_aspectype ),         &   !
          spechygro( maxd_aspectype ),            &   !
          rhcrystal_amode( maxd_amode ),          &   !
          rhdeliques_amode( maxd_amode )
! following not needed
!         sigmaglo_amode( maxd_amode ),           &   !
!         sigmaghi_amode( maxd_amode ),           &   !
!         alnsglo_amode( maxd_amode ),            &   !
!         alnsghi_amode( maxd_amode ),            &   !
!         voltosfc_amode( maxd_amode ),           &   !
!         voltosfclo_amode( maxd_amode ),         &   !
!         voltosfchi_amode( maxd_amode ),         &   !
!         alnv2s_amode( maxd_amode ),             &   !
!         alnv2slo_amode( maxd_amode ),           &   !
!         alnv2shi_amode( maxd_amode ),           &   !

      complex                                     &   !
          specrefndxsw( nswbands, maxd_aspectype ),           &   !
          specrefndxlw( nlwbands, maxd_aspectype )

      character(len=10) :: specname_amode( maxd_aspectype )
      character(len=16) :: modename_amode( maxd_amode )
      character(len=16) :: cnst_name_cw( pcnst )

      character(len=8) :: aodvisname(maxd_amode ),       &
                          ssavisname(maxd_amode )
      character(len=48) :: aodvislongname(maxd_amode ),  &
                           ssavislongname(maxd_amode )

      character(len=8) :: fnactname(maxd_amode ),   &
                          fmactname(maxd_amode ),   &
                          nactname(maxd_amode )
      character(len=48) :: fnactlongname(maxd_amode ),   &
                           fmactlongname(maxd_amode ),   &
                           nactlongname(maxd_amode )

      integer                                       &   !
          lptr_so4_a_amode(maxd_amode),  lptr_so4_cw_amode(maxd_amode), &   !
          lptr_msa_a_amode(maxd_amode),  lptr_msa_cw_amode(maxd_amode), &   !
          lptr_nh4_a_amode(maxd_amode),  lptr_nh4_cw_amode(maxd_amode), &   !
          lptr_no3_a_amode(maxd_amode),  lptr_no3_cw_amode(maxd_amode), &   !
          lptr_pom_a_amode(maxd_amode),  lptr_pom_cw_amode(maxd_amode), &   !
          lptr_soa_a_amode(maxd_amode),  lptr_soa_cw_amode(maxd_amode), &   !
          lptr_bc_a_amode(maxd_amode),   lptr_bc_cw_amode(maxd_amode),  &   !
          lptr_nacl_a_amode(maxd_amode), lptr_nacl_cw_amode(maxd_amode),&   !
          lptr_dust_a_amode(maxd_amode), lptr_dust_cw_amode(maxd_amode),&   !
          modeptr_accum,  modeptr_aitken,                               &   !
          modeptr_ufine,  modeptr_coarse,                               &   !
          modeptr_pcarbon,                                              &   !
          modeptr_finedust,  modeptr_fineseas,                          &   !
          modeptr_coardust,  modeptr_coarseas

      real(r8) ::             &
          specmw_so4_amode,     specdens_so4_amode,       &
          specmw_nh4_amode,     specdens_nh4_amode,       &
          specmw_no3_amode,     specdens_no3_amode,       &
          specmw_pom_amode,     specdens_pom_amode,       &
          specmw_soa_amode,     specdens_soa_amode,       &
          specmw_bc_amode,      specdens_bc_amode,        &
          specmw_dust_amode,    specdens_dust_amode,      &
          specmw_seasalt_amode, specdens_seasalt_amode

	integer species_class(pcnst)	! indicates species class (
				!     cldphysics, aerosol, gas )

	integer     spec_class_undefined
	parameter ( spec_class_undefined = 0 )
	integer     spec_class_cldphysics
	parameter ( spec_class_cldphysics = 1 )
	integer     spec_class_aerosol
	parameter ( spec_class_aerosol = 2 )
	integer     spec_class_gas
	parameter ( spec_class_gas = 3 )
	integer     spec_class_other
	parameter ( spec_class_other = 4 )


!   threshold for reporting negatives from subr qneg3
      real(r8) :: qneg3_worst_thresh_amode(pcnst)

      end module modal_aero_data

!----------------------------------------------------------------
!
!   maxd_amode = maximum allowable number of aerosol modes
!   maxd_aspectype = maximum allowable number of chemical species
!       in each aerosol mode
!
!   ntot_amode = number of aerosol modes
!   ( ntot_amode_gchm = number of aerosol modes in gchm
!     ntot_amode_ccm2 = number of aerosol modes to be made known to ccm2
!       These are temporary until multi-mode activation scavenging is going.
!       Until then, ntot_amode is set to either ntot_amode_gchm or
!       ntot_amode_ccm2 depending on which code is active )
!
!   msectional - if positive, moving-center sectional code is utilized,
!       and each mode is actually a section.
!   msectional_concinit - if positive, special code is used to initialize
!       the mixing ratios of all the sections.
!
!   nspec_amode(m) = number of chemical species in aerosol mode m
!   nspec_amode_ccm2(m) = . . .  while in ccm2 code
!   nspec_amode_gchm(m) = . . .  while in gchm code
!   nspec_amode_nontracer(m) = number of "non-tracer" chemical
!       species while in gchm code
!   lspectype_amode(l,m) = species type/i.d. for chemical species l
!       in aerosol mode m.  (1=sulfate, others to be defined)
!   lmassptr_amode(l,m) = gchm r-array index for the mixing ratio
!       (moles-x/mole-air) for chemical species l in aerosol mode m
!       that is in clear air or interstitial air (but not in cloud water)
!   lmassptrcw_amode(l,m) = gchm r-array index for the mixing ratio
!       (moles-x/mole-air) for chemical species l in aerosol mode m
!       that is currently bound/dissolved in cloud water
!   lwaterptr_amode(m) = gchm r-array index for the mixing ratio
!       (moles-water/mole-air) for water associated with aerosol mode m
!       that is in clear air or interstitial air
!   lkohlercptr_amode(m) = gchm r-array index for the kohler "c" parameter
!       for aerosol mode m.  This is defined on a per-dry-particle-mass basis:
!           c = r(i,j,k,lkohlercptr_amode) * [rhodry * (4*pi/3) * rdry^3]
!   numptr_amode(m) = gchm r-array index for the number mixing ratio
!       (particles/mole-air) for aerosol mode m that is in clear air or
!       interstitial are (but not in cloud water).  If zero or negative,
!       then number is not being simulated.
!   ( numptr_amode_gchm(m) = same thing but for within gchm
!     numptr_amode_ccm2(m) = same thing but for within ccm2
!       These are temporary, to allow testing number in gchm before ccm2 )
!   numptrcw_amode(m) = gchm r-array index for the number mixing ratio
!       (particles/mole-air) for aerosol mode m
!       that is currently bound/dissolved in cloud water
!   lsfcptr_amode(m) = gchm r-array index for the surface area mixing ratio
!       (cm^2/mole-air) for aerosol mode m that is in clear air or
!       interstitial are (but not in cloud water).  If zero or negative,
!       then surface area is not being simulated.
!   lsfcptrcw_amode(m) = gchm r-array index for the surface area mixing ratio
!       (cm^2/mole-air) for aerosol mode m that is currently
!       bound/dissolved in cloud water.
!   lsigptr_amode(m) = gchm r-array index for sigmag for aerosol mode m
!       that is in clear air or interstitial are (but not in cloud water).
!       If zero or negative, then the constant sigmag_amode(m) is used.
!   lsigptrcw_amode(m) = gchm r-array index for sigmag for aerosol mode m
!       that is currently bound/dissolved in cloud water.
!       If zero or negative, then the constant sigmag_amode(m) is used.
!   lsigptrac_amode(m) = gchm r-array index for sigmag for aerosol mode m
!       for combined clear-air/interstial plus bound/dissolved in cloud water.
!       If zero or negative, then the constant sigmag_amode(m) is used.
!
!   dgnum_amode(m) = geometric dry mean diameter (m) of the number
!       distribution for aerosol mode m.
!       (Only used when numptr_amode(m) is zero or negative.)
!   dgnumlo_amode(m), dgnumhi_amode(m) = lower and upper limits on the
!       geometric dry mean diameter (m) of the number distribution
!       (Used when mprognum_amode>0, to limit dgnum to reasonable values)
!   sigmag_amode(m) = geometric standard deviation for aerosol mode m
!   sigmaglo_amode(m), sigmaghi_amode(m) = lower and upper limits on the
!       geometric standard deviation of the number distribution
!       (Used when mprogsfc_amode>0, to limit sigmag to reasonable values)
!   alnsg_amode(m) = alog( sigmag_amode(m) )
!   alnsglo_amode(m), alnsghi_amode(m) = alog( sigmaglo/hi_amode(m) )
!   voltonumb_amode(m) = ratio of number to volume for mode m
!   voltonumblo_amode(m), voltonumbhi_amode(m) = ratio of number to volume
!       when dgnum = dgnumlo_amode or dgnumhi_amode, respectively
!   voltosfc_amode(m), voltosfclo_amode(m), voltosfchi_amode(m) - ratio of
!       surface to volume for mode m (like the voltonumb_amode's)
!   alnv2n_amode(m), alnv2nlo_amode(m), alnv2nhi_amode(m) -
!       alnv2n_amode(m) = alog( voltonumblo_amode(m) ), ...
!   alnv2s_amode(m), alnv2slo_amode(m), alnv2shi_amode(m) -
!       alnv2s_amode(m) = alog( voltosfclo_amode(m) ), ...
!   rhcrystal_amode(m) = crystalization r.h. for mode m
!   rhdeliques_amode(m) = deliquescence r.h. for mode m
!   (*** these r.h. values are 0-1 fractions, not 0-100 percentages)
!
!   mcalcwater_amode(m) - if positive, water content for mode m will be
!       calculated and stored in rclm(k,lwaterptr_amode(m)).  Otherwise, no.
!   mprognum_amode(m) - if positive, number mixing-ratio for mode m will
!       be prognosed.  Otherwise, no.
!   mdiagnum_amode(m) - if positive, number mixing-ratio for mode m will
!       be diagnosed and put into rclm(k,numptr_amode(m)).  Otherwise, no.
!   mprogsfc_amode(m) - if positive, surface area mixing-ratio for mode m will
!       be prognosed, and sigmag will vary temporally and spatially.
!       Otherwise, sigmag is constant.
!       *** currently surface area is not prognosed when msectional>0 ***
!
!   ntot_aspectype = overall number of aerosol chemical species defined (over all modes)
!   specdens_amode(l) = dry density (kg/m^3) of aerosol chemical species type l
!   specmw_amode(l) = molecular weight (kg/kmol) of aerosol chemical species type l
!   specname_amode(l) = name of aerosol chemical species type l
!   specrefndxsw(l) = complex refractive index (visible wavelengths)
!                   of aerosol chemical species type l
!   specrefndxlw(l) = complex refractive index (infrared wavelengths)
!                   of aerosol chemical species type l
!   spechygro(l) = hygroscopicity of aerosol chemical species type l
!
!   lptr_so4_a_amode(m), lptr_so4_cw_amode(m) = gchm r-array index for the
!       mixing ratio for sulfate associated with aerosol mode m
!       ("a" and "cw" phases)
!   (similar for msa, oc, bc, nacl, dust)
!
!   modename_amode(m) = character-variable name for mode m,
!       read from mirage2.inp
!   modeptr_accum - mode index for the main accumulation mode
!       if modeptr_accum = 1, then mode 1 is the main accumulation mode,
!       and modename_amode(1) = "accum"
!   modeptr_aitken - mode index for the main aitken mode
!       if modeptr_aitken = 2, then mode 2 is the main aitken mode,
!       and modename_amode(2) = "aitken"
!   modeptr_ufine - mode index for the ultrafine mode
!       if modeptr_ufine = 3, then mode 3 is the ultrafine mode,
!       and modename_amode(3) = "ufine"
!   modeptr_coarseas - mode index for the coarse sea-salt mode
!       if modeptr_coarseas = 4, then mode 4 is the coarse sea-salt mode,
!       and modename_amode(4) = "coarse seasalt"
!   modeptr_coardust - mode index for the coarse dust mode
!       if modeptr_coardust = 5, then mode 5 is the coarse dust mode,
!       and modename_amode(5) = "coarse dust"
!
!   specdens_XX_amode = dry density (kg/m^3) of aerosol chemical species type XX
!       where XX is so4, om, bc, dust, seasalt
!       contains same values as the specdens_amode array
!       allows values to be referenced differently
!   specmw_XX_amode = molecular weight (kg/kmol) of aerosol chemical species type XX
!       contains same values as the specmw_amode array
!
!-----------------------------------------------------------------------


!--------------------------------------------------------------
!
! ... aerosol size information for the current chunk
!
!--------------------------------------------------------------
!
!  dgncur = current geometric mean diameters (cm) for number distributions
!  dgncur_a - for unactivated particles, dry
!             (in physics buffer as DGNUM)
!  dgncur_awet - for unactivated particles, wet at grid-cell ambient RH
!             (in physics buffer as DGNUMWET)
!
!  the dgncur are computed from current mass and number
!  mixing ratios in the grid cell, BUT are then adjusted to be within
!  the bounds defined by dgnumlo/hi_amode
!
!  v2ncur = current (number/volume) ratio based on dgncur and sgcur
!              (volume in cm^3/whatever, number in particles/whatever)
!         == 1.0 / ( pi/6 * dgncur**3 * exp(4.5*((log(sgcur))**2)) )
!  v2ncur_a - for unactivated particles
!             (currently just defined locally)
!

