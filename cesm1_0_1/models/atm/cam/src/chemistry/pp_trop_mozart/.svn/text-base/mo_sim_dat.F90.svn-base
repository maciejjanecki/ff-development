
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,   only  : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass
      use chem_mods,   only : diag_map
      use chem_mods,   only  : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,   only  : pht_alias_lst, pht_alias_mult
      use chem_mods,   only  : het_lst, extfrc_lst, inv_lst, slvd_lst
      use abortutils,  only  : endrun
      use mo_tracname, only  : solsym
      use chem_mods,   only : frc_from_dataset
      use shr_kind_mod,only : r8 => shr_kind_r8
      use cam_logfile, only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      clscnt(:) = (/  6 , 0 , 0 , 92 , 0 /)

      cls_rxt_cnt(:,1) = (/  29 , 9 , 0 , 1  /)
      cls_rxt_cnt(:,4) = (/  2 , 61 , 139 , 32  /)

      solsym(: 98) = (/ 'O3      ','O       ','O1D     ','N2O     ','NO      ', &
                        'NO2     ','NO3     ','HNO3    ','HO2NO2  ','N2O5    ', &
                        'H2      ','OH      ','HO2     ','H2O2    ','CH4     ', &
                        'CO      ','CH3O2   ','CH3OOH  ','CH2O    ','CH3OH   ', &
                        'C2H5OH  ','C2H4    ','EO      ','EO2     ','CH3COOH ', &
                        'GLYALD  ','C2H6    ','C2H5O2  ','C2H5OOH ','CH3CHO  ', &
                        'CH3CO3  ','CH3COOOH','C3H6    ','C3H8    ','C3H7O2  ', &
                        'C3H7OOH ','PO2     ','POOH    ','CH3COCH3','RO2     ', &
                        'ROOH    ','BIGENE  ','ENEO2   ','MEK     ','MEKO2   ', &
                        'MEKOOH  ','BIGALK  ','ALKO2   ','ALKOOH  ','ISOP    ', &
                        'ISOPO2  ','ISOPOOH ','MVK     ','MACR    ','MACRO2  ', &
                        'MACROOH ','MCO3    ','HYDRALD ','HYAC    ','CH3COCHO', &
                        'XO2     ','XOOH    ','C10H16  ','TERPO2  ','TERPOOH ', &
                        'TOLUENE ','CRESOL  ','TOLO2   ','TOLOOH  ','XOH     ', &
                        'BIGALD  ','GLYOXAL ','PAN     ','ONIT    ','MPAN    ', &
                        'ISOPNO3 ','ONITR   ','CB1     ','CB2     ','OC1     ', &
                        'OC2     ','SOA     ','SO2     ','SO4     ','DMS     ', &
                        'NH3     ','NH4     ','NH4NO3  ','SSLT01  ','SSLT02  ', &
                        'SSLT03  ','SSLT04  ','DST01   ','DST02   ','DST03   ', &
                        'DST04   ','Rn      ','Pb      ' /)

      adv_mass(: 98) = (/ 47.99819946_r8, 15.99940014_r8, 15.99940014_r8, 44.01287842_r8, 30.00613976_r8, &
                          46.00553894_r8, 62.00494003_r8, 63.01234055_r8, 79.01174164_r8, 108.0104828_r8, &
                          2.014800072_r8, 17.00679970_r8, 33.00619888_r8, 34.01359940_r8, 16.04059982_r8, &
                          28.01040077_r8, 47.03200150_r8, 48.03939819_r8, 30.02519989_r8, 32.04000092_r8, &
                          46.06579971_r8, 28.05159950_r8, 61.05780029_r8, 77.05719757_r8, 60.05039978_r8, &
                          60.05039978_r8, 30.06640053_r8, 61.05780029_r8, 62.06520081_r8, 44.05099869_r8, &
                          75.04239655_r8, 76.04979706_r8, 42.07740021_r8, 44.09220123_r8, 75.08360291_r8, &
                          76.09100342_r8, 91.08300018_r8, 92.09040070_r8, 58.07680130_r8, 89.06819916_r8, &
                          90.07559967_r8, 56.10319901_r8, 105.1088028_r8, 72.10260010_r8, 103.0940018_r8, &
                          104.1014023_r8, 72.14379883_r8, 103.1352005_r8, 104.1426010_r8, 68.11419678_r8, &
                          117.1197968_r8, 118.1271973_r8, 70.08779907_r8, 70.08779907_r8, 119.0933990_r8, &
                          120.1007996_r8, 101.0792007_r8, 100.1129990_r8, 74.07620239_r8, 72.06140137_r8, &
                          133.1192017_r8, 134.1266022_r8, 136.2283936_r8, 185.2339935_r8, 186.2413940_r8, &
                          92.13619995_r8, 108.1355972_r8, 141.1417999_r8, 142.1492004_r8, 158.1486053_r8, &
                          98.09819794_r8, 58.03559875_r8, 121.0479431_r8, 119.0743408_r8, 147.0847473_r8, &
                          162.1179352_r8, 147.1259460_r8, 12.01099968_r8, 12.01099968_r8, 12.01099968_r8, &
                          12.01099968_r8, 144.1320038_r8, 64.06479645_r8, 96.06359863_r8, 62.13240051_r8, &
                          17.02894020_r8, 18.03634071_r8, 80.04128265_r8, 58.44246674_r8, 58.44246674_r8, &
                          58.44246674_r8, 58.44246674_r8, 135.0640411_r8, 135.0640411_r8, 135.0640411_r8, &
                          135.0640411_r8, 222.0000000_r8, 207.1999969_r8 /)

      fix_mass(:  4) = (/ 0.00000000_r8, 28.0134792_r8, 31.9988003_r8, 18.0142002_r8 /)

      clsmap(:  6,1) = (/   15,   4,  16,  97,  98,  11 /)
      clsmap(: 92,4) = (/    1,   3,   2,   5,   6,   7,   8,   9,  10,  12, &
                            13,  14,  17,  18,  19,  20,  21,  22,  23,  24, &
                            25,  26,  27,  28,  29,  30,  31,  32,  33,  34, &
                            35,  36,  37,  38,  39,  40,  41,  42,  43,  47, &
                            48,  49,  44,  45,  46,  50,  51,  52,  53,  54, &
                            55,  56,  57,  58,  59,  60,  61,  62,  63,  64, &
                            65,  66,  67,  68,  69,  70,  71,  72,  73,  74, &
                            75,  76,  77,  78,  79,  83,  84,  85,  86,  87, &
                            88,  80,  81,  89,  90,  91,  92,  82,  93,  94, &
                            95,  96 /)

      permute(: 92,4) = (/   90,  29,  56,  88,  92,  84,  55,  40,  30,  86, &
                             87,  25,  89,  42,  77,  59,  32,  37,  26,  45, &
                             49,  60,  18,  68,  41,  73,  85,  48,  72,  19, &
                             69,  39,  65,  54,  64,  75,  33,  20,  34,  21, &
                             67,  62,  47,  63,  35,  70,  80,  57,  82,  81, &
                             83,  36,  91,  43,  76,  78,  79,  27,  61,  74, &
                             50,  23,  24,  51,  44,  28,  52,  38,  53,  46, &
                             58,  66,  71,   1,   2,  22,   3,  31,  17,   4, &
                              5,   6,   7,   8,   9,  10,  11,  12,  13,  14, &
                             15,  16 /)

      diag_map(: 92) = (/    1,   3,   4,   5,   6,   7,   9,  10,  11,  12, &
                            13,  14,  15,  16,  17,  18,  19,  21,  24,  27, &
                            30,  34,  36,  41,  44,  47,  51,  54,  58,  65, &
                            70,  75,  79,  84,  91,  96, 101, 108, 111, 116, &
                           121, 126, 131, 135, 142, 148, 152, 157, 162, 165, &
                           174, 182, 187, 194, 201, 205, 211, 219, 227, 232, &
                           236, 246, 257, 264, 270, 278, 292, 307, 317, 329, &
                           345, 355, 368, 377, 389, 399, 406, 412, 424, 438, &
                           450, 464, 482, 514, 537, 611, 659, 698, 727, 758, &
                           774, 795 /)

      het_lst(: 33) = (/ 'H2O2    ','HNO3    ','CH2O    ','CH3OOH  ','POOH    ', &
                         'CH3COOOH','HO2NO2  ','ONIT    ','MVK     ','MACR    ', &
                         'C2H5OOH ','C3H7OOH ','ROOH    ','CH3COCHO','Pb      ', &
                         'MACROOH ','XOOH    ','ONITR   ','ISOPOOH ','CH3OH   ', &
                         'C2H5OH  ','GLYALD  ','HYAC    ','HYDRALD ','CH3CHO  ', &
                         'ISOPNO3 ','ALKOOH  ','MEKOOH  ','TOLOOH  ','TERPOOH ', &
                         'CH3COOH ','SO2     ','NH3     ' /)

      extfrc_lst(:  2) = (/ 'NO      ','CO      ' /)

      frc_from_dataset(:  2) = (/ .false., .false. /)

      inv_lst(:  4) = (/ 'M       ', 'N2      ', 'O2      ', 'H2O     ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2             ', 'jo1d            ', 'jo3p            ', 'jn2o            ', &
                                     'jno2            ', 'jn2o5           ', 'jhno3           ', 'jno3_a          ', &
                                     'jno3_b          ', 'jho2no2_a       ', 'jho2no2_b       ', 'jch3ooh         ', &
                                     'jch2o_a         ', 'jch2o_b         ', 'jh2o2           ', 'jch3cho         ', &
                                     'jpooh           ', 'jch3co3h        ', 'jpan            ', 'jmpan           ', &
                                     'jmacr_a         ', 'jmacr_b         ', 'jmvk            ', 'jc2h5ooh        ', &
                                     'jc3h7ooh        ', 'jrooh           ', 'jacet           ', 'jmgly           ', &
                                     'jxooh           ', 'jonitr          ', 'jisopooh        ', 'jhyac           ', &
                                     'jglyald         ', 'jmek            ', 'jbigald         ', 'jglyoxal        ', &
                                     'jalkooh         ', 'jmekooh         ', 'jtolooh         ', 'jterpooh        ', &
                                     'usr_O_O2        ', 'o1d_n2          ', 'o1d_o2          ', 'ox_l1           ', &
                                     'ox_l2           ', 'ox_l3           ', 'usr_HO2_HO2     ', 'ox_p1           ', &
                                     'tag_NO2_NO3     ', 'usr_N2O5_M      ', 'tag_NO2_OH      ', 'usr_HNO3_OH     ', &
                                     'tag_NO2_HO2     ', 'usr_HO2NO2_M    ', 'usr_N2O5_aer    ', 'usr_NO3_aer     ', &
                                     'usr_NO2_aer     ', 'ox_p2           ', 'usr_CO_OH_a     ', 'tag_C2H4_OH     ', &
                                     'ox_l6           ', 'ox_p16          ', 'ox_p5           ', 'tag_C3H6_OH     ', &
                                     'ox_l4           ', 'ox_p3           ', 'ox_p4           ', 'tag_CH3CO3_NO2  ', &
                                     'usr_PAN_M       ', 'ox_p9           ', 'usr_CH3COCH3_OH ', 'ox_p10          ', &
                                     'ox_p15          ', 'soa5            ', 'ox_p14          ', 'ox_p17          ', &
                                     'soa4            ', 'ox_p12          ', 'ox_l5           ', 'ox_p6           ', &
                                     'ox_l7           ', 'ox_l8           ', 'ox_p7           ', 'ox_p8           ', &
                                     'usr_MCO3_NO2    ', 'usr_MPAN_M      ', 'soa2            ', 'soa1            ', &
                                     'soa3            ', 'ox_p13          ', 'ox_p11          ', 'usr_XOOH_OH     ', &
                                     'usr_SO2_OH      ', 'usr_DMS_OH      ', 'usr_HO2_aer     ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                                       11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                                       21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                                       31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                                       41,  43,  44,  45,  50,  51,  52,  59,  64,  65, &
                                       66,  67,  69,  71,  72,  73,  74,  77,  84,  85, &
                                       86,  87,  91,  96,  97,  99, 104, 105, 109, 112, &
                                      116, 117, 122, 123, 124, 129, 132, 135, 140, 141, &
                                      148, 150, 151, 158, 164, 165, 166, 167, 168, 169, &
                                      182, 188, 197, 199, 203 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ' /)
      pht_alias_lst(:,2) = (/ 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', '                ', &
                              '                ', 'jn2o5_a         ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'jch3ooh         ', 'jh2o2           ', '                ', 'jpan            ', &
                              '                ', '                ', '                ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3ooh         ', '                ', '                ', &
                              'jch3ooh         ', 'jch3cho         ', 'jch3ooh         ', '                ', &
                              '                ', 'jacet           ', 'jno2            ', 'jmgly           ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ', 'jch3ooh         ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, .28_r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, .2_r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
