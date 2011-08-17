!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: Includes_After_Run.h
!
! !DESCRIPTION: This include file contains the array assignments that need
!  to be made AFTER the call to the Run method of the 
!  GEOSCHEMchem\_GridCompMod.F90 code.  These array assignments take the 
!  the locally-defined arrays (containing output data from the GEOS-Chem
!  column chemistry code) and saves them back into the internal and export
!  states.
!\\
!\\
!  These assignments were placed into this separate include file to avoid 
!  bogging down the GEOSCHEMchem\_GridCompMod.F90 module.
!\\
!\\
! !REVISION HISTORY: 
!  16 Apr 2010 - R. Yantosca - Initial version
!  19 Apr 2010 - R. Yantosca - Put data in export-state fields
!  06 May 2010 - R. Yantosca - Comment out dead species of CSPEC_1d
!  27 May 2010 - R. Yantosca - Bug fix for H2O2s and SO2s
!  02 Jun 2010 - R. Yantosca - Now deallocate Met%DQ*DTMST fields
!  02 Jun 2010 - R. Yantosca - Now deallocate Met%TAUCL* fields
!  01 Jul 2010 - R. Yantosca - Add references to D_OH_MASS and D_AIR_MASS
!  02 Jul 2010 - R. Yantosca - Use OH_SCALE to avoid overflow for mean OH
!  07 Jul 2010 - R. Yantosca - Fixed typo in assignment statement for H2O2s
!EOP
!------------------------------------------------------------------------------
!BOC
       !-------------------------------------------------------------------
       ! Put modified tracers back into the internal state
       !-------------------------------------------------------------------
       TRC_NOx      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%NOx  )
       TRC_Ox       ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%Ox   ) 
       TRC_PAN      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%PAN  ) 
       TRC_CO       ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%CO   ) 
       TRC_ALK4     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%ALK4 ) 
       TRC_ISOP     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%ISOP ) 
       TRC_HNO3     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%HNO3 ) 
       TRC_H2O2     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%H2O2 ) 
       TRC_ACET     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%ACET ) 
       TRC_MEK      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%MEK  ) 
       TRC_ALD2     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%ALD2 ) 
       TRC_RCHO     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%RCHO ) 
       TRC_MVK      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%MVK  ) 
       TRC_MACR     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%MACR ) 
       TRC_PMN      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%PMN  ) 
       TRC_PPN      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%PPN  ) 
       TRC_R4N2     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%R4N2 ) 
       TRC_PRPE     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%PRPE ) 
       TRC_C3H8     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%C3H8 ) 
       TRC_CH2O     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%CH2O ) 
       TRC_C2H6     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%C2H6 ) 
       TRC_N2O5     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%N2O5 ) 
       TRC_HNO4     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%HNO4 ) 
       TRC_MP       ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%MP   ) 
       TRC_DMS      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%DMS  ) 
       TRC_SO2      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SO2  ) 
       TRC_SO4      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SO4  ) 
       TRC_SO4s     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SO4s ) 
       TRC_MSA      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%MSA  ) 
       TRC_NH3      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%NH3  ) 
       TRC_NH4      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%NH4  ) 
       TRC_NIT      ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%NIT  ) 
       TRC_NITs     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%NITs ) 
       TRC_BCPI     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%BCPI ) 
       TRC_OCPI     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%OCPI ) 
       TRC_BCPO     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%BCPO ) 
       TRC_OCPO     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%OCPO ) 
       TRC_ALPH     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%ALPH ) 
       TRC_LIMO     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%LIMO ) 
       TRC_ALCO     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%ALCO ) 
       TRC_SOG1     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOG1 ) 
       TRC_SOG2     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOG2 ) 
       TRC_SOG3     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOG3 ) 
       TRC_SOG4     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOG4 ) 
       TRC_SOA1     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOA1 ) 
       TRC_SOA2     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOA2 ) 
       TRC_SOA3     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOA3 ) 
       TRC_SOA4     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SOA4 ) 
       TRC_DST1     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%DST1 ) 
       TRC_DST2     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%DST2 ) 
       TRC_DST3     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%DST3 ) 
       TRC_DST4     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%DST4 ) 
       TRC_SALA     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SALA ) 
       TRC_SALC     ( I, J, LM:1:-1 )  =  tracer_1d( 1:LM, IdTracers%SALC ) 

       !-------------------------------------------------------------------
       ! Put modified species back into the internal state
       !-------------------------------------------------------------------

       !%%% Active species %%%
       A3O2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%A3O2     )
       ACET         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ACET     )
       ALD2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ALD2     )
       ALK4         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ALK4     )
       ATO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ATO2     )
       B3O2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%B3O2     )
       C2H6         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%C2H6     )
       C3H8         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%C3H8     )
       CH2O         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%CH2O     )
       CO           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%CO       )
       DMS          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DMS      )
       DRYCH2O      ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYCH2O  )
       DRYH2O2      ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYH2O2  )
       DRYHNO3      ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYHNO3  )
       DRYN2O5      ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYN2O5  )
       DRYNO2       ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYNO2   )
       DRYO3        ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYO3    )
       DRYPAN       ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYPAN   )
       DRYPMN       ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYPMN   )
       DRYPPN       ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYPPN   )
       DRYR4N2      ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%DRYR4N2  )
       ETO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ETO2     )
       ETP          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ETP      )
       GCO3         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GCO3     )
       GLYC         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GLYC     )
       GP           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GP       )
       GPAN         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GPAN     )
       H2O2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%H2O2     )
       HAC          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%HAC      )
       HNO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%HNO2     )
       HNO3         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%HNO3     )
       HNO4         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%HNO4     )
       HO2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%HO2      )
       IALD         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%IALD     )
       IAO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%IAO2     )
       IAP          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%IAP      )
       INO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%INO2     )
       INPN         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%INPN     )
       ISN1         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ISN1     )
       ISNP         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ISNP     )
       ISOP         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ISOP     )
       KO2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%KO2      )
       LISOPOH      ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%LISOPOH  )
       MACR         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MACR     )
       MAN2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MAN2     )
       MAO3         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MAO3     )
       MAOP         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MAOP     )
       MAP          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MAP      )
       MCO3         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MCO3     )
       MEK          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MEK      )
       MGLY         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MGLY     )
       MO2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MO2      )
       MP           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MP       )
       MRO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MRO2     )
       MRP          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MRP      )
       SO4          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MSA      )
       MVK          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MVK      )
       MVN2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MVN2     )
       N2O5         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%N2O5     )
       NO           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%NO       )
       NO2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%NO2      )
       NO3          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%NO3      )
       O3           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%O3       )
       OH           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%OH       )
       PAN          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%PAN      )
       PMN          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%PMN      )
       PO2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%PO2      )
       PP           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%PP       )
       PRN1         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%PRN1     )
       PRPE         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%PRPE     )
       PRPN         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%PRPN     )
       R4N1         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%R4N1     )
       R4N2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%R4N2     )
       R4O2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%R4O2     )
       R4P          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%R4P      )
       RA3P         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RA3P     )
       RB3P         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RB3P     )
       RCHO         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RCHO     )
       RCO3         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RCO3     )
       RIO1         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RIO1     )
       RIO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RIO2     )
       RIP          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RIP      )
       RP           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%RP       )
       SO2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%SO2      )
       SO2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%SO4      )
       VRO2         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%VRO2     )
       VRP          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%VRP      )

       !%%% Inactive species in mechanism %%%
       ACTA         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ACTA     )
       CH4          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%CH4      )
       EOH          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%EOH      )
       EMISSION     ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%EMISSION )
       GLCO3        ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GLCO3    )
       GLP          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GLP      )
       GLPAN        ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GLPAN    )
       GLYX         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%GLYX     )
       H            ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%H        )
       H2           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%H2       )
       H2O          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%H2O      )
       HCOOH        ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%HCOOH    )
       ISNO3        ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ISNO3    )
       M            ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%M        )
       MOH          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MOH      )
       MNO3         ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%MNO3     )
       N2           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%N2       )
       NH2          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%NH2      )
       NH3          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%NH3      )
       O            ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%O        )
       O2           ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%O2       )
       O2CH2OH      ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%O2CH2OH  )
       ROH          ( I, J, LM:1:-1 )  =  cspec_1d( 1:LM, IdSpecies%ROH      )
     
       !%%% Dead species in mechanism, comment out (bmy, 4/19/10) %%%
      !CO2          ( I, J, LM:1:-1 )  = cspec_1d( 1:LM, IdSpecies%CO2       ) 
      !DRYDEP       ( I, J, LM:1:-1 )  = cspec_1d( 1:LM, IdSpecies%DRYDEP    )
      !N2O          ( I, J, LM:1:-1 )  = cspec_1d( 1:LM, IdSpecies%N2O       )
      !O1D          ( I, J, LM:1:-1 )  = cspec_1d( 1:LM, IdSpecies%O1D       )
      !PPN          ( I, J, LM:1:-1 )  = cspec_1d( 1:LM, IdSpecies%PPN       )

       !--------------------------------------------------------------------
       ! Put wetdep parameters back into the internal state
       !--------------------------------------------------------------------
       H2O2s        ( I, J, LM:1:-1 )  =  h2o2s_1d  ( 1:LM    )
       SO2s         ( I, J, LM:1:-1 )  =  so2s_1d   ( 1:LM    )

       !--------------------------------------------------------------------
       ! Put seasalt parameters back into the internal state
       !--------------------------------------------------------------------
       ALK_EMIS_SALA( I, J, LM:1:-1 )  =  alkEmis_1d( 1:LM, 1 )
       ALK_EMIS_SALC( I, J, LM:1:-1 )  =  alkEmis_1d( 1:LM, 2 ) 
       N_DENS_SALA  ( I, J, LM:1:-1 )  =  nDens_1d  ( 1:LM, 1 )
       N_DENS_SALC  ( I, J, LM:1:-1 )  =  nDens_1d  ( 1:LM, 2 )

       !-------------------------------------------------------------------
       ! Put SOA parameters back into the internal state
       !-------------------------------------------------------------------

       ! TERP and SESQ emissions
       ORVC_TERP    ( I, J, LM:1:-1 )  =  terp_1d( 1:LM )
       ORVC_SESQ    ( I, J, LM:1:-1 )  =  sesq_1d( 1:LM )
                                          
       ! GPROD parameters                 
       GSOAP11      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 1, 1 )
       GSOAP21      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 2, 1 )
       GSOAP31      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 3, 1 )
       GSOAP12      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 1, 2 )
       GSOAP22      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 2, 2 )
       GSOAP32      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 3, 2 )
       GSOAP13      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 1, 3 )
       GSOAP23      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 2, 3 )
       GSOAP33      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 3, 3 )
       GSOAP14      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 1, 4 )
       GSOAP24      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 2, 4 )
       GSOAP34      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 3, 4 )
       GSOAP15      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 1, 5 )
       GSOAP25      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 2, 5 )
       GSOAP35      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 3, 5 )
       GSOAP16      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 1, 6 )
       GSOAP26      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 2, 6 )
       GSOAP36      ( I, J, LM:1:-1 )  =  gProd_1d( 1:LM, 3, 6 )
                                          
       ! GPROD parameters                 
       ASOAP11      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 1, 1 )
       ASOAP21      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 2, 1 )
       ASOAP31      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 3, 1 )
       ASOAP12      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 1, 2 )
       ASOAP22      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 2, 2 )
       ASOAP32      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 3, 2 )
       ASOAP13      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 1, 3 )
       ASOAP23      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 2, 3 )
       ASOAP33      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 3, 3 )
       ASOAP14      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 1, 4 )
       ASOAP24      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 2, 4 )
       ASOAP34      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 3, 4 )
       ASOAP15      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 1, 5 )
       ASOAP25      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 2, 5 )
       ASOAP35      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 3, 5 )
       ASOAP16      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 1, 6 )
       ASOAP26      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 2, 6 )
       ASOAP36      ( I, J, LM:1:-1 )  =  aProd_1d( 1:LM, 3, 6 )

       !-------------------------------------------------------------------
       ! Put mean OH lifetime parameters back into the internal state
       !
       ! NOTE: D_AIR_MASS and D_OH_MASS are REAL*4, so we need to make
       ! sure that we do not exceed 1e38, which is the max allowable value.
       ! Therefore, D_AIR_MASS and D_OH_MASS will store the values divided
       ! by 1d20 to prevent this overflow situation.
       !-------------------------------------------------------------------
       D_AIR_MASS   ( I, J, LM:1:-1 )  =  airMassDiag( 1:LM ) / OH_SCALE
       D_OH_MASS    ( I, J, LM:1:-1 )  =  ohMassDiag ( 1:LM ) / OH_SCALE

!       !-------------------------------------------------------------------
!       ! Put wet deposition losses into the export state
!       !-------------------------------------------------------------------
!       WD_HNO3( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%HNO3 )    
!       WD_H2O2( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%H2O2 )   
!       WD_CH2O( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%CH2O )   
!       WD_MP  ( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%MP   )   
!       WD_SO2 ( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SO2  )   
!       WD_SO4 ( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SO4  )   
!       WD_SO4s( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SO4s )   
!       WD_MSA ( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%MSA  )   
!       WD_NH3 ( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%NH3  )   
!       WD_NH4 ( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%NH4  )   
!       WD_NIT ( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%NIT  )   
!       WD_NITs( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%NITs )   
!       WD_BCPI( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%BCPI )   
!       WD_OCPI( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%OCPI )   
!       WD_BCPO( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%BCPO )   
!       WD_OCPO( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%OCPO )   
!       WD_ALPH( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%ALPH )   
!       WD_LIMO( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%LIMO )   
!       WD_ALCO( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%ALCO )   
!       WD_SOG1( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOG1 )   
!       WD_SOG2( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOG2 )   
!       WD_SOG3( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOG3 )   
!       WD_SOG4( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOG4 )   
!       WD_SOA1( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOA1 )   
!       WD_SOA2( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOA2 )   
!       WD_SOA3( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOA3 )   
!       WD_SOA4( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SOA4 )   
!       WD_DST1( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%DST1 )   
!       WD_DST2( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%DST2 )   
!       WD_DST3( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%DST3 )   
!       WD_DST4( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%DST4 )   
!       WD_SALA( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SALA )   
!       WD_SALC( I, J, LM:1:-1 )       = wdLossDiag( 1:LM, IdTracers%SALC )   
!
!       !-------------------------------------------------------------------
!       ! Put dry deposition fluxes into the export state
!       !-------------------------------------------------------------------
!       DD_FX_NOx ( I, J )             = ddFluxDiag( IdTracers%NOx  )     
!       DD_FX_O3  ( I, J )             = ddFluxDiag( IdTracers%Ox   ) 
!       DD_FX_PAN ( I, J )             = ddFluxDiag( IdTracers%PAN  ) 
!       DD_FX_HNO3( I, J )             = ddFluxDiag( IdTracers%HNO3 ) 
!       DD_FX_H2O2( I, J )             = ddFluxDiag( IdTracers%H2O2 ) 
!       DD_FX_N2O5( I, J )             = ddFluxDiag( IdTracers%N2O5 ) 
!       DD_FX_PPN ( I, J )             = ddFluxDiag( IdTracers%PPN  ) 
!       DD_FX_R4N2( I, J )             = ddFluxDiag( IdTracers%R4N2 ) 
!       DD_FX_CH2O( I, J )             = ddFluxDiag( IdTracers%CH2O ) 
!       DD_FX_SO2 ( I, J )             = ddFluxDiag( IdTracers%SO2  ) 
!       DD_FX_SO4 ( I, J )             = ddFluxDiag( IdTracers%SO4  ) 
!       DD_FX_SO4S( I, J )             = ddFluxDiag( IdTracers%SO4S ) 
!       DD_FX_MSA ( I, J )             = ddFluxDiag( IdTracers%MSA  ) 
!       DD_FX_NH3 ( I, J )             = ddFluxDiag( IdTracers%NH3  ) 
!       DD_FX_NH4 ( I, J )             = ddFluxDiag( IdTracers%NH4  ) 
!       DD_FX_NIT ( I, J )             = ddFluxDiag( IdTracers%NIT  ) 
!       DD_FX_NITS( I, J )             = ddFluxDiag( IdTracers%NITS ) 
!       DD_FX_BCPI( I, J )             = ddFluxDiag( IdTracers%BCPI ) 
!       DD_FX_OCPI( I, J )             = ddFluxDiag( IdTracers%BCPO ) 
!       DD_FX_BCPO( I, J )             = ddFluxDiag( IdTracers%OCPI ) 
!       DD_FX_OCPO( I, J )             = ddFluxDiag( IdTracers%OCPO ) 
!       DD_FX_ALPH( I, J )             = ddFluxDiag( IdTracers%ALPH ) 
!       DD_FX_LIMO( I, J )             = ddFluxDiag( IdTracers%LIMO ) 
!       DD_FX_ALCO( I, J )             = ddFluxDiag( IdTracers%ALCO ) 
!       DD_FX_SOG1( I, J )             = ddFluxDiag( IdTracers%SOG1 ) 
!       DD_FX_SOG2( I, J )             = ddFluxDiag( IdTracers%SOG2 ) 
!       DD_FX_SOG3( I, J )             = ddFluxDiag( IdTracers%SOG3 ) 
!       DD_FX_SOG4( I, J )             = ddFluxDiag( IdTracers%SOG4 ) 
!       DD_FX_SOA1( I, J )             = ddFluxDiag( IdTracers%SOA1 ) 
!       DD_FX_SOA2( I, J )             = ddFluxDiag( IdTracers%SOA2 ) 
!       DD_FX_SOA3( I, J )             = ddFluxDiag( IdTracers%SOA3 ) 
!       DD_FX_SOA4( I, J )             = ddFluxDiag( IdTracers%SOA4 ) 
!       DD_FX_DST1( I, J )             = ddFluxDiag( IdTracers%DST1 ) 
!       DD_FX_DST2( I, J )             = ddFluxDiag( IdTracers%DST2 ) 
!       DD_FX_DST3( I, J )             = ddFluxDiag( IdTracers%DST3 ) 
!       DD_FX_DST4( I, J )             = ddFluxDiag( IdTracers%DST4 ) 
!       DD_FX_SALA( I, J )             = ddFluxDiag( IdTracers%SALA ) 
!       DD_FX_SALC( I, J )             = ddFluxDiag( IdTracers%SALC ) 
!       
!       !-------------------------------------------------------------------
!       ! Put dry deposition frequencies into the export state
!       !-------------------------------------------------------------------
!       DD_FQ_NOx ( I, J )             = ddFreqDiag( IdTracers%NOx  )     
!       DD_FQ_O3  ( I, J )             = ddFreqDiag( IdTracers%Ox   ) 
!       DD_FQ_PAN ( I, J )             = ddFreqDiag( IdTracers%PAN  ) 
!       DD_FQ_HNO3( I, J )             = ddFreqDiag( IdTracers%HNO3 ) 
!       DD_FQ_H2O2( I, J )             = ddFreqDiag( IdTracers%H2O2 ) 
!       DD_FQ_N2O5( I, J )             = ddFreqDiag( IdTracers%N2O5 ) 
!       DD_FQ_PPN ( I, J )             = ddFreqDiag( IdTracers%PPN  ) 
!       DD_FQ_R4N2( I, J )             = ddFreqDiag( IdTracers%R4N2 ) 
!       DD_FQ_CH2O( I, J )             = ddFreqDiag( IdTracers%CH2O ) 
!       DD_FQ_SO2 ( I, J )             = ddFreqDiag( IdTracers%SO2  ) 
!       DD_FQ_SO4 ( I, J )             = ddFreqDiag( IdTracers%SO4  ) 
!       DD_FQ_SO4S( I, J )             = ddFreqDiag( IdTracers%SO4S ) 
!       DD_FQ_MSA ( I, J )             = ddFreqDiag( IdTracers%MSA  ) 
!       DD_FQ_NH3 ( I, J )             = ddFreqDiag( IdTracers%NH3  ) 
!       DD_FQ_NH4 ( I, J )             = ddFreqDiag( IdTracers%NH4  ) 
!       DD_FQ_NIT ( I, J )             = ddFreqDiag( IdTracers%NIT  ) 
!       DD_FQ_NITS( I, J )             = ddFreqDiag( IdTracers%NITS ) 
!       DD_FQ_BCPI( I, J )             = ddFreqDiag( IdTracers%BCPI ) 
!       DD_FQ_OCPI( I, J )             = ddFreqDiag( IdTracers%BCPO ) 
!       DD_FQ_BCPO( I, J )             = ddFreqDiag( IdTracers%OCPI ) 
!       DD_FQ_OCPO( I, J )             = ddFreqDiag( IdTracers%OCPO ) 
!       DD_FQ_ALPH( I, J )             = ddFreqDiag( IdTracers%ALPH ) 
!       DD_FQ_LIMO( I, J )             = ddFreqDiag( IdTracers%LIMO ) 
!       DD_FQ_ALCO( I, J )             = ddFreqDiag( IdTracers%ALCO ) 
!       DD_FQ_SOG1( I, J )             = ddFreqDiag( IdTracers%SOG1 ) 
!       DD_FQ_SOG2( I, J )             = ddFreqDiag( IdTracers%SOG2 ) 
!       DD_FQ_SOG3( I, J )             = ddFreqDiag( IdTracers%SOG3 ) 
!       DD_FQ_SOG4( I, J )             = ddFreqDiag( IdTracers%SOG4 ) 
!       DD_FQ_SOA1( I, J )             = ddFreqDiag( IdTracers%SOA1 ) 
!       DD_FQ_SOA2( I, J )             = ddFreqDiag( IdTracers%SOA2 ) 
!       DD_FQ_SOA3( I, J )             = ddFreqDiag( IdTracers%SOA3 ) 
!       DD_FQ_SOA4( I, J )             = ddFreqDiag( IdTracers%SOA4 ) 
!       DD_FQ_DST1( I, J )             = ddFreqDiag( IdTracers%DST1 ) 
!       DD_FQ_DST2( I, J )             = ddFreqDiag( IdTracers%DST2 ) 
!       DD_FQ_DST3( I, J )             = ddFreqDiag( IdTracers%DST3 ) 
!       DD_FQ_DST4( I, J )             = ddFreqDiag( IdTracers%DST4 ) 
!       DD_FQ_SALA( I, J )             = ddFreqDiag( IdTracers%SALA ) 
!       DD_FQ_SALC( I, J )             = ddFreqDiag( IdTracers%SALC )
!
!       !-------------------------------------------------------------------
!       ! Put dry deposition velocities into the export state
!       !-------------------------------------------------------------------
!       DD_V_NOx ( I, J )              = ddVelDiag( IdTracers%NOx  )     
!       DD_V_O3  ( I, J )              = ddVelDiag( IdTracers%Ox   ) 
!       DD_V_PAN ( I, J )              = ddVelDiag( IdTracers%PAN  ) 
!       DD_V_HNO3( I, J )              = ddVelDiag( IdTracers%HNO3 ) 
!       DD_V_H2O2( I, J )              = ddVelDiag( IdTracers%H2O2 ) 
!       DD_V_N2O5( I, J )              = ddVelDiag( IdTracers%N2O5 ) 
!       DD_V_PPN ( I, J )              = ddVelDiag( IdTracers%PPN  ) 
!       DD_V_R4N2( I, J )              = ddVelDiag( IdTracers%R4N2 ) 
!       DD_V_CH2O( I, J )              = ddVelDiag( IdTracers%CH2O ) 
!       DD_V_SO2 ( I, J )              = ddVelDiag( IdTracers%SO2  ) 
!       DD_V_SO4 ( I, J )              = ddVelDiag( IdTracers%SO4  ) 
!       DD_V_SO4S( I, J )              = ddVelDiag( IdTracers%SO4S ) 
!       DD_V_MSA ( I, J )              = ddVelDiag( IdTracers%MSA  ) 
!       DD_V_NH3 ( I, J )              = ddVelDiag( IdTracers%NH3  ) 
!       DD_V_NH4 ( I, J )              = ddVelDiag( IdTracers%NH4  ) 
!       DD_V_NIT ( I, J )              = ddVelDiag( IdTracers%NIT  ) 
!       DD_V_NITS( I, J )              = ddVelDiag( IdTracers%NITS ) 
!       DD_V_BCPI( I, J )              = ddVelDiag( IdTracers%BCPI ) 
!       DD_V_OCPI( I, J )              = ddVelDiag( IdTracers%BCPO ) 
!       DD_V_BCPO( I, J )              = ddVelDiag( IdTracers%OCPI ) 
!       DD_V_OCPO( I, J )              = ddVelDiag( IdTracers%OCPO ) 
!       DD_V_ALPH( I, J )              = ddVelDiag( IdTracers%ALPH ) 
!       DD_V_LIMO( I, J )              = ddVelDiag( IdTracers%LIMO ) 
!       DD_V_ALCO( I, J )              = ddVelDiag( IdTracers%ALCO ) 
!       DD_V_SOG1( I, J )              = ddVelDiag( IdTracers%SOG1 ) 
!       DD_V_SOG2( I, J )              = ddVelDiag( IdTracers%SOG2 ) 
!       DD_V_SOG3( I, J )              = ddVelDiag( IdTracers%SOG3 ) 
!       DD_V_SOG4( I, J )              = ddVelDiag( IdTracers%SOG4 ) 
!       DD_V_SOA1( I, J )              = ddVelDiag( IdTracers%SOA1 ) 
!       DD_V_SOA2( I, J )              = ddVelDiag( IdTracers%SOA2 ) 
!       DD_V_SOA3( I, J )              = ddVelDiag( IdTracers%SOA3 ) 
!       DD_V_SOA4( I, J )              = ddVelDiag( IdTracers%SOA4 ) 
!       DD_V_DST1( I, J )              = ddVelDiag( IdTracers%DST1 ) 
!       DD_V_DST2( I, J )              = ddVelDiag( IdTracers%DST2 ) 
!       DD_V_DST3( I, J )              = ddVelDiag( IdTracers%DST3 ) 
!       DD_V_DST4( I, J )              = ddVelDiag( IdTracers%DST4 ) 
!       DD_V_SALA( I, J )              = ddVelDiag( IdTracers%SALA ) 
!       DD_V_SALC( I, J )              = ddVelDiag( IdTracers%SALC )
!
       !-------------------------------------------------------------------
       ! Free all column pointer fields
       !-------------------------------------------------------------------
       IF ( ASSOCIATED( Met%AD       ) ) DEALLOCATE( Met%AD       )
       IF ( ASSOCIATED( Met%AIRDENS  ) ) DEALLOCATE( Met%AIRDENS  )
       IF ( ASSOCIATED( Met%AIRVOL   ) ) DEALLOCATE( Met%AIRVOL   )
       IF ( ASSOCIATED( Met%BXHEIGHT ) ) DEALLOCATE( Met%BXHEIGHT )
       IF ( ASSOCIATED( Met%CLDF     ) ) DEALLOCATE( Met%CLDF     )
       IF ( ASSOCIATED( Met%CMFMC    ) ) DEALLOCATE( Met%CMFMC    )
       IF ( ASSOCIATED( Met%DQIDTMST ) ) DEALLOCATE( Met%DQIDTMST )
       IF ( ASSOCIATED( Met%DQLDTMST ) ) DEALLOCATE( Met%DQLDTMST )
       IF ( ASSOCIATED( Met%DQVDTMST ) ) DEALLOCATE( Met%DQVDTMST )
       IF ( ASSOCIATED( Met%DTRAIN   ) ) DEALLOCATE( Met%DTRAIN   )
       IF ( ASSOCIATED( Met%MOISTQ   ) ) DEALLOCATE( Met%MOISTQ   )
       IF ( ASSOCIATED( Met%OPTD     ) ) DEALLOCATE( Met%OPTD     )
       IF ( ASSOCIATED( Met%PEDGE    ) ) DEALLOCATE( Met%PEDGE    )
       IF ( ASSOCIATED( Met%PMID     ) ) DEALLOCATE( Met%PMID     )
       IF ( ASSOCIATED( Met%RH       ) ) DEALLOCATE( Met%RH       )
       IF ( ASSOCIATED( Met%SPHU     ) ) DEALLOCATE( Met%SPHU     )
       IF ( ASSOCIATED( Met%T        ) ) DEALLOCATE( Met%T        )
       IF ( ASSOCIATED( Met%TAUCLI   ) ) DEALLOCATE( Met%TAUCLI   )
       IF ( ASSOCIATED( Met%TAUCLW   ) ) DEALLOCATE( Met%TAUCLW   )
       IF ( ASSOCIATED( Schem%OH     ) ) DEALLOCATE( Schem%OH     )
       IF ( ASSOCIATED( Schem%PCO    ) ) DEALLOCATE( Schem%PCO    )
       IF ( ASSOCIATED( Schem%LCO    ) ) DEALLOCATE( Schem%LCO    )
       IF ( ASSOCIATED( Schem%JVALUE ) ) DEALLOCATE( Schem%JVALUE )

