/**
 *  \file ScalarUnit.cpp
 *
 *  \details
 *
 *  \project:   EDRSSim
 *  \author(s): Rody Oldenhuis (email: oldenhuis@luxspace.lu)
 *  \version:   0.1
 *  \date:      07-03-2013 10:54:25 W. Europe Standard Time
 *
 *  \copyright 2013 Luxspace Sarl
 *
 *  You may not use this software, in whole or in part, in support
 *  of any product without the express consent of the author.
 *  The intellectual property rights of the algorithms used here
 *  reside with the LuxSpace Sarl company.
 *
 *  \par Change history, in reverse chronological order
 *
 */

#include "ScalarUnit.hpp"


/**
 *  Define the printable SI base unit for the type of physical unit mentioned.
 */


#define DEFINE_SI_BASE_UNIT(unit, str) \
    const std::string ScalarUnitTraits<unit>::SI_unit = str;


// FIXME: (Rody Oldenhuis)
// Displaying characters like ·, ², ³, Ω, etc. is difficult on Windows; it requires the user to change the active
// code page in the console to 1252 (default is the old OEM850, while the rest of windows uses 1252, just for your
// convenience.)

DEFINE_SI_BASE_UNIT(Length               , "m"    )
DEFINE_SI_BASE_UNIT(Mass                 , "kg"   )
DEFINE_SI_BASE_UNIT(Duration             , "s"    )
DEFINE_SI_BASE_UNIT(Current              , "A"    )
DEFINE_SI_BASE_UNIT(Temperature          , "K"    )
DEFINE_SI_BASE_UNIT(Speed                , "m/s"  )
DEFINE_SI_BASE_UNIT(AmountOfSubstance    , "mol"  )
DEFINE_SI_BASE_UNIT(AccelerationComponent, "m/s²" )
DEFINE_SI_BASE_UNIT(JerkComponent        , "m/s³" )
DEFINE_SI_BASE_UNIT(Area                 , "m²"   )
DEFINE_SI_BASE_UNIT(Volume               , "m³"   )
DEFINE_SI_BASE_UNIT(Density              , "kg/m³")
DEFINE_SI_BASE_UNIT(Energy               , "J"    )
DEFINE_SI_BASE_UNIT(SpecificEnergy       , "J/kg" )
DEFINE_SI_BASE_UNIT(Frequency            , "Hz"   )
DEFINE_SI_BASE_UNIT(ForceComponent       , "N"    )
DEFINE_SI_BASE_UNIT(MomentumComponent    , "N·s"  )
DEFINE_SI_BASE_UNIT(Pressure             , "Pa"   )
DEFINE_SI_BASE_UNIT(Charge               , "C"    )
DEFINE_SI_BASE_UNIT(Potential            , "V"    )
DEFINE_SI_BASE_UNIT(Power                , "W"    )
DEFINE_SI_BASE_UNIT(Resistance           , "Ω"    )
DEFINE_SI_BASE_UNIT(Conductance          , "S"    )
DEFINE_SI_BASE_UNIT(Capacitance          , "F"    )
DEFINE_SI_BASE_UNIT(Inductance           , "H"    )
DEFINE_SI_BASE_UNIT(MagneticFluxDensity  , "T"    )
DEFINE_SI_BASE_UNIT(TorqueComponent      , "N·m"  )


#undef DEFINE_SI_BASE_UNIT

