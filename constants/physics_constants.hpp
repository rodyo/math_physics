/**
 *  @file          physics_constants.hpp
 *
 *  @par Description
 *                 TODO multiline description goes here
 *                 and continues here
 *
 *  @par Project
 *                 Gems
 *  @par Author(s)
 *                 Rody Oldenhuis (email: oldenhuis@luxspace.lu)
 *
 *  @date          30 Jan 2013
 *
 *  @version       0.2
 *
 *  @par Copyright
 *                 (C) 2013 Luxspace Sarl. <br><br>
 *                 You may not use this software, in whole or in part, in support
 *                 of any product without the express consent of the author.
 *                 The intellectual property rights of the algorithms used
 *                 here reside with the LuxSpace Sarl company.
 *
 *  @par Changes   Change history, in reverse chronological order
 *       @li       20 Mar 2013 - Mike de Dood - new templates
 *                 Class adapted to the new code templates
 *       @li       30 Jan 2013 - Rody Oldenhuis
 *                 File Creation
 *
 */

#ifndef __PHYS_CONSTANTS_HPP
#define __PHYS_CONSTANTS_HPP

#include "ScalarUnit.hpp"


/// Physics constants
namespace physics
{
    /**
     * Fundamental
     */

    const Speed    c = 299792458_meters_per_second; // speed of light
    const Radius  re = 2.817940285e-15_m;           // Classical electron radius
    const Charge   q = 1.60217646e-19_C;            // Proton charge

    const auto     k = 1.3806488e-23_J / 1_K;       // Boltzmann constant
    const auto     h = 6.62606957e-34_J * 1_s;      // Planck's constant
    const auto  hbar = 1.054571726e-34_J * 1_s;     // Reduced Plankc's constant/Dirac's constant
    const auto    NA = 6.02214129e+23L / 1_mol;     // Avogadro's constant
    const auto sigma = 5.670373e-8_W / 1_m2 / (1_K).pow<4>();  // Stefan-Boltzmann constant


    /**
     * Other
     */
    const Distance au = 149597870700_m;             // Average Sun-Earth distance (roughly)

} // namespace physics


#endif


