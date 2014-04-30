/**
 *  \file ScalarUnit.hpp
 *
 *  \details
 *
 *  \project:   EDRSSim
 *  \author(s): Rody Oldenhuis (email: oldenhuis@luxspace.lu)
 *  \version:   0.1
 *  \date:      07-03-2013 10:41:09 W. Europe Standard Time
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


#ifndef __SCALARUNIT_HPP
#define __SCALARUNIT_HPP


#include <iostream>
#include <assert.h>
#include <cmath>
#include <string>
#include <sstream>

#include "math_constants.hpp"

// FIXME: HACK HACK HACK HACK
// See
// - http://stackoverflow.com/questions/23365562/
// - http://stackoverflow.com/questions/12721584/
#undef _N  // Newtons
#undef _C  // Coulombs


// NOTE: This is the only convenient, readable, manageable way to define all
// units; see the unit definitions towards the end for an explanation.
//
// See also   http://stackoverflow.com/questions/15196955/



#define DECLARE_UNIT_TRAITS(U) \
    template<> struct ScalarUnitTraits<U> { static const std::string SI_unit; };


#define DEFINE_LITERALS(unit, suffix, multiplier) \
    constexpr unit operator"" _##suffix (long double V)        { return unit(V*multiplier); } \
    constexpr unit operator"" _##suffix (unsigned long long V) { return unit(static_cast<long double>(V)*multiplier); }

#define DEFINE_SI_MULTIPLIER(unit, suffix, power, prefix, multiplier) \
    constexpr unit operator"" _##prefix##suffix(long double V)        { return operator"" _##suffix( pow(multiplier, power) * V  ); } \
    constexpr unit operator"" _##prefix##suffix(unsigned long long V) { return operator"" _##suffix( pow(multiplier, power) * static_cast<long double>(V) ); }

#define DEFINE_SHORT_SI_MULTIPLIERS(unit, suffix, power) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     y, 1.0e-24L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     Y, 1.0e+24L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     z, 1.0e-21L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     Z, 1.0e+21L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     a, 1.0e-18L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     E, 1.0e+18L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     f, 1.0e-15L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     P, 1.0e+15L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     p, 1.0e-12L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     T, 1.0e+12L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     n, 1.0e-09L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     G, 1.0e+09L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     u, 1.0e-06L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     M, 1.0e+06L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     m, 1.0e-03L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     k, 1.0e+03L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     c, 1.0e-02L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,     h, 1.0e+02L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     d, 1.0e-01L)       DEFINE_SI_MULTIPLIER (unit,suffix,power,    da, 1.0e+01L)

#define DEFINE_LONG_SI_MULTIPLIERS(unit, suffix, power) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, yocto, 1.0e-24L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, yotta, 1.0e+24L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, zepto, 1.0e-21L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, zetta, 1.0e+21L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, atto , 1.0e-18L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, exa  , 1.0e+18L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, femto, 1.0e-15L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, peta , 1.0e+15L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, pico , 1.0e-12L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, tera , 1.0e+12L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, nano , 1.0e-09L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, giga , 1.0e+09L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, micro, 1.0e-06L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, mega , 1.0e+06L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, milli, 1.0e-03L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, kilo , 1.0e+03L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, centi, 1.0e-02L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, hecto, 1.0e+02L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, deci , 1.0e-01L)       DEFINE_SI_MULTIPLIER (unit,suffix,power, deca , 1.0e+01L)

#define DEFINE_UNIT_SHORT(unit, suffix, multiplier, power) \
    DEFINE_LITERALS(unit, suffix, multiplier) \
    DEFINE_SHORT_SI_MULTIPLIERS(unit, suffix, power)

#define DEFINE_UNIT_LONG(unit, suffix, multiplier, power) \
    DEFINE_LITERALS(unit, suffix, multiplier) \
    DEFINE_LONG_SI_MULTIPLIERS(unit, suffix, power)


// This one should be kept defined
#define convertTo(y) convert(1.0_##y)


// Small traits class to get the SI unit name in (string template argmuents are not supported)
// TODO: Think more about this one...
template<class S>
struct ScalarUnitTraits {
    static const std::string SI_unit;
};


// All units are a permutation of the base units
// TODO: how about things like square-root-Hertz etc.? -> need doubles!
template<
    int L,  //  Length
    int M,  //  Mass
    int t,  //  Time
    int C,  //  Current
    int T,  //  Temperature
    int I,  //  Luminous intensity
    int N,  //  amount of Substance
    unsigned char i=0u // Free counter to resolve unit ambiguities
>
class ScalarUnit
{
  private:

    // the actual value
    long double value;

    // generic pretty-printer
    constexpr std::string toString(const std::string &unit) const {
        std::stringstream output;
        output << value << " " << unit;
        return output.str();
    }


  public:

    // some boilerplate
    constexpr ScalarUnit(const ScalarUnit &x) : value(x.value) {}

    // constructor
    constexpr explicit ScalarUnit(long double value = 0.0L) : value(value) {}


    // TODO: (Rody Oldenhuis) need to optimize; implementations are not NVRO friendly

    // Same-unit operator overloading
    constexpr ScalarUnit& operator+= (const ScalarUnit &rhs) noexcept { value += rhs.value; return *this; }
    constexpr ScalarUnit  operator+  (const ScalarUnit &rhs) const    { return ScalarUnit(*this) += rhs;  }
    constexpr ScalarUnit& operator-= (const ScalarUnit &rhs) noexcept { value -= rhs.value; return *this; }
    constexpr ScalarUnit  operator-  (const ScalarUnit &rhs) const    { return ScalarUnit(*this) -= rhs;  }


    // multiply/divide by double
    constexpr ScalarUnit& operator/= (long double rhs) noexcept { value /= rhs; return *this;      }
    constexpr ScalarUnit  operator/  (long double rhs) const    { return ScalarUnit(*this) /= rhs; }
    constexpr ScalarUnit& operator*= (long double rhs) noexcept { value *= rhs; return *this;      }
    constexpr ScalarUnit  operator*  (long double rhs) const    { return ScalarUnit(*this) *= rhs; }

    // multiply/divide by anything else
    template<typename D> constexpr ScalarUnit& operator/= (D rhs)       { return *this /= static_cast<long double>(rhs); }
    template<typename D> constexpr ScalarUnit  operator/  (D rhs) const { return *this /  static_cast<long double>(rhs); }
    template<typename D> constexpr ScalarUnit& operator*= (D rhs)       { return *this *= static_cast<long double>(rhs); }
    template<typename D> constexpr ScalarUnit  operator*  (D rhs) const { return *this *  static_cast<long double>(rhs); }

    // raise to (integer) power
    template<int P>
    constexpr ScalarUnit<L*P,M*P,t*P,C*P,T*P,I*P,N*P,i> pow() {
        return ScalarUnit<L*P,M*P,t*P,C*P,T*P,I*P,N*P,i>( std::pow(value, P) );
    }


    // multiply/divide by ScalarUnit (with proper typedefs, allows for valid converions between units)
    // FIXME: how to handle the i-counter properly?

    template <int L2, int M2, int t2, int C2, int T2, int I2, int N2>
    constexpr ScalarUnit<L+L2,M+M2,t+t2,C+C2,T+T2,I+I2,N+N2>
    operator* (const ScalarUnit<L2,M2,t2,C2,T2,I2,N2> &rhs) const {
        return ScalarUnit<L+L2,M+M2,t+t2,C+C2,T+T2,I+I2,N+N2> ( value*rhs.getValue() );
    }

    template <int L2, int M2, int t2, int C2, int T2, int I2, int N2>
    constexpr ScalarUnit<L-L2,M-M2,t-t2,C-C2,T-T2,I-I2,N-N2>
    operator/ (const ScalarUnit<L2,M2,t2,C2,T2,I2,N2> &rhs) const {
        return ScalarUnit<L-L2,M-M2,t-t2,C-C2,T-T2,I-I2,N-N2> ( value/rhs.getValue() );
    }
    // (special case)
    //long double operator/ (const ScalarUnit &rhs) const noexcept { return value/rhs.value; }


    // compare ScalarUnits
    constexpr bool operator==(const ScalarUnit &rhs) const noexcept { return value == rhs.value; }
    constexpr bool operator< (const ScalarUnit &rhs) const noexcept { return value <  rhs.value; }
    constexpr bool operator> (const ScalarUnit &rhs) const noexcept { return value >  rhs.value; }
    constexpr bool operator<=(const ScalarUnit &rhs) const noexcept { return value <= rhs.value; }
    constexpr bool operator>=(const ScalarUnit &rhs) const noexcept { return value >= rhs.value; }


    // Convert units
    constexpr long double convert(const ScalarUnit &rhs) const noexcept { return value/rhs.value; }

    // getter for value
    constexpr long double getValue() const noexcept { return value; }

    // DISALLOW implicit conversions; leave unimplemented
    //operator double() const;
   // operator long double() const;

    // toString (for lexical casts)
    constexpr std::string toString() const {
        return toString(ScalarUnitTraits<ScalarUnit<L,M,t,C,T,I,N,i>>::SI_unit);
    }

    //

};


// need opererators "outside" of class to allow these:

template <int L, int M, int t, int C, int T, int I, int N, unsigned char i>
constexpr ScalarUnit<-L,-M,-t,-C,-T,-I,-N,i> operator/ (const long double &lhs, const ScalarUnit<L,M,t,C,T,I,N,i> &rhs) {
    return ScalarUnit<-L,-M,-t,-C,-T,-I,-N,i>( lhs/rhs.getValue() );
}


template <int L, int M, int t, int C, int T, int I, int N, unsigned char i>
constexpr ScalarUnit<L,M,t,C,T,I,N,i> operator*(const long double &lhs, const ScalarUnit<L,M,t,C,T,I,N,i> &rhs) {
    return rhs*lhs;
}

template <int L, int M, int t, int C, int T, int I, int N, unsigned char i>
std::ostream& operator<<(std::ostream &target, const ScalarUnit<L,M,t,C,T,I,N,i> &rhs) {
    target << rhs.toString();
    return target;
}




/**
 *  Define the actual physical units
 */

// L: Length
// M: Mass
// t: Time
// C: Current
// T: Temperature
// I: Luminous intensity
// N: Amount of substance
// i: (free counter)


// Base units        L   M   t   C   T   I   N   i
typedef ScalarUnit <+0, +0, +0, +0, +0, +0, +0, +0u>  Angle; // == dimensionless
typedef ScalarUnit <+1, +0, +0, +0, +0, +0, +0, +0u>  Length;
typedef ScalarUnit <+0, +1, +0, +0, +0, +0, +0, +0u>  Mass;
typedef ScalarUnit <+0, +0, +1, +0, +0, +0, +0, +0u>  Duration;
typedef ScalarUnit <+0, +0, +0, +1, +0, +0, +0, +0u>  Current;
typedef ScalarUnit <+0, +0, +0, +0, +1, +0, +0, +0u>  Temperature;
typedef ScalarUnit <+0, +0, +0, +0, +0, +1, +0, +0u>  LuminousIntensity;
typedef ScalarUnit <+0, +0, +0, +0, +0, +0, +1, +0u>  AmountOfSubstance;

// Derived units     L   M   t   C   T   I   N   i
typedef ScalarUnit <+0, +0, -1, +0, +0, +0, +0, +0u>  AngularSpeed;
typedef ScalarUnit <+0, +0, -2, +0, +0, +0, +0, +0u>  AngularAcceleration;
typedef ScalarUnit <+2, +1, -1, +0, +0, +0, +0, +0u>  AngularMomentum;
typedef ScalarUnit <+2, +0, -1, +0, +0, +0, +0, +0u>  SpecificAngularMomentum;

//                   L   M   t   C   T   I   N   i
typedef ScalarUnit <+1, +0, -1, +0, +0, +0, +0, +0u>  Speed;
typedef ScalarUnit <+1, +0, -2, +0, +0, +0, +0, +0u>  Acceleration;
typedef ScalarUnit <+1, +0, -3, +0, +0, +0, +0, +0u>  Jerk;

//                   L   M   t   C   T   I   N   i
typedef ScalarUnit <+2, +0, +0, +0, +0, +0, +0, +0u>  Area;
typedef ScalarUnit <+3, +0, +0, +0, +0, +0, +0, +0u>  Volume;
typedef ScalarUnit <-3, +1, +0, +0, +0, +0, +0, +0u>  Density;

//                   L   M   t   C   T   I   N   i
typedef ScalarUnit <+1, +1, -2, +0, +0, +0, +0, +0u>  Force;
typedef ScalarUnit <+2, +1, -2, +0, +0, +0, +0, +0u>  Torque;
typedef ScalarUnit <+1, +1, -1, +0, +0, +0, +0, +0u>  Momentum;
typedef ScalarUnit <+1, +0, -1, +0, +0, +0, +0, +0u>  SpecificMomentum;
typedef ScalarUnit <+2, +1, -2, +0, +0, +0, +0, +1u>  Energy;    // NOTE: conflicts with Torque
typedef ScalarUnit <+2, +0, -2, +0, +0, +0, +0, +0u>  SpecificEnergy;
typedef ScalarUnit <+0, +0, -1, +0, +0, +0, +0, +1u>  Frequency; // NOTE: conflicts with AngularSpeed
typedef ScalarUnit <-1, +1, -2, +0, +0, +0, +0, +0u>  Pressure;
typedef ScalarUnit <+2, +1, +0, +0, +0, +0, +0, +0u>  MomentOfInertia;

//                   L   M   t   C   T   I   N   i
typedef ScalarUnit <+0, +0, +1, +1, +0, +0, +0, +0u>  Charge;
typedef ScalarUnit <+2, +1, -3, -1, +0, +0, +0, +0u>  Potential;
typedef ScalarUnit <+2, +1, -3, +0, +0, +0, +0, +0u>  Power;
typedef ScalarUnit <+2, +1, -3, -2, +0, +0, +0, +0u>  Resistance;
typedef ScalarUnit <-2, -1, +3, +2, +0, +0, +0, +0u>  Conductance;
typedef ScalarUnit <-2, -1, +4, +2, +0, +0, +0, +0u>  Capacitance;
typedef ScalarUnit <+2, +1, -2, -2, +0, +0, +0, +0u>  Inductance;
typedef ScalarUnit <+2, +1, -2, -1, +0, +0, +0, +0u>  MagneticFlux;
typedef ScalarUnit <+0, +1, -2, -1, +0, +0, +0, +0u>  MagneticFluxDensity;




// aliases
typedef Angle         DimensionLess;

typedef Length        Displacement;
typedef Length        Distance;
typedef Length        ArcLength;
typedef Length        Radius;

typedef AngularSpeed  AngularRate;

typedef Duration      Time;

typedef Energy        Work;
typedef Energy        KineticEnergy;
typedef Energy        PotentialEnergy;
typedef Energy        BindingEnergy;

typedef Torque        Moment;




/**
 *  User-defined literals
 */


// free literals for multiplers of pi
// FIXME: (Rody Oldenhuis) move to constants
constexpr long double operator"" _pi(long double x)        { return x * ::math::pi; }
constexpr long double operator"" _pi(unsigned long long x) { return x * ::math::pi; }


// Angle
// ------------------------------------
DECLARE_UNIT_TRAITS(Angle)

DEFINE_UNIT_SHORT(Angle, rad, 1,1)

DEFINE_LITERALS(Angle, deg, ::math::pi/180.0)



// free literals for Angles that are multiples of pi radians
constexpr Angle operator"" _pi_rad(long double x)        { return Angle(x * ::math::pi); }
constexpr Angle operator"" _pi_rad(unsigned long long x) { return Angle(x * ::math::pi); }



// Length
// ------------------------------------
DECLARE_UNIT_TRAITS(Length)

DEFINE_UNIT_SHORT(Length, m, 1,1)

DEFINE_UNIT_LONG(Length,  meter, 1,1)
DEFINE_UNIT_LONG(Length, meters, 1,1)


// Mass
// ------------------------------------
DECLARE_UNIT_TRAITS(Mass)

// NOTE: SI base unit is kg (NOT g)...we have to do some workarounds:
DEFINE_UNIT_SHORT(Mass, g, 1.0e-3, 1)

DEFINE_UNIT_LONG(Mass,    gram, 1.0e-3, 1)
DEFINE_UNIT_LONG(Mass,   grams, 1.0e-3, 1)
DEFINE_UNIT_LONG(Mass,  gramme, 1.0e-3, 1)
DEFINE_UNIT_LONG(Mass, grammes, 1.0e-3, 1)



// Duration
// ------------------------------------
DECLARE_UNIT_TRAITS(Duration)

DEFINE_UNIT_SHORT(Duration, s, 1,1)

DEFINE_UNIT_LONG(Duration,  second, 1,1)
DEFINE_UNIT_LONG(Duration, seconds, 1,1)

DEFINE_LITERALS(Duration,   sec, 1)
DEFINE_LITERALS(Duration,   min, 60)
DEFINE_LITERALS(Duration,     h, 3600)
DEFINE_LITERALS(Duration,  hour, 3600)
DEFINE_LITERALS(Duration, hours, 3600)
DEFINE_LITERALS(Duration,   day, 86400)
DEFINE_LITERALS(Duration,  days, 86400)
DEFINE_LITERALS(Duration,     y, 365.25*86400)
DEFINE_LITERALS(Duration,  year, 365.25*86400)
DEFINE_LITERALS(Duration, years, 365.25*86400)


// Current
// ------------------------------------
DECLARE_UNIT_TRAITS(Current)

DEFINE_UNIT_SHORT(Current, A, 1,1)

DEFINE_UNIT_LONG(Current,  ampere, 1,1)
DEFINE_UNIT_LONG(Current, amperes, 1,1)
DEFINE_UNIT_LONG(Current,  Ampere, 1,1)
DEFINE_UNIT_LONG(Current, Amperes, 1,1)



// Temperature
// ------------------------------------
DECLARE_UNIT_TRAITS(Temperature)

DEFINE_UNIT_SHORT(Temperature, K, 1,1)

DEFINE_UNIT_LONG(Temperature,  kelvin, 1,1)
DEFINE_UNIT_LONG(Temperature, kelvins, 1,1)
DEFINE_UNIT_LONG(Temperature,  Kelvin, 1,1)
DEFINE_UNIT_LONG(Temperature, Kelvins, 1,1)



// LuminousIntensity
// ------------------------------------
// TODO: placeholder; future work



// AmountOfSubstance
// ------------------------------------
DECLARE_UNIT_TRAITS(AmountOfSubstance)

DEFINE_UNIT_SHORT(AmountOfSubstance, mol, 1,1)



// Speed
// ------------------------------------
DECLARE_UNIT_TRAITS(Speed)

DEFINE_UNIT_SHORT(Speed, mps, 1,1)

DEFINE_UNIT_LONG(Speed,  meter_per_second, 1,1)
DEFINE_UNIT_LONG(Speed, meters_per_second, 1,1)


// Acceleration
// ------------------------------------
DECLARE_UNIT_TRAITS(Acceleration)

DEFINE_UNIT_SHORT(Acceleration, mps2, 1,1)


// Jerk
// ------------------------------------
DECLARE_UNIT_TRAITS(Jerk)

DEFINE_UNIT_SHORT(Jerk, mps3, 1,1)




// Area
// ------------------------------------
DECLARE_UNIT_TRAITS(Area)

DEFINE_UNIT_SHORT(Area, m2, 1, 2)

DEFINE_UNIT_LONG(Area,  meter_squared, 1, 2)
DEFINE_UNIT_LONG(Area, meters_squared, 1, 2)



// Volume
// ------------------------------------
DECLARE_UNIT_TRAITS(Volume)

DEFINE_UNIT_SHORT(Volume, m3, 1, 3)

DEFINE_UNIT_LONG(Volume,  meter_cubed, 1, 3)
DEFINE_UNIT_LONG(Volume, meters_cubed, 1, 3)



// Density
// ------------------------------------
DECLARE_UNIT_TRAITS(Density)

DEFINE_UNIT_SHORT(Density, gpm3, 1.0e-9, 3)

DEFINE_UNIT_LONG(Density,    gram_per_meter_cubed, 1.0e-9, 3)
DEFINE_UNIT_LONG(Density,   grams_per_meter_cubed, 1.0e-9, 3)
DEFINE_UNIT_LONG(Density,  gramme_per_meter_cubed, 1.0e-9, 3)
DEFINE_UNIT_LONG(Density, grammes_per_meter_cubed, 1.0e-9, 3)




// Energy
// ------------------------------------
DECLARE_UNIT_TRAITS(Energy)

DEFINE_UNIT_SHORT(Energy, J, 1,1)

DEFINE_UNIT_LONG(Energy,  joule, 1,1)
DEFINE_UNIT_LONG(Energy, joules, 1,1)
DEFINE_UNIT_LONG(Energy,  Joule, 1,1)
DEFINE_UNIT_LONG(Energy, Joules, 1,1)


// SpecificEnergy
// ------------------------------------
DECLARE_UNIT_TRAITS(SpecificEnergy)

DEFINE_UNIT_SHORT(SpecificEnergy, Jpkg, 1,1)

DEFINE_UNIT_LONG(SpecificEnergy,  joule_per_kilogram, 1,1)
DEFINE_UNIT_LONG(SpecificEnergy, joules_per_kilogram, 1,1)
DEFINE_UNIT_LONG(SpecificEnergy,  Joule_per_kilogram, 1,1)
DEFINE_UNIT_LONG(SpecificEnergy, Joules_per_kilogram, 1,1)


// Frequency
// ------------------------------------
DECLARE_UNIT_TRAITS(Frequency)

DEFINE_UNIT_SHORT(Frequency, Hz, 1,1)

DEFINE_UNIT_LONG(Frequency, hertz, 1,1)
DEFINE_UNIT_LONG(Frequency, Hertz, 1,1)




// Force
// ------------------------------------
DECLARE_UNIT_TRAITS(Force)

DEFINE_UNIT_SHORT(Force, N, 1,1)

DEFINE_UNIT_LONG(Force,  newton, 1,1)
DEFINE_UNIT_LONG(Force, newtons, 1,1)
DEFINE_UNIT_LONG(Force,  Newton, 1,1)
DEFINE_UNIT_LONG(Force, Newtons, 1,1)


// Momentum
// ------------------------------------
DECLARE_UNIT_TRAITS(Momentum)

DEFINE_UNIT_SHORT(Momentum, Ns, 1,1)

DEFINE_UNIT_LONG(Momentum, newton_seconds, 1,1)
DEFINE_UNIT_LONG(Momentum,  newton_second, 1,1)
DEFINE_UNIT_LONG(Momentum, Newton_seconds, 1,1)
DEFINE_UNIT_LONG(Momentum,  Newton_second, 1,1)



// SpecificMomentum
// ------------------------------------


// Pressure
// ------------------------------------
DECLARE_UNIT_TRAITS(Pressure)

DEFINE_UNIT_SHORT(Pressure, Pa, 1,1)

DEFINE_UNIT_LONG(Pressure, pascals, 1,1)
DEFINE_UNIT_LONG(Pressure,  pascal, 1,1)
DEFINE_UNIT_LONG(Pressure, Pascals, 1,1)
DEFINE_UNIT_LONG(Pressure,  Pascal, 1,1)



// Torque
// ------------------------------------
// NOTE: == energy
DECLARE_UNIT_TRAITS(Torque)

DEFINE_UNIT_SHORT(Torque, Nm, 1,1)

DEFINE_UNIT_LONG(Torque,  newton_meter, 1,1)
DEFINE_UNIT_LONG(Torque, newton_meters, 1,1)
DEFINE_UNIT_LONG(Torque,  Newton_meter, 1,1)
DEFINE_UNIT_LONG(Torque, Newton_meters, 1,1)



// MomentOfInertia
// ------------------------------------


// Charge
// ------------------------------------
DECLARE_UNIT_TRAITS(Charge)

DEFINE_LITERALS(Charge, C, 1)
DEFINE_SHORT_SI_MULTIPLIERS(Charge, C, 1)

DEFINE_LITERALS(Charge, coulombs, 1)
DEFINE_LONG_SI_MULTIPLIERS(Charge, coulombs, 1)

DEFINE_LITERALS(Charge, coulomb, 1)
DEFINE_LONG_SI_MULTIPLIERS(Charge, coulomb, 1)

DEFINE_LITERALS(Charge, Coulombs, 1)
DEFINE_LONG_SI_MULTIPLIERS(Charge, Coulombs, 1)

DEFINE_LITERALS(Charge, Coulomb, 1)
DEFINE_LONG_SI_MULTIPLIERS(Charge, Coulomb, 1)


// Potential
// ------------------------------------
DECLARE_UNIT_TRAITS(Potential)

DEFINE_LITERALS(Potential, V, 1)
DEFINE_SHORT_SI_MULTIPLIERS(Potential, V, 1)

DEFINE_LITERALS(Potential, volts, 1)
DEFINE_LONG_SI_MULTIPLIERS(Potential, volts, 1)

DEFINE_LITERALS(Potential, volt, 1)
DEFINE_LONG_SI_MULTIPLIERS(Potential, volt, 1)

DEFINE_LITERALS(Potential, Volts, 1)
DEFINE_LONG_SI_MULTIPLIERS(Potential, Volts, 1)

DEFINE_LITERALS(Potential, Volt, 1)
DEFINE_LONG_SI_MULTIPLIERS(Potential, Volt, 1)


// Power
// ------------------------------------
DECLARE_UNIT_TRAITS(Power)

DEFINE_UNIT_SHORT(Power, W, 1,1)

DEFINE_UNIT_LONG(Power,  watt, 1,1)
DEFINE_UNIT_LONG(Power, watts, 1,1)
DEFINE_UNIT_LONG(Power,  Watt, 1,1)
DEFINE_UNIT_LONG(Power, Watts, 1,1)



// Resistance
// ------------------------------------
DECLARE_UNIT_TRAITS(Resistance)

DEFINE_LITERALS(Resistance, ohms, 1)
DEFINE_LONG_SI_MULTIPLIERS(Resistance, ohms, 1)

DEFINE_LITERALS(Resistance, ohm, 1)
DEFINE_LONG_SI_MULTIPLIERS(Resistance, ohm, 1)

DEFINE_LITERALS(Resistance, Ohms, 1)
DEFINE_LONG_SI_MULTIPLIERS(Resistance, Ohms, 1)

DEFINE_LITERALS(Resistance, Ohm, 1)
DEFINE_LONG_SI_MULTIPLIERS(Resistance, Ohm, 1)


// Conductance
// ------------------------------------
DECLARE_UNIT_TRAITS(Conductance)

//DEFINE_LITERALS(Conductance, S, 1)
//DEFINE_SHORT_SI_MULTIPLIERS(Conductance, S, 1)

DEFINE_LITERALS(Conductance, siemens, 1)
DEFINE_LONG_SI_MULTIPLIERS(Conductance, siemens, 1)

DEFINE_LITERALS(Conductance, Siemens, 1)
DEFINE_LONG_SI_MULTIPLIERS(Conductance, Siemens, 1)


// Capacitance
// ------------------------------------
DECLARE_UNIT_TRAITS(Capacitance)

DEFINE_LITERALS(Capacitance, F, 1)
DEFINE_SHORT_SI_MULTIPLIERS(Capacitance, F, 1)

DEFINE_LITERALS(Capacitance, farads, 1)
DEFINE_LONG_SI_MULTIPLIERS(Capacitance, farads, 1)

DEFINE_LITERALS(Capacitance, farad, 1)
DEFINE_LONG_SI_MULTIPLIERS(Capacitance, farad, 1)

DEFINE_LITERALS(Capacitance, Farads, 1)
DEFINE_LONG_SI_MULTIPLIERS(Capacitance, Farads, 1)

DEFINE_LITERALS(Capacitance, Farad, 1)
DEFINE_LONG_SI_MULTIPLIERS(Capacitance, Farad, 1)



// Inductance
// ------------------------------------
DECLARE_UNIT_TRAITS(Inductance)

DEFINE_LITERALS(Inductance, H, 1)
DEFINE_SHORT_SI_MULTIPLIERS(Inductance, H, 1)

DEFINE_LITERALS(Inductance, henries, 1)
DEFINE_LONG_SI_MULTIPLIERS(Inductance, henries, 1)

DEFINE_LITERALS(Inductance, henry, 1)
DEFINE_LONG_SI_MULTIPLIERS(Inductance, henry, 1)

DEFINE_LITERALS(Inductance, Henries, 1)
DEFINE_LONG_SI_MULTIPLIERS(Inductance, Henries, 1)

DEFINE_LITERALS(Inductance, Henry, 1)
DEFINE_LONG_SI_MULTIPLIERS(Inductance, Henry, 1)



// MagneticFlux
// ------------------------------------


// MagneticFluxDensity
// ------------------------------------
DECLARE_UNIT_TRAITS(MagneticFluxDensity)

DEFINE_LITERALS(MagneticFluxDensity, T, 1)
DEFINE_SHORT_SI_MULTIPLIERS(MagneticFluxDensity, T, 1)

DEFINE_LITERALS(MagneticFluxDensity, teslas, 1)
DEFINE_LONG_SI_MULTIPLIERS(MagneticFluxDensity, teslas, 1)

DEFINE_LITERALS(MagneticFluxDensity, tesla, 1)
DEFINE_LONG_SI_MULTIPLIERS(MagneticFluxDensity, tesla, 1)

DEFINE_LITERALS(MagneticFluxDensity, Teslas, 1)
DEFINE_LONG_SI_MULTIPLIERS(MagneticFluxDensity, Teslas, 1)

DEFINE_LITERALS(MagneticFluxDensity, Tesla, 1)
DEFINE_LONG_SI_MULTIPLIERS(MagneticFluxDensity, Tesla, 1)




// We're done with the macros:
#undef DECLARE_UNIT_TRAITS
#undef DEFINE_LITERALS
#undef DEFINE_SI_MULTIPLIER
#undef DEFINE_SHORT_SI_MULTIPLIERS
#undef DEFINE_LONG_SI_MULTIPLIERS
#undef DEFINE_UNIT_SHORT
#undef DEFINE_UNIT_LONG

#endif // include guard


