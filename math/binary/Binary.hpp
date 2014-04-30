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


#ifndef __BINARY_HPP
#define __BINARY_HPP


#include <iostream>
#include <assert.h>
#include <cmath>
#include <string>
#include <sstream>


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
    constexpr unit operator"" _##prefix##suffix(long double V)        { return operator"" _##suffix( pow(multiplier,power) * V  ); } \
    constexpr unit operator"" _##prefix##suffix(unsigned long long V) { return operator"" _##suffix( pow(multiplier,power) * static_cast<long double>(V) ); }

#define DEFINE_SHORT_SI_MULTIPLIERS(unit, suffix, power) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     Y, 1.0e+24L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     Z, 1.0e+21L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     E, 1.0e+18L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     P, 1.0e+15L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     T, 1.0e+12L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     G, 1.0e+09L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     M, 1.0e+06L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     k, 1.0e+03L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,     h, 1.0e+02L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    da, 1.0e+01L)

#define DEFINE_LONG_SI_MULTIPLIERS(unit, suffix, power) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, yotta, 1.0e+24L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, zetta, 1.0e+21L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, exa  , 1.0e+18L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, peta , 1.0e+15L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, tera , 1.0e+12L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, giga , 1.0e+09L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, mega , 1.0e+06L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, kilo , 1.0e+03L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, hecto, 1.0e+02L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power, deca , 1.0e+01L)

#define DEFINE_SHORT_BINARY_MULTIPLIERS(unit, suffix, power) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    Yi, 1.0e+24L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    Zi, 1.0e+21L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    Ei, 1.0e+18L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    Pi, 1.0e+15L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    Ti, 1.0e+12L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    Gi, 1.0e+09L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    Mi, 1.0e+06L) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,    ki, 1.0e+03L)

#define DEFINE_LONG_BINARY_MULTIPLIERS(unit, suffix, power) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  yobi, pow(1024L,8)) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  zebi, pow(1024L,7)) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  exbi, pow(1024L,6)) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  pebi, pow(1024L,5)) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  tebi, pow(1024L,4)) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  gibi, pow(1024L,3)) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  mebi, pow(1024L,2)) \
    DEFINE_SI_MULTIPLIER (unit,suffix,power,  kibi, pow(1024L,1))


#define DEFINE_UNIT_SHORT(unit, suffix, multiplier, power) \
    DEFINE_LITERALS(unit, suffix, multiplier) \
    DEFINE_SHORT_SI_MULTIPLIERS(unit, suffix, power)

#define DEFINE_UNIT_LONG(unit, suffix, multiplier, power) \
    DEFINE_LITERALS(unit, suffix, multiplier) \
    DEFINE_LONG_SI_MULTIPLIERS(unit, suffix, power)


// Small traits class to get the SI unit name in (string template argmuents are not supported)
// TODO: Think more about this one...
template<class S>
struct ScalarUnitTraits {
    static const std::string SI_unit;
};



class Bits
{
  private:

    // the actual number of bits
    unsigned long long int number;

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

    // toString (for lexical casts)
    constexpr std::string toString() const {
        std::stringstream output;
        output << value << " " << " bits";
        return output.str();
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



// We're done with the macros:
#undef DECLARE_UNIT_TRAITS
#undef DEFINE_LITERALS
#undef DEFINE_SI_MULTIPLIER
#undef DEFINE_SHORT_SI_MULTIPLIERS
#undef DEFINE_LONG_SI_MULTIPLIERS
#undef DEFINE_UNIT_SHORT
#undef DEFINE_UNIT_LONG

#endif // include guard


