/**
 *  @file          constants.hpp
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

/*
 * Macros other than the include guards and includes are discouraged
 */
#ifndef _GEMS_CONSTANTS_HPP_
#define _GEMS_CONSTANTS_HPP_

/*
 *-------------------------------------------------------------------------
 * Includes
 * All includes go here
 *-------------------------------------------------------------------------
 */
#include <limits>
#include <map>
#include <vector>
#include <cmath>
#include <string>

/// SI multiplier table
namespace si
{
    const double

    // names          prefixes
    yotta = 1e24,     Y  = yotta,
    zetta = 1e21,     Z  = zetta,
    exa   = 1e18,     E  = exa  ,
    peta  = 1e15,     P  = peta ,
    tera  = 1e12,     T  = tera ,
    giga  = 1e9,      G  = giga ,
    mega  = 1e6,      M  = mega ,
    kilo  = 1e3,      k  = kilo ,
    hecto = 1e2,      h  = hecto,
    deca  = 1e1,      da = deca ,

    deci  = 1e-1,     d  = deci ,
    centi = 1e-2,     c  = centi,
    milli = 1e-3,     m  = milli,
    micro = 1e-6,     mu = micro,
    nano  = 1e-9,     n  = nano ,
    pico  = 1e-12,    p  = pico ,
    femto = 1e-15,    f  = femto,
    atto  = 1e-18,    a  = atto ,
    zepto = 1e-21,    z  = zepto,
    yocto = 1e-24,    y  = yocto;
} /* namespace si */


/// Binary prefixes
namespace binary
{
    const double

    // names          prefixes
    kibi = 1.*1024,   ki = kibi,
    mebi = ki*1024,   Mi = mebi,
    gibi = Mi*1024,   Gi = gibi,
    tebi = Gi*1024,   Ti = tebi,
    pebi = Ti*1024,   Pi = pebi,
    exbi = Pi*1024,   Ei = exbi,
    zebi = Ei*1024,   Zi = zebi,
    yobi = Zi*1024,   Yi = yobi;

} // namespace binary


/// Physics constants
namespace physics
{
    /**
     * Fundamental constants of nature (and a few not-so-fundamental ones)
     *
     *
     */
    const double

    c     = 299792458,       // speed of light [m/s]
    k     = 1.3806488e-23,   // Boltzmann constant [J K^-1]
    sigma = 5.670373e-8,     // Stefan-Boltzmann constant [W m^-2 K^-4]
    NA    = 6.02214129e+23,  // Avogadro's constant [1/mol]
    h     = 6.62606957e-34,  // Planck's constant [Js]
    hbar  = 1.054571726e-34, // Reduced Plankc's constant, or Dirac's constant [Js]
    re    = 2.817940285e-15, // Classical electron radius [m]
    q     = 1.60217646e-19,  // Proton charge [C]

    au    = 149597870700;    // Average Sun-Earth distance (roughly speaking) [m]

} /* namespace physics */


/// Mathematical constants
namespace math
{
    // TODO: use __float128? -> GCC-only for the moment...Leave the digits in though
    const long double

    // computed with GNU bc arbitrary precision calculator
    pi         = 3.14159265358979323846264338327950288419716939937510L, // π
    pio2       = 1.57079632679489661923132169163975144209858469968755L, // π/2
    pio4       = 0.78539816339744830961566084581987572104929234984377L, // π/4
    pio6       = 0.52359877559829887307710723054658381403286156656251L, // π/6
    threepio2  = 4.71238898038468985769396507491925432629575409906265L, // 3π/2
    twopi      = 6.28318530717958647692528676655900576839433879875020L, // 2π
    tau        = twopi, // AKA, the TRUE circle constant!!
    twotau     = 12.5663706143591729538505735331180115367886775975004L,
    oneopi     = 0.31830988618379067153776752674502872406891929148091L, // 1/π
    twoopi     = 0.63661977236758134307553505349005744813783858296182L, // 2/π
    oneosqrtpi = 0.56418958354775628694807945156077258584405062932900L, // 1/sqrt(π)
    twoosqrtpi = 1.12837916709551257389615890312154517168810125865800L, // 2/sqrt(π)

    sqrt2      = 1.41421356237309504880168872420969807856967187537694L, // sqrt(2)
    sqrt3      = 1.73205080756887729352744634150587236694280525381038L, // sqrt(3)
    sqrt1o2    = 0.70710678118654752440084436210484903928483593768847L, // sqrt(½)
    sqrt1o3    = 0.57735026918962576450914878050195745564760175127012L, // sqrt(⅓)

    epsilon    = std::numeric_limits<double>::epsilon(),
    dx         = epsilon,  // alias
    eps        = epsilon,  // alias
    INF        = std::numeric_limits<double>::infinity(),
    inf        = INF,      // alias
    NaN        = std::numeric_limits<double>::quiet_NaN(),
    nan        = NaN,      // alias

    epsilonf   = std::numeric_limits<float>::epsilon(),
    dxf        = epsilonf, // alias
    epsf       = epsilonf, // alias
    INFf       = std::numeric_limits<float>::infinity(),
    inff       = INFf,     // alias
    NaNf       = std::numeric_limits<float>::quiet_NaN(),
    nanf       = NaNf,     // alias

    // Angular conversions
    DEG2RAD    = pi/180,
    RAD2DEG    = 180/pi,
    // 1° = 60"
    ARCMIN2RAD = DEG2RAD / 60,
    RAD2ARCMIN = RAD2DEG * 60,
    // 1" = 60'
    ARCSEC2RAD = ARCMIN2RAD / 60,
    RAD2ARCSEC = RAD2ARCMIN * 60;

    // handy dandy utilities
    // ---------------------

    // Check for NaN, inf
    // FIXME: (Rody Oldenhuis) is this needed...?
    template<typename T>
    bool isFinite(const T &x){
        return std::isfinite(x);
    }

    // Check for almost-equality of floats (operator==() is unreliable)
    // see
    //     http://www.boost.org/doc/libs/1_36_0/libs/test/doc/html/utf/testing-tools/floating_point_comparison.html
    //     http://stackoverflow.com/questions/17333/most-effective-way-for-float-and-double-comparison
    bool isAlmostEqual(double, double, double tol = 1.0 );
    bool isAlmostEqual(float , float , float  tol = 1.0f);

    // signum: -1 for anything <0, +1 for anything >0, 0 for anything ==0.
    template <typename T>
    char signum(const T &in) {
        return (in > static_cast<T>(0)) - (in < static_cast<T>(0));
    }

    // round towards zero
    template<typename T>
    T fix(const T &x) {
        if (static_cast<T>(0.0) > x)
            return std::ceil(x);
        else
            return std::floor(x);
    }

    // modulus after division
    // NOTE: NOT the same as % or fmod()
    template <typename T>
    T mod(const T &x, const T &y) {
        return (static_cast<T>(0)==y) ? x : ( x-y*std::floor(x/y) );
    }

    // remainder after division
    template <typename T>
    T rem(const T &x, const T &y) {
        return (static_cast<T>(0)==y) ? x : ( x-y*fix(x/y) );
    }

    // all elements are non-zero
    template<typename T>
    bool all(const T &V) {
        for (typename T::const_iterator i=V.begin(); i<V.end(); ++i)
            if (!(*i))
                return false;
        return true;
    }

    // any element in iterable container is non-zero
    template<typename T>
    bool any(const T &V) {
        for (typename T::const_iterator i=V.begin(); i<V.end(); ++i)
            if (!(*i))
                return true;
        return false;
    }

    /**
     * Wrap values in a given interval
     */

    /// Circular wrapping: the amount by which x exceeds one of the limits is added recursively to the other limit
    template<typename T>
    T wrap_circular(const T &x, const T &min_x, const T &max_x) {
        return min_x + math::mod(x-min_x, max_x-min_x);
    }

    /// Saturation wrapping: if x exceeds any of the limits, the limit itself is returned.
    template<typename T>
    T wrap_saturate(const T &x, const T &min_x, const T &max_x) {
        return std::min( max_x, std::max(min_x,x) );
    }
    template<typename T> // alias for saturation wrap
    T ensure_interval(const T &x, const T &min_x, const T &max_x) {
        return wrap_saturate(x, min_x, max_x);
    }

} /* namespace math */

#endif


