/*
 * Rody Oldenhuis
 * cosine Science & Computing BV
 * roldenhuis@cosine.nl
 *
 * PhysUnits.hh
 * Created: 31.08.2012 16:20:53 CEST
 */

#ifndef _PHYSUNITS_HH
#define _PHYSUNITS_HH

#include <thread>
#include <boost/date_time/posix_time/posix_time.hpp>

// global helper variables
namespace
{
    union endian_t {
        long l;
        char c[sizeof(long)];
    } __endian__t;
}


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
}


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

}


/// Physics constants
namespace physics
{
    const double

    c     = 299792458,       // speed of light [m/s]
    au    = 1.49597870e11,   // unit distance [m]
    KSun  = 1370,            // Solar constant [W/m^2]

    // FIXME: (Rody Oldenhuis) DEPRACATED
    //    Using Solar mass (or masses in general) introduces terrible numerical
    //    inaccuracies everywhere due to current uncertainty in the value of
    //    Newton's G. Use standard gravitational parameters (=GM); these have
    //    up to 6 orders of magnitude better accuracy.
    G     = 6.67384e-11,     // Gravitational constant [m^3/kg/s^2] (CODATA 2010 recommended value)
    MSun  = 1.98855e30,      // Solar mass [kg]
    // -----------------------------------------

    k     = 1.3806488e-23,   // Boltzmann constant [J K^-1]
    sigma = 5.670373e-8,     // Stefan-Boltzmann constant [W m^-2 K^-4]

    // FIXME: (Rody Oldenhuis) Not really a "constant", more a unit...
    // Use class Length for these things.
    A     = 1.0e-10,         // Angstrom
    // -----------------------------------------

    NA    = 6.02214129e+23,  // Avogadro's constant [1/mol]
    h     = 6.62606957e-34,  // Planck's constant [Js]
    hbar  = 1.054571726e-34, // Reduced Plankc's constant, or Dirac's constant [Js]
    re    = 2.817940285e-15, // Classical electron radius [m]
    q     = 1.60217646e-19,  // Proton charge [C]

    // FIXME: (Rody Oldenhuis) DEPRECATED
    //   Only included because full conversion takes too much time.
    //   Find files in need of conversion with "grep 'physics::day' *.cc *.hh"
    //   and convert everything to class Time.
    s        = 1,
    minute   = 60*s,
    hour     = 60*minute,
    day      = 24*hour,
    jyear    = 365.25*day,
    jcentury = 100*jyear;
    // -----------------------------------------

}


/// Mathematical constants
namespace math
{
    const double

    pi         = 3.14159265358979323846264338327950288419716939937510,
    pio2       = 1.57079632679489661923132169163975144209858469968755,
    pio4       = 0.78539816339744830961566084581987572104929234984377,
    pio6       = 0.52359877559829887307710723054658381403286156656251,
    threepio2  = 4.71238898038468985769396507491925432629575409906265,
    twopi      = 6.28318530717958647692528676655900576839433879875020,
    tau        = twopi, // the TRUE circle constant!!
    twotau     = 12.5663706143591729538505735331180115367886775975004,
    oneopi     = 0.31830988618379067153776752674502872406891929148091, // 1/pi
    twoopi     = 0.63661977236758134307553505349005744813783858296182, // 2/pi
    oneosqrtpi = 0.56418958354775628694807945156077258584405062932900, // 1/sqrt(pi)
    twoosqrtpi = 1.12837916709551257389615890312154517168810125865800, // 2/sqrt(pi)

    sqrt2      = 1.41421356237309504880168872420969807856967187537694,
    sqrt3      = 1.73205080756887729352744634150587236694280525381038,
    sqrt1o2    = 0.70710678118654752440084436210484903928483593768847, // sqrt(1/2)
    sqrt1o3    = 0.57735026918962576450914878050195745564760175127012, // sqrt(1/3)

    epsilon    = std::numeric_limits<double>::epsilon(),
    dx         = epsilon,
    eps        = epsilon,
    INF        = std::numeric_limits<double>::infinity(),
    inf        = INF,
    NaN        = std::numeric_limits<double>::quiet_NaN(),
    nan        = NaN;


    // some things needfloat versions of the above
    const float

    epsilonf = std::numeric_limits<float>::epsilon(),
    dxf      = epsilonf,
    epsf     = epsilonf;


    // handy utilities
    // ----------------

    // signum: -1 for anything <0, +1 for anything >0, 0 for anything ==0.
    template <typename T>
    int signum(T in) {
        return (in>0)-(in<0);
    }

    // negative-value-correct modulus
    // NOTE: NOT the same as % or fmod()
    template <typename T>
    T mod(T x, T y) {
        return (0==y) ? x : x-y*std::floor(x/y);
    }

    // FIXME: (Rody Oldenhuis) DEPRECATED
    //   Only included because full conversion takes too much time.
    //   Find files in need of conversion with "grep 'math::deg' *.cc *.hh"
    //   and convert everything to class Angle.
    const double
        rad = 1.0,
        deg = pi/180.0;
    // -----------------------------------------
}


/// Constants specific to the system HIPSIM's being run on
namespace machine
{
    const unsigned int

    // get number of physical cores present in the system
    // platform-independent, > C++11
    numCores = std::thread::hardware_concurrency();


    // This machine's endianness
    const bool

    machine_is_msb = (__endian__t.l = 1, __endian__t.c[sizeof(long)-1] == 1),
    machine_is_lsb = (__endian__t.l = 1, __endian__t.c[sizeof(long)-1] != 1);


    // handy dandy utilities

    template<typename T>
    T change_endianness(T x)
    {
        char *b = (char *)(&x);
        switch (sizeof(T))
        {
            case 1:
                break;

            case 2:
                std::swap(b[0], b[1]);
                break;

            case 4:
                std::swap(b[0], b[3]);
                std::swap(b[1], b[2]);
                break;

            case 8:
                std::swap(b[0], b[7]);
                std::swap(b[1], b[6]);
                std::swap(b[2], b[5]);
                std::swap(b[3], b[4]);
                break;

            default:
                static const bool cannot_swap = true;
                assert(cannot_swap);
        }
        return x;
    }

    template<typename T>
    void write_msb(std::ostream & os, T x) {
        if (machine_is_lsb)
            x = change_endianness(x);
        os.write(reinterpret_cast<char *>(&x), sizeof(T));
    }

    template<typename T>
    void write_lsb(std::ostream & os, T x) {
        if (machine_is_msb)
            x = change_endianness(x);
        os.write(reinterpret_cast<char *>(&x), sizeof(T));
    }

    template<typename T>
    void read_msb(std::istream & is, T & x) {
        is.read(reinterpret_cast<char *>(&x), sizeof(T));
        if (machine_is_lsb)
        x = change_endianness(x);
    }

    template<typename T>
    void read_lsb(std::istream & is, T & x) {
        is.read(reinterpret_cast<char *>(&x), sizeof(T));
        if (machine_is_msb)
        x = change_endianness(x);
    }

}


// base classes for all units
// --------------------------

// fwd declare everything
namespace _ignore{
    class Unit;
    template<typename T>class UnitOperators;
}


// fwd declare everything
class Angle;      class Energy;
class Length;     class Force;
class Area;       class Temperature;
class Volume;     class Mass;
class Time;       class Speed;
class Density;    class Pressure;


// Math operators
template <typename T> T  abs(const _ignore::UnitOperators<T>& t);
template <typename T> T fabs(const _ignore::UnitOperators<T>& t);


// definition of unit base classes
namespace _ignore
{
    // Base class for all units
    class Unit
    {
    protected:
        double value;
        Unit(double value) : value(value) {}

    public:
        virtual ~Unit(){}

        Unit() : value(0.0) {}
        Unit(const Unit& obj) : value(obj.value) {}

        // enforce string representation
        virtual std::string toString() const = 0;

    };

    // RO: use CRTP pattern to have the definition of the operators all in
    // one place, while guaranteeing that only same-type operators work. See
    // http://stackoverflow.com/questions/12742728/
    template <typename T>
    class UnitOperators
        : public Unit
    {
        friend T  abs<T>(const UnitOperators<T>& t);
        friend T fabs<T>(const UnitOperators<T>& t);

    protected:
        UnitOperators(double value) : Unit(value){}

    public:
        virtual ~UnitOperators(){}

        UnitOperators()             : Unit()    {}
        UnitOperators(const T& obj) : Unit(obj) {}

        // RO : use static_cast, and make this class a friend class of all
        // classes subclassing it, as per
        // http://stackoverflow.com/questions/12910883/
        T&     operator= (const T& rhs)       { if (this!=&rhs) value = rhs.value; return static_cast<T&>(*this); }

        bool   operator> (const T& rhs) const { return value >  rhs.value; }
        bool   operator< (const T& rhs) const { return value <  rhs.value; }
        bool   operator>=(const T& rhs) const { return value >= rhs.value; }
        bool   operator<=(const T& rhs) const { return value <= rhs.value; }
        bool   operator==(const T& rhs) const { return value == rhs.value; }
        bool   operator!=(const T& rhs) const { return value != rhs.value; }

        T&     operator+=(const T& rhs)       { value += rhs.value; return static_cast<T&>(*this); }
        T&     operator-=(const T& rhs)       { value -= rhs.value; return static_cast<T&>(*this); }
        T&     operator*=(double F)           { value *= F; return static_cast<T&>(*this); }
        T&     operator/=(double F)           { value /= F; return static_cast<T&>(*this); }

        T      operator+ (const T& rhs) const { return T(*this) += rhs; }
        T      operator- (const T& rhs) const { return T(*this) -= rhs; }
        T      operator- ()             const { return T(-value);}
        T      operator+ ()             const { return T(+value);}

        T      operator* (double F)     const { return T(*this) *= F; }
        T      operator/ (double F)     const { return T(*this) /= F; }
        double operator* (const T& rhs) const { return value*rhs.value; }
        double operator/ (const T& rhs) const { return value/rhs.value; }

        // some nice constant instances
        static const T zero        ;
        static const T dx          ;
        static const T infinitesmal;
        static const T small       ;
        static const T NaN         ;
        static const T nan         ;
        static const T inf         ;
        static const T infinite    ;
        static const T infinity    ;

    };

}



/**
 * @brief Object with units [angle]. Aimed at avoiding confusion between degrees, radians,...
 */

/**
 * Angle
 * avoid confusion between degrees, radians, gradians
 *
 * @see Length @see Time @see Mass @see Temperature @see Energy
 * @see Area @see Force @see Speed
 *
 */


class Angle final
    : public _ignore::UnitOperators<Angle>
{
    typedef _ignore::UnitOperators<Angle> base;
    friend class _ignore::UnitOperators<Angle>;

private:

    // value is always in radians
    explicit Angle(double value) : base(value) {};

public:
   ~Angle(){}

    // default value: 0
    Angle() : base(){};
    // copy constructor
    Angle(const Angle& a2) : base(a2){};
    Angle(const base&  a2) : base(a2){};

    /// Angle is instantiated through named constructors
    static Angle
        radians    (double rad),
        degrees    (double deg),
        arcminutes (double moa),
        arcseconds (double soa),
        GRADIANS   (double grad); // NOTE: capitalized to avoid
                                  // confusion w/ radians
    /// Angle's different return types
    double
        radians    () const,
        degrees    () const,
        arcminutes () const,
        arcseconds () const,
        GRADIANS   () const; // NOTE: capitalized to avoid
                             // confusion w/ radians
    /// alternatively, one can use various pre-defined angles
    static const Angle
        tau ,   _360,
        pi  ,   _180,
        pio2,   _90 ,
        pio4,   _45 ,
        pio6,   _30 ;

    /// wrap value into [-pi +pi) / [-180 +180)
    Angle& wrap_posneg();
    friend Angle wrap_posneg(const Angle& a); // free function

    /// wrap value into [0 2pi) / [0 360)
    Angle& wrap_positive();
    friend Angle wrap_positive(const Angle& a); // free function

    // NOTE: angles have no units -- it's OK for these
    // operators to return Angles
    using _ignore::UnitOperators<Angle>::operator/;
    using _ignore::UnitOperators<Angle>::operator*;
    Angle
        operator/ (const Angle& a2) const,
        operator* (const Angle& a2) const;

    /// human-readable string representation
    std::string toString() const;
};


/**
 * @brief Object with units of Length. Makes metric/imperial less prone to error.
 */

/**
 * Lengths
 * make metric <-> imperial a bit less prone to error
 *
 * @see Angle @see Time @see Mass @see Temperature @see Energy
 * @see Area @see Force @see Speed
 */
class Length final
    : public _ignore::UnitOperators<Length>
{
    typedef _ignore::UnitOperators<Length> base;
    friend class _ignore::UnitOperators<Length>;

private:
    // value is always in [m]
    explicit Length(double value) : base(value) {}

public:
    ~Length(){}

    // default
    Length() : base() {}
    // Copy constructor
    Length (const Length& L2) : base(L2) {}
    Length (const base&   L2) : base(L2) {}

    /// Length is instantiated through named constructors
    static Length
        meters           (double m)  ,     m (double   m),
        kilometers       (double km) ,    km (double  km),
        nanometers       (double nm) ,    nm (double  nm),
        microns          (double mu) ,    um (double  um),
        micrometers      (double mu) ,
        inches           (double in) ,    in (double  in),
        feet             (double ft) ,    ft (double  ft),
        foot             (double ft) ,
        yards            (double yrd),
        miles            (double mi) ,    mi (double  mi),
        nauticalMiles    (double nmi),   nmi (double nmi),
        astronomicalUnits(double au ),    AU (double  au),
        angstroms        (double a)  ,     A (double   a);


    /// Length can return various units of length
    double
        meters            () const,   m() const,
        kilometers        () const,  km() const,
        nanometers        () const,  nm() const,
        microns           () const,  um() const,
        micrometers       () const,
        inches            () const,  in() const,
        feet              () const,  ft() const,
        foot              () const,
        yards             () const,
        miles             () const,  mi() const,
        nauticalMiles     () const, nmi() const,
        astronomicalUnits () const,  AU() const,
        angstroms         () const,   A() const;

    /// alternatively, one can use various constant lengths
    static const Length

        // SI
        one_km, kilometer,          // 1 kilometer
        one_m,  meter,              // 1 meter
        one_mm, millimeter,         // 1 millimeter
        one_um, micrometer, micron, // 1 micron
        one_nm, nanometer,          // 1 nanometer

        // NON-SI
        one_AU, astronomicalUnit,
        inch, one_ft, yard, mile, nauticalMile, angstrom;

    ///cast operators
    operator Energy() const; // wavelength to energy

    /// Human-readable string representation
    std::string toString() const;
};


/**
 * @brief
 */

/**
 * Times
 * make hours, minute, seconds, etc. readable
 *
 * @see Angle @see Length @see Mass @see Temperature @see Energy
 * @see Area @see Force @see Speed
 */
class Time final
    : public _ignore::UnitOperators<Time>
{
    typedef _ignore::UnitOperators<Time> base;
    friend class _ignore::UnitOperators<Time>;

private:
    /// value is always given in seconds
    explicit Time(double value) : _ignore::UnitOperators<Time>(value){}

public:
   ~Time(){}

    // default value
    Time() : base(){}
    // Copy constructor
    Time (const Time& t2) : base(t2){}
    Time (const base& t2) : base(t2){}

    /// Time is instantiated through named constructors
    static Time
        microseconds(double ms),
        milliseconds(double ms),
        seconds     (double s),
        minutes     (double m),
        hours       (double h),
        days        (double d),
        years       (double y),
        centuries   (double c);

    /// Length can return various units of time
    double
        microseconds() const,
        milliseconds() const,
        seconds     () const,
        minutes     () const,
        hours       () const,
        days        () const,
        years       () const,
        centuries   () const;

    /// alternatively, one can use various pre-defined lengths
    static const Time
        microsecond,
        millisecond,
        second,
        minute,
        hour,
        day,
        year,
        century;

    // current unix time (seconds+microseconds)
    static Time now ();
    // boost version
    static boost::posix_time::ptime now(void *);

    /// Human-readable string representation
    std::string toString() const;
};



/**
 * @brief
 */

/**
 * Mass
 * make metric <-> imperial a bit less prone to error
 *
 * @see Angle @see Length @see Temperature @see Energy
 * @see Area @see Force @see Speed @see Time
 */
class Mass final
    : public _ignore::UnitOperators<Mass>
{
    typedef _ignore::UnitOperators<Mass> base;
    friend class _ignore::UnitOperators<Mass>;

private:

    // value is always in kg
    explicit Mass(double value) : base(value){}

public:
    ~Mass(){}

    // default value
    Mass() : base(){}
    // Copy constructor
    Mass (const Mass& m2) : base(m2){}
    Mass (const base& m2) : base(m2){}

    // nice named constructors
    static Mass
        kilogrammes (double kg ),  kg  (double kg ),  kilograms (double kg ),
        pounds      (double lbs),  lbs (double lbs),
        tonnes      (double T  ),
        longTonnes  (double LT ),
        shortTonnes (double sT ),
        stones      (double st ),
        ounces      (double st );

    // various return values
    double
        kilogrammes() const,    kg  () const,  kilograms () const,
        pounds     () const,    lbs () const,
        tonnes     () const,
        longTonnes () const,
        shortTonnes() const,
        stones     () const,
        ounces     () const;

    // constants
    static const Mass
        kilo,
        pound,
        tonne,
        longTonne,
        shortTonne,
        stone,
        ounce;

    // Human-readable string representation
    std::string toString() const;
};



/**
 * @brief
 */

/**
 * Temperature
 *
 * @see Angle @see Length @see Mass @see Energy
 * @see Area @see Force @see Speed @see Time
 */
class Temperature final
    : public _ignore::UnitOperators<Temperature>
{
    typedef _ignore::UnitOperators<Temperature> base;
    friend class _ignore::UnitOperators<Temperature>;

private:

    // value is always in K
    // NOTE: ensure temperature is POSITIVE
    explicit Temperature(double value) : base(value<0.0?0.0:value) {}

public:
   ~Temperature(){}

    // default value
    Temperature() : base(){}
    // copy constructor
    Temperature (const Temperature& T2) : base(T2){}
    Temperature (const base&        T2) : base(T2){}

    // nice named constructors
    static Temperature
        kelvin     (double K),
        celsius    (double C),
        fahrenheit (double F);

    // various return values
    double
        kelvin     () const,
        celsius    () const,
        fahrenheit () const;

    // constants
    static const Temperature
        K, // ZERO kelvin
        C, // ZERO celcius
        F; // ZERO fahrenheit

    // NOTE: negative temperatures do not exist; exclude unary minus operator
    Temperature& operator-();

    // cast operators
    operator Energy() const;

    // Human-readable string representation
    std::string toString() const;
};



/// Combined unit classes / operations


/// Area = Length * Length
class Area final
    : public _ignore::UnitOperators<Area>
{
    typedef _ignore::UnitOperators<Area> base;
    friend class _ignore::UnitOperators<Area>;

private:
    // value is always in m²
    explicit Area (double value) : base(value){}

public:
   ~Area(){}

    // default value
    Area() : base(){}
    // copy constructor
    Area (const Area& E2) : base(E2){}
    Area (const base& E2) : base(E2){}

    // nice named constructors
    static Area
        squareMeters     (double m2  ),    m2 (double m2),
        squareKilometers (double m2  ),   km2 (double km2),
        squareMillimeters(double m2  ),   mm2 (double mm2),
        ares             (double are ),
        acres            (double acre),
        hectares         (double ha  );

    // various return values
    double
        squareMeters     () const,    m2 () const,
        squareKilometers () const,   km2 () const,
        squareMillimeters() const,   mm2 () const,
        ares             () const,
        acres            () const,
        hectares         () const;

    // constants
    const static Area
        squareMeter,
        squareKilometer,
        squareMillimeter,
        are,
        acre,
        hectare;

    // Human-readable string representation
    std::string toString() const;

};


/// Volume = Length * Length * Length
class Volume final
    : public _ignore::UnitOperators<Volume>
{
    typedef _ignore::UnitOperators<Volume> base;
    friend class _ignore::UnitOperators<Volume>;

private:
    // value is always in m³
    explicit Volume (double value) : base(value){}

public:
   ~Volume(){}

    // default value
    Volume() : base(){}
    // copy constructors
    Volume (const Volume& V2) : base(V2){}
    Volume (const base&   V2) : base(V2){}

    // nice named constructors
    static Volume
        cubicKilometers (double km3  ),
        cubicMeters     (double  m3  ),
        cubicCentimeters(double cm3  ),
        cubicMillimeters(double mm3  ),
        cubicFeet       (double ft3  ),
        cubicInches     (double in3  ),
        litres          (double l    ),
        gallons         (double g    ), // imperial gallon
        USgallons       (double g    ); // US gallon

    // various return values
    double
        cubicKilometers  () const,
        cubicMeters      () const,
        cubicCentimeters () const,
        cubicMillimeters () const,
        cubicFeet        () const,
        cubicInches      () const,
        litres           () const,
        gallons          () const,
        USgallons        () const;

    // constants
    const static Volume
        cubicKilometer,    km3,
        cubicMeter,         m3,
        cubicCentimeter,   cm3,
        cubicMillimeter,   mm3,
        cubicFoot,         ft3,
        cubicInch,         in3,
        litre,               l,
        gallon,              g,
        USgallon,           Ug;

    // Human-readable string representation
    std::string toString() const;

};

/// Energy = Mass * Length / Time
class Energy final
    : public _ignore::UnitOperators<Energy>
{
    typedef _ignore::UnitOperators<Energy> base;
    friend class _ignore::UnitOperators<Energy>;

private:
    // value is always in Joules
    explicit Energy (double value) : base(value){}

public:
   ~Energy(){}

    // default value
    Energy() : base(){}
    // copy constructor
    Energy (const Energy& E2) : base(E2){}
    Energy (const base&   E2) : base(E2){}

    // nice named constructors
    static Energy
        joules        (double J),
        electronVolts (double eV);

    // various return values
    double
        joules       () const,
        electronVolts() const;

    // handy constants
    static const Energy
        joule       ,   J,
        electronVolt,   eV;

    // cast operators
    operator Length()       const; // to wavelength
    operator Temperature () const; // to temperature

    // Human-readable string representation
    std::string toString() const;
};



/// Force = Mass * Length / Time / Time
class Force final
    : public _ignore::UnitOperators<Force>
{
    typedef _ignore::UnitOperators<Force> base;
    friend class _ignore::UnitOperators<Force>;

private:

    // value is always in Newtons
    explicit Force(double value) : base(value){}

public:
   ~Force(){}

    // default value
    Force() : base(){}
    // copy constructor
    Force (const Force& F2) : base(F2){}
    Force (const base&  F2) : base(F2){}

    // nice named constructors
    static Force
        newtons (double N),  N (double n);

    // various return values
    double
        newtons() const,   N () const;

    // constants
    static const Force
        newton;

    /// Human-readable string representation
    std::string toString() const;
};


/// Speed = Length / Time
class Speed final
    : public _ignore::UnitOperators<Speed>
{
    typedef _ignore::UnitOperators<Speed> base;
    friend class _ignore::UnitOperators<Speed>;

private:

    // value is always given in m/s
    explicit Speed (double value) : base(value){}

public:
    ~Speed(){}

    // default value
    Speed() : base(){}
    // copy  constructor
    Speed (const Speed& s2) : base(s2){}
    Speed (const base&  s2) : base(s2){}

    // nice named constructors
    static Speed
        metersPerSecond     (double mps),    mps  (double mps),
        knots               (double knt),    kt   (double knt),
        kilometersPerHour   (double kph),    kph  (double kph),
        kilometersPerSecond (double kps),    kps  (double kps),
        milesPerHour        (double mph),    miph (double mph),
        milesPerSecond      (double mps),    mips (double mps);

    // various return values
    double
        metersPerSecond     () const,    mps  () const,
        knots               () const,    kt   () const,
        kilometersPerHour   () const,    kph  () const,
        kilometersPerSecond () const,    kps  () const,
        milesPerHour        () const,    miph () const,
        milesPerSecond      () const,    mips () const;

    // constants
    static const Speed
        meterPerSecond    ,
        knot              ,
        kilometerPerHour  ,
        kilometerPerSecond,
        milePerHour       ,
        milePerSecond     ,

        of_light          , c;

    // NOTE: in principle, speed is a magnitude and thus always positive.
    // However, negative speeds are used, so include unary minus operator.

    // Human-readable string representation
    std::string toString() const;
};


/// Pressure = Force / Area
class Pressure final
    : public _ignore::UnitOperators<Pressure>
{
    typedef _ignore::UnitOperators<Pressure> base;
    friend class _ignore::UnitOperators<Pressure>;

private:

     // value is always given in N/m²
    explicit Pressure (double value) : base(value){}


public:
    ~Pressure(){}

    // default value
    Pressure() : base(){}
    // copy  constructors
    Pressure (const Pressure& p2) : base(p2){}
    Pressure (const base&     p2) : base(p2){}

    // nice named constructors
    static Pressure
        pascals     (double Pa ),    Pa (double Pa ),
        kiloPascals (double kPa),   kPa (double kPa),
        bars        (double bar),
        atmospheres (double atm),   atm(double atm);

    // various return values

    // constants

    /// Human-readable string representation
    std::string toString() const;
};

/// Density = Mass / Volume
class Density final
    : public _ignore::UnitOperators<Density>
{
    typedef _ignore::UnitOperators<Density> base;
    friend class _ignore::UnitOperators<Density>;

private:

     // value is always given in kg/m³
    explicit Density (double value) : base(value){}

public:
    ~Density(){}

    // default value
    Density() : base(){}
    // copy  constructors
    Density (const Density& d2) : base(d2){}
    Density (const base&    d2) : base(d2){}

    // nice named constructors
    static Density
        kilogramsPerCubicMeter  (double kgm3),
        gramsPerCubicCentimeter (double gcm3),
        tonnesPerCubicMeter     (double Mgm3),
        kilogramsPerLiter       (double kgl ),
        poundsPerCubicInch      (double lbscuin),
        poundsPerCubicFoot      (double lbscuft);

    // various return values
    double
        kilogramsPerCubicMeter  (),
        gramsPerCubicCentimeter (),
        tonnesPerCubicMeter     (),
        kilogramsPerLiter       (),
        poundsPerCubicInch      (),
        poundsPerCubicFoot      ();


    // constants
    static const Density
        kgm3,    kilogramPerCubicMeter,
        gcm3,    gramPerCubicCentimeter,
        Mgm3,    tonPerCubicMeter,         Tm3,
        kgl ,    kilogramPerLiter,
        lbscuin, poundPerCubicInch,        lbsin,
        lbscuft, poundPerCubicFoot,        lbsft;

    /// Human-readable string representation
    std::string toString() const;

};


/// overload various relevant math operators

// Angle
// ---------------------------------

// wrap value into [-pi +pi) / [-180 +180)
Angle wrap_posneg(const Angle& a);

// wrap value into [0 2pi) / [0 360)
Angle wrap_positive(const Angle& a);

// regular trig functions
double sin  (const Angle& a);
double cos  (const Angle& a);
double tan  (const Angle& a);
// FIXME: namespace conflict:
//double csc  (const Angle& a);
double sec  (const Angle& a);
double cot  (const Angle& a);

// hyperbolic trig functions
double sinh (const Angle& a);
double cosh (const Angle& a);
double tanh (const Angle& a);
double csch (const Angle& a);
double sech (const Angle& a);
double coth (const Angle& a);

// Inverse trig functions need to be overloaded by return type,
// which is "obviously" not possible. Therefore, define a few
// leaf helper classes:

namespace _ignore
{
    class Invtrig : public UnitOperators<Angle> {
    protected:
        double value;
    public:
        operator double();
        operator Angle ();
        std::string toString() const; // NOTE: leave unimplemented
    };
}

// regular inverse functions
class Asin   final : public _ignore::Invtrig { public: Asin (double value); };
class Acos   final : public _ignore::Invtrig { public: Acos (double value); };
class Atan   final : public _ignore::Invtrig { public: Atan (double value); };
class Atan2  final : public _ignore::Invtrig { public: Atan2(double y, double x); };
class Acsc   final : public _ignore::Invtrig { public: Acsc (double value); };
class Asec   final : public _ignore::Invtrig { public: Asec (double value); };
class Acot   final : public _ignore::Invtrig { public: Acot (double value); };

// hyperbolic inverse functions
class Asinh  final : public _ignore::Invtrig { public: Asinh(double value); };
class Acosh  final : public _ignore::Invtrig { public: Acosh(double value); };
class Atanh  final : public _ignore::Invtrig { public: Atanh(double value); };
class Acsch  final : public _ignore::Invtrig { public: Acsch(double value); };
class Asech  final : public _ignore::Invtrig { public: Asech(double value); };
class Acoth  final : public _ignore::Invtrig { public: Acoth(double value); };


/// Make  std::cout << {class}  work
std::ostream& operator<<(std::ostream& target, const _ignore::Unit& U);


/// Make some cross-unit transformations valid

// sqrt(Area) = Length
Length sqrt(const Area& A);

// Length * Length = Area
Area   operator*(const Length& L1, const Length& L2);

// Length * Area = Volume
Volume operator*(const Length&  L, const Area&   A);
Volume operator*(const Area&    A, const Length& L);

// Mass / Volume = Density
Density operator/(const Mass&   M, const Volume& V);

// Length / Time = Speed
Speed operator/ (const Length&  L, const Time&   t);

// Speed * Time = Length
Length operator*(const Speed&   s, const Time&   t);
Length operator*(const Time&    t, const Speed&  s);

// Length * Angle = (arc)Length
Length operator*(const Length&  L, const Angle&  a);
Length operator*(const Angle&   a, const Length& L);

// Force * Length = Energy
Energy operator*(const Force&   F, const Length& L);
Energy operator*(const Length&  L, const Force&  F);

// Force / Area = Pressure
Pressure operator/(const Force& F, const Area&   A);


// Handy dandy physics utilities
// (defined HERE because otherwise the unit classes are incomplete)
namespace physics
{
    // FIXME: (Rody Oldenhuis) DEPRECATED
    // replaced by operator Length() in class Energy
    inline double energy2wavelength(const double e) {
        BOOST_ASSERT(e > 0);
        return h*c/e;
    }
    inline Length energy2wavelength(const Energy& E){
        BOOST_ASSERT(E > Energy::zero);
        return Length::meters(h*c/E.joules());
    }
    // -----------------------------------------

    // FIXME: (Rody Oldenhuis) DEPRECATED
    // replaced by operator Energy() in class Length
    inline double wavelength2energy(const double w) {
        BOOST_ASSERT(w > 0);
        return h*c/w;
    }
    inline Energy wavelength2energy(const Length& w) {
        BOOST_ASSERT(w > Length::zero);
        return Energy::joules(h*c/w.meters());
    }
    // -----------------------------------------

    // FIXME: (Rody Oldenhuis) DEPRECATED
    // replaced by operator Temperature() in class Energy
    inline double energy2temperature(const double e) {
        BOOST_ASSERT(e >= 0);
        return e/k;
    }
    inline Temperature energy2temperature(const Energy& E) {
        BOOST_ASSERT(E >= Energy::zero);
        return Temperature::kelvin(E.joules()/k);
    }
    // -----------------------------------------

    // FIXME: (Rody Oldenhuis) DEPRECATED
    // replaced by operator Energy() in class Temperature
    inline double temperature2energy(const double T) {
        BOOST_ASSERT(T >= 0);
        return T*k;
    }
    inline Energy temperature2energy(const Temperature& T) {
        BOOST_ASSERT(T >= Temperature::zero);
        return Energy::joules(T.kelvin()*k);
    }
    // -----------------------------------------
}




#endif



