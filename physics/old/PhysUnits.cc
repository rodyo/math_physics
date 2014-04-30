/*
 * Rody Oldenhuis
 * cosine measurement systems
 * roldenhuis@cosine.nl
 *
 * PhysUnits.cc
 * Created: 05.10.2012 10:03:41 CEST
 */

#include "PhysUnits.hh"
#include "Exception.hh"


namespace _ignore
{
    template <typename T> const T UnitOperators<T>::zero         =  T();
    template <typename T> const T UnitOperators<T>::dx           =  T(math::dx);
    template <typename T> const T UnitOperators<T>::infinitesmal =  T(math::dx);
    template <typename T> const T UnitOperators<T>::small        =  T(math::dx);
    template <typename T> const T UnitOperators<T>::NaN          =  T(math::NaN);
    template <typename T> const T UnitOperators<T>::nan          =  T(math::NaN);
    template <typename T> const T UnitOperators<T>::inf          =  T(math::inf);
    template <typename T> const T UnitOperators<T>::infinite     =  T(math::inf);
    template <typename T> const T UnitOperators<T>::infinity     =  T(math::inf);
}

template <typename T> T  abs(const _ignore::UnitOperators<T>& t) { return T(std::fabs(t.value)); }
template <typename T> T fabs(const _ignore::UnitOperators<T>& t) { return T(std::fabs(t.value)); }


// Angle
// --------------------------------------------------------------

// named constructors
Angle Angle::radians    (double rad) { return Angle( rad);                       }
Angle Angle::degrees    (double deg) { return Angle( deg*math::pi/180.0);        }
Angle Angle::arcminutes (double moa) { return Angle( moa*math::pi/180.0/60.0);   }
Angle Angle::arcseconds (double soa) { return Angle( soa*math::pi/180.0/3600.0); }
Angle Angle::GRADIANS   (double grad){ return Angle(grad*math::pi/200.0);        }

Angle Angle::operator/ (const Angle& a2) const { return Angle(value/a2.value); }
Angle Angle::operator* (const Angle& a2) const { return Angle(value*a2.value); }

// constant angles
const Angle
    Angle::tau  = Angle(math::tau),
    Angle::_360 = Angle(math::tau),

    Angle::pi   = Angle(math::pi),
    Angle::_180 = Angle(math::pi),

    Angle::pio2 = Angle(math::pio2),
    Angle::_90  = Angle(math::pio2),

    Angle::pio4 = Angle(math::pio4),
    Angle::_45  = Angle(math::pio4),

    Angle::pio6 = Angle(math::pio6),
    Angle::_30  = Angle(math::pio6);

// different return types
double Angle::radians    () const { return value;                       }
double Angle::degrees    () const { return value*180.0/math::pi;        }
double Angle::arcminutes () const { return value*180.0*60.0/math::pi;   }
double Angle::arcseconds () const { return value*180.0*3600.0/math::pi; }
double Angle::GRADIANS   () const { return value*200.0/math::pi;        }

// wrap value into [-pi +pi) / [-180 +180)
Angle& Angle::wrap_posneg() {
    value = math::mod(value + math::pi, math::tau) - math::pi;
    return *this;
}

// wrap value into [0 2pi) / [0 360)
Angle& Angle::wrap_positive(){
    value = math::mod(value, math::tau);
    return *this;
}

// Human-readable string representation
std::string Angle::toString() const {
    std::stringstream output;
    output << degrees() << "°";
    return output.str();
}

// wrap value into [-pi +pi) / [-180 +180)
Angle wrap_posneg(const Angle& a){
    return Angle::radians(math::mod(a.value + math::pi, math::tau) - math::pi);
}

// wrap value into [0 2pi) / [0 360)
Angle wrap_positive(const Angle& a){
    return Angle::radians(math::mod(a.value, math::tau));
}


// overload various relevant math operators
//-----------------------------------------

Angle fabs  (const Angle& a) { return Angle::radians( fabs(a.radians()) );}
Angle  abs  (const Angle& a) { return Angle::radians( fabs(a.radians()) );}

// regular trig functions
double sin  (const Angle& a) { return std::sin(a.radians()); }
double cos  (const Angle& a) { return std::cos(a.radians()); }
double tan  (const Angle& a) { return std::tan(a.radians()); }
// FIXME: namespace conflict:
//double csc  (const Angle& a) { return 1.0/std::sin(a.radians()); }
double sec  (const Angle& a) { return 1.0/std::cos(a.radians()); }
double cot  (const Angle& a) { return 1.0/std::tan(a.radians()); }

// hyperbolic trig functions
double sinh (const Angle& a) { return std::sinh(a.radians()); }
double cosh (const Angle& a) { return std::cosh(a.radians()); }
double tanh (const Angle& a) { return std::tanh(a.radians()); }

double csch (const Angle& a) { return 1.0/std::sinh(a.radians()); }
double sech (const Angle& a) { return 1.0/std::cosh(a.radians()); }
double coth (const Angle& a) { return 1.0/std::tanh(a.radians()); }

// inverse trig functions need to be overloaded by return type,
// which is not possible in C. Therefore:

namespace _ignore
{
    std::string Invtrig::toString() const {
        std::stringstream output;
        output << Angle::radians(value).degrees() << "°";
        return output.str();
    }
    Invtrig::operator double() { return value;                 }
    Invtrig::operator Angle () { return Angle::radians(value); }

}

// regular trig

std::ostream& operator<<(std::ostream& target, const _ignore::Invtrig& T){
    target << T.toString();
    return target;
}

namespace
{   // Constrain arguments to valid values (on real plane)

    // regular trig
    double safe_acos (double v){ return std::acos ( (v<-1.0) ? -1.0 : (v>+1.0?+1.0:v) ); }
    double safe_asin (double v){ return std::asin ( (v<-1.0) ? -1.0 : (v>+1.0?+1.0:v) ); }

    // hyperbolic trig
    double safe_atanh(double v){ return std::atanh( (v<-1.0) ? -1.0 : (v>+1.0?+1.0:v) ); }
    double safe_acosh(double v){ return std::acosh( (v<+1.0) ? +1.0 :  v              ); }

}

Asin::Asin(double v) { value = safe_asin(v); }
Acos::Acos(double v) { value = safe_acos(v); }
Atan::Atan(double v) { value = std::atan(v); }
Acsc::Acsc(double v) { value = safe_asin(1.0/v); }
Asec::Asec(double v) { value = safe_acos(1.0/v); }
Acot::Acot(double v) { value = std::atan(1.0/v); }
Atan2::Atan2(double y, double x) { value = std::atan2(y,x); }

// hyperbolic inverse trig functions
Asinh::Asinh(double v) { value = std::asinh(v); }
Acosh::Acosh(double v) { value = safe_acosh(v); }
Atanh::Atanh(double v) { value = safe_atanh(v); }
Acsch::Acsch(double v) { value = std::asinh(1.0/v); }
Asech::Asech(double v) { value = safe_acosh(1.0/v); }
Acoth::Acoth(double v) { value = safe_atanh(1.0/v); }




// Length
// --------------------------------------------------------------

// various named constructors                                                     // and some aliases
Length Length::meters           (double m)   { return Length(m);             }    Length Length::m   (double   m) { return meters            (  m); }
Length Length::kilometers       (double km)  { return Length(km *si::k);     }    Length Length::km  (double  km) { return kilometers        ( km); }
Length Length::nanometers       (double nm)  { return Length(nm *si::n);     }    Length Length::nm  (double  nm) { return nanometers        ( nm); }
Length Length::microns          (double mu)  { return Length(mu *si::mu);    }    Length Length::um  (double  um) { return micrometers       ( um); }
Length Length::micrometers      (double mu)  { return Length(mu *si::mu);    }
Length Length::inches           (double in)  { return Length(in *0.0254);    }    Length Length::in  (double  in) { return inches            ( in); }
Length Length::feet             (double ft)  { return Length(ft *0.304800);  }    Length Length::ft  (double  ft) { return feet              ( ft); }
Length Length::foot             (double ft)  { return feet(ft);              }
Length Length::yards            (double yrd) { return Length(yrd*0.914400);  }
Length Length::miles            (double mi)  { return Length(mi *1609.344);  }    Length Length::mi  (double  mi) { return miles             ( mi); }
Length Length::nauticalMiles    (double nmi) { return Length(nmi*1852.000);  }    Length Length::nmi (double nmi) { return nauticalMiles     (nmi); }
Length Length::astronomicalUnits(double au)  { return Length(au*physics::au);}    Length Length::AU  (double  au) { return astronomicalUnits ( au); }
Length Length::angstroms        (double A)   { return Length(A*physics::A);  }    Length Length::A   (double   a) { return angstroms         (  a); }

// various return values                                                  // and some aliases
double Length::meters           () const { return value;             }    double Length::m  () const { return meters            (); }
double Length::kilometers       () const { return value/si::k;       }    double Length::km () const { return kilometers        (); }
double Length::nanometers       () const { return value/si::n;       }    double Length::nm () const { return nanometers        (); }
double Length::microns          () const { return value/si::mu;      }    double Length::um () const { return microns           (); }
double Length::micrometers      () const { return value/si::mu;      }
double Length::inches           () const { return value/0.0254;      }    double Length::in () const { return inches            (); }
double Length::feet             () const { return value*3.28084;     }    double Length::ft () const { return feet              (); } double Length::foot() const { return feet(); }

double Length::yards            () const { return value*1.09361;     }
double Length::miles            () const { return value*0.000621371; }    double Length::mi () const { return miles             (); }
double Length::nauticalMiles    () const { return value*0.000539957; }    double Length::nmi() const { return nauticalMiles     (); }
double Length::astronomicalUnits() const { return value/physics::au; }    double Length::AU () const { return astronomicalUnits (); }
double Length::angstroms        () const { return value/physics::A;  }    double Length::A  () const { return angstroms         (); }


// cast operators
Length::operator Energy() const { return value >= 0.0 ? Energy::joules(physics::h*physics::c/value) : Energy::NaN; }

// various constants
const Length
    Length::one_km = Length( si::k),         Length::kilometer       = Length::one_km,
    Length::one_m  = Length(   1.0),         Length::meter           = Length::one_m ,
    Length::one_mm = Length( si::m),         Length::millimeter      = Length::one_mm,
    Length::one_um = Length(si::mu),         Length::micrometer      = Length::one_um,
    Length::micron = Length(si::mu),
    Length::one_nm = Length( si::n),         Length::nanometer       = Length::one_nm,
    Length::one_AU = Length(physics::au),    Length::astronomicalUnit= Length::one_AU;

// Human-readable string representation
std::string Length::toString() const {
    std::stringstream output;
    output << value << " m";
    return output.str();
}


// Time
// --------------------------------------------------------------

// nice named constructors
Time Time::microseconds(double ms){ return Time(ms*si::mu); }
Time Time::milliseconds(double ms){ return Time(ms*si::m);  }
Time Time::seconds     (double s) { return Time(s);         }
Time Time::minutes     (double m) { return Time(m*60.0);    }
Time Time::hours       (double h) { return Time(h*3600.0);  }
Time Time::days        (double d) { return Time(d*86400.0); }
Time Time::years       (double y) { return Time(y*365.2563835*86400); }
Time Time::centuries   (double c) { return Time(c*36525.63835*86400); }



// various return values
double Time::microseconds() const { return value/si::mu;  }
double Time::milliseconds() const { return value/si::m;   }
double Time::seconds     () const { return value;         }
double Time::minutes     () const { return value/60.0;    }
double Time::hours       () const { return value/3600.0;  }
double Time::days        () const { return value/86400.0; }
double Time::years       () const { return value/365.2563835/86400; }
double Time::centuries   () const { return value/36525.63835/86400; }

// Get current time in seconds+microseconds past 1/1/1970 00:00
Time Time::now(){
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return Time( (double)tv.tv_sec + ((double)tv.tv_usec*1.0e-6) );
}
boost::posix_time::ptime Time::now(void *P){
    return boost::posix_time::second_clock::universal_time();
}

// constants
const Time
    Time::microsecond = microseconds(1.0),
    Time::millisecond = milliseconds(1.0),
    Time::second      = seconds     (1.0),
    Time::minute      = minutes     (1.0),
    Time::hour        = hours       (1.0),
    Time::day         = days        (1.0),
    Time::year        = years       (1.0),
    Time::century     = centuries   (1.0);

// Human-readable string representation
std::string Time::toString() const {
    std::stringstream output;
    output << value << " s";
    return output.str();
}




// Mass
// --------------------------------------------------------------

// nice named constructors
Mass Mass::kg         (double kg ) { return Mass(kg);              }
Mass Mass::kilogrammes(double kg ) { return Mass(kg);              }
Mass Mass::kilograms  (double kg ) { return Mass(kg);              }
Mass Mass::lbs        (double pnd) { return Mass(pnd * 0.4535920); }
Mass Mass::pounds     (double pnd) { return lbs(pnd);              }
Mass Mass::tonnes     (double T)   { return Mass(T   * si::k);     }
Mass Mass::longTonnes (double LT)  { return Mass(LT  * 1016.0500); }
Mass Mass::shortTonnes(double sT)  { return Mass(sT  * 907.18500); }
Mass Mass::stones     (double st)  { return Mass(st  * 6.3502900); }
Mass Mass::ounces     (double st)  { return Mass(st  * 0.0283495); }

// various return values
double Mass::kg         () const { return value;           }
double Mass::kilogrammes() const { return value;           }
double Mass::kilograms  () const { return value;           }
double Mass::lbs        () const { return value/0.4535920; }
double Mass::pounds     () const { return lbs();           }
double Mass::tonnes     () const { return value/1000.0000; }
double Mass::longTonnes () const { return value/1016.0500; }
double Mass::shortTonnes() const { return value/907.18500; }
double Mass::stones     () const { return value/6.3502900; }
double Mass::ounces     () const { return value/0.0283495; }

// handy dandy constant values
const Mass
    Mass::kilo       = kg         (1.0),
    Mass::pound      = lbs        (1.0),
    Mass::tonne      = tonnes     (1.0),
    Mass::longTonne  = longTonnes (1.0),
    Mass::shortTonne = shortTonnes(1.0),
    Mass::stone      = stones     (1.0),
    Mass::ounce      = ounces     (1.0);

// Human-readable string representation
std::string Mass::toString() const {
    std::stringstream output;
    output << value << " kg";
    return output.str();
}




// Temperature
// --------------------------------------------------------------

// nice named constructors
Temperature Temperature::kelvin     (double K) { return Temperature(K);                         }
Temperature Temperature::celsius    (double C) { return Temperature(C+273.15);                  }
Temperature Temperature::fahrenheit (double F) { return Temperature((F-32.0)*5.0/9.0 + 273.15); }

// various return values
double Temperature::kelvin     () const { return  value                        ; }
double Temperature::celsius    () const { return  value-273.15                 ; }
double Temperature::fahrenheit () const { return (value-273.15)*9.0/5.0 + 32.0 ; }

// NOTE: negative temperatures do not exist; excelude unary minus operator
Temperature& Temperature::operator-() {
    throw csc::ex::NotImplemented("On the absolute temperature scale, negative temperatures do not exist.");
    return *this;
}

// cast operators
Temperature::operator Energy() const { return value >= 0.0 ? Energy::joules(physics::k*value) : Energy::NaN; }

// constants
const Temperature
    Temperature::K = kelvin(0),
    Temperature::C = celsius(0),
    Temperature::F = fahrenheit(0);

// Human-readable string representation
std::string Temperature::toString() const {
    std::stringstream output;
    output << value << " K";
    return output.str();
}



// Area
// --------------------------------------------------------------

// named constructors
Area Area::squareMeters     (double m2  ) { return Area(m2);              }
Area Area::squareKilometers (double km2 ) { return Area(km2*si::k*si::k); }
Area Area::squareMillimeters(double mm2 ) { return Area(mm2*si::m*si::m); }
Area Area::ares             (double are ) { return Area(are*si::h);       }
Area Area::acres            (double acre) { return Area(acre*4046.86);    }
Area Area::hectares         (double ha  ) { return Area(ha*si::h*si::h);  }

// and some of their aliases
Area Area::m2  (double m2  ) { return squareMeters      (m2); }
Area Area::km2 (double km2 ) { return squareKilometers (km2); }
Area Area::mm2 (double mm2 ) { return squareMillimeters(mm2); }

// various return values
double Area::squareMeters     () const { return value;             }
double Area::squareKilometers () const { return value/si::k/si::k; }
double Area::squareMillimeters() const { return value/si::m/si::m; }
double Area::ares             () const { return value/si::h;       }
double Area::acres            () const { return value/4046.86;     }
double Area::hectares         () const { return value/si::h/si::h; }

// and some of their aliases
double Area::m2 () const { return squareMeters     (); }
double Area::km2() const { return squareKilometers (); }
double Area::mm2() const { return squareMillimeters(); }

// constants
const Area
    Area::squareMeter      = m2       (1.0),
    Area::squareKilometer  = km2      (1.0),
    Area::squareMillimeter = mm2      (1.0),
    Area::are              = ares     (1.0),
    Area::acre             = acres    (1.0),
    Area::hectare          = hectares (1.0);

// Human-readable string representation
std::string Area::toString() const {
    std::stringstream output;
    output << value << " m²";
    return output.str();
}



// Volume
// --------------------------------------------------------------


// nice named constructors

Volume Volume::cubicKilometers (double km3) { return Volume(km3*si::k*si::k*si::k); }
Volume Volume::cubicMeters     (double  m3) { return Volume(m3                   ); }
Volume Volume::cubicCentimeters(double cm3) { return Volume(cm3*si::c*si::c*si::c); }
Volume Volume::cubicMillimeters(double mm3) { return Volume(mm3*si::m*si::m*si::m); }
Volume Volume::cubicFeet       (double ft3) { return Volume(ft3*0.0283168        ); }
Volume Volume::cubicInches     (double in3) { return Volume(in3*1.6387e-5        ); }
Volume Volume::litres          (double l  ) { return Volume(l  *si::d*si::d*si::d); }
Volume Volume::gallons         (double g  ) { return Volume(g  *0.00454609       ); } // imperial gallon
Volume Volume::USgallons       (double g  ) { return Volume(g  *0.00378541       ); } // US gallon

// various return values
double Volume::cubicKilometers  () const { return value/si::k/si::k/si::k;}
double Volume::cubicMeters      () const { return value;                  }
double Volume::cubicCentimeters () const { return value/si::c/si::c/si::c;}
double Volume::cubicMillimeters () const { return value/si::m/si::m/si::m;}
double Volume::cubicFeet        () const { return value/0.0283168;        }
double Volume::cubicInches      () const { return value/1.6387e-5;        }
double Volume::litres           () const { return value/si::d/si::d/si::d;}
double Volume::gallons          () const { return value/0.00454609;       }
double Volume::USgallons        () const { return value/0.00378541;       }

// constants
const Volume
    Volume::cubicKilometer  = Volume::cubicKilometers (1.0),   Volume::km3 = cubicKilometer,
    Volume::cubicMeter      = Volume::cubicMeters     (1.0),   Volume::m3  = cubicMeter,
    Volume::cubicCentimeter = Volume::cubicCentimeters(1.0),   Volume::cm3 = cubicCentimeter,
    Volume::cubicMillimeter = Volume::cubicMillimeters(1.0),   Volume::mm3 = cubicMillimeter,
    Volume::cubicFoot       = Volume::cubicFeet       (1.0),   Volume::ft3 = cubicFoot,
    Volume::cubicInch       = Volume::cubicInches     (1.0),   Volume::in3 = cubicInch,
    Volume::litre           = Volume::litres          (1.0),   Volume::  l = litre,
    Volume::gallon          = Volume::gallons         (1.0),   Volume::  g = gallon,
    Volume::USgallon        = Volume::USgallons       (1.0),   Volume:: Ug = USgallon;


// Human-readable string representation
std::string Volume::toString() const {
    std::stringstream output;
    output << value << " m³";
    return output.str();
}



// Energy
// --------------------------------------------------------------

// nice named constructors
Energy Energy::joules        (double  J) { return Energy(J);             }
Energy Energy::electronVolts (double eV) { return Energy(eV*physics::q); }

// various return values
double Energy::joules()        const { return value;            }
double Energy::electronVolts() const { return value/physics::q; }

// cast operators

// Energy -> wavelength
Energy::operator Length()       const { return value >= 0.0 ? Length::meters(physics::h*physics::c/value) : Length::NaN;      }
Energy::operator Temperature () const { return value >= 0.0 ? Temperature::kelvin(value/physics::k)       : Temperature::NaN; }

// handy constant
const Energy
    Energy::joule        = joules(1.0),        Energy::J  = joule,
    Energy::electronVolt = electronVolts(1.0), Energy::eV = electronVolt;

// Human-readable string representation
std::string Energy::toString() const {
    std::stringstream output;
    output << value << " J";
    return output.str();
}



// Force
// --------------------------------------------------------------

// nice named constructors
Force Force::newtons (double N) { return Force(N); }

Force Force::N(double n) { return Force(n); }

// various return values
double Force::newtons   () const { return value; }

double Force::N()          const { return newtons(); }

// constants
const Force
    Force::newton = Force(1.0);

// Human-readable string representation
std::string Force::toString() const {
    std::stringstream output;
    output << value << " N";
    return output.str();
}



// Speed
// --------------------------------------------------------------

// nice named constructors
Speed Speed::metersPerSecond     (double mps) { return Speed ( mps                                                      ); }
Speed Speed::knots               (double knt) { return Speed ( Length::nauticalMiles(knt).meters()/Time::hour.seconds() ); }
Speed Speed::kilometersPerHour   (double kph) { return Speed ( Length::kilometers(kph).meters()   /Time::hour.seconds() ); }
Speed Speed::kilometersPerSecond (double kps) { return Speed ( Length::kilometers(kps).meters()                         ); }
Speed Speed::milesPerHour        (double mph) { return Speed ( Length::miles(mph).meters()        /Time::hour.seconds() ); }
Speed Speed::milesPerSecond      (double mps) { return Speed ( Length::miles(mps).meters()                              ); }

Speed Speed::mps  (double mps) { return metersPerSecond     (mps);}
Speed Speed::kt   (double knt) { return knots               (knt);}
Speed Speed::kph  (double kph) { return kilometersPerHour   (kph);}
Speed Speed::kps  (double kps) { return kilometersPerSecond (kps);}
Speed Speed::miph (double mph) { return milesPerHour        (mph);}
Speed Speed::mips (double mps) { return milesPerSecond      (mps);}

// various return values
double Speed::metersPerSecond     () const { return value;                                                      }
double Speed::knots               () const { return Length::meters(value).nauticalMiles()/Time::second.hours(); }
double Speed::kilometersPerHour   () const { return Length::meters(value).kilometers()   /Time::second.hours(); }
double Speed::kilometersPerSecond () const { return Length::meters(value).kilometers();                         }
double Speed::milesPerHour        () const { return Length::meters(value).miles()        /Time::second.hours(); }
double Speed::milesPerSecond      () const { return Length::meters(value).miles();                              }

double Speed::mps  () const { return metersPerSecond     ();}
double Speed::kt   () const { return knots               ();}
double Speed::kph  () const { return kilometersPerHour   ();}
double Speed::kps  () const { return kilometersPerSecond ();}
double Speed::miph () const { return milesPerHour        ();}
double Speed::mips () const { return milesPerSecond      ();}

// some constant speeds
const Speed
    Speed::meterPerSecond     = metersPerSecond     (1.0),
    Speed::knot               = knots               (1.0),
    Speed::kilometerPerHour   = kilometersPerHour   (1.0),
    Speed::kilometerPerSecond = kilometersPerSecond (1.0),
    Speed::milePerHour        = milesPerHour        (1.0),
    Speed::milePerSecond      = milesPerSecond      (1.0),

    Speed::of_light = Speed(physics::c),
    Speed::c        = of_light;

// Human-readable string representation
std::string Speed::toString() const {
    std::stringstream output;
    output << value << " m/s";
    return output.str();
}





// Density
// --------------------------------------------------------------

// named constructors
Density Density::kilogramsPerCubicMeter  (double kgm3)    { return Density(kgm3);       }
Density Density::gramsPerCubicCentimeter (double gcm3)    { return Density(gcm3*si::k); }
Density Density::tonnesPerCubicMeter     (double Mgm3)    { return Density(Mgm3*si::k); }
Density Density::kilogramsPerLiter       (double kgl )    { return Density(kgl*si::k);  }
Density Density::poundsPerCubicInch      (double lbscuin) { return Density(lbscuin*27679.9047); }
Density Density::poundsPerCubicFoot      (double lbscuft) { return Density(lbscuft*16.0184634); }

// various return values
double Density::kilogramsPerCubicMeter  () { return value;      }
double Density::gramsPerCubicCentimeter () { return value/si::k;}
double Density::tonnesPerCubicMeter     () { return value/si::k;}
double Density::kilogramsPerLiter       () { return value/si::k;}
double Density::poundsPerCubicInch      () { return value/27679.9047;}
double Density::poundsPerCubicFoot      () { return value/16.0184634;}


// constants
const Density
    Density::kgm3     = kilogramsPerCubicMeter  (1.0), Density::kilogramPerCubicMeter  = kgm3   ,
    Density::gcm3     = gramsPerCubicCentimeter (1.0), Density::gramPerCubicCentimeter = gcm3   ,
    Density::Mgm3     = tonnesPerCubicMeter     (1.0), Density::tonPerCubicMeter       = Mgm3   ,  Density::Tm3   = tonPerCubicMeter,
    Density::kgl      = kilogramsPerLiter       (1.0), Density::kilogramPerLiter       = kgl    ,
    Density::lbscuin  = poundsPerCubicInch      (1.0), Density::poundPerCubicInch      = lbscuin,  Density::lbsin = poundPerCubicInch,
    Density::lbscuft  = poundsPerCubicFoot      (1.0), Density::poundPerCubicFoot      = lbscuft,  Density::lbsft = poundPerCubicFoot;

/// Human-readable string representation
std::string Density::toString() const {
    std::stringstream output;
    output << value << " kg/m³";
    return output.str();
}




// Cross-unit transformations
// -=-=-=-=-=-=-=-=-=-=-=-=-=

// sqrt(Area) = Length
Length sqrt(const Area& A) { return Length::meters(std::sqrt(A.squareMeters())); }

// Length * Length = Area
Area operator*(const Length& L1, const Length& L2) { return Area::squareMeters(L1.m()*L2.m()); }

// Length * Area = Volume
Volume operator*(const Length& L1, const Area&   A2) { return Volume::cubicMeters(L1.m()*A2.squareMeters()); }
Volume operator*(const Area&   A1, const Length& L2) { return L2*A1; }

// Mass / Volume = Density
Density operator/(const Mass&  M, const Volume& V) { return Density::kilogramsPerCubicMeter( M.kilograms()/V.cubicMeters() ); }

// Length / Time = Speed
Speed operator/ (const Length& L, const Time&   t) { return Speed::metersPerSecond( L.m()/t.seconds() ); }

// Speed * Time = Length
Length operator*(const Speed&  s, const Time&   t) { return Length::meters(s.metersPerSecond()*t.seconds()); }
Length operator*(const Time&   t, const Speed&  s) { return s*t; }

// Length * Angle = (arc)Length
Length operator*(const Length& L, const Angle&  a) { return Length::meters( L.m() * a.radians() ); }
Length operator*(const Angle&  a, const Length& L) { return L*a; }

// Force * Length = Energy
Energy operator*(const Force&  F, const Length& L) { return Energy::joules( F.N() * L.m() ); }
Energy operator*(const Length& L, const Force&  F) { return F*L; }


// Force / Area = Pressure
//Pressure operator/(const Force& F, const Area&   A) { return Pressure::Pa(F.N()/A.m2()); }


// Make  std::cout << {class}  work
std::ostream& operator<<(std::ostream& target, const _ignore::Unit& U) { target << U.toString(); return target; }

