// useful definitions


// standard includes
#include "math.h"
#include <limits>

namespace Math{

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	// useful global constants
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	const double inf = std::numeric_limits<double>::infinity(); // infinity (e.g., 1/0)
	const double NaN = std::numeric_limits<double>::quiet_NaN();// Not-a-Number (e.g., 0/0)
	const double  pi = 3.1415926535897;  // pi of course...	
	const double   c = 299792.458;       // speed of light [km/s]
	const double  AU = 149597870.700;    // astronomical unit [km]
	const double   G = 6.67300e-20;      // astronomical unit [km3/kg/s2]
	const double  g0 = 9.80665e-3;       // avg. gravitational acceleration at sealevel [km/s2]

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	// functions
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

	/* define abs() */
	template <class genType>
	   inline genType abs(const genType x){return ((x>0.0)?x:-x);}

	/* signum function */
	template <class genType>
	inline int sign(const genType x){return ((x>0.0)-(x<0.0));}

	/* NAND */
	inline bool nand(bool a, bool b){return (!(a&&b));}
	
	/* min, max */
	template <class genType>
	   inline genType max(const genType a, const genType b ){return (b<a)?a:b;}
	template <class genType>
	   inline genType min(const genType a, const genType b ){return (b<a)?b:a;}
	   
	/* "safe" acos/asin */
	template <class genType>
		inline genType acos(const genType x){return (std::acos(min(1.0, max(-1.0,x))));}
	template <class genType>
		inline genType asin(const genType x){return (std::asin(min(1.0, max(-1.0,x))));}

	/* csc,sec,cot, and csch,sech,coth */
	template <class genType>
	   inline genType csc (const genType x){return (1.0/std::sin(x));}
	template <class genType>
	   inline genType sec (const genType x){return (1.0/std::cos(x));}
	template <class genType>
	   inline genType cot (const genType x){return (1.0/std::tan(x));}
	template <class genType>
	   inline genType csch(const genType x){return (1.0/std::sinh(x));}
	template <class genType>
	   inline genType sech(const genType x){return (1.0/std::cosh(x));}
	template <class genType>
	   inline genType coth(const genType x){return (1.0/std::tanh(x));}

	/* acsc,asec,acot and acsch,asech,acoth */
	template <class genType>
	   inline genType acsc (const genType x){return (std::asin(1/x));}
	template <class genType>
	   inline genType asec (const genType x){return (std::acos(1/x));}
	template <class genType>
	   inline genType acot (const genType x){
		  if (x>0) return (std::atan(1/x));
		  else return (pi+std::atan(1/x));
	}
	//inline double acot2(double x);
	//inline double acsch(double x);
	//inline double asech(double x);
	//inline double acoth(double x);
	
	/* cross and dot products */
	template <class genType>
		inline genType* cross(const genType a[], const genType b[]){
			genType c[] = {a[1]*b[2]-a[2]*b[1], /
			               a[2]*b[0]-a[0]*b[2], /
						   a[0]*b[1]-a[1]*b[0]};
			return &c;		  
		}
	template <class genType>
		inline genType dot(const genType a[], const genType b[]){
			return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);		
		}
		

	/* hypot = sqrt(a*a + b*b) (preventing overflow) */
	template <class genType>
	   inline genType hypot(const genType a, const genType b){return (a/std::sin(std::atan(a/b)));

	/* round towards zero */
	template <class genType>
	   inline genType fix(const genType x){return ( (x<0)?std::ceil(x):std::floor(x) );}
	   
   
};