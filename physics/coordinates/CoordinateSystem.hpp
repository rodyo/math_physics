#ifndef __COORDINATESYSTEM_HPP
#define __COORDINATESYSTEM_HPP

#include <string>

#include "Quaternion.hpp"
#include "ScalarUnit.hpp"

using namespace math;

// NOTE: THIS IS A SINGLETON!
class BaseSystem final
{
public:
	BaseSystem* getInstance();
	~BaseSystem();
	
	constexpr bool operator== (BaseSystem &other);
	
private:
	BaseSystem();	
	
};


// NOTE: vectorunit instances need a coordinate system in which the vector is defined
//       that wayt, when doing operations (cross, dot, etc.), automatically the correct
//       coordinate system is used. So: rethink the design!
class CoordinateSystem final
{
public:
	Direction  linearOffset;
	Quaternion rotationalOffset;
	
	std::string name;

public:
	CoordinateSystem(BaseSystem &base, 
					 Direction  &linearOffset,
					 Quaternion &rotationalOffset);
					 
	CoordinateSystem(CoordinateSystem &parent, 
					 Direction        &linearOffset = Direction::zero,
					 Quaternion       &rotationalOffset);

	template<class T>
	VectorUnit<T>
	inParent(T x, T y, T z) { 
	}
	
	template<class T>
	VectorUnit<T>
	inBase(Tx, Ty, Tz) {
	}

private:
	BaseSystem*       base;
	CoordinateSystem* parent;
};

#endif
