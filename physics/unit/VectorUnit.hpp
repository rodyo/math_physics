#ifndef __VECTORUNIT_HPP
#define __VECTORUNIT_HPP

#include <iostream>

#include "Matrix.hpp"
#include "ScalarUnit.hpp"
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
class VectorUnit
        : public math::Vector<ScalarUnit<L,M,r,C,T,I,N,i>>
{
    // TODO: also do named constructurs (useful for things like ForceComponent::zero;
    // TODO: Vectorunits have an associated coorddinate system. All operations will
    //       take place only in the SAME coordinate system

    // TODO: Adding forces might give rise to torques (forces not at the same
    // point of origin), etc.

};

// Aliases
typedef VectorUnit<Length>                           Direction;

typedef VectorUnit<ForceComponent>                   Force;
typedef VectorUnit<TorqueComponent>                  Torque;
typedef Torque                                       Moment;

typedef VectorUnit<Speed>                            Velocity;
typedef VectorUnit<AccelerationComponent>            Acceleration;
typedef VectorUnit<JerkComponent>                    Jerk;
typedef VectorUnit<MomentumComponent>                Momentum;
typedef VectorUnit<SpecificMomentumComponent>        SpecificMomentum;

typedef VectorUnit<AngularSpeed>                     AngularVelocity;
typedef VectorUnit<AngularAccelerationComponent>     AngularAcceleration;
typedef VectorUnit<AngularJerkComponent>             AngularJerk;
typedef VectorUnit<AngularMomentumComponent>         AngularMomentum;
typedef VectorUnit<SpecificAngularMomentumComponent> SpecificAngularMomentum;


#endif
