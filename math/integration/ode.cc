#include <string>
#include <vector>

#include "ode.hh"

namespace ODE
{

    // Call Runge-Kutta type integrator
    void
    integrate(
        (void)(*dydt) (
            const double t, const Vec &y,
                                  Vec &dy),
        const Vec &y0,
        const Vec &tspan,
        const odeOptions &opts,
              Vec &t,
              Vec &y)
    {
        switch (opts.integrator)
        {
            // Error
            case (RKN1210):
                break;

            //
            case (ODE86):
            case (ODE45):
                break;

            default:
                break;
        }




    }


    // Call Runge-Kutta-Nystrom type integrator
    void
    integrate(
        (void)(*d2ydt2) (
            const double t, const Vec &y, const Vec &dy,
                                  Vec &dy,      Vec &ddy),
        const Vec &y0,
        const Vec &dy0,
        const Vec &tspan,
        const odeOptions &opts,
              Vec &t,
              Vec &y,
              Vec &dy)
    {
        switch (opts.integrator)
        {
            //
            case (RKN1210):
                rkn1210(d2ydt2, y0,dy0, tspan, opts, t,y,dy);
                break;

            // Error
            case (ODE86):
            case (ODE45):
                break;

            default:
                break;
        }
    }

}
