/*
 * Rody Oldenhuis
 * oldenhuis@gmail.com
 *
 * rkn1210.hh
 * Created: 29.11.2012 21:11:28 CET
 */

#ifndef _ODE_RKN1210_HH
#define _ODE_RKN1210_HH

#include "ode.hh"

namespace ODE
{
    /**
     * RKN1210       12th/10th order Runge-Kutta-Nystrom integrator
     *
     * RKN1210() is a 12th/10th order numerical integrator for ordinary
     * differential equations of the special form
     *
     * y'' = f(t, y)                     (1)
     *
     * with initial conditions
     *
     * y(t0) = y0,  y'(t0) = yp0         (2)
     *
     * This second-order differential equation is integrated with a
     * Runge-Kutta-Nystrom method, with 17 function evaluations per step. The
     * RKN-class of integrators is especially suited for this purpose, since
     * compared to a classic Runge-Kutta integration scheme the same accuracy
     * can be obtained with less function evaluations.
     *
     * This RKN12(10) method is a very high-order method, to be used in problems
     * with *extremely* stringent error tolerances. As the name implies, the
     * error should be less than O(h^13). In verious studies, it has been shown
     * that this particular integration technique is overall more efficient for
     * ODE's of the form (1) than multi-step or extrapolation methods that give
     * the same accuracy.
     *
     * RKN1210()'s behavior is very similar MATLAB's ODE-integrator suite:
     *
     * USAGE:
     * ----------------
     *
     * [t, y, yp] = RKN1210(funfcn, tspan, y0, yp0)
     * [t, y, yp] = RKN1210(funfcn, tspan, y0, yp0, options)
     *
     * [t, y, yp, exitflag, output] = RKN1210(...)
     * [t, y, yp, TE, YE, YPE, IE, exitflag, output] = RKN1210(...)
     *
     * INPUT ARGUMENTS:
     * ----------------
     *
     * funfcn  - definition of the second-derivative function f(t, y)
     *             (See (1)). It should accept a scalar value [t] and a
     *             column vector [y] which has the same number of elements
     *             as the initial values [y0] and [dy0] provided.
     *
     *     tspan - time interval over which to integrate. It can be a
     *             two-element vector, in which case it will be interpreted
     *             as an interval. In case [tspan] has more than two
     *             elements, the integration is carried out to all times in
     *             [tspan]. Only the values for those times are then
     *             returned in [y] and [yp].
     *
     * y0, yp0 - initial values as in (2). Both should contain the same number
     *             of elements.
     *
     * options - options structure, created with ODESET(). Used options are
     *             MaxStep, InitialStep, AbsTol, Stats, Event, OutputFcn,
     *             OutputSel, and Refine. See the help for ODESET() for more
     *             information.
     *
     * How to use Event and/or Output functions is described in the documentation
     * on ODESET(). There is one difference: RKN1210() now also passes the first
     * derivative [yp] at each step as an argument:
     *
     * status = outputFcn(t, y, yp, flag)
     * [value, isterminal, direction] = event(t, y, yp)
     *
     * where [t] is scalar, and [y] and [yp] are column vectors, as with f(t,y).
     *
     *
     * OUTPUT ARGUMENTS:
     * ----------------
     *
     * t, y, yp - The approximate solutions for [y] and [y'] at times [t].
     *             All are concatenated row-wise, that is
     *
     *             t  = N-by-1
     *             y  = N-by-numel(y0)
     *             y' = N-by-numel(y0)
     *
     *             with N the number of sucessful steps taken during the
     *             integration.
     *
     * exitflag - A scalar value, indicating the termination conditions
     *             of the integration:
     *
     *             -2: a non-finite function value was encountered during the
     *                 integration (INF of NaN); the integration was stopped.
     *             -1: the step size [h] fell below  the minimum acceptable
     *                 value at some time(s) [t]; results may be inaccurate.
     *             0: nothing was done; initial state.
     *             +1: sucessful integration, normal exit.
     *             +2: integration was stopped by one of the output
     *                 functions.
     *             +3: One or more events were detected, and their
     *                 corresponding [isterminal] condition also evaluated to
     *                 [true].
     *
     * TE,YE,    - These arguments are only returned when one or more event
     * YPE,IE     functions are used. [TE] contains the times at which events
     *             were detected. [YE] and [YPE] lists the corresponding values
     *             of the solution [y] and the first derivative [yp] at these
     *             times. [IE] contains indices to the event-functions with
     *             which these events were detected. Use a smaller value for
     *             AbsTol (in [options]) to increase the accuracy of these
     *             roots when required.
     *
     *     output - structure containing additional information about the
     *             integration. It has the fields:
     *
     *             output.h              step size (sucesful steps only) at
     *                                 each time [tn]
     *             output.rejected       amount of rejected steps
     *             output.accepted       amount of accepted steps
     *             output.delta          estimate of the largest possible
     *                                 error at each time [tn]
     *             output.message        Short message describing the
     *                                 termination conditions
     *
     *             Note that these fields contain the information of ALL
     *             steps taken, even for cases where [tspan] contains
     *             more than 2 elements.
     *
     *
     * See also ODE45, ODE86, RKN86.
     *
     *
     * Based on the codes for ODE86 and RKN86, also available on the MATLAB
     * FileExchange.

     *
     *
     * ELEMENTARY EXAMPLE
     * (2-body gravitational interaction / circular orbit):
     *
     * f = @(t, y) [-y(1)/sqrt(y(1)^2+y(2)^2)^3;-y(2)/sqrt(y(1)^2+y(2)^2)^3];
     * [t, y] = rkn1210(f, [0, 50*pi], [1; 0], [0; 1]);
     * plot(y(:,1), y(:,2), '-k');  // orbit should be exactly circular
     * figure, plot( abs(y(:,1) - cos(t)), 'k'); // plot the x-error
     * hold on, plot( abs(y(:,2)- sin(t)), 'r'); // plot the y-error
     */

    // the return structure
    struct rkn1210_result
    {
        Vec
            t, y, dy,         // time, y-solution, yprime-solution
            TE, YE, YPE;      // Event solutions
        VecUI
            IE;               // Index to event function

        ODE::Stats
            output;

        int
            exitflag;

        // assign default values
        rkn1210_result()
            : t  {}, y  {}, dy  {}
            , TE {}, YE {}, YPE {}, IE {}
            , output {} , exitflag {}
        {}
    };

    /* The integrator
     *
     * The construction of RKN12(10) is described in
     * High-Order Embedded Runge-Kutta-Nystrom Formulae
     * J. R. DORMAND, M. E. A. EL-MIKKAWY, AND P. J. PRINCE
     * IMA Journal of Numerical Analysis (1987) 7, 423-430
     *
     * Coefficients obtained from
     * http:// www.tampa.phys.ucl.ac.uk/rmat/test/rknint.f
     * These are also available in any format on request to these authors.
     *
     *
     * Author:
     * Name    : Rody P.S. Oldenhuis
     * E-mail  : oldenhuis@gmail.com
     */
    rkn1210_result
    rkn1210(
        Vector<> (*d2ydt2) (
            const double /*t*/, const Vector<>& /*&y*/),
        const Vector<> &y0, const Vector<> &dy0, const Vector<> &tspan,
        const ODE::Options &opts);



    // complete the result
    rkn1210_result
    finalize(rkn1210_result &results, const ODE::Options &opts, double hmin);

}


#endif
