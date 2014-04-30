/*
 * Rody Oldenhuis
 * oldenhuis@gmail.com
 *
 * ode.hh
 * Created: 16.11.2012 16:29:30 CET
 */

#ifndef _ODE_HH
#define _ODE_HH

#include "Matrix.hh"

namespace ODE
{
    // defines options for the integrator
    struct Options
    {
        double
            abstol,
            reltol,
            maxstep,
            initial_step;

        bool
            show_stats;

        // Which integrator will be used?
        enum integrator_t
        {
            // Runge-Kutta-Nystrom types (RKN)
            RKN1210,

            // Runge-Kutta types (RK)
            ODE86,
            ODE45

        } integrator;


        // ------------------------------------
        // Event functions
        // ------------------------------------

        // options for event functions
        struct eventOpts_t
        {
            std::vector<unsigned int>
                eventSel;

            int8_t
                direction;

            bool
                isTerminal;
        };
        eventOpts_t eventOpt;

        // vector of pair of event function/options
        std::vector<
            std::pair<
                        // value,  isterminal, direction
                std::tuple<double, bool,       int8_t>   (*)(
                    const double &t, const Vector<> &y, // const Vector<> &dy,
                    const eventOpts_t &opts),
                eventOpts_t>
            > eventFcns;


        // ------------------------------------
        // Output functions
        // ------------------------------------

        // options for output functions
        struct outputOpts_t
        {
            std::vector<unsigned int>
                outputSel;

            double
                refine;
        };
        outputOpts_t outputOpt;

        // vector of pair of output function/options
        std::vector<
            std::pair<
                // halt
                bool    (*)(
                    const double &t, const Vector<> &y, // const Vector<> &dy,
                    const std::string &type,
                    const outputOpts_t &opts),
                //
                outputOpts_t>
            > outputFcns;


        // ------------------------------------
        // Default options
        // ------------------------------------

        // constructor sets default options
        Options()
            : abstol         (1e-12)
            , reltol         (1e-6)

            , maxstep        (0) // 0 means full interval
            , initial_step   (0) // 0 means use default

            , show_stats     (false)

            , integrator     (RKN1210)

            , eventOpt       (eventOpts_t())
            , eventFcns      {}

            , outputOpt      (outputOpts_t())
            , outputFcns     {}

        {}

    };


    // defines an output structure containing
    // performance statistics of the integration
    struct Stats
    {
        unsigned int
            rejected_steps,
            accepted_steps,
            fevals;

        Vec
            step_size,
            delta;

        std::string
            message;

        // constructor sets defaults
        Stats()
            : rejected_steps(0)
            , accepted_steps(0)
            , fevals        (0)
            , step_size     (Vec())
            , delta         (Vec())
            , message       ("Integration not yet started.")
        {}

    };



    // Call any type integrator
    void
    integrate(
        Vec (*dydt) (
            const double /*t*/, const Vec &/*y*/),
        const Vec &y0,
        const Vec &tspan,
        const Options &opts,
        Vec &t, Vec &y, double &exitflag, Stats &output);

}


#endif
