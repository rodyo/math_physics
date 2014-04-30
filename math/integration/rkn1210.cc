#include <cmath>
#include <sstream>
#include <iostream>

#include "ode.hh"
#include "rkn1210.hh"

namespace ODE
{

    rkn1210_result
    rkn1210(
        Vector<> (*d2ydt2) (
            const double /*t*/, const Vector<>& /*&y*/),
        const Vector<> &y0, const Vector<> &dy0, const Vector<> &tspan,
        const Options &opts)
    {

        // Check input
        // TODO

        // load the coefficients
        // ------------------------------------------------------------
#if 1

        // c
        static const Vec c_tmp = {
            0.0e0,
            2.0e-2,
            4.0e-2,
            1.0e-1,
            1.33333333333333333333333333333e-1,
            1.6e-1,
            5.0e-2,
            2.0e-1,
            2.5e-1,
            3.33333333333333333333333333333e-1,
            5.0e-1,
            5.55555555555555555555555555556e-1,
            7.5e-1,
            8.57142857142857142857142857143e-1,
            9.45216222272014340129957427739e-1,
            1.0e0,
            1.0e0};
        static const Vector<> c(c_tmp);

        // A
        static const Vec _A_tmp = {
                                0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            2.000000000000000e-04,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            2.666666666666667e-04,     5.333333333333334e-04,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            2.916666666666667e-03,    -4.166666666666667e-03,     6.250000000000000e-03,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            1.646090534979424e-03,                         0,     5.486968449931412e-03,     1.755829903978052e-03,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            1.945600000000000e-03,                         0,     7.151746031746032e-03,     2.912711111111111e-03,     7.899428571428571e-04,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            5.664062500000000e-04,                         0,     8.809730489417989e-04,    -4.369212962962963e-04,     3.390066964285714e-04,    -9.946469907407407e-05,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            3.083333333333333e-03,                         0,                         0,     1.777777777777778e-03,     2.700000000000000e-03,     1.578282828282828e-03,     1.086060606060606e-02,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            3.651839374801130e-03,                         0,     3.965171714072343e-03,     3.197258262930628e-03,     8.221467306855435e-03,    -1.313092695957238e-03,     9.771586968064868e-03,     3.755769069232834e-03,                         0,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            3.707241068718501e-03,                         0,     5.082045854555286e-03,     1.174708002175412e-03,    -2.114762991512699e-02,     6.010463698107881e-02,     2.010573476850619e-02,    -2.835075012293358e-02,     1.487956891858193e-02,                         0,                         0,                         0,                         0,                         0,                         0,                         0,
            3.512537656073344e-02,                         0,    -8.615749195138479e-03,    -5.791448051007917e-03,     1.945554823782616e+00,    -3.435123867456514e+00,    -1.093070110747522e-01,     2.349638311899517e+00,    -7.560094086870229e-01,     1.095289722215693e-01,                         0,                         0,                         0,                         0,                         0,                         0,
            2.052779253748250e-02,                         0,    -7.286446764480180e-03,    -2.115355607961840e-03,     9.275807968723522e-01,    -1.652282484425737e+00,    -2.107956300568657e-02,     1.206536432620787e+00,    -4.137144770010661e-01,     9.079873982809654e-02,     5.355552600533985e-03,                         0,                         0,                         0,                         0,                         0,
            -1.432407887554552e-01,                         0,     1.252870377309182e-02,     6.826019163969827e-03,    -4.799555395574387e+00,     5.698625043951941e+00,     7.553430369523645e-01,    -1.275548785828108e-01,    -1.960592605111738e+00,     9.185609056635262e-01,    -2.388008550528443e-01,     1.591108135723422e-01,                         0,                         0,                         0,                         0,
            8.045019205520489e-01,                         0,    -1.665852706701125e-02,    -2.141583404262973e-02,     1.682723592896247e+01,    -1.117283535717610e+01,    -3.377159297226324e+00,    -1.524332665536085e+01,     1.717983573821542e+01,    -5.437719239823995e+00,     1.387867161836466e+00,    -5.925827732652812e-01,     2.960387317129735e-02,                         0,                         0,                         0,
            -9.132967666973580e-01,                         0,     2.411272575780518e-03,     1.765812269386174e-02,    -1.485164977972038e+01,     2.158970867004576e+00,     3.997915583117880e+00,     2.843415180023223e+01,    -2.525936435494160e+01,     7.733878542362238e+00,    -1.891302894847867e+00,     1.001484507022472e+00,     4.641199599109052e-03,     1.121875502214896e-02,                         0,                         0,
            -2.751962972055940e-01,                         0,     3.661188877915492e-02,     9.789519688231562e-03,    -1.229306234588621e+01,     1.420722645393790e+01,     1.586647690678954e+00,     2.457773532759595e+00,    -8.935193694403273e+00,     4.373672731613407e+00,    -1.834718176544949e+00,     1.159208528906149e+00,    -1.729025316538392e-02,     1.932597790446077e-02,     5.204442937554993e-03,                         0,
            1.307639184740406e+00,                         0,     1.736410918974584e-02,    -1.854445645426580e-02,     1.481152203286773e+01,     9.383176308482470e+00,    -5.228426199944542e+00,    -4.895128052584765e+01,     3.829709603433793e+01,    -1.058738133697598e+01,     2.433230437622627e+00,    -1.045340604257544e+00,     7.177320950867259e-02,     2.162210970808278e-03,     7.009595759602514e-03,                         0};
        static const Matrix<> A(17, 16, _A_tmp);

        // Bhat (high-order b)
        static const Vec _Bhat_tmp = {
            +1.21278685171854149768890395495e-2,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +8.62974625156887444363792274411e-2,
            +2.52546958118714719432343449316e-1,
            -1.97418679932682303358307954886e-1,
            +2.03186919078972590809261561009e-1,
            -2.07758080777149166121933554691e-2,
            +1.09678048745020136250111237823e-1,
            +3.80651325264665057344878719105e-2,
            +1.16340688043242296440927709215e-2,
            +4.65802970402487868693615238455e-3,
            +0.0e0,
            +0.0e0};
        static const Vector<> Bhat(_Bhat_tmp);

        // BprimeHat (high-order b-prime)
        static const Vec _Bphat_tmp  = {
            +1.21278685171854149768890395495e-2,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +9.08394342270407836172412920433e-2,
            +3.15683697648393399290429311645e-1,
            -2.63224906576909737811077273181e-1,
            +3.04780378618458886213892341513e-1,
            -4.15516161554298332243867109382e-2,
            +2.46775609676295306562750285101e-1,
            +1.52260530105866022937951487642e-1,
            +8.14384816302696075086493964505e-2,
            +8.50257119389081128008018326881e-2,
            -9.15518963007796287314100251351e-3,
            +2.5e-2};
        static const Vector<> Bphat(_Bphat_tmp);

        // B (low-order b)
        static const Vec _B_tmp = {
            +1.70087019070069917527544646189e-2,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +7.22593359308314069488600038463e-2,
            +3.72026177326753045388210502067e-1,
            -4.01821145009303521439340233863e-1,
            +3.35455068301351666696584034896e-1,
            -1.31306501075331808430281840783e-1,
            +1.89431906616048652722659836455e-1,
            +2.68408020400290479053691655806e-2,
            +1.63056656059179238935180933102e-2,
            +3.79998835669659456166597387323e-3,
            +0.0e0,
            +0.0e0};
        static const Vector<> B(_B_tmp);

        // Bprime (low-order bprime)
        static const Vec _Bp_tmp = {
            +1.70087019070069917527544646189e-2,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +0.0e0,
            +7.60624588745593757356421093119e-2,
            +4.65032721658441306735263127583e-1,
            -5.35761526679071361919120311817e-1,
            +5.03182602452027500044876052344e-1,
            -2.62613002150663616860563681567e-1,
            +4.26221789886109468625984632024e-1,
            +1.07363208160116191621476662322e-1,
            +1.14139659241425467254626653171e-1,
            +6.93633866500486770090602920091e-2,
            +2.0e-2,
            +0.0e0};
        static const Vector<> Bp(_Bp_tmp);

#endif


        // Initialize
        // ------------------------------------

        bool
            have_eventFcn  = (opts.eventFcns.size()  != 0),
            have_outputFcn = (opts.outputFcns.size() != 0);

        rkn1210_result
            results = {};

        Vector<>
            yt  = y0,
            dyt = dy0;

        Matrix<>
            f = y0*zeros<>(1u,17u);

        double
            tt = tspan(0),
            t0 = tt,
            tfinal = tspan(tspan.numel-1),
            power = 1/12,
            hmin = abs(tfinal-tt)/1e12,
            hmax = abs(tfinal-tt),
            h = 0,
            h2 = 0;

        // times might be given in reverse
        int
            direction = 1 - 2*(tspan.numel == 2 && tfinal < t0);


        // initial step
        try {
            f(ALL,1) = d2ydt2(tt, yt);
            results.output.fevals++;
        }
        catch (...)
        {
            //ME2 = MException('rkn1210:incorrect_funfcnoutput', sprintf(...
                //'Derivative function should return a %3.0f-element column vector.', numel(y0)));
            //rethrow(addCause(ME,ME2));
        }

        if (!opts.initial_step) {
            // default initial step
            h = std::pow(opts.abstol,power) / max( max( max(abs(dyt)), max(f(ALL,1)) ), 1e-4 );
            h = h>hmin ? (h<hmax?h:hmax) : hmin;
        }
        else
            // user provided initial step
            h = opts.initial_step;

        h *= direction; // take care of direction
        h2 = h*h; // pre-compute the square
        results.output.step_size.push_back(h);


        // main loop
        while ( abs(tt-tfinal) > 0 )
        {
            // take care of final step
            if ((tt + h) > tfinal){
                h  = (tfinal - tt);
                h2 = h*h;
            }

            // Compute the second-derivative
            for (unsigned int j=0; j<17; ++j){
                f(ALL, j) = d2ydt2( tt + c(j)*h, yt + c(j)*h*dyt + h2*f*A(ALL,j) );
                results.output.fevals++;
            }

            // check for inf or NaN
            if ( any(!isfinite(f(ALL))) ) {
                results.exitflag = -2;
                // use warning (not error) to preserve output thus far
                //warning('rkn1210:nonfinite_values',...
                //     ['INF or NAN value encountered during the integration.\n',...
                //      'Terminating integration...']);
                return finalize(results, opts, hmin);
            }

            // pre-compute the sums of the products with the coefficients
            Matrix<>
                fBphat = f*Bphat,
                fBhat  = f*Bhat;

            // Estimate the error and the acceptable error
            double
                delta1 = max(abs(h2*(fBhat - f*B))), // error ~ |Y - y|
                delta2 = max(abs(h*(fBphat - f*Bp))),// error ~ |dot{Y} - dot{y}|
                delta  = max(delta1, delta2);        // use the maximum of these errors

            // update the solution only if the error is acceptable
            if ( (delta <= opts.abstol) &&
                (delta / max(norm(yt+h*dyt+h2*fBhat),norm(dyt+h*fBphat)) <= opts.reltol) )
            {

                // update the new solution
                tt  += h;
                yt  += h*dyt + h2*fBhat;
                dyt += h*fBphat;

                // This construction is WAY better than growing the arrays on
                // EACH iteration; especially for "cheap" integrands, this
                // construction causes a lot less overhead.
                results.t.push_back(tt);
                results.y.push_back(yt);
                results.dy.push_back(dyt);
                results.output.step_size.push_back(h);
                results.output.delta.push_back(delta);
                results.output.accepted_steps++;

                // event functions
                if (have_eventFcn)
                {
                    /*

                    // number of functions provided
                    num_events = numel(Event);
                    // initialize TE (event times), YE (event solutions) YPE (event
                    // derivs) and IE (indices to corresponding event function). Check
                    // user-provided event functions at the same time
                    previous_event_values = zeros(num_events,1);
                    for k = 1:num_events
                        try
                            previous_event_values(k) = feval(Event{k}, t, y, dy);
                        catch ME
                            ME2 = MException('rkn1210:eventFcn_dont_evaluate',...
                                sprintf('Event function #%1d failed to evaluate on initial call.', k));
                            rethrow(addCause(ME,ME2));
                        end
                    end
                    TE = []; YE = []; YPE = []; IE = [];

                    */
                }

                // output functions
                if (have_outputFcn)
                {
                    /*

                    // WHICH elements should be passed to the output function?
                    OutputSel = odeget(options, 'OutputSel', 1:numel(y0));
                    *
                    // adjust the number of points passed to the output functions by this factor
                    Refine = odeget(options, 'Refine', 1);
                    *
                    // number of functions provided
                    num_outputFcn = numel(OutputFcn);
                    // call each output function with 'init' flag. Also check whether
                    // the user-provided output function evaluates
                    for k = 1:num_outputFcn
                        try
                            feval(OutputFcn{k}, t, y(OutputSel), dy(OutputSel), 'init');
                        catch ME
                            ME2 = MException('rkn1210:OutputFcn_dont_evaluate',...
                                sprintf('Output function #%1d failed to evaluate on initial call.', k));
                            rethrow(addCause(ME,ME2));
                        end
                    end
                    *

                */

                    // evaluate all output-functions
                    for (auto k : opts.outputFcns)
                    {
                        // evaluate output function
                        try {
                            bool halt = k.first( t, yt(k.second.outputSel), dyt(k.second.outputSel), "iter" );
                        }
                        catch (...) {
                            // ME2 = MException('rkn1210:OutputFcn_failure_integration',...
                            //     sprintf('Output function #%1d failed to evaluate during integration.', k));
                            // rethrow(addCause(ME,ME2));
                        }

                        // halt integration when requested
                        if (halt) {
                            results.exitflag = 2;
                            return finalize(results, opts, hmin);
                        }
                    }
                }

            }

            // rejected step: just increase its counter
            else
                results.output.rejected_steps++;




            // adjust the step size
            if (delta)
            {
                h  = sign(h)*min(hmax, 0.9*abs(h)*std::pow(opts.abstol/delta,power));

                // Use [Refine]-option when output functions are present
                if (have_outputFcn)
                    h /= opts.outputOpt.refine;

                // pre-compute the square
                h2 = h*h;

                // check the new stepsize
                if (abs(h) < hmin) {
                    results.exitflag = -1;
                    // use warning to preserve results thus far
                    //warning('rkn1210:stepsize_too_small', ...
                    //    ['Failure at time t = %6.6e: \n',...
                    //    'Step size fell below the minimum acceptible value of %6.6d.\n',...
                    //    'A singularity is likely.'], t, hmin);
                }

            } // adjust step-size

        } // main loop

        // if the algorithm ends up here, all was ok
        return finalize(results, opts, hmin);

#if 0

        /* Different case: numel(tspan) > 2


        // do the calculation recursively if [tspan] has
        // more than two elements
        if (numel(tspan) > 2)
            // the output times are already known
            tout = tspan(:);
            // call this function as many times as there are times in [tspan]
            for i = (1:numel(tspan)-1)
                // new initial values
                tspanI = tspan(i:i+1);
                y0     = yout(end, :);
                yp0    = dyout(end, :);
                // call the integrator
                if have_events
                    [toutI, youtI, dyoutI, TEI, YEI, DYEI, IEI, exitflag, outputI] = ...
                        rkn1210(funfcn, tspanI, y0, yp0, options);
                else
                    [toutI, youtI, dyoutI, exitflag, outputI] = ...
                        rkn1210(funfcn, tspanI, y0, yp0, options);
                end
                // append the solutions
                yout  = [yout;  youtI(end, :)];  %#ok
                dyout = [dyout; dyoutI(end, :)]; %#ok
                if have_events
                    TE   = [TE; TEI];  %#ok
                    YEI  = [YE; YEI];  %#ok
                    DYEI = [YPE; DYEI];%#ok
                    IEI  = [IE; IEI];  %#ok
                end
                // process the output
                output.h        = [output.h; outputI.h];
                output.fevals   = output.fevals + outputI.fevals;
                output.rejected = output.rejected + outputI.rejected;
                output.accepted = output.accepted + outputI.accepted;
                output.delta    = [output.delta; outputI.delta];
                // evaluate any output functions at each [t] in [tspan]
                if have_outputFcn
                    // evaluate all functions
                    for k = 1:num_outputFcn
                        try
                            halt = feval(OutputFcn{k}, toutI, youtI(OutputSel), dyoutI(OutputSel), []);
                        catch ME
                            ME2 = MException('rkn1210:OutputFcn_failure_integration',...
                                sprintf('Output function #%1d failed to evaluate during integration.', k));
                            rethrow(addCause(ME,ME2));
                        end
                    end
                    // halt integration when requested
                    if (halt) {
                        results.exitflag = 2;
                        return return finalize(results, opts, hmin);
                    }
                end
                // should we quit?
                if exitflag == -2 || exitflag == 3, break, end
            end
            // we're done.
            return finalize(results, opts, hmin);
        end

        */


        // evaluate event-funtions
        if have_events
            // evaluate all event-functions
            for k = 1:num_events
                // evaluate event function, and check if any have changed sign
                %
                // NOTE: although not really necessary (event functions have been
                // checked upon initialization), use TRY-CATCH block to produce
                // more useful errors in case something does go wrong.
                terminate = false;
                try
                    // evaluate function
                    [value, isterminal, zerodirection] = ...
                        feval(Event{k}, t, y, dy);

                    // look for sign change
                    if (previous_event_values(k)*value < 0)
                        // ZERODIRECTION:
                        // 0: detect all zeros (default
                        // +1: detect only INcreasing zeros
                        // -1: detect only DEcreasing zeros
                        if (zerodirection == 0) ||...
                        (sign(value) == sign(zerodirection))
                            // terminate?
                            terminate = terminate || isterminal;
                            // detect the precise location of the zero
                            // NOTE: try-catch is necessary to prevent things like
                            // discontinuous event-functions from resulting in
                            // unintelligible error messages
                            try
                                detect_Events(k, tp, previous_event_values(k), t, value);
                            catch ME
                                ME2 = MException('rkn1210:eventFcn_failure_zero',...
                                    sprintf('Failed to locate a zero for event function #%1d.', k));
                                throwAsCaller(addCause(ME,ME2));
                            end
                        end
                    end

                    // save new value
                    previous_event_values(k) = value;

                catch ME
                    ME2 = MException('rkn1210:eventFcn_failure_integration',...
                        sprintf('Event function #%1d failed to evaluate during integration.', k));
                    rethrow(addCause(ME,ME2));
                end

                // do we need to terminate?
                if (terminate) {
                    results.exitflag = 3;
                    return finalize(results, opts, hmin);
                }

            end
        end // Event functions

        // evaluate output functions
        if  have_outputFcn
            // evaluate all event-functions
            for k = 1:num_outputFcn
                // evaluate output function
                try
                    halt = feval(OutputFcn{k}, t, y(OutputSel), dy(OutputSel), []);
                catch ME
                    ME2 = MException('rkn1210:OutputFcn_failure_integration',...
                        sprintf('Output function #%1d failed to evaluate during integration.', k));
                    rethrow(addCause(ME,ME2));
                end
                // halt integration when requested
                if (halt) {
                    results.exitflag = 2;
                    return finalize(results, opts, hmin);
                }
            end
        end // Output functions


    #endif

    }

    #if 0

    // Detect events
    // NOTE: simple Regula-Falsi method to detect the zero
    detect_Events(which_event, ta, fa, tb, fb)

        // initialize
        iterations = 0;
        maxiterations = 1e4;
        opts = options;
        y0 = yp;
        dy0 = dyp;
        tt = ta;

        // prune unnessesary options
        opts = odeset(opts,...
            'Event'    , [],...
            'OutputFcn', [],...
            'Stats'    , 'off');

        // start iteration
        while ( min(abs(fa),abs(fb)) > opts.abstol )
        {
            // increase no. of iterations
            iterations = iterations + 1;
            // make step
            ttp = tt;
            tt  = (0.5*fb*ta - fa*tb) / (0.5*fb-fa);

            // safety precaution
            if (ttp == tt), break, end

            // evaluating the event-function at this new trial location is
            // somewhat complicated. We need to recursively call this
            // RKN1210-routine to get appropriate values for [y] and [dy] at
            // the new time [tt] into the event function:
            [~, Zyout, Zdyout, ~, Zoutput] = ...
                rkn1210(funfcn, [ttp, tt], y0, dy0, opts);

            // save old values for next iteration
            y0 =  Zyout(end,:).';    dy0 = Zdyout(end,:).';
            // NOW evaluate event-function with these values
            fval = feval(Event{which_event}, tt, y0, dy0);
            // keep track of number of function evaluations
            output.fevals = output.fevals + Zoutput.fevals;

            // compute new step
            if (fb*fval>0)
                tb = tt; fb = fval;
            else
                ta = tt; fa = fval;
            end

            // check no. of iterations
            if (iterations > maxiterations) {
                error('rkn1210:rootfinder_exceeded_max_iterations',...
                    'Root could not be located within %d iterations. Exiting...',...
                    maxiterations);
            }
        }

        // The zero has been found; insert values into proper arrays
        TE  = [TE; tt];   YPE = [YPE; dy0];
        YE  = [YE; y0];   IE  = [IE; which_event];

        // if new values won't fit, first grow the arrays
        if index > size(yout,1)
            grow_arrays; end

        // the integrand first overshoots the zero; that's how it's
        // detected. We want the zero to be in the final arrays, but we also
        // want them in chronological orderd. So move the overshoot one
        // down, and insert the zero in its place:
        index = index + 1;
        yout (index,:) = yout (index-1,:);     yout (index-1,:) = y0.';
        dyout(index,:) = dyout(index-1,:);     dyout(index-1,:) = dy0.';
        tout (index,:) = tout (index-1,:);     tout (index-1,:) = tt;
        output.h(index,:) = tt-tout(index-1,:);
        output.h(index,:) = tt-tout(index-1,:);

    end // find zeros of Event-functions


    #endif

    // clean up and finalize
    rkn1210_result
    finalize(rkn1210_result &results, const Options &opts, double hmin)
    {
        // final call to the output functions
        if (opts.outputFcns.size() != 0) {
            for (auto kk : opts.outputFcns) {
                try {
                    ////feval(OutputFcn{kk}, t, y(OutputSel), dy(OutputSel), 'done');
                }
                catch (...){
                    ////ME2 = MException('rkn1210:OutputFcn_failure_finalization',...
                    ////    sprintf('Output function #%1d failed to evaluate during finalization call.', k));
                    ////rethrow(addCause(ME,ME2));
                }
            }
        }

        // add message to output structure
        switch (results.exitflag)
        {
            case +0:
                results.output.message = "Integration completed sucessfully.";
                break;

            case +2:
                results.output.message = "Integration terminated by one of the output functions.";
                break;

            case +3:
                results.output.message = "Integration terminated by one of the event functions.";
                break;

            case -1: {
                std::stringstream msg;
                msg
                    << "Integration terminated successfully, but the step size\n"
                    << "fell below the minimum allowable value of " << hmin <<  "\n"
                    << "for one or more steps. Results may be inaccurate.\n";
                results.output.message = msg.str();
                break;
            }

            case -2:
                results.output.message = "Integration unsuccessful; second derivative function \n \
                    returned a NaN or INF value.";
                break;
        }

        // display stats
        if (opts.show_stats){
            std::cout
                << std::endl << std::endl
                << "Number of successful steps     : " << results.output.accepted_steps << std::endl
                << "Number of rejected steps       : " << results.output.rejected_steps << std::endl
                << "Number of evaluations of f(t,y): " << results.output.fevals << std::endl << std::endl
                << results.output.message << std::endl;
        }

    } // finalize the integration


} // ODE namespace










