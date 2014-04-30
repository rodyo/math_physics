


// Spherical Geometry
namespace spherespace
{
    // Primary use case
    Coordinate::Coordinate(const Angle& lat, const Angle& lon)
        : latitude  (lat)
        , longitude (lon)
    {
        latitude.wrap_posneg();
        longitude.wrap_posneg();

        // check bounds
        if (latitude.degrees() > +90.0)
            latitude = Angle::degrees(+180.0) - latitude;
        if (latitude.degrees() < -90.0)
            latitude = Angle::degrees(-180.0) - latitude;
    }

    // Secondary use case
    Coordinate::Coordinate(const csc::Vector3& V)
    {
        // TODO: shouldn't this be in Geometry.cc?

        // NOTE: vector is *not* necessarily a unit vector, so
        // don't check for it. We should check for zero-length
        // vectors though:
        if (V.isZeroLength())
            throw csc::ex::Condition(
                "Coordinate: spherical coordinates not \
                defined for zero-length vector.");

        // Cartesian -> Spherical. Bounds are automatically correct.
        double hypotxy = hypot(V[0],V[1]);
        latitude  = Angle::radians( atan2(V[2], hypotxy) );
        longitude = Angle::radians( atan2(V[1],V[0]) );
    }

    // convert to vector/direction
    csc::Vector3 Coordinate::toVector(){
        return csc::make_direction_equatorial(
            latitude.radians(), longitude.radians());
    }

    // displayable string
    std::string Coordinate::toString() const {
        std::stringstream output;
        output
            << "Spherical coordinate:\n"
            << "Latitude : " << latitude << std::endl
            << "Longitude: " << longitude;
        return output.str();
    }

    // It's also convenient to have a toVector() as a free function
    csc::Vector3 toVector(const Coordinate& C){
        return C.toVector();
    }

}


namespace ellipsoidspace
{
    Coordinate(const Coordinate& p) : spherespace::Coordinate(p) {}


    // Primary use case
    Coordinate(const Angle& lat, const Angle& lon,
        double flattening, const Length& equatiorialRadius)
        : spherespace::Coordinate(lat,lon)
        , flattening       (flattening)
        , altitude         (Length::NaN)
        , equatiorialRadius(equatiorialRadius)
    {}

    // Secondary use case
    Coordinate(const csc::Vector3& V,
        double flattening, const Length& equatiorialRadius)
        : spherespace::Coordinate(V)
        , flattening       (flattening)
        , equatiorialRadius(equatiorialRadius)
    {
        // Closed-form solution is a quite involved. Here, an simpler,
        // iterative solution is used:

        double
            RN, h=0, cP, sP,
            phiP = math::inf,
            hP   = math::inf,
            R    = equatiorialRadius.meters(),
            phi  = latitude.radians(),
            e2   = (2.0-flattening)*flattening,
            p    = std::sqrt(V[0]*V[0]+V[1]+V[1]);

        while ( fabs(h-hp)>math::eps || fabs(h-hp)>math::eps )
        {
            phiP = phi;
            hp   = h;

            RN = R / std::sqrt(1-e2*sin(phi)*sin(phi));
            cP = cos(phi);

            if (fabs(cP) < math::eps)
                sP = sin(phi),
                h  = (V[2]+e2*sP)/sP - RN;
            else
                h = p/cP - RN;

            phi = atan2( V[2], p*(1-e2*(RN/(RN+h))) );
        }

        altitude = Length::meters(h);
        latitude = Angle::radians(phi);

    }

    // convert to vector
    csc::Vector3 Coordinate::toVector(){
        // closed form solution:
        double
            e2  = (2.0-flattening)*flattening,
            RN = equatiorialRadius.meters() /
                 std::sqrt(1-e2*sin(latitude)*sin(latitude));

        return csc::Vector3(
            (RN+altitude.meters())*cos(latitude)*cos(longitude),
            (RN+altitude.meters())*cos(latitude)*sin(longitude),
            ((1-e2)*RN+altitude.meters())*sin(latitude));
    }

    // displayable string
    std::string Coordinate::toString() const {
        std::stringstream output;
        output
            << "Ellipsoidal coordinate:\n"
            << "Latitude : " << latitude
            << "Longitude: " << longitude
            << "Altitude : " << altitude;
        return output.str();
    }

    // It's also convenient to have a toVector() as a free function
    csc::Vector3 toVector(const Coordinate& C){
        return C.toVector();
    }

}

// functions to make  std::cout << {class}  work
std::ostream& operator<<(std::ostream& target, const spherespace::Coordinate&    c){ target << c.toString(); return target; }
std::ostream& operator<<(std::ostream& target, const ellipsoidspace::Coordinate& c){ target << c.toString(); return target; }


// All things spherical trig
namespace shapes
{
    // spherical polygon
    // ------------------------------------

    SphericalPolygon::SphericalPolygon(){}


    /// human readable string representaiton
    std::string SphericalPolygon::toString() const {
        std::stringstream output;
        output << "SphericalPolygon: \n";
        // TODO
        return output.str();
    }


    // spherical triangle
    // ------------------------------------

    /**
     * 4-quadrant arccosine function
     *
     * Computes the four-quadrant arccosine of the value G. The
     * resulting value is not uniquely determined by G, nor by
     * the lengths or order of the sides of the triangle. Therefore,
     * an angle beta is required to recover this information.
     *
     * The information lies in the Hemisphere function @ref H(beta).
     * If H(beta) < pi/2, the acos() of the small angle
     * (0 <= alpha <= pi/2) will be returned. If H(beta) > pi/2, the
     * acos() of the large angle (pi/2 < alpha < pi) will be returned.
     *
     * @param[in] G The value to take the arccosine of.
     * @param[in] beta The angle deciding whether to use the acos() or its complement
     * @return the four-quadrant arccosine of alpha.
     */

    // Hemisphere function
    double SphericalTriangle::H(const Angle& beta){
        return (wrap_positive(beta) < Angle::pi) ? 1.0 : -1.0;
    }
    // 4-quadrant acos
    Angle SphericalTriangle::acos2(double G, const Angle& beta ){
        // nominal calculation
        double
            Hval = H(beta),
            retval = Hval * acos(G);
        // handle special case of acos(G) = 0
        if (retval==0)
            retval = Hval<0 ? math::pi : 0;

        return Angle::radians(retval);
    }


    Angle SphericalTriangle::mal(
        const Angle& A, const Angle& B,
        const Angle& a, const Angle& b)
    {
        // sine & cosine of side C
        double
            sinA = sin(A), cosA = cos(A),
            sinB = sin(B), cosB = cos(B),
            cosa = cos(a), cosb = cos(b),

            sinC =  sinA*cosB*cosb + sinB*cosA*cosa,
            cosC = -cosA*cosB      + sinB*sinA*cosa*cosb;

        // internal angle C:
        return Angle::radians(atan2(sinc,cosc)).wrap_positive();
    }

    /**
     * Solutions to the angle-angle-angle problem:
     *    "given all three angles, compute the sides"
     * This problem has zero or two solutions.
     *
     * @param[in] A first angle
     * @param[in] B second angle
     * @param[in] C third angle
     */
    void SphericalTriangle::aaa(const Angle& A, const Angle& B, const Angle& C)
    {
        // first assign everything obvious
        Angle
            Aa  = wrap_positive(A),
            Bb  = wrap_positive(B),
            Cc  = wrap_positive(C),
            tau = Angle::tau,
            pi  = Angle::pi,

            // FIXME: we've just wrapped them, does this make sense?
            fA  = fabs(Aa),
            fB  = fabs(Bb),
            fC  = fabs(Cc);

        A1 = Aa;   A2 = tau - Aa;
        B1 = Bb;   B2 = tau - Bb;
        C1 = Cc;   C2 = tau - Cc;

        double
            cosA = cos(Aa), sinA = sin(Aa),
            cosB = cos(Bb), sinB = sin(Bb),
            cosC = cos(Cc), sinC = sin(Cc);

        // check constraint
        if (fabs(pi-fA-fB) <= fC && fC <= pi-fabs(fA-fB))
           throw csc::ex::Condition("SphericalTriangle: AAA problem has no solution.");

        // first solution
        a1 = acos2( (cosA + cosB*cosC)/sinB/sinC, A);
        b1 = acos2( (cosB + cosA*cosC)/sinA/sinC, B);
        c1 = acos2( (cosC + cosA*cosB)/sinA/sinB, C);

        // second solution
        a2 = tau - a1;
        b2 = tau - b1;
        c2 = tau - c1;
    }

    /**
     * Solutions to the side-side-side problem:
     *    "given all three sides, compute the angles"
     * This problem has zero or two solutions.
     *
     * @param[in] a first side
     * @param[in] b second side
     * @param[in] c third side
     */
    void SphericalTriangle::sss(const Angle& a, const Angle& b, const Angle& c)
    {
        // first assign everything obvious
        Angle
            aa   = wrap_positive(a),
            bb   = wrap_positive(a),
            cc   = wrap_positive(a),
            tau  = Angle::tau,
            pi   = Angle::pi,
            fpma = fabs(pi-aa),
            fpmb = fabs(pi-bb),
            fpmc = fabs(pi-cc);

        a1 = aa;   a2 = tau - aa;
        b1 = bb;   b2 = tau - bb;
        c1 = cc;   c2 = tau - cc;

        double
            cosa = cos(aa), sina = sin(aa),
            cosb = cos(bb), sinb = sin(bb),
            cosc = cos(cc), sinc = sin(cc);

        // check condition
        if (fpma-fpmb <= fpmc && fpmc <= fpma+fpmb)
            throw csc::ex::Condition("SphericalTriangle: SSS problem has no solution.");

        // first solution
        A1 = acos2( (cosa - cosb*cosc)/sinb/sinc, a);
        B1 = acos2( (cosb - cosa*cosc)/sina/sinc, b);
        C1 = acos2( (cosc - cosa*cosb)/sina/sinb, c);

        // second solution
        A2 = tau - A1;
        B2 = tau - B1;
        C2 = tau - C1;
    }

    /**
     * Solutions to the angle-angle-side problem:
     *    "given two angles and one side opposite one of the angles,
     *     compute the two other sides and other angle."
     * This problem has zero or two solutions.
     *
     * @param[in] A first angle
     * @param[in] B second angle
     * @param[in] a side opposite A
     */
    void SphericalTriangle::aas(const Angle& A, const Angle& B, const Angle& a)
    {
        // first assign everything obvious
        Angle
            Aa  = wrap_positive(A),
            Bb  = wrap_positive(B),
            aa  = wrap_positive(a),
            pi  = Angle::pi,
            tau = Angle::tau;

        A1 = Aa;   A2 = tau - Aa;
        B1 = Bb;   B2 = tau - Bb;
        a1 = aa;   a2 = tau - aa;

        double
            sinA = sin(Aa), sinB = sin(Bb), sina = sin(aa),
            sinBsina = sinB*sina;

        // check condition
        if (fabs(sinBsina) > fabs(sinA))
            throw csc::ex::Condition("SphericalTriangle: AAS problem has no solution.");

        // first solution
        b1 = Angle::radians(asin( sinBsina/sinA )).wrap_positive();
        c1 = msl(aa, b1, Aa, Bb);
        C1 = mal(Aa, Bb, aa, b1);

        // second solution
        b2 = (pi - b1).wrap_positive();
        c2 = msl(aa, b2, Aa, Bb);
        C2 = mal(Aa, Bb, aa, b2);
    }

    /**
     * Solutions to the angle-side-angle problem:
     *    "given two angles and the side between their associated
     *     vertices, compute the two other sides and other angle."
     * This problem has two solutions.
     *
     * @param[in] A first angle
     * @param[in] B second angle
     * @param[in] c side not opposite A or B
     */
    void SphericalTriangle::asa(const Angle& A, const Angle& B, const Angle& c)
    {
        // first assign everything obvious
        Angle
            Aa  = wrap_positive(A),
            Bb  = wrap_positive(B),
            cc  = wrap_positive(c),
            pi  = Angle::pi(),
            tau = Angle::tau();

        A1 = Aa;   A2 = tau - Aa;
        B1 = Bb;   B2 = tau - Bb;
        c1 = cc;   c2 = tau - cc;

        double
            cosA = cos(Aa),  sinA = sin(Aa),
            cosB = cos(Bb),  sinB = sin(Bb),
            cosC          ,  sinC;

        // first solution
        // NOTE: normal acos (in stead of acos2) is indeed correct.
        C1 = Acos(-cosA*cosB + sinA*sinB*cos(cc));
        cosC = cos(C1);
        sinC = sin(C1);
        a1 = Acos( (cosA + cosB*cosC)/sinB/sinC );
        b1 = Acos( (cosB + cosA*cosC)/sinA/sinC );

        // second solution
        C2 = tau - C1;
        a2 = (a1 + pi).wrap_positive();
        b2 = (b1 + pi).wrap_positive();
    }

    /**
     * Solutions to the side-angle-side problem:
     *    "given two sides and the angle between them, compute
     *     the other side and the two other angles."
     * This problem has two solutions.
     *
     * @param[in] a first side
     * @param[in] C angle between sides a and b
     * @param[in] b seconds side
     */
    void SphericalTriangle::sas(const Angle& a, const Angle& C, const Angle& b)
    {
        // first assign everything obvious
        Angle
            aa  = wrap_positive(a),
            Cc  = wrap_positive(C),
            bb  = wrap_positive(b),
            pi  = Angle::pi(),
            tau = Angle::tau();

        a1 = aa;   a2 = tau - aa;
        C1 = Cc;   C2 = tau - Cc;
        b1 = bb;   b2 = tau - bb;

        double
            cosa = cos(aa),  sina = sin(aa),
            cosb = cos(bb),  sinb = sin(bb),
            cosC = cos(Cc),  cosc, sinc;

        // first solution
        c1 = acos2( (sina + cosa*cosb)*sinb*cosC, Cc);
        cosc = cos(c1);
        sinc = sin(c1);
        A1 = acos2( (cosa - cosb*cosc)/sinb/sinc, aa);
        B1 = acos2( (cosb - cosa*cosc)/sina/sinc, bb);

        // second solution
        c2 = tau - c1;
        A2 = (A1 + pi).wrap_positive();
        B2 = (B1 + pi).wrap_positive();

    }

    /**
     * Solutions to the side-side-angle problem:
     *    "given two sides and the angle opposite them, compute
     *     the other side and the two other angles."
     * This problem has zero or two solutions.
     *
     * @param[in] a first side
     * @param[in] b seconds side
     * @param[in] A angle opposite sides a and b
     */
    void SphericalTriangle::ssa(const Angle& a, const Angle& b, const Angle& A)
    {
                    // first assign everything obvious
        Angle
            aa  = wrap_positive(a),
            bb  = wrap_positive(b),
            Aa  = wrap_positive(A),
            pi  = Angle::pi(),
            tau = Angle::tau();

        a1 = aa;   a2 = tau - aa;
        b1 = bb;   b2 = tau - bb;
        A1 = Aa;   A2 = tau - Aa;

        double
            sina = sin(aa), sinb = sin(bb), sinA = sin(Aa),
            sinbsinA = sinb*sinA;

        // check condition
        if ( fabs(sinbsinA) <= fabs(sina) )
            throw csc::ex::Condition("SphericalTriangle: SSA problem has no solution.");

        // first solution
        B1 = Asin(sinbsinA/sina).wrap_positive();
        C1 = mal(Aa, B1, aa, bb);
        c1 = msl(aa, bb, Aa, B1);

        // second solution
        B2 = (pi - B1).wrap_positive();
        C2 = mal(Aa, B2, aa, bb);
        c2 = msl(aa, bb, Aa, B2);
    }


    // Constructor from three Cartesian vectors
    SphericalTriangle::SphericalTriangle(
        const csc::Vector3& v1,
        const csc::Vector3& v2,
        const csc::Vector3& v3)
    {
        // check equality of radii
        if (fabs(norm_2(v1)-norm_2(v2)) > 2.0*math::dx ||
            fabs(norm_2(v1)-norm_2(v3)) > 2.0*math::dx )
        {
            throw csc::ex::Condition("SphericalTriangle: vectors must have the same length.");
        }

        // find sides
        Angle
            s1 = Angle::radians(csc::angle_between(v1,v2)),
            s2 = Angle::radians(csc::angle_between(v2,v3)),
            s3 = Angle::radians(csc::angle_between(v3,v1));

        // TODO: make sure they're counterclockwise
        // TODO: make sure the triangle is not degenerate

        // complete triangle
        sss(s1, s2, s3);
    }

    // Constructor from three coordinates on a Sphere
    SphericalTriangle::SphericalTriangle(
        const sphere::Coordinate& c1,
        const sphere::Coordinate& c2,
        const sphere::Coordinate& c3)
    {
        auto
            v1 = toVector(c1),
            v2 = toVector(c2),
            v3 = toVector(c3);
        SphericalTriangle(v1,v2,v3);
    }


    // Constructor from three angles
    SphericalTriangle::SphericalTriangle(
        const Angle& R1,
        const Angle& R2,
        const Angle& R3,
        const std::string& triangle_type)
    {
        std::string lc(triangle_type);
        boost::algorithm::to_lower(lc);

        if (lc.compare("aaa") == 0 || lc.compare("angle-angle-angle") == 0)
            aaa(R1,R2,R3);

        else if (lc.compare("sss") == 0 || lc.compare("side-side-side") == 0)
            sss(R1,R2,R3);

        else if (lc.compare("aas") == 0 || lc.compare("angle-angle-side") == 0)
            aas(R1,R2,R3);

        else if (lc.compare("asa") == 0 || lc.compare("angle-side-angle") == 0)
            asa(R1,R2,R3);

        else if (lc.compare("sas") == 0 || lc.compare("side-angle-side") == 0)
            sas(R1,R2,R3);

        else if (lc.compare("ssa") == 0 || lc.compare("side-side-angle") == 0)
            ssa(R1,R2,R3);

        else
            throw csc::ex::Argument("SphericalTriangle: invalid problem type specified.");
    }

    // Area of the triangle
    double SphericalTriangle::area(
        unsigned int solution,
        double R)
    {
        if (solution == 0)
            return R*R*( (A1+B1+C1)-Angle::pi );
        else
            return R*R*( (A2+B2+C2)-Angle::pi );
    }

    // human-readable string representation
    std::string SphericalTriangle::toString() const
    {
        std::stringstream output;
        output
            << "\nSphericalTriangle:          \n"
            << "          .        .          \n"
            << "           .      .           \n"
            << "           \\    /            \n"
            << "         a2--\\  /--b2         \n"
            << "              C2              \n"
            << "              /\\              \n"
            << "             /C1\\             \n"
            << "            /    \\            \n"
            << "           /      \\           \n"
            << "          /-b1  a1-\\          \n"
            << "         /          \\         \n"
            << "   c2   /      c1    \\    c2  \n"
            << "  _|__ /A1_____|____B1\\ __|__ \n"
            << "      A2               B2     \n"
            << "     /                  \\     \n"
            << " b2-/                    \\-a2 \n"
            << "   .                      .   \n"
            << "  .                        .  \n"
            << "angle A1: " << A1 << ", angle A2: " << A2 << "\n"
            << "angle B1: " << B1 << ", angle B2: " << B2 << "\n"
            << "angle C1: " << C1 << ", angle C2: " << C2 << "\n"
            << "side  a1: " << a1 << ", side  a2: " << a2 << "\n"
            << "side  b1: " << b1 << ", side  b2: " << b2 << "\n"
            << "side  c1: " << c1 << ", side  c2: " << c2 << "\n\n";
        return output.str();
    }


    // Make  std::cout << {class}  work
    std::ostream& operator<<(std::ostream& target, const SphericalTriangle& st){
        target << st.toString();
        return target;
    }

    std::ostream& operator<<(std::ostream& target, const SphericalPolygon& sp){
        target << sp.toString();
        return target;
    }


} // namespace sphericaltrig



