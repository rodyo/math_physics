#ifndef _SPHERICALTRIG_H
#define _SPHERICALTRIG_H




// Spherical Geometry
namespace spherespace
{
    struct Coordinate;

    /**
     * std::pair-like thing for latitude/longitude on the sphere.
     *
     * Creates a named pair latitude/longitude on the sphere. Provides
     * a means to convert these coordinates from the spherical 2-space
     * to the Cartesian, intertial 3-space.
     */
    struct Coordinate
    {
        Angle
            latitude, longitude;

        // Primary use case
        Coordinate(const Angle& lat, const Angle& lon);

        // Secondary use case
        Coordinate(const csc::Vector3& V);


        // copy constructor
        Coordinate(const Coordinate& p)
            : latitude (p.latitude), longitude(p.longitude)
        {}

        // convert to vector/direciton
        csc::Vector3 toVector();

        // displayable string
        std::string toString() const;
    };

    /// It's also convenient to have a toVector() as a non-member
    csc::Vector3 toVector(const Coordinate& C);

}


// Ellipsoidal geometry
// NOTE: actually, the "ellipsoid" in common use in Geodesy, is
// a spheroid. See WGS84.
namespace ellipsoidspace
{
    struct Coordinate;

    /**
     * std::pair-like thing for geodetic latitude/longitude
     *
     * Creates a named pair latitude/longitude on the sphere. Provides
     * a means to convert these coordinates from the spherical 2-space
     * to the Cartesian, intertial 3-space.
     */
    class Coordinate : public spherespace::Coordinate
    {
    private:
        double
            flattening;
        Length
            equatiorialRadius;

    public:
        Angle
            latitude, longitude;
        Length
            altitude;

        // Primary use case
        Coordinate(const Angle& lat, const Angle& lon,
            double flattening = csc::SB_EarthData.flattening,
            const Length& equatiorialRadius = csc::SB_EarthData.eradius);

        // Secondary use case: compute lat/lon for given flattening
        Coordinate(const csc::Vector3& V,
            double flattening  = csc::SB_EarthData.flattening,
            const Length& equatiorialRadius = csc::SB_EarthData.eradius);


        // copy constructor
        Coordinate(const Coordinate& p);

        // convert to vector/direciton
        csc::Vector3 toVector();

        // displayable string
        std::string toString() const;
    };

    /// It's also convenient to have a toVector() as a non-member
    csc::Vector3 toVector(const Coordinate& C);

}

/// function to make  std::cout << {class}  work
std::ostream& operator<<(std::ostream& target, const ellipsoidspace::Coordinate& c);
std::ostream& operator<<(std::ostream& target, const spherespace::Coordinate&    c);


/**
 * All things spherical trig
 */
namespace shapes
{
    class  SphericalPolygon;
    class  SphericalTriangle;

    /** Spherical polygon
    // TODO: well, all of it.
    */
    class SphericalPolygon
    {
    private:

    protected:
        std::vector<Angle>
            angles,   // internal angles
            sides;    // angular sides

    public:
        ~SphericalPolygon(){}

        SphericalPolygon();


        /// human readable string representaiton
        std::string toString() const;
    };


    /**
     * Spherical triangle: a special kind of spherical polygon
     *
     * This implementation covers all possible ways to construct
     * a spherical triangle. Methodology can be found in
     * "Mission Geometry and Constellation Design", Wertz
     *
     */
    class SphericalTriangle
        : public sphericaltrig::SphericalPolygon
    {
    private:

        // Hemisphere function
        double H(const Angle& beta);
        // 4-quadrant arccosine function
        Angle acos2(double G, const Angle& beta );

        // Middle Angle Law
        Angle mal(
            const Angle& A, const Angle& B,
            const Angle& a, const Angle& b);

        // Middle Side Law
        Angle msl(
            const Angle& a, const Angle& b,
            const Angle& A, const Angle& B);

        // Angle-angle-angle problem
        // "given all three angles, compute the sides"
        void aaa(const Angle& A, const Angle& B, const Angle& C);

        // side-side-side problem:
        // "given all three sides, compute the angles"
        void sss(const Angle& a, const Angle& b, const Angle& c);

        // angle-angle-side problem:
        // "given two angles and one side opposite one of the angles,
        //  compute the two other sides and other angle."
        void aas(const Angle& A, const Angle& B, const Angle& a);


        // Angle-side-angle problem:
        // "given two angles and the side between their associated
        //  vertices, compute the two other sides and other angle."
        void asa(const Angle& A, const Angle& B, const Angle& c);

        // side-angle-side problem:
        // "given two sides and the angle between them, compute
        //  the other side and the two other angles."
        void sas(const Angle& a, const Angle& C, const Angle& b);

        // Side-side-angle problem:
        // "given two sides and the angle opposite them, compute
        //  the other side and the two other angles."
        void ssa(const Angle& a, const Angle& b, const Angle& A);

    public:

        /*           .        .
                      .      .
                       \    /
                    a2--\  /--b2
                         C2
                         /\
                        /C1\
                       /    \
                      /      \
                     /-b1  a1-\
                    /          \
              c2   /      c1    \    c2
             _|__ /A1_____|____B1\ __|__
                 A2               B2
                /                  \
            b2-/                    \-a2
              .                      .
             .                        .     */

        Angle
            a1,b1,c1, // "short" sides
            a2,b2,c2, // "long" sides

            A1,B1,C1, // "internal" angles
            A2,B2,C2; // "external" angles

        /// Constructor from three Cartesian vectors
        SphericalTriangle(
            const csc::Vector3& v1,
            const csc::Vector3& v2,
            const csc::Vector3& v3);

        /// Constructor from three coordinates on a Sphere
        SphericalTriangle(
            const sphere::Coordinate& c1,
            const sphere::Coordinate& c2,
            const sphere::Coordinate& c3);



        /// Constructor from three angles
        SphericalTriangle(
            const Angle& R1,
            const Angle& R2,
            const Angle& R3,
            const std::string& triangle_type);

        /// Area of the triangle
        double area(
            unsigned int solution = 0,
            double R = 1.0);

        /// human-readable string representation
        std::string toString() const;
    };

    /// Make  std::cout << {class}  work
    std::ostream& operator<<(std::ostream& target, const SphericalTriangle& st);
    std::ostream& operator<<(std::ostream& target, const SphericalPolygon&  sp);

} // namespace sphericaltrig



#endif
