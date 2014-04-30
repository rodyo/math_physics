/*
 * Rody Oldenhuis
 * cosine measurement systems
 * roldenhuis@cosine.nl
 *
 * Shapes.hh
 * Created: 05.10.2012 09:34:14 CEST
 */

#ifndef _SHAPES_HH
#define _SHAPES_HH

// FIXME: (Rody Oldenhuis) Should not be defined here, but....well, where?
//#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
//#define BOOST_MPL_LIMIT_LIST_SIZE 50
//#define BOOST_MPL_LIMIT_VECTOR_SIZE 50



/**
 * Basic shapes
 * Useful for doing intersections ect.
 */
namespace shapes
{
    // Prerequisites
    // ------------------------------------

    namespace _ignore {
        class Shape;           // Generic shape (NOTE: baseclass, don't use)
        template <typename T>
        class ShapeOps;        // helper class for overloading operators
    }

    class ReferenceFrame;      // Generic Reference frame, useful to connect different shapes

    // Basics
    class Nil;                 // abscence of shape -- used for ill-defined shapes
    class Point;               // A 0D point

    class Plane;               // your basic plane
    class Line;                // these could be defined as degernerate cases
                               // of conics, but other definitions are better
    class Polygon;             // polygon (wrapper of polygon2d)
    class Triangle;            // your basic flat-triangle
    class Rectangle;           // your basic flat-rectangle
    class Square;              // a special recangle

    // Conics
    namespace _ignore {
        class Conic;           // Generic quadric line = conic (NOTE: baseclass, don't use)
    }
    class Hyperbola;           //
    class Parabola;            //
    class Ellipse;             //
    class Circle;              // = special case of ellipse

    // Quadrics
    namespace _ignore {
        class Quadric;         // Generic quadric surface (NOTE: baseclass, don't use)
    }
    class Ellipsoid;           //
    class Spheroid;            // = special case of ellipsoid
    class Sphere;              // = special case of spheroid
    class EllipticCylinder;    //
    class Cylinder;            // = circular, special case of elliptic cylinder
    class HyperbolicCylinder;  //
    class EllipticCone;        //
    class Cone;                // = circular, special case of elliptic cone
    class EllipticParaboloid;  //
    class Paraboloid;          //
    class EllipticHyperboloid; //
    class Hyperboloid;         //

    class Pyramid;             // A pyramid is a collection of planes
    class Prism;               // Like a cirular cylinder, but with arbitrary cross secion

    class Composite;           // composite object

    // define a few aliases for consistency
    typedef Cylinder    CircularCylinder;
    typedef Cone        CircularCone;
    typedef Paraboloid  CircularParaboloid;
    typedef Hyperboloid CircularHyperboloid;

    // FIXME: (Rody Oldenhuis) See FIXME at the top
    // Container for ALL shapes (used in Composite)
    //using allShapes = boost::variant<
        //Nil,                   Point,         Plane,                Line,
        //Polygon,               Triangle,      Rectangle,            Square,
        //Hyperbola,             Parabola,      Ellipse,              Circle,
        //Ellipsoid,             Spheroid,      Sphere,
        //EllipticCylinder,      Cylinder,      HyperbolicCylinder,
        //EllipticCone,          Cone,
        //EllipticParaboloid,    Paraboloid,
        //EllipticHyperboloid,   Hyperboloid,
        //Pyramid,               Prism
    //>;
    using allShapes = boost::variant<
        Nil,                   Point,         Plane,                Line,
        Polygon,               Triangle,      Rectangle,            Square,
        Hyperbola,             Parabola,      Ellipse,              Circle,
        Ellipsoid,             Spheroid,      Sphere,
                               Cylinder,
                               Cone,
        Pyramid,               Prism
    >;


    /**
     * CoordinateSystem
     *
     *
     */
    class ReferenceFrame
    {
    private:

    protected:

        csc::Vector3
            i,j,k,    // directions of the x,y,z axes
            offset;   // offset of the origin
                      // (both w.r.t. the master frame)
        std::string
            _name;    // reference frame can have a human-readable name

    public:
        const std::string& name;

        ~ReferenceFrame(){}

        /// default frame: [1 0 0] [0 1 0] [0 0 1] (like the master)
        ReferenceFrame();

        /// define frame from scratch
        ReferenceFrame(
            csc::Vector3 i,
            csc::Vector3 j,
            csc::Vector3 k,
            csc::Vector3 offset = csc::Vector3());

        /// define frame as rotated/translated version of existing frame
        ReferenceFrame(
            const ReferenceFrame& F,
            const csc::RotationMatrix3x3& R,
            csc::Vector3 offset = csc::Vector3());

        /**
         * Convert coordinates given in this frame to their
         * counterparts in the master frame
         */
        Point        inMaster(Point& P);
        csc::Vector3 inMaster(csc::Vector3& V);

        /**
         * Convert coordinates given in the master frame to
         * their counterparts in this frame
         */
        Point        inThis  (Point& P);
        csc::Vector3 inThis  (csc::Vector3& V);

        /// Human-readable string representation
        std::string toString() const;
    };



    // Basic shapes
    // ------------------------------------

    /**
     * Generic shape
     *
     *
     */
    namespace _ignore
    {
        class Shape
        {
        protected:

            /// human-readable shape name
            std::string
                _name;

            /**
             * Most shapes can become degenerate or ill-defined.
             * Also, the user might intend to use the "inverse" of the shape
             * (i.e., everything that's NOT the shape)
             */
            bool
                _degenerate,
                _ill_defined,
                _inverse;

            /**
             * coordinates of a point on the shape for some combination
             * of the shape's parameters
             */
            csc::Vector3
                _current_coordinate;

            /**
             * compute coordinates for some combination of the shape's
             * parameters.
             *
             * To be Used in isMember(Vector3) below.
             */
            virtual void evaluate() = 0;

            /**
             * Check shape for degeneracy and validity.
             *
             * To be used in shape constructors.
             */
            virtual void checkConsistency() = 0;

        public:
            virtual ~Shape(){}

            const std::string
                &name;
            const bool
                &degenerate,
                &ill_defined,
                &inverse;

            Shape()
                : _name        ("(no name)")
                , _degenerate  (false)
                , _ill_defined (false)
                , _inverse     (false)
                , _current_coordinate {}

                , name         (_name)
                , degenerate   (_degenerate)
                , ill_defined  (_ill_defined)
                , inverse      (_inverse)
            {}

            /// test if some point (or other shape) lies on the shape
            virtual bool isMember(const csc::Point& coord) = 0;
            bool isMember(const Shape& S) {
                throw csc::ex::NotImplemented("Membership tests of shapes in shapes has not yet been implemented.");
            }

            /// area, circumference, ect.
            virtual Area   area()          const = 0;  // alias
            virtual Length circumference() const = 0;  Length perimeter(){ return circumference(); }
            virtual Volume volume()        const = 0;

            /// setter for the name
            void setName(const char* nm) { _name = std::string(nm); }
            void setName(const std::string& nm) { _name = nm; }

            /// Human-readable string representation
            virtual std::string toString() const;

        };

        /**
         * ShapeOps
         *
         * Helper class to easily define overloaded operators for
         * all shapes.
         *
         */

        // Use virtual inheritance to avoid diamond problem later on. See
        // http://stackoverflow.com/questions/13031327/
        template <typename T>
        class ShapeOps
            : virtual public Shape
        {
        public:
            virtual ~ShapeOps(){}
            ShapeOps() : Shape() {}

            /// use inverse of the shape
            T& invert    () { _inverse = true; return static_cast<T&>(*this);}
            T& operator~ () { return invert();}
        };
    }


    /**
     * Nil: the abscence of shape, the empty set, vacuum, nothing, nada
     *
     * NOTE: the complement of Nil is a degenerate form, meaning
     * "all of space".
     *
     */
    class Nil
        : public _ignore::ShapeOps<Nil>
    {
    private:
        typedef _ignore::ShapeOps<Nil> bc;

        void checkConsistency(){}
        void evaluate(){}

    public:
        ~Nil(){}
        Nil();

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Point:
     *
     *
     */
    class Point
        : public _ignore::ShapeOps<Point>
    {
    private:
        typedef _ignore::ShapeOps<Point> bc;

        void evaluate();
        void checkConsistency();

    public:
        ~Point();

        csc::Vector3
            coordinates;

        Point();
        Point(const csc::Vector3& P);
        Point(const Point& P);

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        operator Nil();
        operator csc::Vector3();

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Plane
     *
     *
     */
    class Plane
        : public _ignore::ShapeOps<Plane>
    {
    private:
        typedef _ignore::ShapeOps<Plane> bc;
        friend class Line;
        friend class Sphere;
        friend class Ellipsoid;

    protected:

        csc::Vector3
            normal_vector,  // vector normal to the plane
            offset,         // offset from origin
            u,v;            // plane P is also defined by
                            // P = offset + su + tv

        double
            a, b, c, d;     // scalar equation

        void
            evaluate(),
            checkConsistency();

    public:
        ~Plane();

        /// normal vector defines plane
        Plane(const csc::Vector3& N);

        /// two vectors define plane
        Plane(const csc::Vector3& v1,
              const csc::Vector3& v2);

        /// three points define plane
        Plane(const csc::Point& p1,
              const csc::Point& p2,
              const csc::Point& p3);

        /// Line and point define plane
        Plane(const Line& L,
              const csc::Point& p);

        /// scalar equation of plane (ax + by + cz = d)
        Plane(double  a, double  b, double  c, double  d);
        Plane(double& a, double& b, double& c, double& d);

        /// get coordinates for some value of s, t
        Point evaluate(double s, double t) const;
        /// get missing coordinate for some combination of x,y | x,z | y,z
        void evaluate(double& x, double  y, double  z) const;
        void evaluate(double  x, double& y, double  z) const;
        void evaluate(double  x, double  y, double& z) const;

        /**
         * Membership test:
         * - does the given point lie in the plane?
         *
         * @param[in]: P the point to test
         * @return true if point lies in the plane.
         */
        bool isMember(const csc::Point& P);

        /**
         * Membership test:
         * - does the given point lie in the plane?
         *
         * @param[in]: P the point to test
         * @return true if point lies in the plane.
         */
        bool isMember(const Point& P);

        /**
         * Membership test:
         * - does the given line lie in the plane?
         *
         * @param[in]: L the line to test
         * @return true if line indeed lies in the plane.
         */
        bool isMember(const Line& L);

        /// Area, circumference
        Area   area()          const;
        Length circumference() const;
        Volume volume()        const;

        /// Intersections
        using line_intersection      = boost::variant<Nil, Point, Line>;
        using plane_intersection     = boost::variant<Nil, Line , Plane>;
        using sphere_intersection    = boost::variant<Nil, Point, Circle>;
        using ellipsoid_intersection = boost::variant<Nil, Point, Circle, Ellipse>;

        auto intersect(const Line&      L) -> line_intersection;
        auto intersect(const Plane&     P) -> plane_intersection;
        auto intersect(const Sphere&    S) -> sphere_intersection;
        auto intersect(const Ellipsoid& E) -> ellipsoid_intersection;

        /// Distances
        Length distance(const Line&   L);
        Length distance(const Plane&  P);
        Length distance(const Sphere& S);

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Line:
     *
     *
     */
    class Line
        : public _ignore::ShapeOps<Line>
    {
    private:
        typedef _ignore::ShapeOps<Line> bc;

    protected:

        csc::Vector3
            direction,      // direction of the line
            support;        // support vector so that the line L is defined by
                            // L = support + t*direction
        double
            t,              // value of the parameter
            _t_min, _t_max; // limits on the parameter

        void evaluate();
        void checkConsistency();

        friend class Plane;

    public:
        ~Line(){}

        const double
            &t_min, &t_max;

        /// line parallel to X-axis
        Line();
        /// line defined by two vectors
        Line(const csc::Vector3& v1, const csc::Vector3& v2);
        /// line defined by two Points
        Line(const Point& v1, const Point& v2);
        /// Cartesian equation: (x-x0)/c1 = (y-y0)/c2 = (z-z0)/c3
        Line(
            double x0, double c1,
            double y0, double c2,
            double z0, double c3);
        /// line through origin, defined by direction cosines
        Line(double cx, double cy, double cz);
        /// line through origin, defined by Angles w.r.t. coordinate axes
        Line(const Angle& ax, const Angle& ay, const Angle& az);

        /// check if a point lies on the line
        bool isMember(const Point& pt);
        bool isMember(const csc::Point& pt);

        /// length of the line segment
        Length length() const;

        /// set lower/upper limits
        void setLimits(double tmin, double tmax);
        void setLimits(const csc::Vector3& tmin, const csc::Vector3& tmax);
        void setLimits(const Point& tmin, const Point& tmax);

        /// get coordinates of the line at some value of the parameter t
        Point operator[] (double tt);
        Point at(double tt);

        /// get value for the parameter for some coordinates on the line
        double operator[](const csc::Point& p);
        double operator[](const Point& p);
        double at(const Point& p);
        double at(const csc::Point& p);

        /// Area, circumference
        Area   area()          const;
        Length circumference() const;
        Volume volume()        const;

        /// allow type casting to simpler types for degenerate lines
        operator Nil();
        operator Point();

        // Intersections
        // -------------
        typedef boost::variant<Nil, Point, Line>         line_intersection;
        typedef boost::variant<Nil, Point, Line>         plane_intersection;
        typedef boost::variant<Nil, std::vector<Point>>  circle_intersection;
        typedef boost::variant<Nil, std::vector<Point>>  ellipse_intersection;

        auto intersect(const    Line& L) -> line_intersection;
        auto intersect(const   Plane& P) -> plane_intersection;
        auto intersect(const  Circle& C) -> circle_intersection;
        auto intersect(const Ellipse& E) -> ellipse_intersection;

        // Distances
        // ---------

        Length distance(const  Line& L);
        Length distance(const Plane& P);


        // Human-readable string representation
        std::string toString() const;
    };


    /**
     * Polygon
     * NOTE: wrapper around polygon defined in polygon2d.*.
     */
    class Polygon
        : public csc::Polygon2D
        , virtual public _ignore::ShapeOps<Polygon>
    {
    private:
        typedef _ignore::ShapeOps<Polygon> bc;

    protected:

        csc::Vector3
            u,v,
            support;

        double
            t;

        unsigned int
            _num_sides;           // how many sides/vertices?
        bool
            _is_convex,           // all angles <= 180 and !_is_star?
            _is_star,             // at least two connecting line segments intersect
            _is_counterclockwise; // are they ordered counter-clockwise?

        std::vector<Angle>
            getAngles();

        void
            sortCounterClockwise(),
            evaluate(),
            checkConsistency();

        bool
            isSorted();

    public:
        virtual ~Polygon();

        const bool
            &is_convex,
            &is_concave,
            &is_star,
            &self_intersects;
        unsigned int
            &num_sides,
            &num_vertices;

        Polygon();
        Polygon(const Polygon& P);
        Polygon(const std::vector<Point>& pts);
        Polygon(const std::vector<csc::Vector3>& pts);

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Find centroid
        Point centroid() const;

        /// Get angle at given vertex
        Angle angle(unsigned int i) const;

        /// Get basic info on the polygon
        bool
            isConvex() const,
            isConcave() const,
            isStar() const,
            selfIntersects() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Triangle
     *
     *
     */
    class Triangle
        : public Polygon
        , virtual public _ignore::ShapeOps<Triangle>
    {
    private:
        typedef _ignore::ShapeOps<Triangle> Bc;

        csc::Vector3
            a,b,c;

        void evaluate();
        void checkConsistency();

    public:
        ~Triangle(){}

        Triangle(const Point& a, const Point& b, const Point& c);
        Triangle(const csc::Vector3& a, const csc::Vector3& b, const csc::Vector3& c);

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Return angles in the triangle
        Angle ab();
        Angle bc();
        Angle ca();

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Rectangle
     *
     *
     */
    class Rectangle
        : public Polygon
        , virtual public _ignore::ShapeOps<Rectangle>
    {
    private:
        typedef _ignore::ShapeOps<Rectangle> bc;

        void evaluate(){}
        void checkConsistency(){}

    public:
        virtual ~Rectangle();
        Rectangle();

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Square
     *
     *
     */
    class Square
        : public Rectangle
        , virtual public _ignore::ShapeOps<Square>
    {
    private:
        typedef _ignore::ShapeOps<Square> bc;

        void evaluate(){}
        void checkConsistency(){}

    public:
        ~Square();
        Square();

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    // Conics
    // -------------------------------------

    /**
     * Generic conic (=1D quadric)
     *
     *
     */
    namespace _ignore
    {
        class Conic
            : public _ignore::ShapeOps<Conic>
        {
        private:
            typedef _ignore::ShapeOps<Conic> bc;

        protected:
            double     // conics have equations of this form:
                a, b;  //     ± x²/a² ± y²/b² = ±1 v 0

        public:
            virtual ~Conic(){}
            Conic();

            Area   area();
            Volume volume() const;
        };
    }


    /**
     * Ellipse
     *
     *
     */
    class Ellipse
        : public _ignore::Conic
        , virtual public _ignore::ShapeOps<Ellipse>
    {
    private:
        typedef _ignore::ShapeOps<Ellipse> bc;

        Length
            _a, _b;        // semi-major, semi-minor axes, respectfully
        double
            _e;            // eccentricity

    protected:
                           // needed for segments of the ellipse
        Angle              // these define min/max angles in the
            _t_min, _t_max,// parametric equation  C+
            t;             // value of the parameter

        csc::Vector3
            u,v, C;        // Vectors pointing from center C to
                           // where semi-major axis crosses the ellipse (u) and
                           // where semi-minor axis crosses the ellipse (v)
        void
            evaluate(),
            checkConsistency();

    public:
        virtual ~Ellipse();

        const Length
            &a, &b;
        const double
            &e;
        const Angle
             &t_min, &t_max;

        Ellipse();
        Ellipse(const Length& a, const Length& b);

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /**
         * Evaluate the ellipse for some value of the parameter
         *
         */
        double at(double tt);
        double operator[](double(tt));

        double at(const Point& P);
        double operator[](const Point& P);

        /// Allow conversion to vector (since it's basically the same anyway)
        operator csc::Vector3();

        /// allow type casting to simpler types for degenerate Ellipses
        operator Nil();
        operator Point();
        operator Line();
        operator Circle();

        operator Hyperbola();
        operator Parabola();

        // intersections
        typedef boost::variant<Nil, Point, std::vector<Point>, Ellipse>  plane_intersection;
        typedef boost::variant<Nil, Point, std::vector<Point>>          line_intersection;

        auto intersect(const Plane& P) const -> plane_intersection;
        auto intersect(const Line&  L) const -> line_intersection;


        // distances
        Length distance(const Point&   P) const;
        Length distance(const Plane&   P) const;
        Length distance(const Ellipse& E) const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Circle: a special Ellipse
     *
     *
     */
    class Circle
        : public Ellipse
        , virtual public _ignore::ShapeOps<Circle>
    {
    private:
        typedef _ignore::ShapeOps<Circle> bc;

        Length
            _radius;

        void
            evaluate(),
            checkConsistency();

    public:
        ~Circle(){}

        const Length
            &radius;

        Circle();
        Circle(const Length& L, const csc::Vector3& C, const csc::Vector3& u, const csc::Vector3& v);

        /// Area, circumference
        bool   isMember (const csc::Point& coord);
        Area   area()          const;
        Length circumference() const;
        Volume volume()        const;

        /// allow type casting for degenerate circles
        operator Nil();
        operator Point();

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Parabola
     *
     *
     */
    class Parabola
        : public _ignore::Conic
        , virtual public _ignore::ShapeOps<Parabola>
    {
    private:
        typedef _ignore::ShapeOps<Parabola> bc;

        void
            evaluate(),
            checkConsistency();

    public:
        ~Parabola(){}
        Parabola();

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area()  const;
        Length circumference() const;
        Volume volume() const;

        ///
        operator Nil();
        operator Point();
        operator Ellipse();
        operator Line();

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Hyperbola
     *
     *
     */
    class Hyperbola
        : public _ignore::Conic
        , virtual public _ignore::ShapeOps<Hyperbola>
    {
    private:
        typedef _ignore::ShapeOps<Hyperbola> bc;

        void
            evaluate(),
            checkConsistency();

    public:
        ~Hyperbola(){}
        Hyperbola();

        /// Area, circumference
        bool isMember(const csc::Point& coord);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        ///
        operator Nil();
        operator Point();
        operator Ellipse();
        operator Line();

        /// Human-readable string representation
        std::string toString() const;
    };


    // Quadrics
    // -------------------------------------

    /**
     * Generic quadric surface
     *
     *
     */
    namespace _ignore
    {
        class Quadric
            : public _ignore::ShapeOps<Quadric>
        {
        private:
            typedef _ignore::ShapeOps<Quadric> bc;

        protected:
            double        // quadrics have equations of this form:
                a, b, c;  //     ± x²/a² ± y²/b² ± z²/c² = ±1 v 0

        public:
            virtual ~Quadric(){}
        };
    }


    /**
     * Ellipsoid
     *
     *
     */
    class Ellipsoid
        : public _ignore::Quadric
        , virtual public _ignore::ShapeOps<Ellipsoid>
    {
    private:
        typedef _ignore::ShapeOps<Ellipsoid> bc;

        void checkConsistency(){}
        void evaluate(){}

    public:
        virtual ~Ellipsoid(){}
        Ellipsoid();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Spheroid
     *
     *
     */
    class Spheroid
        : public Ellipsoid
        , virtual public _ignore::ShapeOps<Spheroid>
    {
    private:
        typedef _ignore::ShapeOps<Spheroid> bc;

        void checkConsistency(){}
        void evaluate(){}

    public:
        virtual ~Spheroid(){}
        Spheroid();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;


        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Sphere
     *
     *
     */
    class Sphere
        : public Spheroid
        , virtual public _ignore::ShapeOps<Sphere>
    {
    private:
        typedef _ignore::ShapeOps<Sphere> bc;

        Length
            radius;
        csc::Vector3
            center;

        void checkConsistency(){}
        void evaluate(){}

    public:
        ~Sphere(){}

        /// construct unit sphere
        Sphere();
        /// construct arbitrary other sphere
        Sphere(Length radius);
        Sphere(csc::Vector3 center);
        Sphere(Length radius, csc::Vector3 center);

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area   area()          const;
        Length circumference() const;
        Volume volume()        const;

        // intersections
        // ------------------------------

        using line_intersection   = boost::variant<Nil, Point, std::vector<Point>>;
        using plane_intersection  = boost::variant<Nil, Point, Circle>;
        using sphere_intersection = boost::variant<Nil, Point, Circle, Sphere>;

        /**
         *
         *
         */
        auto intersect(const Line&   L) const -> line_intersection;
        auto intersect(const Plane&  P) const -> plane_intersection;
        auto intersect(const Sphere& P) const -> sphere_intersection;

        // distances
        // ------------------------------

        /**
         * Plane-sphere distance
         *
         * Distance is always positive, except when the sphere intersects with the
         * plane. In that case, the negative value is the distance between the plane
         * and the centre of the sphere.
         *
         * @param[in] P plane
         * @return the distance (Length) between the sphere and the plane.
         */
        Length distance(const Plane& P)  const;


        Length distance(const Sphere& S) const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Hyperbolic cylinder
     *
     *
     */
    class HyperbolicCylinder
        : public _ignore::Quadric
        , virtual public _ignore::ShapeOps<HyperbolicCylinder>
    {
    private:
        typedef _ignore::ShapeOps<HyperbolicCylinder> bc;

    protected:

        void checkConsistency(){}
        void evaluate(){}

    public:
        ~HyperbolicCylinder();
        HyperbolicCylinder();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Elliptic cylinder
     *
     *
     */
    class EllipticCylinder
        : public _ignore::Quadric
        , virtual public _ignore::ShapeOps<EllipticCylinder>
    {
    private:
        typedef _ignore::ShapeOps<EllipticCylinder> bc;

    protected:

        void checkConsistency(){}
        void evaluate(){}

    public:
        virtual ~EllipticCylinder();
        EllipticCylinder();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Cylinder
     *
     *
     */
    class Cylinder
        : public EllipticCylinder
        , virtual public _ignore::ShapeOps<Cylinder>
    {
    private:
        typedef _ignore::ShapeOps<Cylinder> bc;

        Length
            radius,
            height;
        csc::Vector3
            offset,
            orientation;

        void evaluate();
        void checkConsistency();

    public:
        Cylinder();
        ~Cylinder();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Elliptic cone
     *
     *
     */
    class EllipticCone
        : public _ignore::Quadric
        , virtual public _ignore::ShapeOps<EllipticCone>
    {
    private:
        typedef _ignore::ShapeOps<EllipticCone> bc;

    protected:

        void evaluate(){}
        void checkConsistency(){}

    public:
        virtual ~EllipticCone(){};
        EllipticCone();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Cone
     *
     *
     */
    class Cone
        : public EllipticCone
        , public _ignore::ShapeOps<Cone>
    {
    private:
        typedef _ignore::ShapeOps<Cone> bc;

        void checkConsistency(){}
        void evaluate(){}

    protected:

    public:
        Cone();
        ~Cone();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Elliptic paraboloid
     *
     *
     */
    class EllipticParaboloid
        : public _ignore::Quadric
        , virtual public _ignore::ShapeOps<EllipticParaboloid>
    {
    private:
        typedef _ignore::ShapeOps<EllipticParaboloid> bc;

    protected:

        void checkConsistency(){}
        void evaluate(){}

    public:
        virtual ~EllipticParaboloid();
        EllipticParaboloid();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Paraboloid
     *
     *
     */
    class Paraboloid
        : public EllipticParaboloid
        , virtual public _ignore::ShapeOps<Paraboloid>
    {
    private:
        typedef _ignore::ShapeOps<Paraboloid> bc;

    protected:

        void checkConsistency(){}
        void evaluate(){}

    public:
        ~Paraboloid();
        Paraboloid();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Elliptic hyperboloid
     *
     *
     */
    class EllipticHyperboloid
        : public _ignore::Quadric
        , virtual public _ignore::ShapeOps<EllipticHyperboloid>
    {
    private:
        typedef _ignore::ShapeOps<EllipticHyperboloid> bc;

    protected:

        void checkConsistency(){}
        void evaluate(){}

    public:
        virtual ~EllipticHyperboloid();
        EllipticHyperboloid();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Hyperboloid
     *
     *
     */
    class Hyperboloid
        : public EllipticHyperboloid
        , virtual public _ignore::ShapeOps<Hyperboloid>
    {
    private:
        typedef _ignore::ShapeOps<Hyperboloid> bc;

    protected:

        void checkConsistency(){}
        void evaluate(){}

    public:
        ~Hyperboloid();
        Hyperboloid();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    // Compound shapes
    // --------------

    /**
     * @brief Pyramid of arbitrary cross section. Both finite and infinite.
     */

    /**
     * Pyramid
     *
     * A Pyramid is a shape of which all outer surfaces are triangular and
     * converge to a single point at the top. The base (or cross-section) is
     * a Polygon.
     *
     * Instantiations of Pyramid implemented here have two interpretations:
     *
     *  1. The "ordinary pyramid": the finite volume contained by the
     *     polygon base area, and all of the pyramid's sides.
     *
     *  2. The "infinite pyramid": the infinite volume that has the
     *     ordinary pyramid at its top. It is the natural extension of the
     *     pyramid ad infinitum. The polygon base area is placed infinitely
     *     far away, with a corresponding area increase.
     *
     * Both interpretations have their own specific methods, that are
     * available with every instance.
     *
     * @see Cone @see Prism
     */
    class Pyramid
        : public _ignore::ShapeOps<Pyramid>
    {
    private:
        typedef _ignore::ShapeOps<Pyramid> bc;

        csc::Vector3
            tip,
            direction;
        Polygon
            crossSection;

        std::vector<Plane>
            planes;

        void
            checkConsistency(),
            evaluate();

    public:
        ~Pyramid(){}

        Pyramid(const Point& p, const Polygon& P, const csc::Vector3& D);

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    /**
     * Prism
     *
     *
     */
    class Prism
        : public _ignore::ShapeOps<Prism>
    {
    private:
        typedef _ignore::ShapeOps<Prism> bc;

        void checkConsistency(){}
        void evaluate(){}

    public:
        ~Prism();
        Prism();

        /// Area, circumference
        bool isMember(const csc::Point& pt);
        Area area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };


    // Composite objects
    // ---------------------------------

    // This is to do Constructive solid modelling in RSS. See
    // http://en.wikipedia.org/wiki/Constructive_solid_geometry

    /**
     * @brief Shape formed by combining other shapes
     */

    /**
     * Composite shape
     *
     * A composite shape ("composite") is a shape made by combining an arbitrary
     * number of base shapes defined above. Combining can be accomplished by any
     * operation valid on sets:
     *
     *  - intersect
     *  - union
     *  - subtract
     *  - complement
     *
     * Any of the following operations can be performed on composites:
     *
     *  - area, volume
     *  - minimum distance from another shape to composite
     *  -
     *
     */
    class Composite
        : public _ignore::ShapeOps<Composite>
    {
    private:
        typedef _ignore::ShapeOps<Composite> bc;

    protected:

        std::vector<allShapes>
            shapeSet;

        void
            evaluate(),
            checkConsistency();

    public:
        ~Composite();
        Composite();

        /// Constructor operations: union
        Composite
            Union     (const _ignore::Shape& S, ...),
            operator+ (const _ignore::Shape& S),
            operator| (const _ignore::Shape& S);
        Composite&
            operator+= (const _ignore::Shape& S),
            operator|= (const _ignore::Shape& S);

        /// intersect
        Composite
            Intersect (const _ignore::Shape& S, ...),
            operator/ (const _ignore::Shape& S),
            operator& (const _ignore::Shape& S);
        Composite&
            operator/= (const _ignore::Shape& S),
            operator&= (const _ignore::Shape& S);

        /// complement
        Composite
            Complement(const _ignore::Shape& S, ...),
            operator- (const _ignore::Shape& S),
            operator% (const _ignore::Shape& S);
        Composite&
            operator-= (const _ignore::Shape& S);

        /// additional useful operators
        Composite
            operator^ (const _ignore::Shape& S);  // (A ⋃ B) \ (A ⋂ B)


        /// The basics
        bool   isMember(const csc::Point& pt);
        Area   area() const;
        Length circumference() const;
        Volume volume() const;

        /// Human-readable string representation
        std::string toString() const;
    };



    // Forward output of toString() method when streaming
    // to an ostream (like std::cout)

    /// functions to make  std::cout << {shape::class}  etc. work
    std::ostream& operator<<(std::ostream& target, const _ignore::Shape& S);

}




#endif



