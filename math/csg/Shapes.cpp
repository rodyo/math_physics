/*
 * Rody Oldenhuis
 * cosine measurement systems
 * roldenhuis@cosine.nl
 *
 * Shapes.cc
 * Created: 05.10.2012 09:33:47 CEST
 */





namespace shapes
{

    // Reference frame
    // ----------------------------




    // Nil
    // ========================================
    Nil::Nil()
        : bc()
    {
        _degenerate  = true;
        _ill_defined = true;
    }

    bool Nil::isMember(const csc::Point& coord) { return false; }
    Area Nil::area()            const { return Area::NaN; }
    Length Nil::circumference() const { return Length::NaN; }
    Volume Nil::volume()        const { return Volume::NaN; }




    // Point
    // ========================================

    void Point::checkConsistency() {
    }

    Point::Point()
        : bc()
        , coordinates{}
    {
        checkConsistency();
    }

    Point::Point(const csc::Vector3& P)
        : bc()
        , coordinates (P)
    {
        checkConsistency();
    }

    Point::Point(const Point& P)
        : bc()
        , coordinates (P.coordinates) {
        _name = P.name;
    }

    Area Point::area()            const { return {}; }
    Length Point::circumference() const { return {}; }
    Volume Point::volume()        const { return {}; }
    bool Point::isMember(const csc::Point& coord) {
        return ((coord-coordinates).isZeroLength());
    }

    Point::operator csc::Vector3(){
        if (ill_defined || degenerate)
            throw csc::ex::Condition("Typecast from Point to Vector3 is only possible for well-defined Points.");
        return coordinates;
    }

    Point::operator Nil() {
        if (ill_defined || degenerate)
            return {};
        else
            throw csc::ex::Condition("Typecast from Point to Nil is only possible for ill-defined or degenerate Points.");
    }




    // Plane
    // ========================================

    // normal vector defines plane through origin
    Plane::Plane(const csc::Vector3& N)
        : bc()
        , normal_vector (N)
        , offset        {}
        , a (normal_vector[0])
        , b (normal_vector[1])
        , c (normal_vector[2])
        , d (offset.length() )
    {}

    // two vectors define plane through origin
    Plane::Plane(const csc::Vector3& v1, const csc::Vector3& v2)
        : bc()
        , normal_vector ( csc::normalize(csc::vector_product(v1,v2)) )
        , offset        {}
        , a (normal_vector[0])
        , b (normal_vector[1])
        , c (normal_vector[2])
        , d (offset.length() )
    {}

    // three points define plane
    Plane::Plane(const csc::Point& p1, const csc::Point& p2, const csc::Point& p3)
        : bc()
        , normal_vector ( csc::normalize(csc::vector_product(p2-p1,p3-p2))    )
        , offset        ( normal_vector*csc::scalar_product(normal_vector,p1) )
        , a (normal_vector[0])
        , b (normal_vector[1])
        , c (normal_vector[2])
        , d (offset.length() )
    {}

    // Line and point define plane
    Plane::Plane(const Line& L, const csc::Point& p)
    {
        // FIXME: (Rody Oldenhuis) this is a stub
        normal_vector =
        u =
        v =
        offset = {};
    }

    // scalar equation of plane
    //  ax + by + cz = d
    Plane::Plane(double a, double b, double c, double d)
        : bc()
        , normal_vector (csc::Vector3(a,b,c))
        , a(a) , b(b)
        , c(c) , d(d)
    {
        double L = normal_vector.length();
        normal_vector /= L;
        offset = normal_vector*d/L;
    }

    // membership tests
    bool Plane::isMember(const csc::Point& coord) {
        double angle = csc::angle_between(offset-coord, normal_vector);
        return (fabs(angle-math::pio2) < math::epsilon);
    }

    bool Plane::isMember(const Line& L) {
        // FIXME: (Rody Oldenhuis)
        return false;
    }

    // get coordinates for some value of s, t
    Point Plane::evaluate(double s, double t) const {
        return {offset + u*s + v*t};
    }

    // get missing coordinate for some combination of x,y | x,z | y,z
    void Plane::evaluate(double& x, double y, double z) const {
        x = (offset.length() - normal_vector[1]*y - normal_vector[2]*z)/normal_vector[0];
    }
    void Plane::evaluate(double x, double& y, double z) const {
        y = (offset.length() - normal_vector[0]*x - normal_vector[2]*z)/normal_vector[1];
    }
    void Plane::evaluate(double x, double y, double& z) const {
        z = (offset.length() - normal_vector[0]*x - normal_vector[1]*y)/normal_vector[2];
    }

    // Human-readable string representation
    std::string Plane::toString() const
    {
        std::stringstream output;
        output
            << "Plane: \"" << _name << "\"\n"
            << "(ax + by + cz = d:)\n"
            << "a: " << normal_vector[0]
            << "b: " << normal_vector[1]
            << "c: " << normal_vector[2]
            << "d: " << offset.length();
        return output.str();
    }

    // Intersections
    // ------------------------------------

    auto Plane::intersect(const Line& L) -> line_intersection
    {
        // Line lies inside plane: intersection is the line
        if (isMember(L))
            return L;

        // line is parallel to plane: Nil
        if (fabs(csc::angle_between(normal_vector, L.direction)-math::pio2) < math::epsilon)
            return Nil();

        // TODO: point


        return Nil();
    }

    auto Plane::intersect(const Plane& P) -> plane_intersection
    {
        // plane is parallel. fully contained?
        if (csc::angle_between(normal_vector, P.normal_vector) < math::epsilon) {
            if (offset==P.offset)
                return *this; // yes
            else
                return Nil(); // no
        }

        // Offset can be found by solving the equation
        //     f1 + c*f1hat = f2 + d*f2hat = offset
        csc::Vector3
            f1    = offset,
            f2    = P.offset,
            // shorthand
            phi   = f2-f1,
            // vector perpendicular to f1
            f1hat = csc::vector_product(csc::vector_product(f1,f2),f1),
            // vector perpendicular to f2
            f2hat = csc::vector_product(csc::vector_product(f2,f1),f2);

        // there are 3 equations, 2 unknowns. Compute all solutions
        // and take the average

        double d = (
            (phi[1]*f1hat[0]-phi[0]*f1hat[1]) / (f2hat[1]*f1hat[0] - f2hat[0]*f1hat[1]) +
            (phi[2]*f1hat[1]-phi[1]*f1hat[2]) / (f2hat[2]*f1hat[1] - f2hat[1]*f1hat[2]) +
            (phi[0]*f1hat[2]-phi[2]*f1hat[0]) / (f2hat[0]*f1hat[2] - f2hat[2]*f1hat[0]) ) / 3.0;

        double c = (
            (phi[1]*f2hat[0]-phi[0]*f2hat[1]) / (f1hat[1]*f2hat[0] - f1hat[0]*f2hat[1]) +
            (phi[2]*f2hat[1]-phi[1]*f2hat[2]) / (f1hat[2]*f2hat[1] - f1hat[1]*f2hat[2]) +
            (phi[0]*f2hat[2]-phi[2]*f2hat[0]) / (f1hat[0]*f2hat[2] - f1hat[2]*f2hat[0]) ) / 3.0;

        // finally
        return Line(
            ( (f1+f1hat*c) + (f2+f2hat*d) )/2.0,
            // direction vector is trivial:
            csc::vector_product(normal_vector, P.normal_vector));
    }

    auto Plane::intersect(const Sphere& S) -> sphere_intersection {
        return S.intersect(*this);
    }

    // distances
    // ------------------------------------

    Length Plane::distance(const Line& L){
        // TODO: (Rody Oldenhuis)
        return Length();
    }

    Length Plane::distance(const Plane& P){
        // TODO: (Rody Oldenhuis)
        return Length();
    }

    Length Plane::distance(const Sphere& S){ return S.distance(*this); }


    // Line
    // ========================================

    void Line::checkConsistency(){

        // TODO: (Rody Oldenhuis) should this counts as _ill_defined?
        if (_t_max < _t_min)
            std::swap(_t_max, _t_min);

        _degenerate = _t_max==_t_min || direction.isZeroLength();
    }

    // compute current coordinate
    void Line::evaluate(){
        _current_coordinate = support + direction*t;
    }

    Line::Line()
        : bc()
        , direction (csc::Vector3(1,0,0))
        , support   (csc::Vector3())
        , t         (0)
        , _t_min    (-math::inf)
        , _t_max    (+math::inf)
        , t_min     (_t_min)
        , t_max     (_t_max) {
        _name = "Line";
        checkConsistency();
    }

    // line defined by two vectors
    Line::Line(const csc::Vector3& v1, const csc::Vector3& v2)
        : bc()
        , direction (v2-v1)
        , support   (v1)
        , t         (0)
        , _t_min    (-math::inf)
        , _t_max    (+math::inf)
        , t_min     (_t_min)
        , t_max     (_t_max) {
        _name = "Line";
        checkConsistency();
    }

    // line defined by two Points
    Line::Line(const Point& v1, const Point& v2)
        : bc()
        , t_min (_t_min)
        , t_max (_t_max) {
        Line(v1.coordinates, v2.coordinates);
    }

    // scalar equation of the line
    //  (x-x0)/c1 = (y-y0)/c2 = (z-z0)/c3
    Line::Line(
        double x0, double c1,
        double y0, double c2,
        double z0, double c3)
        : bc()
        , direction (csc::Vector3(c1, c2, c3))
        , support   (csc::Vector3(x0, y0, z0))
        , t         (0)
        , _t_min    (-math::inf)
        , _t_max    (+math::inf)
        , t_min     (_t_min)
        , t_max     (_t_max) {
        _name = "Line";
        checkConsistency();
    }

    // line through origin, defined by direction cosines
    Line::Line(double cx, double cy, double cz)
        : bc()
        , direction (csc::Vector3(cx,cy,cz) )
        , support   (csc::Vector3() )
        , t         (0)
        , _t_min    (-math::inf)
        , _t_max    (+math::inf)
        , t_min     (_t_min)
        , t_max     (_t_max) {
        _name = "Line";
        checkConsistency();
    }

    // line through origin, defined by Angles w.r.t. coordinate axes
    Line::Line(const Angle& ax, const Angle& ay, const Angle& az)
        : bc()
        , t_min (_t_min)
        , t_max (_t_max) {
        Line(cos(ax), cos(ay), cos(az));
    }


    // Length of the line segment
    Length Line::length() const {
        if (!_degenerate && !_ill_defined &&
            !isinf(_t_max) && !isinf(_t_min))
            return Length::meters( ((direction*_t_min)-(direction*_t_min)).length() );
        else
            return Length::infinite;
    }

    // Area, circumference are trivial
    Area Line::area()            const { return {}; }
    Length Line::circumference() const { return length(); }
    Volume Line::volume()        const { return {}; }


    // check if a point lies on the line
    bool Line::isMember(const csc::Point& pt){
        csc::Point p2 = (pt-support)/direction;
        t = (p2[0]+p2[1]+p2[2])/3.0;
        return (t <= _t_max && t >= _t_min);
    }
    bool Line::isMember(const Point& pt){
        return isMember(pt.coordinates);
    }

    // set upper/lower limits
    void Line::setLimits(double tmin, double tmax){
        _t_min = tmin;
        _t_max = tmax;
        checkConsistency();
    }
    void Line::setLimits(const csc::Point& tmin, const csc::Point& tmax) {
        _t_min = at(tmin);
        _t_max = at(tmax);
        checkConsistency();
    }
    void Line::setLimits(const Point& tmin, const Point& tmax){
        setLimits(tmin.coordinates, tmax.coordinates);
    }


    // Get coordinates of the line at some value for the parameter t
    Point Line::at (double tt){
        t = (tt < t_min ? t_min : (tt>t_max?t_max:tt));
        evaluate();
        return Point(_current_coordinate);
    }
    Point Line::operator[] (double tt){
        return at(tt);
    }

    // get value for the parameter for some coordinates on the line
    double Line::at(const csc::Point& p) {
        if (isMember(p)) // isMember sets t
            return t;
        else
            throw csc::ex::Condition("Point has to lie on the line before it defines the parameter.");
    }
    double Line::at(const Point& p) {
        return at(p.coordinates);
    }
    double Line::operator[](const Point& p) {
        return at(p);
    }
    double Line::operator[](const csc::Point& p) {
        return at(p);
    }

    // allow casting to simpler types in case Line is degenerate
    Line::operator Point() {
        if (degenerate && !ill_defined)
            return at(t);
        else
            throw csc::ex::NotImplemented("Cast from Line to Point only possible for degenerate yet well-defined Lines.");
    }
    Line::operator Nil() {
        if (ill_defined)
            return Nil();
        else
            throw csc::ex::NotImplemented("Cast from Line to Nil only possible for ill-defined Lines.");
    }

    /// Human-readable string representation
    std::string Line::toString() const {
        std::stringstream output;
        output
            << "Line: \"" << _name << "\"\n"
            << "name      : " << name << "\n"
            << "direction : " << direction << "\n"
            << "support   : " << support << "\n";
        return output.str();
    }





    // Polygon
    // ========================================


    void Polygon::evaluate(){}

    void Polygon::checkConsistency(){

        _num_sides = get_edges().size();

        //_is_counterclockwise

       // list is empty

         for (unsigned int i=0; i<_num_sides; ++i)
            if (angle(i) > Angle::_180 )
                _is_convex =  false;
        //_is_convex &= !is_star();


    }

    // Get a vector containing all internal angles
    std::vector<Angle> Polygon::getAngles()
    {
        auto elist = get_edges();
        std::vector<Angle> angles;
        csc::Vector2D v1, v2;

        for (auto i : elist)
        {
            v1 = i.first();
            v2 = i.last();

            // Angle is pi - atan2(cross,dot), which in 2D reduces to:
            _is_counterclockwise ?
                angles.push_back(Angle::pi - Atan2(
                    v1.x*v2.y - v1.y*v2.x,
                    v1.x*v2.x + v1.y*v2.y ))
            :
                angles.push_back(Angle::pi + Atan2(
                        v1.x*v2.y - v1.y*v2.x,
                        v1.x*v2.x + v1.y*v2.y ));
        }

        return angles;
    }

    // check whether vertices are sorted
    bool Polygon::isSorted(){
        auto a = getAngles();

        // TODO: (Rody Oldenhuis)
        //return std::is_sorted(a.cbegin(), a.cend());
        return true;
    }


    //// Sort points counterclockwise
    //void Polygon::sortCounterClockwise()
    //{
        //if (isSorted()) {
            //_is_counterclockwise = true;
            //return;
        //}

        //auto elist = get_edges();
        //Point C = centroid();

        //auto angles = getAngles();

        //using iAngle = std::pair<Angle,unsigned int>;
        //std::vector<iAngle> angles_with_inds;
        //for (auto i : angles)
            //angles_with_inds.push_back(iAngle(i, std::distance(angles.cbegin(),i)));

        //std::sort(
            //angles_with_inds.begin(), angles_with_inds.end(),
            //[](const iAngle& l, const iAngle& r){ return l.first < r.first; });

        //auto new_elist = elist;
        //for (auto i : angles_with_inds)
            //new_elist[std::distance(angles_with_inds.cbegin(),i)] =
                //elist[i.second];

        //elist.swap(new_elist);
        //_is_counterclockwise = true;
    //}

    // empty polygon - return square on XY plane, with Z axis crossing
    // at centroid, sides of length 1
    Polygon::Polygon()
        : is_convex      (!_is_convex)
        , is_concave     (_is_convex)
        , is_star        (_is_star)
        , self_intersects(_is_star)
        , num_sides      (_num_sides)
        , num_vertices   (_num_sides)
    {

    }

    // Given a set of points
    Polygon::Polygon(const std::vector<Point>& pts)
        : Polygon()
    {
        std::vector<csc::Vector3> vec;
        for (auto i : pts) vec.push_back(i.coordinates);
        // TODO: (Rody Oldenhuis)
        //Polygon(vec);
    }

    // Given a set of vectors
    Polygon::Polygon(const std::vector<csc::Vector3>& pts)
        : Polygon()
    {
        // FIXME: (Rody Oldenhuis) we have to translate the 3D
        // representation to a list of 2D vectors, a coordinate
        // basis u,v and an offset w.r.t. the origin.

        // First, some checks
        auto elist = get_edges();
        for (auto i = pts.begin(); i < pts.end()-1; ++i){
            // TODO: (Rody Oldenhuis)
        }

    }

    // return *internal* angle at vertex i
    Angle Polygon::angle(unsigned int i) const
    {
        // TODO: (Rody Oldenhuis)
        return {};

        //std::list<csc::Segment2D>::const_iterator p1,p2;

        //if (i==0) {
            //p1 = end();
            //p2 = begin();
        //}
        //else {
            //p1 = begin()+i-1;
            //p2 = p1+1;
        //}


    }

    // centroid
    Point Polygon::centroid() const {
        csc::Vector3 C();
        for (auto i : get_edges())
            C += {i};
        C /= num_vertices;
        return C;
    }

    //// wrapper for area
    //Area Polygon::area() const {
        //return Area::squareMeters(csc::Polygon2D::area());
    //}

    //// compute circumference
    //Length Polygon::circumference() const {
        //Length len();
        //auto lst = get_edges();
        //for (auto j=i=lst.begin(); i<lst.end(); j=i,++i)
            //len += {(i-j).length()};
        //len += {(i-lst[0]).length()};
        //return len;
    //}

    //// Human-readable
    //std::string Polygon::toString() const;




    // Triangle
    // ========================================

    void Triangle::evaluate(){}
    void Triangle::checkConsistency(){}

    Triangle::Triangle(const Point& a, const Point& b, const Point& c)
        : Bc()
        , a(a.coordinates)
        , b(b.coordinates)
        , c(c.coordinates)
    {}

    Triangle::Triangle(const csc::Vector3& a, const csc::Vector3& b, const csc::Vector3& c)
        : Bc()
        , a(a), b(b), c(c)
    {}


    // Area
    Area Triangle::area() const {
        return Area::squareMeters( csc::vector_product(b-a,c-b).length()/2.0 );
    }

    // circumference
    Length Triangle::circumference() const {
        return Length::meters(
            (b-a).length() +
            (c-b).length() +
            (a-c).length() );
    }

    Volume Triangle::volume() const { return {}; }

    /// Return angles in the triangle
    Angle Triangle::ab() { return angle(0); }
    Angle Triangle::bc() { return angle(1); }
    Angle Triangle::ca() { return angle(2); }

    /// Human-readable string representation
    std::string Triangle::toString() const {
        std::stringstream output;
        output
            << "Triangle \"" << _name << "\"\n"
            << "a: " << a << "\n"
            << "b: " << b << "\n"
            << "c: " << c;
        return output.str();
    }





    // Rectangle
    // ========================================




    // Square
    // ========================================



    // Conics
    // ========================================




    // Generic Conic
    // ========================================




    // Hyperbola
    // ========================================



    // Parabola
    // ========================================




    // Ellipse
    // ========================================
    void Ellipse::checkConsistency()
    {
        double small = math::epsilon;
        Length zero  = {};

        _t_max.wrap_positive();
        _t_min.wrap_positive();

        if (_t_max < _t_min)
            std::swap(_t_max, _t_min);

        if (a < b){
            std::swap(_a,_b);
            std::swap(u,v);
        }

        _e = std::sqrt(1.0-((b/a)*(b/a)));

        // Really hopeless situations:
        if ( e<0.0 || (a<Length::small && b<Length::small) ||
            (u.isZeroLength() && v.isZeroLength() ) )
        {
            _ill_defined = true;
            return;
        }

        // LOTS of degenerate cases:
        if (fabs(_t_max-_t_min) < Angle::small || // == Point
            abs(a-b) < Length::small || (e<small && e > 0.0) || // == Circle
            ( (a < Length::zero) ^ (b < Length::zero) ) || // == Hyperbola
            (e == 1.0 && b!=Length::zero) ||  // == Parabola
            ((a<Length::small && a>=Length::zero) ^ (b<Length::small && b>=Length::zero)) || (u.isZeroLength() ^ v.isZeroLength()) ) // == Line or Point
            _degenerate = true;
    }

    void Ellipse::evaluate(){
        _current_coordinate = C + u*cos(t) + v*sin(t);
    }


    // default ellipse: a = 1, b = 0.5
    Ellipse::Ellipse()
        : bc()

        , _a(Length::meter)
        , _b(Length::meters(0.5))
        , _e(cos(math::pio6)) // Ha!

        , _t_min (Angle::zero)
        , _t_max (Angle::tau )
        , t      (Angle::zero)

        , u (csc::Vector3(1,0,0))
        , v (csc::Vector3(0,1,0))
        , C (csc::Vector3(     ))

        , a(_a)
        , b(_b)
        , e(_e)
        , t_min(_t_min)
        , t_max(_t_max)
    {
        _name = "(default ellipse)";
        checkConsistency();
    }

    Ellipse::Ellipse(const Length& a, const Length& b)
        : bc()

        , _a(Length::meter)
        , _b(Length::meters(0.5))
        , _e(cos(math::pio6)) // Ha!

        , _t_min (Angle::zero)
        , _t_max (Angle::tau )
        , t      (Angle::zero)

        , u (csc::Vector3(1,0,0))
        , v (csc::Vector3(0,1,0))
        , C (csc::Vector3(     ))

        , a(_a)
        , b(_b)
        , e(_e)
        , t_min(_t_min)
        , t_max(_t_max)
    {
        _name = "(no name)";
        checkConsistency();
    }

    Ellipse::operator Nil() {
        if (_ill_defined) return Nil();
        throw csc::ex::Condition("Only ill-defined ellipses are Nil.");
    }

    Ellipse::operator Point() {
        if (_degenerate &&
           (((a<Length::small && a>=Length::zero) ^ (b<Length::small && b>=Length::zero)) ||
            (u.isZeroLength() ^ v.isZeroLength())  ||
            abs(_t_max-_t_min) < Angle::small) )
        {
           evaluate();
           return {_current_coordinate};
        }

        throw csc::ex::Condition("Ellipse does not degenerate into a point.");
    }


    Ellipse::operator Circle(){
        if (_degenerate &&
            abs(a-b) < Length::small &&
            e > 0 && e < math::epsilon)
            return Circle( (a+b)/2.0, C,u,v);

        throw csc::ex::Condition("Semi-major and minor axes must be equal for an ellipse to equal a circle.");
    }

    Ellipse::operator Hyperbola() {
        // TODO
        return {};
    }

    Ellipse::operator Parabola() {
        // TODO
        return {};
    }

    // Intersections
    // --------------------------------

    auto Ellipse::intersect(const Line& L) const -> line_intersection
    {
        // TODO: (Rody Oldenhuis)
        return Nil();
    }

    auto Ellipse::intersect(const Plane& L) const -> plane_intersection
    {
        // TODO: (Rody Oldenhuis)
        return Nil();
    }

    // Distances
    // --------------------------------

    Length Ellipse::distance(const Point&   P) const
    {
        // TODO: (Rody Oldenhuis)
        return Length::zero;
    }

    Length Ellipse::distance(const Plane&   P) const
    {
        // TODO: (Rody Oldenhuis)
        return Length::zero;
    }

    Length Ellipse::distance(const Ellipse& E) const
    {
        // TODO: (Rody Oldenhuis)
        return Length::zero;
    }

    // Human-readable string representation
    std::string Ellipse::toString() const
    {
        std::stringstream output;
        output
            << "Ellipse \"" << _name << "\"\n"
            << "\n"
            << "\n";
        return output.str();
    }



    // Circle
    // ========================================

    void Circle::evaluate(){
        _current_coordinate = (u*cos(t) + v*sin(t))*_radius.meters() + C;
    }

    void Circle::checkConsistency() {

    }

    // unit circle
    Circle::Circle()
        : Ellipse()
        , _radius(Length::meter)
        , radius(_radius)
    {
        u = csc::Vector3(1, 0, 0);
        v = csc::Vector3(0, 1, 0);
        C = {};
        _name = "Unit circle";
    }

    // general constructor
    Circle::Circle(const Length& r, const csc::Vector3& u, const csc::Vector3& v, const csc::Vector3& C)
        : Ellipse()
        , _radius(r)
        , radius(_radius)
    {
        this->u = u;
        this->v = v;
        this->C = C;
        _name = "(no name)";
    }

    // type casting for degenerate circles
    Circle::operator Point() {
        if (degenerate){
            evaluate();
            return _current_coordinate;
        }
        else
            throw csc::ex::NotImplemented("Circle cast to Point is only implemented for degenerate circles.");
    }


    // Human-readable string representation
    std::string Circle::toString() const
    {
        std::stringstream output;
        output
            << "Circle \"" << _name << "\"\n"
            << "radius : " << radius << "\n"
            << "center : " << C << "\n"
            << "axis u : " << u << "\n"
            << "axis v : " << v << "\n";
        return output.str();
    }




    // Quadrics
    // ========================================

    // Generic quadric

    namespace _ignore
    {

    }

    // Ellipsoid
    // ========================================


    // Spheroid: a special ellipsoid
    // ========================================


    // Sphere : a special spheroid
    // ========================================
    Sphere::Sphere()
        : bc()
        , radius  (Length::meter)
        , center  {}
    {
        _name = "Unit sphere";
    }

    // construct arbitrary other sphere
    Sphere::Sphere(Length radius, csc::Vector3 center)
        : radius(radius)
        , center(center)
    {
        _name = "(no name)";
    }

    Sphere::Sphere(Length radius)       { Sphere(radius, {}); }
    Sphere::Sphere(csc::Vector3 center) { Sphere(Length::meter, center); }

    // the basics
    Area   Sphere::area()          const { return radius*radius * 4.0*math::pi; }
    Volume Sphere::volume()        const { return area()*radius/3.0; }
    Length Sphere::circumference() const { return radius*2.0*math::pi; }

    // intersections
    // ----------------------------

    auto Sphere::intersect(const Plane& P) const -> plane_intersection
    {
        // centre-plane distance
        Length dist = Length::meters(
            csc::scalar_product(P.offset-center, csc::normalize(P.offset)) );

        // plane just touches: Point
        if (dist == radius) return Point(P.normal_vector*dist.meters());

        // no intersection: Nil
        if (dist > radius) return Nil();

        // plane is inside the sphere: circle
        return Circle(
            sqrt(radius*radius - dist*dist), // radius
            center + csc::normalize(P.offset)*dist.meters(), // center
            P.u, P.v); // plane is parallel to circle, obviously
    }

    auto Sphere::intersect(const Sphere& P) const -> sphere_intersection
    {
        // TODO
        return Circle();

    }

    // distances
    // ----------------------------

    // Plane-sphere distance
    // Distance is always positive, except when the sphere intersects with the
    // plane. In that case, the negative value is the distance between the plane
    // and the centre of the sphere.
    Length Sphere::distance(const Plane& P) const {
        // centre-plane distance
        Length d = abs(Length::meters(
            csc::scalar_product(P.offset-center, csc::normalize(P.offset)) ));
        // subtract radius, or negate upon intersect
        return d<radius ? -d : d-radius;
    }

    // Sphere-sphere distance
    // distance is always positive, exceptwhen the two spheres intersect.
    // In that case, the negative value willequal the distance between the centres
    // of the spheres.
    Length Sphere::distance(const Sphere& S) const
    {
        Length
            // distance between centers
            d = Length::meters((S.center-center).length()),
            // sum of radii
            sumR = radius+S.radius;

        // sumtract sum when appropriate, negate upon intersect
        return d<sumR ? -d : d-sumR;
    }

    // Human-readable string representation
    std::string Sphere::toString() const
    {
        std::stringstream output;
        output
            << "Sphere \"" << _name << "\"\n"
            << "radius: " << radius
            << "center: " << center;
        return output.str();
    }






    // Hyperbolic cylinder
    // ========================================





    // Elliptic cylinder
    // ========================================




    // Cylinder
    // ========================================

    Cylinder::Cylinder()
    {
    }

    // Human-readable string representation
    std::string Cylinder::toString() const
    {
        std::stringstream output;
        output
            << "Cylinder: \"" << _name << "\"\n"
            << "radius : \n"
            << "height : \n";
        return output.str();
    }


    // Elliptic cone
    // ----------------------------

    // Cone
    // ----------------------------
    Cone::Cone(){}

    // Human-readable string representation
    std::string Cone::toString() const
    {
        std::stringstream output;
        output
            << "Cone: \"" << _name << "\"\n"
            << "\n"
            << "\n";
        return output.str();
    }




    // Elliptic paraboloid
    // ========================================




    // Paraboloid
    // ========================================



    // Elliptic hyperboloid
    // ========================================




    // Hyperboloid
    // ========================================




    // Pyramid
    // ========================================

    void Pyramid::checkConsistency(){}
    void Pyramid::evaluate(){}

    Pyramid::Pyramid(const Point& P, const Polygon& poly, const csc::Vector3& D)
        : bc()
        , tip          (P.coordinates)
        , direction    (D)
        , crossSection (poly)
    {


    }

    // Area
    Area Pyramid::area() const {
        // FIXME: (Rody Oldenhuis)
        return Area();

        //Area A();
        //A += crossSection.area();
        //auto pts = crossSection.get_edges();
        //for (auto i=pts.cbegin()+1,j=i-1; i<pts.cend(); j=i,++i)
            //A += Triangle(tip,i,j).area();
    }

    // Volume
    Volume Pyramid::volume() const {
        throw csc::ex::NotImplemented("Volume of pyramid not yet fully implemented.");
        return Volume::cubicMeters(
            // TODO: multiply by height of pyramid.
            crossSection.area().squareMeters()/3.0);
    }

    // Human-readable string representation
    std::string Pyramid::toString() const {
        std::stringstream output;
        output
            << "Pyramid \"" << _name << "\"\n"
            << "\n"
            << "\n";
        return output.str();
    }





    // Shape >> ostream defs
    // ========================================

    std::ostream& operator<<(std::ostream& target, const _ignore::Shape& S) {
        target << S.toString();
        return target;
    }

}







