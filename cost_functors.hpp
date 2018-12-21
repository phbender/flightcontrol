#include <boost/geometry.hpp>


template<typename POINT_T, typename POLYGON_T>
struct PolygonCostFunctor {
    const POLYGON_T poly;
    PolygonCostFunctor(const POLYGON_T& poly): poly{poly} {
        
    }

    bool operator()(const double* p0, double* residual) const {
        const POINT_T pt(p0[0], p0[1]);
        bool is_inside = boost::geometry::within(pt, this->poly);
        residual[0] = is_inside ? 0 : boost::geometry::distance(pt, this->poly);
        return true;
    }
};

struct VelCostFunctor {

    const double v_des;
    const double dt;

    VelCostFunctor(const double& v_des, const double& dt) : v_des{v_des}, dt{dt}
    { }

    /// Template signature <2, 2, 2>
    template<typename T>
    bool operator()(const T* const p0, const T* const p1, T* residual) const
    {
        T vx = (p1[0] - p0[0]) / T(dt);
        T vy = (p1[1] - p0[1]) / T(dt);

        residual[0] = (v_des - vx);
        residual[1] =  vy;

        return true;
    }
};


struct AccCostFunctor {



    AccCostFunctor() {}

    template<typename T>
    bool operator()(const T* const p0, const T* const p1, const T* const p2, T* residual) const {
        for( int i = 0; i < 2; ++i )
        {
            residual[i] = p0[i] - T(2) * p1[i] + p2[i];
        }
        return true;
    }
};

/// This cost functor takes four parameter blocks and returns the two residuals jerk_x and jerk_y.
struct JerkCostFunctor {

    JerkCostFunctor() {}

    /// Template signature <2, 2, 2, 2, 2>
    template<typename T>
    bool operator()(const T* const p0, const T* const p1,
                    const T* const p2, const T* const p3,
                    T* residual) const {
        for( int i = 0; i < 2; ++i )
        {
            residual[i] = (-p0[i] + T(3) * p1[i] - T(3) * p2[i] + p3[i]) ;
        }
        return true;
    }
};
