#include <boost/geometry.hpp>

// puts cost on leaving a polygon
template <typename POINT_T, typename POLYGON_T>
struct PolygonCostFunctor
{
    const POLYGON_T poly;
    PolygonCostFunctor(const POLYGON_T &poly) : poly{poly}
    {
    }

    bool operator()(const double *p0, double *residual) const
    {
        const POINT_T pt(p0[0], p0[1]);
        bool is_inside = boost::geometry::within(pt, this->poly);
        residual[0] = is_inside ? 0 : boost::geometry::distance(pt, this->poly);
        return true;
    }
};

// puts cost on deviation of desired flight height
struct ElevationCostFunctor
{
    const double elevation_des;

    ElevationCostFunctor(const double &elevation_des) : elevation_des{elevation_des}
    {
    }

    template <typename T>
    bool operator()(const T *const p0, T *residual) const
    {
        residual[0] = (T(elevation_des) - p0[1]);
        return true;
    }
};

// puts cost on acceleration
struct AccCostFunctor
{
    AccCostFunctor() {}

    template <typename T>
    bool operator()(const T *const p0, const T *const p1, const T *const p2, T *residual) const
    {
        for (int i = 0; i < 2; ++i)
        {
            residual[i] = p0[i] - T(2) * p1[i] + p2[i];
        }
        return true;
    }
};

// puts cost on jerk
struct JerkCostFunctor
{
    JerkCostFunctor() {}

    template <typename T>
    bool operator()(const T *const p0, const T *const p1,
                    const T *const p2, const T *const p3,
                    T *residual) const
    {
        for (int i = 0; i < 2; ++i)
        {
            residual[i] = (-p0[i] + T(3) * p1[i] - T(3) * p2[i] + p3[i]);
        }
        return true;
    }
};