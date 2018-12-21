#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <json/json.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <ceres/ceres.h>

#include "cost_functors.hpp"
#include "prettyprint.hpp"

using point_t = boost::geometry::model::d2::point_xy<double>;
using polygon_t = boost::geometry::model::polygon<point_t>;

struct sample_t
{
    double x;
    double y;
};

std::ostream &operator<<(std::ostream &out, const sample_t &sample)
{
    out << "<" << sample.x << ", " << sample.y << ">";
    return out;
}

std::vector<sample_t> resample(const std::vector<sample_t> &points, const size_t &N)
{
    std::vector<double> offsets(points.size());
    offsets[0] = 0;
    for (size_t i = 0; i < points.size() - 1; ++i)
    {
        offsets[i + 1] = offsets[i] + std::hypot(points[i + 1].x - points[i].x, points[i + 1].y - points[i].y);
    }

    std::cout << "offsets" << offsets << std::endl;

    const double &length = offsets.back();
    const double increment = length / (N - 1);

    std::vector<sample_t> result(N);

    for (size_t i = 0; i < N; ++i)
    {
        double arclen = i * increment;
        const auto it = std::lower_bound(offsets.cbegin(), offsets.cend(), arclen);
        const size_t index = std::distance(offsets.cbegin(), it);
        const double remainder = arclen - offsets[index - 1];
        const double frac = remainder / (offsets[index] - offsets[index - 1]);
        result[i] = sample_t{
            points[index - 1].x + frac * (points[index].x - points[index - 1].x),
            points[index - 1].y + frac * (points[index].y - points[index - 1].y)};
    }

    return result;
}

Json::Value samples_to_json(const std::vector<sample_t>& samples) {
    Json::Value result;
    for(const sample_t& sample: samples) {
        Json::Value s;
        s.append(sample.x);
        s.append(sample.y);
        result.append(s);
    }
    return result;
}

std::vector<sample_t> json_to_samples(const Json::Value& json_samples) {
    std::vector<sample_t> result(json_samples.size());
    for(Json::ArrayIndex i = 0; i < json_samples.size(); ++i) {
        result[i] = sample_t{
            json_samples[i][0].asDouble(),
            json_samples[i][1].asDouble()
        };
    }

    return result;
}

std::vector<sample_t> initialize(const std::vector<sample_t> raw_init, const size_t N_SAMPLES)
{
    const std::vector<sample_t> resampled = resample(raw_init, N_SAMPLES - 4);
    std::vector<sample_t> initialization(N_SAMPLES);
    std::copy(resampled.cbegin(), resampled.cend(), initialization.begin() + 2);
    initialization[0] = initialization[1] = initialization[2];
    initialization[N_SAMPLES - 2] = initialization[N_SAMPLES - 1] = initialization[N_SAMPLES - 3];
    return initialization;
}

int main()
{
    Json::Value world;
    std::ifstream("world.json") >> world;
    
    polygon_t poly;

    for(const sample_t& point: json_to_samples(world["bounds"])) {
        boost::geometry::append(poly.outer(), point_t{
            point.x,
            point.y
        });
    }

    std::vector<sample_t> raw_init = json_to_samples(world["rawInit"]);

    const size_t N_SAMPLES = 200;
    const double w_vel(world["weights"]["vel"].asDouble()), 
        w_acc(world["weights"]["acc"].asDouble()), 
        w_jerk(world["weights"]["jerk"].asDouble()), 
        w_bound(world["weights"]["bound"].asDouble());

    std::vector<sample_t> pts = initialize(raw_init, N_SAMPLES);

    Json::Value json_out;
    json_out["initialization"] = samples_to_json(pts);

    std::cout << pts << std::endl;

    ceres::Problem ceres_problem;
    ceres::CostFunction *cost_function;

    for (int i = 0; i < N_SAMPLES; ++i)
    {
        double *p0 = (double *)&(pts[i]);

        cost_function =
            new ceres::NumericDiffCostFunction<PolygonCostFunctor<point_t, polygon_t>, ceres::FORWARD, 1, 2>(
                new PolygonCostFunctor<point_t, polygon_t>(poly));
        ceres_problem.AddResidualBlock(cost_function, new ceres::ScaledLoss(NULL, w_bound, ceres::TAKE_OWNERSHIP), p0);
    }

    for (int i = 0; i < N_SAMPLES - 2; ++i)
    {
        double *p0 = (double *)&(pts[i]);
        double *p1 = (double *)&(pts[i + 1]);
        double *p2 = (double *)&(pts[i + 2]);

        cost_function =
            new ceres::AutoDiffCostFunction<AccCostFunctor, 2, 2, 2, 2>(
                new AccCostFunctor());
        ceres_problem.AddResidualBlock(cost_function, new ceres::ScaledLoss(NULL, w_acc, ceres::TAKE_OWNERSHIP), p0, p1, p2);
    }

    for (int i = 0; i < N_SAMPLES - 3; ++i)
    {
        double *p0 = (double *)&(pts[i]);
        double *p1 = (double *)&(pts[i + 1]);
        double *p2 = (double *)&(pts[i + 2]);
        double *p3 = (double *)&(pts[i + 3]);

        cost_function =
            new ceres::AutoDiffCostFunction<JerkCostFunctor, 2, 2, 2, 2, 2>(
                new JerkCostFunctor());
        ceres_problem.AddResidualBlock(cost_function, new ceres::ScaledLoss(NULL, w_jerk, ceres::TAKE_OWNERSHIP), p0, p1, p2, p3);
    }

    const std::vector<size_t> constant_blocks = {{0, 1, 2, N_SAMPLES - 3, N_SAMPLES - 2, N_SAMPLES - 1}};

    for (const auto &block : constant_blocks)
    {
        ceres_problem.SetParameterBlockConstant((double *)&(pts[block]));
    }

    ceres::SubsetParameterization *pin_x = new ceres::SubsetParameterization(2, {{0}});

    for (size_t i = 0; i < 20; ++i)
    {
        ceres_problem.SetParameterization((double *)&(pts[i]), pin_x);
        ceres_problem.SetParameterization((double *)&(pts[N_SAMPLES - i - 1]), pin_x);
    }

    ceres::Solver::Options options;
    options.max_num_iterations = 1000;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;

    ceres::Solve(options, &ceres_problem, &summary);
    std::cout << summary.BriefReport() << std::endl;

    json_out["solution"] = samples_to_json(pts);

    std::ofstream("result.json") << json_out;

    return 0;
}
