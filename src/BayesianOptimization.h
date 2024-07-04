#ifndef CORE_H
#define CORE_H


#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_octree.h>
#include <cinolib/drawable_segment_soup.h>

#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_octree.h>
#include <cinolib/drawable_segment_soup.h>

#include <Eigen/Core>
#include <memory>
#include <cmath>
#include <vector>
#include <queue>



#include <ShellThickening.h>
#include <universal.h>
#include <stdio.h>
#include <SelfIntersect.h>
#include <Deformation.h>
#include <cmath>
#include <iostream>
#include <random>
#include <sequential-line-search/acquisition-function.hpp>
#include <sequential-line-search/gaussian-process-regressor.hpp>
#include <sequential-line-search/utils.hpp>

namespace sequential_line_search
{
    class GaussianProcessRegressor;
}

class Core
{
public:
    Core();

    static Core& getInstance()
    {
        static Core core;
        return core;
    }

    std::shared_ptr<sequential_line_search::GaussianProcessRegressor> regressor;

    Eigen::MatrixXd X;
    Eigen::VectorXd y;

    double evaluateObjectiveFunction(const Eigen::VectorXd& x) const;

    // For optimization
    void            proceedOptimization();
    Eigen::VectorXd x_max;
    double          y_max;

    // For regression
    void addData(const Eigen::VectorXd& x, double y);
    void computeRegression();
    void setBasicData(cinolib::DrawableTrimesh<> m, double K_SHAPE = 1.0, double K_DETAIL = 1.0);
    void exportModel(std::string output_path, std::string model_name);
    cinolib::DrawableTrimesh<> MODEL;
    double K_SHAPE ; 
    double K_DETAIL ;
    std::pair<std::vector<double>, std::vector<double>> *RHO_RANGE;
};

#endif // CORE_H
