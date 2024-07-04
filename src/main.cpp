#include <cinolib/meshes/meshes.h>
#include <cinolib/ARAP.h>
#include <cinolib/gl/glcanvas.h>
#include <time.h>
#include <limits.h>
#include <omp.h>
#include <yaml-cpp/yaml.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <filesystem>
#include <algorithm>
#include <Evaluation.h>
#include <SelfIntersect.h>
#include <Deformation.h>
using namespace cinolib;
#include <BayesianOptimization.h>
#include <iostream>
#include <iostream>
#include <string>
#include <filesystem>
#include <direct.h>
#include <ShellThickening.h>
#include <mosek.h>
#include<ModelOperation.h>
void bys_test() {

    using Eigen::MatrixXd;
    using Eigen::VectorXd;

   
        Core core;

        constexpr unsigned n_trials = 6;
        constexpr unsigned n_iterations = 15;

        for (unsigned trial = 0; trial < n_trials; ++trial)
        {
            std::cout << "========================" << std::endl;
            std::cout << "Trial " << trial + 1 << std::endl;
            std::cout << "========================" << std::endl;

            // Iterate optimization steps
            for (unsigned i = 0; i < n_iterations; ++i)
            {
                std::cout << "---- Iteration " << i + 1 << " ----" << std::endl;
                core.proceedOptimization();
            }

            std::cout << std::endl << "Found minimizer: " << core.x_max.transpose() << std::endl;
            std::cout << "Found minimum: " << -core.y_max << std::endl << std::endl;

            // Reset the optimization
            core.X = MatrixXd::Zero(0, 0);
            core.y = VectorXd::Zero(0);
            core.x_max = VectorXd::Zero(0);
            core.y_max = NAN;
            core.computeRegression();
        }

  

}
void bys_opt(cinolib::DrawableTrimesh<> m, double K_SHAPE, double K_DETAIL,std::string output_path, std::string model_name,uint n_trials = 20) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;


    Core core;
    std::cout << "========================" << std::endl;
    std::cout << "Bayesian Optimization " << std::endl;
    std::cout << "========================" << std::endl;
    core.setBasicData(m, K_SHAPE, K_DETAIL);
    // Iterate optimization steps
    for (unsigned i = 0; i < n_trials; ++i)
    {
        std::cout << "---- Iteration " << i + 1 << " ----" << std::endl;
        core.proceedOptimization();
    }

    std::cout << std::endl << "Found maximizer: " << core.x_max.transpose() << std::endl;
    std::cout << "Found maximum: " << -core.y_max << std::endl << std::endl;
    core.exportModel(output_path, model_name);
    
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char* argv[]) {
//    Env::Mosek_test();
    //bys_test();
//    msk::Mosek_test();
    clock_t start, finish;
    start = clock();
    YAML::Node config = YAML::LoadFile("../config.yaml");
    std::string datapath = config["InputPath"].as<std::string>();
    std::string outputpath= config["OutputPath"].as<std::string>();
    std::string filename = config["InputModel"].as<std::string>();
    std::string debugpath=config["DebugPath"].as<std::string>();
    std::string mode= config["Mode"].as<std::string>();
    double scale = config["Scale"].as<double>();
    double K_SHAPE = config["K_shape"].as<double>();
    double K_DETAIL = config["K_detail"].as<double>();

    std::string data = datapath + "\\" + filename;
    cinolib::DrawableTrimesh<> m(data.c_str());
    if (scale != 0) {
        m = Ope::scale(m, scale);
    }

    if (mode == "DEBUG-boundary") {
        //debug-boundary test
        std::vector <std::vector<double>>* thickening_range = Env::get_thickening_range(m);
        Env::output_range(m, *thickening_range, debugpath, filename);
    }

    else if (mode == "DEBUG-thickening") {
        //debug-thickening test
        std::vector <std::vector<double>>* thickening_range = Env::get_thickening_range(m);
        std::pair<std::vector<double>, double> delta_m = Env::shell_THICKENING(m);
        std::string DEBUG_thickening_path = debugpath + "//" + filename +  "thickening_Debug_"  + ".obj";
        Env::output_volume_model(DEBUG_thickening_path, m, delta_m.first);
    }


    
    else {
        Deform::get_deform_database(m);
        //bys_opt(m, K_SHAPE, K_DETAIL, outputpath, filename);
        bys_opt(m, 0.2, 0.8, outputpath, filename);
        bys_opt(m, 0.4, 0.6, outputpath, filename);
        bys_opt(m, 0.5, 0.5, outputpath, filename);
        bys_opt(m, 0.6, 0.4, outputpath, filename);
        bys_opt(m, 0.8, 0.2, outputpath, filename);
        finish = clock();
        std::cout << "Total Time cost: " << (double)(finish - start) / CLOCKS_PER_SEC << std::endl;

    }
    return 0;
}
 