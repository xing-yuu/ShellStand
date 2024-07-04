#include "BayesianOptimization.h"


using namespace sequential_line_search;

using Eigen::MatrixXd;
using Eigen::VectorXd;

Core::Core()
{
    X = MatrixXd::Zero(0, 0);
    y = VectorXd::Zero(0);
    x_max = VectorXd::Zero(0);
    y_max = NAN;

    computeRegression();
}
void Core::setBasicData(cinolib::DrawableTrimesh<> m, double K_SHAPE , double K_DETAIL) {
    this->MODEL = m;
    this->K_SHAPE = K_SHAPE;
    this->K_DETAIL = K_DETAIL;

    std::vector<double> RHO;
    std::cout << "++++++++++++++++++++" << std::endl;
    std::cout << "Start calculate saliency." << std::endl;
    clock_t start = clock();
    //RHO
    std::vector<double>* saliency_score = Env::get_saliency_score(m, unv::FACE_DATA);
    clock_t end = clock();
    std::cout << "End saliency." << std::endl;
    std::cout << "Time cost: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
    std::cout << "++++++++++++++++++++" << std::endl << std::endl;



    std::cout << "++++++++++++++++++++" << std::endl;
    std::cout << "Start calculate visibility." << std::endl;
    clock_t start_2 = clock();
    std::vector <std::vector<double>>* visibility_score = Env::get_visibility_score(m, unv::FACE_DATA, 200, false);
    clock_t end_2 = clock();
    std::cout << "End visibility." << std::endl;
    std::cout << "Time cost: " << (double)(end_2 - start_2) / CLOCKS_PER_SEC << std::endl;
    std::cout << "++++++++++++++++++++" << std::endl << std::endl;

    for (uint i = 0; i < saliency_score->size(); i++) {
        double vi = (saliency_score->at(i) + std::max(visibility_score->at(0)[i], visibility_score->at(1)[i])) / 2;
        if (isnan(vi))
            RHO.push_back(0.5);
        else {
            RHO.push_back(vi);
        }
    }
    std::cout << "SALIENCY Time cost: " << (double)(end_2 - start) / CLOCKS_PER_SEC << std::endl;
    //RANGE
    std::cout << "++++++++++++++++++++" << std::endl;
    std::cout << "Start calculate upper bound." << std::endl;
    clock_t start_3 = clock();
    std::vector <std::vector<double>>* thickening_range = Env::get_thickening_range(m);
    std::vector<double>* face_range_a = unv::point2face(thickening_range->at(0), m);
    std::vector<double>* face_range_b = unv::point2face(thickening_range->at(1), m);
    std::vector<double>range;
    for (uint i = 0; i < face_range_b->size(); i++) {
        range.push_back(std::max(1.2, std::max(face_range_a->at(i), face_range_b->at(i))));
    }
    clock_t end_3 = clock();
    std::cout << "End upper bound.." << std::endl;
    std::cout << "Time cost: " << (double)(end_3 - start_3) / CLOCKS_PER_SEC << std::endl;
    std::cout << "++++++++++++++++++++" << std::endl << std::endl;
    double MAX_rho = INT_MIN, MIN_rho = INT_MAX;
    double MAX_range = INT_MIN, MIN_range = INT_MAX;
    for (uint i = 0; i < m.num_polys(); i++) {
        MAX_rho = std::max(MAX_rho, RHO[i]);
        MIN_rho = std::min(MIN_rho, RHO[i]);
        MAX_range = std::max(MAX_range, range[i]);
        MIN_range = std::min(MIN_range, range[i]);
    }
    std::cout << "MAX rho is \t" << MAX_rho << std::endl;
    std::cout << "MIN rho is \t" << MIN_rho << std::endl;
    std::cout << "MAX range is \t" << MAX_range << std::endl;
    std::cout << "MIN range is \t" << MIN_range << std::endl;

    this->RHO_RANGE = new std::pair<std::vector<double>, std::vector<double>>(RHO, range);


}
void Core::proceedOptimization()
{
    const VectorXd x = (X.cols() == 0) ? utils::GenerateRandomVector(1) : acquisition_func::FindNextPoint(*regressor);
    const double   y = evaluateObjectiveFunction(x);

    std::cout << "x: " << x.transpose() << ((X.cols() == 0) ? " (randomly chosen)" : "") << std::endl;
    std::cout << "y: " << y*(-1.0) << std::endl;

    addData(x, y);
    computeRegression();

    const int num_data_points = X.cols();

    VectorXd f(num_data_points);
    for (int i = 0; i < X.cols(); ++i)
    {
        f(i) = regressor->PredictMu(X.col(i));

        int best_index;
        y_max = f.maxCoeff(&best_index);
        x_max = X.col(best_index);
    }
}
void Core::exportModel(std::string output_path, std::string _model_name) {

    std::string model_name = _model_name.substr(0, _model_name.length() - 4);
    if (_access(output_path.c_str(), 0) == -1)	//如果文件夹不存在
        _mkdir(output_path.c_str());
    std::string  output_path_2 = output_path + "//" + model_name;
    if (_access(output_path_2.c_str(), 0) == -1)	//如果文件夹不存在
        _mkdir(output_path_2.c_str());
    std::string information = "_kd" + std::to_string(this->K_DETAIL) + "_ks" + std::to_string(this->K_SHAPE);
    double alpha = x_max(0);
    std::cout << "*********************************" << std::endl;
    std::cout << "Start output deformation model." << std::endl;
    std::pair<cinolib::DrawableTrimesh<>, double> defor_res = Deform::deform_ARAP_DATABASE(this->MODEL, alpha);
    cinolib::DrawableTrimesh<> deformation_model = defor_res.first;
    std::string deformation = output_path_2 + "//" + "deformation_kd" + std::to_string(this->K_DETAIL) + "_ks" + std::to_string(this->K_SHAPE)+".obj";
    deformation_model.save(deformation.data());
    std::cout << "End ." << std::endl;
    std::cout << "*********************************" << std::endl << std::endl;
    std::cout << "*********************************" << std::endl;
    std::cout << "Start output Thickening model." << std::endl;
    std::pair<std::vector<double>, double> thick_res = Env::shell_THICKENING(defor_res.first, *this->RHO_RANGE);
    std::string thickening = output_path_2 + "//" + "volume_kd" + std::to_string(this->K_DETAIL) + "_ks" + std::to_string(this->K_SHAPE) + "_alpha" + std::to_string(alpha) + ".obj";
    Env::output_volume_model(thickening, defor_res.first, thick_res.first);


    std::cout << "End ." << std::endl;
    std::cout << "*********************************" << std::endl;
    
}


void Core::addData(const VectorXd& x, double y)
{
    if (X.rows() == 0)
    {
        this->X = x;
        this->y = VectorXd::Constant(1, y);
        return;
    }

    const unsigned D = X.rows();
    const unsigned N = X.cols();

    MatrixXd newX(D, N + 1);
    newX.block(0, 0, D, N) = X;
    newX.col(N) = x;
    this->X = newX;

    VectorXd newY(this->y.rows() + 1);
    newY << this->y, y;
    this->y = newY;
}
/*
double Core::evaluateObjectiveFunction(const Eigen::VectorXd& x) const
{
    return 1.0 - 1.5 * x(0) * std::sin(x(0) * 13.0);
}
*/
double Core::evaluateObjectiveFunction(const Eigen::VectorXd& x) const
{
    double alpha = x(0);
    clock_t start = clock();
    std::cout << "========================" << std::endl;
    std::cout << "Start deformation." << std::endl;
    std::pair<cinolib::DrawableTrimesh<>, double> defor_res = Deform::deform_ARAP_DATABASE(this->MODEL, alpha);
    clock_t end = clock();
    std::cout << "End deformation." << std::endl;
    std::cout << "Time cost: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
    std::cout << "========================" << std::endl << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << "Start Thickening." << std::endl;
    std::pair<std::vector<double>, double> thick_res = Env::shell_THICKENING(defor_res.first, *this->RHO_RANGE);
    clock_t end2 = clock();
    std::cout << "End thickening." << std::endl;
    std::cout << "Time cost: " << (double)(end2 - end) / CLOCKS_PER_SEC << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << "Now Alpha is \t" << alpha << std::endl;
    std::cout << "Now F_shape is \t" << defor_res.second << std::endl;
    std::cout << "Now F_detail is \t" << thick_res.second << std::endl;

    return (this->K_DETAIL * thick_res.second + this->K_SHAPE * defor_res.second) * (-1.0);

}

void Core::computeRegression()
{
    regressor = std::make_shared<GaussianProcessRegressor>(X, y);
}
