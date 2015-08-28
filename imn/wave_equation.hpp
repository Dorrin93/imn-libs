/**
 * @file   wave_equation.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   August, 2015
 * @brief  Function solving wave equations.
 *
 * Header for implementation of numerical methods created for solving wave equations.
 */
#ifndef __IMN_WAVE__
#define __IMN_WAVE__
#include <functional>
#include <fstream>
#include <cmath>
#include <vector>
#include "utils.hpp"

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn{

    using func2dOm = std::function<double(double, double, double)>;

    enum class WaveInitType {RIGID, LOOSE, TWIN_CENTER};


    void verlet_scheme(unsigned grid_density, double dx, double t_min, double t_max, double dt, func2d u, func2d v,
                       func2d a, std::ofstream& file, double t_write = -1, func2d u_sol = nullptr);

    void boundary_of_center(unsigned grid_density, double dx, double t_min, double t_max, double dt, initFunc u,
                            std::ofstream& file, WaveInitType condition_type, double density = 1, double boundary_x = 0.5);

    void damped_oscillations(unsigned grid_density, double dx, double t_min, double t_max, double dt, initFunc u,
                             double beta, std::ofstream& file);

    void forced_oscillations(unsigned grid_density, double dx, double t_min, double t_max, double dt, double beta,
                             func2d force, std::ofstream &file);

    void resonances_enegry(unsigned grid_density, double dx, double t_min, double t_max, double dt, double t_res,
                           double beta, double om_min, double om_max, double dom, std::ofstream &file, func2dOm force);


    enum class WaveProblemType {VERLET, BOC_R, BOC_L, BOC_TC, DAMPED, FORCED};

    void wave_solver(std::vector<double> &u_vals, std::vector<double> &v_vals, std::vector<double> &a_vals, double t_min,
                     double t_max, double dt, double dx, WaveProblemType type, std::ofstream &ofile, double t_w, func2d u_sol,
                     double rho, double boundary, double beta, func2d force);

}

/** @} End of group */

#endif
