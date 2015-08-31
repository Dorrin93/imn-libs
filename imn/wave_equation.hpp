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

    /**
     * @brief 2 dimensional function with omega.
     *
     * First arguemnt is position, second time and third angular velocity.
     */
    using func2dOm = std::function<double(double, double, double)>;

    /**
     * @brief Type of initial condition in boundary of center function.
     * @var RIGID Rigid initial conditions.
     * @var LOOSE Loose ends of the string.
     * @var TWIN_CENTER String has two different media.
     */
    enum class WaveInitType {RIGID, LOOSE, TWIN_CENTER};


    /**
     * @brief Velocity Verlet Scheme.
     *
     * Finite difference method used for solving wave equasion.
     *
     * @param grid_density Density of the computational grid.
     * @param dx Grid step.
     * @param t_min Initial time. Have to be greater than 0.
     * @param t_max Stop time. Have to be greater than t_min.
     * @param dt Time step.
     * @param u 2D function of position. First parameter is x, second is t.
     * @param v 2D function of velocity. Same parameters as position.
     * @param a 2D function of acceleration. Same parameters as position.
     * @param file Output file stream we will be writing resuts. May be open or closed.
     * @param t_write Time when we want to write resuts. If negative, results will be written after all iterations.
     * Default do -1.
     * @param u_sol Solution function. If given, result will be written along result from algorithm.
     * Usefull for comparsion. Deafault to nullptr.
     */
    void verlet_scheme(unsigned grid_density, double dx, double t_min, double t_max, double dt, func2d u, func2d v,
                       func2d a, std::ofstream& file, double t_write = -1, func2d u_sol = nullptr);

    /**
     * @brief Velocity Verlet Scheme with initial condition.
     *
     * It's similar to verlet_scheme function, but initial condition is different. Instead of giving specyfic funtions
     * of position, velocity and acceleration, we give only position in t=0 and status of string
     * (given in condition_type variable).
     *
     * @param grid_density Density of the computational grid.
     * @param dx Grid step.
     * @param t_min Initial time. Have to be greater than 0.
     * @param t_max Stop time. Have to be greater than t_min.
     * @param dt Time step.
     * @param u 2D function of position in t=0. Only one argument (x).
     * @param file Output file stream we will be writing resuts. May be open or closed.
     * @param condition_type Type of initial condition, indicating status of string.
     * @param density Used only if condition_type=TWIN_CENTER. Density of second center(first is 1). Default to 1.
     * @param boundary_x Used only if condition_type=TWIN_CENTER. Position of boundary between centers.
     */
    void boundary_of_center(unsigned grid_density, double dx, double t_min, double t_max, double dt, initFunc u,
                            std::ofstream& file, WaveInitType condition_type, double density = 1, double boundary_x = 0.5);

    /**
     * @brief Velocity Verlet Scheme with damping factor.
     *
     * Similar to boundary_of_center with condition_type=RIGID but damping factor is added.
     *
     * @param grid_density Density of the computational grid.
     * @param dx Grid step.
     * @param t_min Initial time. Have to be greater than 0.
     * @param t_max Stop time. Have to be greater than t_min.
     * @param dt Time step.
     * @param u 2D function of position in t=0. Only one argument (x).
     * @param beta Damping factor. If negative, probably oscillations will be increasing.
     * @param file Output file stream we will be writing resuts. May be open or closed.
     */
    void damped_oscillations(unsigned grid_density, double dx, double t_min, double t_max, double dt, initFunc u,
                             double beta, std::ofstream& file);

    /**
     * @brief Velocity Verlet Scheme with damping factor and driving force.
     *
     * Similar to damped_oscillations, but external force added. Also, all three: position, velocity and accelration
     * start from 0.
     *
     * @param grid_density Density of the computational grid.
     * @param dx Grid step.
     * @param t_min Initial time. Have to be greater than 0.
     * @param t_max Stop time. Have to be greater than t_min.
     * @param dt Time step.
     * @param beta Damping factor.
     * @param force External force function. Notice, angular velocity must be constant if included.
     * @param file Output file stream we will be writing resuts. May be open or closed.
     */
    void forced_oscillations(unsigned grid_density, double dx, double t_min, double t_max, double dt, double beta,
                             func2d force, std::ofstream &file);

    /**
     * @brief Velocity Verlet Scheme used for calculating average energy of steady state when resonance occur.
     *
     * To cause resonance, external force with changing angular velocity is needed. Complexity is basically
     * n * forced_oscillations, so calculations takes much more time.
     *
     * @param grid_density Density of the computational grid.
     * @param dx Grid step.
     * @param t_min Initial time. Have to be greater than 0.
     * @param t_max Stop time. Have to be greater than t_min.
     * @param dt Time step.
     * @param t_res Time when resonnace starts to occur. Have to be grater than t_min and less than t_max.
     * @param beta Damping factor.
     * @param om_min Initial angular velocity.
     * @param om_max Stop angular velocity.
     * @param dom Angular velocity step.
     * @param file Output file stream we will be writing resuts. May be open or closed.
     * @param force External force function. Notice, angular velocity must varry, so there are three input variables.
     */
    void resonances_enegry(unsigned grid_density, double dx, double t_min, double t_max, double dt, double t_res,
                           double beta, double om_min, double om_max, double dom, std::ofstream &file, func2dOm force);


    /**
     * @brief Problem type indicator for solver.
     * @var VERLET Most simple Verlet case.
     * @var BOC_R Boundary of center: rigid.
     * @var BOC_L Boundary of center: loose.
     * @var BOC_TC Boundary of center: twin center.
     * @var DAMPED Verlet with damped oscillations.
     * @var FORCED Verlet with damped and forced oscillations.
     */
    enum class WaveProblemType {VERLET, BOC_R, BOC_L, BOC_TC, DAMPED, FORCED};

    /**
     * @brief Solver for Verlet equasions.
     *
     * Not meant to be used by standard user - it's called from different functions.
     *
     * @param u_vals Position vector.
     * @param v_vals Velocity vector.
     * @param a_vals Acceleration vector.
     * @param t_min Initial time.
     * @param t_max Stop time.
     * @param dt Time step.
     * @param dx Grid (position) step.
     * @param type Type of solved problem.
     * @param ofile Open ofstream.
     * @param t_w Write time.
     * @param u_sol Position solution function
     * @param rho Density.
     * @param boundary Boundary position.
     * @param beta Damping factor.
     * @param force Force function.
     */
    void wave_solver(std::vector<double> &u_vals, std::vector<double> &v_vals, std::vector<double> &a_vals, double t_min,
                     double t_max, double dt, double dx, WaveProblemType type, std::ofstream &ofile, double t_w, func2d u_sol,
                     double rho, double boundary, double beta, func2d force);

}

/** @} End of group */

#endif
