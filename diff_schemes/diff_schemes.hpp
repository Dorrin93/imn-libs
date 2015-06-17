#ifndef __DIFF_SCHEMES__
#define __DIFF_SCHEMES__

#include <fstream>
#include <cmath>
#include <string>
#include <iostream>

/**
 * \addtogroup imn
 * \{
 */

//* Engineering Numerical Methods functions
namespace imn{

/**
 * \typedef Input function
 * Differetial equasion we try to solve
 * First argument is time and second is some variable
 * Could also be its own derivative
 */
typedef double(*input)(double, double);
/**
 * \typedef Solution fuction
 * Funtion returning accurate solution of differential equasion
 */
typedef double(*solution)(double);
/**
 * \typedef Function for Runge-Kutta scheme
 * First argument is time, next two are variables
 */
typedef double(*rkfunc)(double, double, double);


/**
 * \brief Explicit Euler method
 *
 * Function writes results to file.
 * If file is open, function just writes results.
 * If not, function opens file, writes result to an "explicit_euler.dat" file
 * and closes file.
 * Format:
 * moment_in_time euler_result
 * If solution function given:
 * moment_in_time euler_result accurate_result difference
 *
 * \param inp_func Function we want to resolve
 * \param icondition Initial condition
 * \param dt Time step
 * \param begin Beginning point of interval
 * \param end Ending point of interval
 * \param file Output file stream to write results
 * \param sol_funct Accurate solution function. Default to 0 (NULL).
 *
 */
void explicit_euler(input inp_func, double icondition, double dt,
        double begin, double end, std::ofstream &file, 
        solution sol_func = 0);

/**
 * \brief Implicit Euler method
 *
 * Function writes results to file.
 * If file is open, function just writes results.
 * If not, function opens file, writes result to an "implicit_euler.dat" file
 * and closes file.
 * Format:
 * moment_in_time euler_result
 * If solution function given:
 * moment_in_time euler_result accurate_result difference
 *
 * \param inp_func Function we want to resolve
 * \param inp_deriv Derivaive of function we want to resolve
 * \param icondition Initial condition
 * \param dt Time step
 * \param begin Beginning point of interval
 * \param end Ending point of interval
 * \param file Output file stream to write results
 * \param sol_funct Accurate solution function. Default to 0 (NULL).
 *
 */
void implicit_euler(input inp_func, input inp_deriv, double icondition,
        double dt, double begin, double end, 
        std::ofstream &file, solution sol_func = 0);

/**
 * \brief Implicit Euler method with Newton iteration
 *
 * Function writes results to file.
 * If file is open, function just writes results.
 * If not, function opens file, writes result to an "newton_euler.dat" file
 * and closes file.
 * Format:
 * moment_in_time iteration_count euler_result
 * If solution function given:
 * moment_in_time iteration_count euler_result accurate_result difference
 *
 * \param inp_func Function we want to resolve
 * \param inp_deriv Derivaive of function we want to resolve
 * \param icondition Initial condition
 * \param dt Time step
 * \param begin Beginning point of interval
 * \param end Ending point of interval
 * \param file Output file stream to write results
 * \param stop Iteration stop condition (accuracy). Default to 1e-6.
 * \param sol_funct Accurate solution function. Default to 0 (NULL).
 *
 */
void newton_euler(input inp_func, input inp_deriv, double icondition,
        double dt, double begin, double end, 
        std::ofstream &file, double stop = 1e-6, solution sol_func = 0);

/**
 * \brief Explicit Euler method with Richardson extrapolation
 *
 * Function writes results to file.
 * If file is open, function just writes results.
 * If not, function opens file, writes result to an "richardson_euler.dat" file
 * and closes file.
 * Format:
 * moment_in_time euler_result
 * If solution function given:
 * moment_in_time euler_result accurate_result difference
 *
 * \param inp_func Function we want to resolve
 * \param icondition Initial condition
 * \param dt Time step
 * \param begin Beginning point of interval
 * \param end Ending point of interval
 * \param file Output file stream to write results
 * \param sol_funct Accurate solution function. Default to 0 (NULL).
 *
 */
void richardson_euler(input inp_func, double icondition, double dt,
        double begin, double end, std::ofstream &file, 
        solution sol_func = 0);

/**
 * \brief Explicit Euler method with Richardson extrapolation and 
 * automatic time step selection algorithm
 *
 * Function writes results to file.
 * If file is open, function just writes results.
 * If not, function opens file, writes result to an "autostep_euler.dat" file
 * and closes file.
 * Format:
 * moment_in_time time_step euler_result
 * If solution function given:
 * moment_in_time time_step euler_result accurate_result difference
 *
 * \param inp_func Function we want to resolve
 * \param icondition Initial condition
 * \param init_dt Initial time step. If too high, file will contain only
 * initial condition value.
 * \param begin Beginning point of interval
 * \param end Ending point of interval
 * \param file Output file stream to write results
 * \param tolerance Truncation error tolerance. Default to 1e-5.
 * \param S Safety parameter. Must be S < 1. Default to 0.75.
 * \param sol_funct Accurate solution function. Default to 0 (NULL).
 *
 */
void autostep_euler(input inp_func, double icondition, double init_dt,
        double begin, double end, std::ofstream &file,
        double tolerance = 1e-5, double S = 0.75, solution sol_func = 0);

/**
 * \brief Explicit Runge-Kutta 4 scheme
 *
 * Function writes results to file.
 * If file is open, function just writes results.
 * If not, function opens file, writes result to an "rk4.dat" file
 * and closes file.
 * Format: 
 * moment_in_time x_value y_value
 *
 * \param equasion1 First function to solve
 * \param equasion2 Second function to solve
 * \param cond_x Initial condition for x (1) variable
 * \param cond_y Initial condition for y (2) variable
 * \param dt Time step
 * \param begin Beginning point of interval
 * \param end Ending point of interval
 * \param file Output file stream to write results
 *
 */
void rk4(rkfunc equasion1, rkfunc equasion2, double cond_x, double cond_y,
        double dt, double begin, double end, std::ofstream &file);

}

/** \} End of group */

#endif

