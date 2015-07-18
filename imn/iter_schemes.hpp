/**
 * @file   diff_schemes.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   June, 2015
 * @brief  Differential equasions numerical methods
 *
 * Header for implementation of numerical methods created for solving partial elliptical differential equasions.       |
 * Methods incudes standard point relaxation and relaxation with grid denser in every iteration.
 */

#ifndef __IMN_ITER_SCHEMES__
#define __IMN_ITER_SCHEMES__
#include <cmath>
#include <fstream>
#include <string>
#include "grid.hpp"
#include "utils.hpp"

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn{

    /**
     * @brief Point relaxation, under-relaxation and over-relaxation of Poisson function.
     *
     * Funtion takes computational grid, relaxation factor and density funtion to calculate grid values for given
     * parameters. Stop condition is an integral, computed from grid values. Function creates default output as integral
     * of function of iteration. If one pases true as w_grid parameter, function will create new file with grid values
     * in gnuplot friendly format. May be used for optimum omega research.
     *
     * @param grd Computational grid
     * @param omega Relaxation factor
     * @param density 2D potential density function
     * @param file Output file stream we will write integral. May be open or closed.
     * If closed, it will be open with default name containing omega factor, and then closed.
     * @param tolerance Integral change tolerance. Lower value means more accurate and longer computations. Default to 1e-8.
     * @param w_grid Indicator telling function to write grid to file. Default to false.
     */
    void poisson_pr(Grid &grd, double omega, func2d density, std::ofstream &file, double tolerance = 1e-8,
                    bool w_grid = false);

    /**
     * @brief Relaxation with dense grid.
     *
     * Function similar to Point relaxation, but algorithm assumes calculating grid every k*dx or k*dy values instead
     * of dx or dy. Then, k is divided by two and grid is calculated once again, k is divided etc. That allows to decrease
     * general interation amount.
     *
     * @param grd Computational grid
     * @param omega Relaxation factor
     * @param k Initial grid step. Preferrable 2^n. It must be less than grid.x_size() and grid.y_size().
     * @param density 2D potential density function
     * @param file Output file stream we will write integral. May be open or closed.
     * If closed, it will be open with default name containing omega factor, and then closed.
     * @param tolerance Integral change tolerance. Lower value means more accurate and longer computations. Default to 1e-8.
     * @param w_grid Indicator telling function to write grid to file. Default to false.
     */
    void poisson_dens_grid(Grid &grd, double omega, unsigned k, func2d density, std::ofstream &file,
                           double tolerance = 1e-8, bool w_grid = false);
}

/** @} End of group */

#endif
