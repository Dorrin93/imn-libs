/**
 * @file   vi_flow.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   June, 2015
 * @brief  Viscous flow numerical methods
 *
 * Header for implementation of numerical methods created for solving viscous flow of liquid.
 * Everything is based on solving Navier-Stokes equasion.
 */

#ifndef __VI_FLUX__
#define __VI_FLUX__

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "ogrid.hpp"

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn {

    /**
     * @brief Navier-Stokes equasion solver
     *
     * It's recommended for user to not use this function, but one of two wrappers,
     * so no exact description is given.
     *
     * @param flux_grid Reference to previously created Obstacle Grid object. It'll hold values of aproximated flux.
     * @param vort_grid Reference to previously created Obstacle Grid object. It'll hold values of aproximated vorticity.
     * @param flux Function defining flux of liquid.
     * @param vort Function defining vorticity of liquid.
     * @param file Reference to output file stream, both closed or opened are accepted.
     * @param prec Precision of calculations. Default to 1e-7.
     * @param min_iter Minimal loop iterations protection. Default to 500.
     * Values if first iterations doesn't change to much, so one can increase precision, or
     * add minimal number of iterations. Second option is usually more favorable if one
     * cares about time.
     * @param obstackle Flag indicating if function is solving flow with obstackle or not.
     */
    void navier_stokes(OGrid& flux_grid, OGrid& vort_grid, func2d flux, func2d vort, std::ofstream &file, double prec,
                       int min_iter, bool ns_obs);

    /**
     * @brief Poiseuille flow solver
     *
     * Solver uses Navier-Stokes equasion to solve simple flow of liquid in pipe without any oblstacle.
     *
     * @param flux_grid Reference to previously created Obstacle Grid object. It must not have obstacle function defined.
     * @param vort_grid Reference to previously created Obstacle Grid object. It must not have obstacle function defined.
     * @param flux Function defining flux of liquid.
     * @param vort Function defining vorticity of liquid.
     * @param file Reference to output file stream, both closed or opened are accepted.
     * @param prec Precision of calculations. Default to 1e-7.
     * @param min_iter Minimal loop iterations protection. Default to 200.
     * Values if first iterations doesn't change to much, so one can increase precision, or
     * add minimal number of iterations. Second option is usually more favorable if one
     * cares about time.
     */
    void poiseuille(OGrid& flux_grid, OGrid& vort_grid, func2d flux, func2d vort, std::ofstream &file,
                    double prec = 1e-7, int min_iter = 200);

    /**
     * @brief Solver of Navier-Stokes equasion with obstackle in computational grid.
     *
     * Solver uses Navier-Stokes equasion to solve complex flow of liquid in pipe with some obstacle. Right now, only
     * one sided (upper wall, lower wall, middle of pipe) obstacle is supported, defining two or three may cause
     * undefined behaviour.
     *
     * @param flux_grid Reference to previously created Obstacle Grid object. It must have obstacle function defined.
     * @param vort_grid Reference to previously created Obstacle Grid object. It must have obstacle function defined.
     * @param flux Function defining flux of liquid. It does not have any data about obstacle.
     * @param vort Function defining vorticity of liquid. It does not have any data about obstacle.
     * @param file Reference to output file stream, both closed or opened are accepted.
     * @param prec Precision of calculations. Default to 1e-7.
     * @param min_iter Minimal loop iterations protection. Default to 500.
     * Values if first iterations doesn't change to much, so one can increase precision, or
     * add minimal number of iterations. Second option is usually more favorable if one
     * cares about time.
     */
    void ns_obstacle(OGrid& flux_grid, OGrid& vort_grid, func2d flux, func2d vort, std::ofstream &file,
                     double prec = 1e-7, int min_iter = 500);

}

/** @} End of group */

#endif

