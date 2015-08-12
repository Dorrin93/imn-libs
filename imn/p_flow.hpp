/**
 * @file   p_flow.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Potential flow numerical methods.
 *
 * Header for implementation of numerical methods created for solving potential flow of liquids (Laplace equasion).
 */
#ifndef __IMN_P_FLOW__
#define __IMN_P_FLOW__
#include "ogrid.hpp"
#include "utils.hpp"
#include <fstream>

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn {

    /**
     * @brief Function for solving potetntial flow using stream function.
     *
     * Function takes obstacle grid object and iterates through it, until difference of integrals is smaller than
     * eps. Integrals are written to file through ofstream. Function also can write grid representation to separate
     * file in gnuplot friendly format.
     *
     * @param grd Computation grid with obstacle.
     * @param u0 Initial x-axis velocity.
     * @param file Output file stream function will be writing to. Could be open or closed, if open user will have to close
     * it by himself.
     * @param eps Integral change tolerance. Lower value means more accurate and longer computations. Default to 1e-6.
     * @param w_grid Writing grid to file indicator. Default to true.
     */
    void stream_lines(OGrid& grd, double u0, std::ofstream& file, double eps = 1e-6, bool w_grid = true);

    /**
     * @brief Function for solving potetntial flow using relaxation by dyscretization.
     *
     * Function takes obstacle grid object and iterates through it, until difference of integrals is smaller than
     * eps. Integrals are written to file through ofstream. Function also can write grid potential lines to separate
     * file in gnuplot friendly format.
     *
     * @param grd Computation grid with obstacle.
     * @param u0 Initial x-axis velocity.
     * @param file Output file stream function will be writing to. Could be open or closed, if open user will have to close
     * it by himself.
     * @param eps Integral change tolerance. Lower value means more accurate and longer computations. Default to 1e-6.
     * @param w_grid Writing grid to file indicator. Default to true.
     */
    void potential_lines(OGrid& grd, double u0, std::ofstream& file, double eps = 1e-6, bool w_grid = true);

    /**
     * @brief Base function, previous two are wrappers to it. User should use one of previous two.
     * @param grd Computation grid with obstacle.
     * @param u0 Initial x-axis velocity.
     * @param file Output file stream function will be writing to.
     * @param eps Integral change tolerance. Lower value means more accurate and longer computations. Default to 1e-6.
     * @param w_grid Writing grid to file indicator. Default to true.
     * @param stream Method usage indicator.
     */
    void base_lines(OGrid& grd, double u0, std::ofstream& file, double eps, bool w_grid, bool stream);

}

/** @} End of group */

#endif
