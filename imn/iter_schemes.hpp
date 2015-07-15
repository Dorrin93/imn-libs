/**
 * \file   diff_schemes.hpp
 * \Author Bartłomiej Meder (bartem93@gmail.com)
 * \date   June, 2015
 * \brief  Differential equasions numerical methods
 *
 * Header for implementation of numerical methods created for solving partial elliptical differential equasions.       |
 * Methods incudes TODO
 */

#ifndef __IMN_ITER_SCHEMES__
#define __IMN_ITER_SCHEMES__
#include <cmath>
#include <fstream>
#include <string>
#include "grid.hpp"
#include "utils.hpp"

/**
 * \addtogroup imn
 * \{
 */

//* Engineering Numerical Methods functions
namespace imn{

    void poisson_pr(Grid &grd, double omega, func2d density, std::ofstream &file, double tolerance = 1e-8,
                    bool w_grid = false);

    void poisson_dens_grid(Grid &grd, double omega, unsigned k, func2d density, std::ofstream &file,
                           double tolerance = 1e-8, bool w_grid = false);
}

/** \} End of group */

#endif
