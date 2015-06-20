/**
 * \file   vi_flow.hpp
 * \Author Bartłomiej Meder (bartem93@gmail.com)
 * \date   June, 2015
 * \brief  Viscous flow numerical methods
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

/**
 * \addtogroup imn
 * \{
 */

//* Engineering Numerical Methods functions
namespace imn {

    /**
     * Flux or vorticity fuction.
     * It's one of two functions we will use to solve the Navier-Stokes equasion.
     */
    typedef double(*func)(double, double);


    /**
     * \brief Simple 2D point representation.
     */
    struct Point {
        double x;
        double y;
    };

    /**
     * \brief Representation of 2D wall.
     *
     * Important info: this struct is used to represent walls of obstackle in
     * Navier-Stokes equasion. Currently only vertical and horizontal walls are
     * suppoerted.
     */
    struct Wall {
        /**
         * \brief Wall type
         * Values indicates, if wall is vertical/horizontal and internal/external.
         */
        enum WallType {
            UP, DOWN, LEFT, RIGHT
        };
        Point A;
        Point B;
        WallType type;
    };


    /**
     * \brief Navier-Stokes equasion solver
     *
     * It's recommended for user to not use this function, but one of two wrappers,
     * so no exact description is given.
     *
     * \param xmin Minimal horizontal value of computational grid.
     * \param xmax Maximal horizontal value of computational grid.
     * \param ymin Minimal vertical value of computational grid.
     * \param ymax Maximal vertical value of computational grid.
     * \param dx Horizontal axis step.
     * \param dy Vertical axis step.
     * \param flux Pointer to flux function.
     * \param vort Pointer to vorticity function.
     * \param file Reference to output file stream, both closed or opened are accepted.
     * \param prec Precision of calculations. Default to 1e-7.
     * \param min_iter Minimal loop iterations protection. Default to 500.
     * Values if first iterations doesn't change to much, so one can increase precision, or
     * add minimal number of iterations. Second option is usually more favorable if one
     * care about time.
     * \param obstackle Flag indicating if function is solving flow with obstackle or not. Default to false.
     * \param walls Pointer to vector of Wall pointers. Needed only if obstackle = true.
     * \param top Flag indicating if obstackle is on top or bottom of computational grid. Default to falne.
     * Currently, both of top and bottom obstackles are not supported. Well, I'll think of
     * addind this option in future.
     */
    void navier_stokes(double xmin, double xmax, double ymin, double ymax, double dx, double dy, func flux, func vort,
                       std::ofstream &file, double prec = 1e-7, int min_iter = 500, bool obstacle = false,
                       const std::vector<Wall *> *walls = 0, bool top = false);

    /**
     * \brief Poiseuille flow solver
     *
     * Poiseuille flow assumes no obstackles in pipe (our computational grid).
     * User specifies pipe dimensions, step and flux and vorticity functions (they have
     * to include parameters like density or viscosity, solver doesn't handle it).
     * In order to obtain results, function uses the finite difference method.
     * Function gives output to specified file in following manner:
     * x_coord y_coord x_velocity y_velocity
     * Next, they can be easily utilized by gnuplot or similar software.
     *
     * \param xmin Minimal horizontal value of computational grid.
     * \param xmax Maximal horizontal value of computational grid.
     * \param ymin Minimal vertical value of computational grid.
     * \param ymax Maximal vertical value of computational grid.
     * \param dx Horizontal axis step.
     * \param dy Vertical axis step.
     * \param flux Pointer to flux function.
     * \param vort Pointer to vorticity function.
     * \param file Reference to output file stream, both closed or opened are accepted.
     * \param prec Precision of calculations. Default to 1e-7.
     * \param min_iter Minimal loop iterations protection. Default to 200.
     * Values if first iterations doesn't change to much, so one can increase precision, or
     * add minimal number of iterations. Second option is usually more favorable if one
     * care about time.
     */
    void poiseuille(double xmin, double xmax, double ymin, double ymax, double dx, double dy, func flux, func vort,
                    std::ofstream &file, double prec = 1e-7, int min_iter = 200);

    /**
     * \brief Solver of Navier-Stokes equasion with obstackle in computational grid.
     *
     * User specifies pipe dimensions, step and flux and vorticity functions (they have
     * to include parameters like density or viscosity, solver doesn't handle it).
     * Also, a vector of Wall pointers representing obstackle must be given.
     * VERY IMPORTANT! When creating walls, A.x must be less or equal B.x, same goes for y.
     * In order to obtain results, function uses the finite difference method.
     * Function gives output to specified file in following manner:
     * x_coord y_coord x_velocity y_velocity
     * Next, they can be easily utilized by gnuplot or similar software.
     *
     * \param xmin Minimal horizontal value of computational grid.
     * \param xmax Maximal horizontal value of computational grid.
     * \param ymin Minimal vertical value of computational grid.
     * \param ymax Maximal vertical value of computational grid.
     * \param dx Horizontal axis step.
     * \param dy Vertical axis step.
     * \param flux Pointer to flux function.
     * \param vort Pointer to vorticity function.
     * \param walls Reference to vector of Wall pointers (obstackle).
     * \param file Reference to output file stream, both closed or opened are accepted.
     * \param prec Precision of calculations. Default to 1e-7.
     * \param min_iter Minimal loop iterations protection. Default to 500.
     * Values if first iterations doesn't change to much, so one can increase precision, or
     * add minimal number of iterations. Second option is usually more favorable if one
     * care about time.
     * \param top Flag indicating if obstackle is on top or bottom of computational grid. Default to falne.
     * Currently, both of top and bottom obstackles are not supported.
     */
    void ns_obstacle(double xmin, double xmax, double ymin, double ymax, double dx, double dy, func flux, func vort,
                     const std::vector<Wall *> &walls, std::ofstream &file, double prec = 1e-7, int min_iter = 500,
                     bool upper = false);

    /**
     * \biref Vorticity in-point value
     *
     * Function returns value of vorticity in given point of computational grid.
     * This value depend on both flux and vorticity grids. Not really useful for
     * standard user.
     *
     * \param f 2D flux computational grid (2D array)
     * \param v 2D vorticity computational grid (2D array)
     * \param x x location of given point in f and v arrays
     * \param y y location of given point in f and v arrays
     * \returns vorticity in v[x][y]
     */
    double vortHelper(double **f, double **v, long x, long y);

    /**
     * \biref Flux in-point value
     *
     * Function returns value of flux in given point of computational grid.
     * This value depend on both flux and vorticity grids, and grid steps.
     * Not really useful for standard user.
     *
     * \param f 2D flux computational grid (2D array).
     * \param v 2D vorticity computational grid (2D array).
     * \param x x location of given point in f and v arrays.
     * \param y y location of given point in f and v arrays.
     * \param dx x axis step.
     * \param dy y axis step.
     * \return Flux in f[x][y].
     */
    double fluxHelper(double **f, double **v, long x, long y, double dx, double dy);

    /**
     * \brief Position indicator
     *
     * Function checks, if given point is inside or outside obstackle, expressed by
     * vector of walls. Not really useful for standard user.
     *
     * \param x x coordinate of point.
     * \param y y coordinate of point.
     * \param dx x axis step. Used for casting.
     * \param dy y axis step. Used for casting.
     */
    bool inside(double x, double y, double dx, double dy, const std::vector<Wall *> *walls);

}

/** \} End of group */

#endif

