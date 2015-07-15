/**
 * \file   iter_schemes.hpp
 * \Author Bartłomiej Meder (bartem93@gmail.com)
 * \date   July, 2015
 * \brief  Differential equasions numerical methods
 *
 * Source for implementation of numerical methods created for solving partial elliptical differential equasions.       |
 * Methods incudes TODO
 */

#include "iter_schemes.hpp"
#include <iostream>

auto tt = [](double n){return n*n;};

auto fdm = [](imn::Grid& grd, int x, int y, unsigned k, imn::func2d dens)
{
    return 0.25 * (grd(x+k, y) + grd(x-k, y) + grd(x, y+k) + grd(x, y-k) +
           k*k*grd.dx()*grd.dy()*dens(grd.x_min() + grd.dx()*x, grd.y_min() + grd.dy()*y));
};

auto integ = [](imn::Grid& grd, int x, int y, unsigned k, imn::func2d dens)
{
    return 0.5 * (tt(grd(x+k, y)-grd(x, y)) + tt(grd(x, y+k)-grd(x,y))) -
           k*k*grd.dx()*grd.dy()*dens(grd.x_min() + x*grd.dx(), grd.y_min() + y*grd.dy())*grd(x, y);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void imn::poisson_pr(Grid &grd, double omega, func2d density, std::ofstream &file, double tolerance, bool w_grid)
{
    auto o_s = double_to_string(omega);
    // checking if file is was open and eventually opening it
    auto was_open = file.is_open();
    if(!was_open)
        file.open("poisson_pr_" + o_s + "_integral.dat");

    file << "# n a" << std::endl;

    auto iter = 1L;

    // integral values
    auto a = 0.;
    auto a_bf = 1.;

    for(; fabs(a - a_bf) >= tolerance; ++iter) {
        a_bf = a; a = 0;

        for(std::size_t i = 1; i < grd.x_size() - 1; ++i){
            for(std::size_t j = 1; j < grd.y_size() - 1; ++j){
                grd(i, j) = (1-omega) * grd(i, j) + omega*fdm(grd, i, j, 1, density);
            }
        }

        for(std::size_t i = 0; i < grd.x_size() - 1; ++i){
            for(std::size_t j = 0; j < grd.y_size() - 1; ++j){
                a += integ(grd, i, j, 1, density);
            }
        }


        file << iter << " " << a << std::endl;
    }

    if(w_grid){
        std::ofstream grd_file;



        grd_file.open("grid_pr_" + o_s + ".dat");

        grd_file << "# x y phi(x, y)" << std::endl;

        grd_file << grd;

        grd_file.close();

    }

    if(!was_open)
        file.close();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void imn::poisson_dens_grid(Grid &grd, double omega, unsigned k, func2d density, std::ofstream &file,
                            double tolerance, bool w_grid)
{
    auto o_s = double_to_string(omega);

    auto was_open = file.is_open();
    if(!was_open)
        file.open("poisson_dens_" + o_s + "_integral.dat");

    file << "# n a" << std::endl;

    auto iter = 1L;

    for(; k > 0; k /= 2){
        auto a = 0.;
        auto a_bf = 1.;

        for(; fabs(a - a_bf) >= tolerance; ++iter) {
            a_bf = a; a = 0;

            for(std::size_t i = k; i < grd.x_size() - k; i += k){
                for(std::size_t j = k; j < grd.y_size() - k; j += k){
                    grd(i, j) = (1-omega) * grd(i, j) + omega*fdm(grd, i, j, k, density);
                }
            }

            for(std::size_t i = 0; i < grd.x_size() - k; i += k){
                for(std::size_t j = 0; j < grd.y_size() - k; j += k){
                    a += integ(grd, i, j, k, density);
                }
            }

            file << iter << " " << a << std::endl;
        }

        if(w_grid){
            std::ofstream grd_file;

            grd_file.open("grid_dens_" + o_s + "_" + std::to_string(k) + ".dat");

            grd_file << "# x y phi(x, y)" << std::endl;

            grd.write_to_file(grd_file, false);

            grd_file.close();

        }

        // grid densification
        for(std::size_t i = k; i < grd.x_size() - k; i += k){
            for(std::size_t j = k/2; j < grd.y_size() - k/2; j += k){
                grd(i, j) = ( grd(i, j-k/2) + grd(i, j+k/2) ) / 2.;
            }
        }
        for(std::size_t i = k/2; i < grd.x_size() - k/2; i += k){
            for(std::size_t j = k; j < grd.y_size() - k; j += k){
                grd(i, j) = ( grd(i-k/2, j) + grd(i+k/2, j) ) / 2.;
            }
        }
        for(std::size_t i = k/2; i < grd.x_size() - k/2; i += k){
            for(std::size_t j = k/2; j < grd.y_size() - k/2; j += k){
                grd(i, j) = 0.25 * ( grd(i-k/2, j-k/2) + grd(i+k/2, j-k/2) + grd(i-k/2, j+k/2) + grd(i+k/2, j+k/2) );
            }
        }

    }

    if(!was_open)
        file.close();
}
