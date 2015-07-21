#include<array>
#include<memory>
#include<iostream>
#define __LAB3__

#ifdef __LAB1__
#include "imn/diff_schemes.hpp"
#endif

#ifdef __LAB2__
#include "imn/iter_schemes.hpp"
#endif

#ifdef __LAB3__
#include "imn/p_flow.hpp"
#endif

int main(){

    std::ofstream happiness;

    #ifdef __LAB1__
    // differential equasion and its solution (for comparsion)
    auto eufunc = [](double t, double u){ return u*cos(t); };
    auto solu = [](double t){ return exp(sin(t)); };

    // if your derivative don't need 'u', you have to do something like that
    auto pack_cos = [](double t, double u){ return cos(t); };

    // rk4 functions with The Beatles reference
    auto warm = [](double t, double paul, double john){ return john; };
    auto gun = [](double t, double paul, double john){ return -paul; };

    imn::explicit_euler(eufunc, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::implicit_euler(eufunc, pack_cos, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::newton_euler(eufunc, pack_cos, 1, M_PI/400., 0, 4.*M_PI, happiness, 1e-6, solu);
    imn::richardson_euler(eufunc, 1, M_PI/400., 0, 4.*M_PI, happiness, solu);
    imn::autostep_euler(eufunc, 1, 5.*M_PI, 0, 4.*M_PI, happiness, 1e-5, 0.75, solu);
    imn::rk4(warm, gun, 0, 1, 0.1, 0, 4.*M_PI, happiness);
    #endif

    #ifdef __LAB2__
    // first we need to create grid
    // let's say it's grounded, so we don't apply potential function
    auto grid = std::make_unique<imn::Grid>(-129, 129, -129, 129, 2., 2.);

    // then we need to define potential density
    // it's prety complicated, thus we will use two lambdas
    const auto rho1 = [](double x, double y){ return (1./(900.*M_PI))*exp(-(x/30.)*(x/30.)-(y/30.)*(y/30.)); };
    const auto pdens = [&rho1](double x, double y) { return rho1(x - 30, y - 30) + rho1(x + 30, y + 30); };

    // we want to find optimal omega value, so we will iterate through some candidates
    std::array<double, 6> omegas{0.75, 0.95, 1.2, 1.5, 1.95, 1.99};

    for(auto o : omegas){
        imn::poisson_pr(*grid, o, pdens, happiness);
        grid->clear();
    }

    imn::poisson_dens_grid(*grid, 1, 16, pdens, happiness, 1e-8, true);
    #endif

    #ifdef __LAB3__
    // creating grid with obstacle possibilty
    auto grid = std::make_unique<imn::OGrid>(0., 200., 0., 100., 1., 1.);

    // setting obstacle - function taking doubles and returning bool, could be also lambda
    grid->set_obstackle([](double x, double y){ return x >= 80 && x <= 120 && y >= 20 && y <= 80; } );

    imn::stream_lines(*grid, 1., happiness);

    #endif

    return 0;
}
