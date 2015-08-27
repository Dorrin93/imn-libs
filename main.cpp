#include<array>
#include<memory>
#include<iostream>
#define __LAB7__

#ifdef __LAB1__
#include "imn/diff_schemes.hpp"
#endif

#ifdef __LAB2__
#include "imn/iter_schemes.hpp"
#endif

#ifdef __LAB3__
#include "imn/p_flow.hpp"
#endif

#ifdef __LAB4__
#include "imn/vi_flow.hpp"
#endif

#ifdef __LAB7__
#include "wave_equation.hpp"
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

    for(const auto& o : omegas)
        imn::poisson_pr(*grid, o, pdens, happiness);

    imn::poisson_dens_grid(*grid, 1, 16, pdens, happiness, 1e-8, true);
    #endif

    #ifdef __LAB3__
    // creating grid with obstacle possibilty
    auto grid = std::make_unique<imn::OGrid>(0., 200., 0., 100., 1., 1.);

    // setting obstacle - function taking doubles and returning bool, could be also lambda
    grid->set_obstackle([](double x, double y){ return x >= 80 && x <= 120 && y >= 20 && y <= 80; } );

    imn::stream_lines(*grid, 1., happiness);

    imn::potential_lines(*grid, 1., happiness);

    #endif

    #ifdef __LAB4__
    // creating two grids - both needs to have same dimensions
    auto fgrid = std::make_unique<imn::OGrid>(-1., 2., -0.5, 0.5, 0.01, 0.01);
    auto vgrid = std::make_unique<imn::OGrid>(-1., 2., -0.5, 0.5, 0.01, 0.01);

    // setting obstacle - function taking doubles and returning bool, could be also lambda, it must be set to both grids
    auto obs = [](double x, double y){ return x >= -0.1 && x <= 0.1 && y >= -0.5 && y <= 0; };
    fgrid->set_obstackle(obs);
    vgrid->set_obstackle(obs);

    // creating functions defining flux and vorticity
    auto ffunc = [](double x, double y){ return -0.5 * ((y*y*y/3) - 0.5*0.5*y ); };
    auto vfunc = [](double x, double y){ return -0.5 * 2*y; };

    imn::poiseuille(*fgrid, *vgrid, ffunc, vfunc, happiness);

    imn::ns_obstacle(*fgrid, *vgrid, ffunc, vfunc, happiness);

    #endif

    #ifdef __LAB7__
    // first we need to define functions of posiztion, velocity and acceleration
    auto u = [](double x, double t){ return cos(M_PI*t) * sin(M_PI*x) - 0.5 * cos(2*M_PI*t) * sin(2*M_PI*x); };
    auto v = [](double x, double t){ return -M_PI * sin(M_PI*t) * sin(M_PI*x) + M_PI * sin(2*M_PI*t) * sin(2*M_PI*x); };
    auto a = [](double x, double t){ return -M_PI*M_PI * cos(M_PI*t)*sin(M_PI*x) + 2*M_PI*M_PI * cos(2*M_PI*t)*sin(2*M_PI*x); };

    // we will check solution in few time stamps
    std::array<double, 9> stamps{0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2};

    for(auto& it: stamps)
        imn::verlet_scheme(101, 0.01, 0, 2, 0.005, u, v, a, happiness, it, u);

    // for next function we need to define initial position function
    auto u0 = [](double x) { return exp(-100 * (x - 0.5) * (x - 0.5)); };

    imn::boundary_of_center(101, 0.01, 0, 2, 0.005, u0, happiness, imn::InitType::RIGID);
    imn::boundary_of_center(101, 0.01, 0, 2, 0.005, u0, happiness, imn::InitType::LOOSE);
    imn::boundary_of_center(101, 0.01, 0, 2, 0.005, u0, happiness, imn::InitType::TWIN_CENTER, 3., 0.75);

    // we will check damped oscillation for three betas
    std::array<double, 3> betas{0.2, 1.0, 3.0};
    for(auto& it : betas)
        imn::damped_oscillations(101, 0.01, 0, 4, 0.005, u0, it, happiness);

    #endif

    return 0;
}
