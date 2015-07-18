/**
 * @file   diff_schemes.cxx
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   June, 2015
 * @brief  Differential equasions numerical methods
 *
 * Source for implementation of numerical methods created for solving differential equasions.
 * Methods incudes explicit and impicit Euler methods and Runge-Kutta 4 scheme.
 */

#include "diff_schemes.hpp"

/////////////////////////////////////////////////////////////////////////////////////
void imn::explicit_euler(input inp_func, double icondition, double dt, double begin, double end,
                         std::ofstream &file, solution sol_func)
{
    // checking difference between solution(begin) and condition,
    // if more then about 0 user probably uses bad input
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){

        std::cerr << "Difference between solution function and initial condition is "
                  << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect." << std::endl;

    }

    // checking if file is was open and eventually opening it
    auto was_open = file.is_open();
    if(!was_open)
        file.open("explicit_euler.dat");

    file << "# t u_n(t)";

    if(sol_func)
        file << " u0(t) u0(t)-u_n(t)";

    file << std::endl;

    file << begin << " " << icondition;
    if(sol_func)
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    auto u_n = icondition;

    for(auto i = begin + dt; i <= end; i += dt){

        u_n = u_n + dt * inp_func(i - dt, u_n);

        file << i << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void imn::implicit_euler(input inp_func, input inp_deriv, double icondition,
        double dt, double begin, double end, std::ofstream &file,
        solution sol_func)
{
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    auto was_open = file.is_open();
    if(!was_open)
        file.open("implicit_euler.dat");

    file << "# t u_n(t)";
    if(sol_func)
        file << " u0(t) u0(t)-u_n(t)";
    file << std::endl;

    file << begin << " " << icondition;
    if(sol_func)
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    auto u_n = icondition;

    for(auto i = begin + dt; i <= end; i += dt){

        u_n = u_n / (1. - dt * inp_deriv(i, u_n));

        file << i << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void imn::newton_euler(input inp_func, input inp_deriv, double icondition,
        double dt, double begin, double end,
        std::ofstream &file, double stop, solution sol_func)
{
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    auto was_open = file.is_open();
    if(!was_open)
        file.open("newton_euler.dat");

    file << "# t iters u_n(t)";
    if(sol_func)
        file << " u0(t) u0(t)-u_n(t)";
    file << std::endl;


    file << begin << " " << 0 << " " << icondition;
    if(sol_func)
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    auto u_n = icondition;
    // value from external iteration before
    auto u_pn = u_n;
    // value from internal iteration before
    auto u_wn = u_n;
    auto iter = 0;

    //TD<decltype(u_wn)> u_wnType; // compiler says it's double :3

    for(auto i = begin + dt; i <= end; i += dt, iter = 0, u_pn = u_n){

        while(stop < fabs(u_n - u_wn) || !iter){
            u_wn = u_n;
            u_n = u_wn - (u_wn - u_pn - dt * inp_func(i, u_wn)) / (1. - dt * inp_deriv(i, u_wn));

            ++iter;
        }

        file << i << " " << iter << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void imn::richardson_euler(input inp_func, double icondition, double dt,
        double begin, double end, std::ofstream &file, solution sol_func)
{
    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    auto was_open = file.is_open();
    if(!was_open)
        file.open("richardson_euler.dat");

    file << "# t u_n(t)";
    if(sol_func)
        file << " u0(t) u0(t)-u_n(t)";
    file << std::endl;


    file << begin << " " << icondition;
    if(sol_func)
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    auto u_n = icondition;

    for(auto i = begin + 2.*dt; i <= end; i += 2.*dt){

        // first short step
        auto u_1 = u_n + dt * inp_func(i - 2. * dt, u_n);
        // second short step
        auto u_2 = u_1 + dt * inp_func(i - dt, u_1);
        // long step
        auto u_l = u_n + 2. * dt * inp_func(i - 2. * dt, u_n);

        u_n = 2. * u_2 - u_l;

        file << i << " " << u_n;
        if(sol_func)
            file << " " << sol_func(i) << " " << sol_func(i) - u_n;
        file << std::endl;

    }

    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void imn::autostep_euler(input inp_func, double icondition, double init_dt,
        double begin, double end, std::ofstream &file,
        double tolerance, double S, solution sol_func)
{
    if(S >= 1.){
        std::cerr << "S parameter must be S < 1 !" << std::endl;
        return;
    }

    if(sol_func && fabs(sol_func(begin) - icondition) > 1e-9){
        std::cerr << "Difference between solution function and initial \
            condition is " << sol_func(begin) - icondition << std::endl;
        std::cerr << "It should be about 0, so input could be incorrect."
            << std::endl;
    }

    auto was_open = file.is_open();
    if(!was_open)
        file.open("autostep_euler.dat");

    file << "# t dt u_n(t)";
    if(sol_func)
        file << " u0(t) u0(t)-u_n(t)";
    file << std::endl;


    file << begin << " " << init_dt << " " << icondition;
    if(sol_func)
        file << " " << sol_func(begin) << " " << sol_func(begin) - icondition;
    file << std::endl;

    auto u_n = icondition;
    auto dt = init_dt;
    auto cnt = 0;

    for(auto i = begin + 2*dt; i <= end || cnt < 5; i += 2*dt, ++cnt){

        auto u_1 = u_n + dt * inp_func(i - 2 * dt, u_n);
        auto u_2 = u_1 + dt * inp_func(i - dt, u_1);
        auto u_l = u_n + 2 * dt * inp_func(i - 2 * dt, u_n);

        if(fabs(u_2 - u_l) < tolerance){
            u_n = u_2;

            file << i << " " << dt << " " << u_n;
            if(sol_func)
                file << " " << sol_func(i) << " " << sol_func(i) - u_n;
            file << std::endl;
        }
        else{
            // getting back in time :3
            i -= 2 * dt;
        }

        dt = sqrt( (S * tolerance) / fabs(u_2 - u_l) ) * dt;

    }

    if(!was_open)
        file.close();

}

/////////////////////////////////////////////////////////////////////////////////////
void imn::rk4(rkfunc equasion1, rkfunc equasion2, double cond_x, double cond_y,
        double dt, double begin, double end, std::ofstream &file)
{
    auto was_open = file.is_open();
    if(!was_open)
        file.open("rk4.dat");

    file << "# t x y(x)" << std::endl;

    file << begin << " " << cond_x << " " << cond_y << std::endl;

    auto x = cond_x;
    auto y = cond_y;


    for(auto i = begin + dt; i <= end; i += dt){

        auto k1x = equasion1(i - dt, x, y);
        auto k1y = equasion2(i - dt, x, y);

        auto k2x = equasion1(i - dt + dt / 2, x + dt * k1x / 2, y + dt * k1y / 2);
        auto k2y = equasion2(i - dt + dt / 2, x + dt * k1x / 2, y + dt * k1y / 2);

        auto k3x = equasion1(i - dt + dt / 2, x + dt * k2x / 2, y + dt * k2y / 2);
        auto k3y = equasion2(i - dt + dt / 2, x + dt * k2x / 2, y + dt * k2y / 2);

        auto k4x = equasion1(i, x + dt * k3x, y + dt * k3y);
        auto k4y = equasion2(i, x + dt * k3x, y + dt * k3y);

        x = x + (dt / 6.) * (k1x + 2 * k2x + 2 * k3x + k4x);
        y = y + (dt / 6.) * (k1y + 2 * k2y + 2 * k3y + k4y);

        file << i << " " << x << " " << y << std::endl;

    }

    if(!was_open)
        file.close();
}
