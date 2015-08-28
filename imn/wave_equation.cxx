/**
 * @file   wave_equation.cxx
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   August, 2015
 * @brief  Function solving wave equations.
 *
 * Source for implementation of numerical methods created for solving wave equations.
 */
#include "wave_equation.hpp"

void imn::verlet_scheme(unsigned grid_density, double dx, double t_min, double t_max, double dt, func2d u, func2d v,
                   func2d a, std::ofstream &file, double t_write, func2d u_sol) {

    auto was_open = file.is_open();
    if(!was_open)
        t_write >= 0 ? file.open("verlet_" + imn::double_to_string(t_write) + ".dat") : file.open("verlet.dat");

    std::vector<double> u_vals(grid_density);
    std::vector<double> v_vals(grid_density);
    std::vector<double> a_vals(grid_density);

    for(auto i = 0u; i < grid_density; ++i){
        u_vals[i] = u(i*dx, 0);
        v_vals[i] = v(i*dx, 0);
        a_vals[i] = a(i*dx, 0);
    }

    wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, WaveProblemType::VERLET, file, t_write,
                u_sol, 1, 0, 0, nullptr);

    if(!was_open)
        file.close();

}


void imn::boundary_of_center(unsigned grid_density, double dx, double t_min, double t_max, double dt, initFunc u,
                        std::ofstream &file, WaveInitType condition_type, double density, double boundary_x)
{
    auto was_open = file.is_open();
    if(!was_open) {
        switch (condition_type){
            case WaveInitType::RIGID:       file.open("boundary_rigid.dat"); break;
            case WaveInitType::LOOSE:       file.open("boundary_loose.dat"); break;
            case WaveInitType::TWIN_CENTER: file.open("boundary_twin.dat");  break;
        }
    }

    std::vector<double> u_vals(grid_density);
    std::vector<double> v_vals(grid_density, 0);
    std::vector<double> a_vals(grid_density, 0);

    for(auto i = 0u; i < grid_density; ++i){
        u_vals[i] = u(i*dx);
    }

    if(condition_type == WaveInitType::RIGID || condition_type == WaveInitType::TWIN_CENTER)
        u_vals[0] = u_vals[grid_density-1] = 0;

    if(condition_type == WaveInitType::LOOSE){
        u_vals[0] = u_vals[1];
        u_vals[grid_density-1] = u_vals[grid_density-2];
    }

    switch (condition_type){
        case WaveInitType::RIGID:
            wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, WaveProblemType::BOC_R, file, -1,
                        nullptr, density, boundary_x, 0, nullptr); break;
        case WaveInitType::LOOSE:
            wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, WaveProblemType::BOC_L, file, -1,
                        nullptr, density, boundary_x, 0, nullptr); break;
        case WaveInitType::TWIN_CENTER:
            wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, WaveProblemType::BOC_TC, file, -1,
                        nullptr, density, boundary_x, 0, nullptr); break;

    }

    if(!was_open)
        file.close();
}


void ::imn::damped_oscillations(unsigned grid_density, double dx, double t_min, double t_max, double dt,
                                imn::initFunc u, double beta, std::ofstream &file)
{
    auto was_open = file.is_open();
    if(!was_open)
        file.open("damped_oscillations_" + imn::double_to_string(beta) + ".dat");

    std::vector<double> u_vals(grid_density);
    std::vector<double> v_vals(grid_density, 0);
    std::vector<double> a_vals(grid_density, 0);

    for(auto i = 0u; i < grid_density; ++i)
        u_vals[i] = u(i*dx);

    u_vals[0] = u_vals[grid_density-1] = 0;

    wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, WaveProblemType::DAMPED, file, -1,
                nullptr, 1, 0, beta, nullptr);

    if(!was_open)
        file.close();

}

void imn::forced_oscillations(unsigned grid_density, double dx, double t_min, double t_max, double dt, double beta,
                         func2d force, std::ofstream &file)
{
    auto was_open = file.is_open();
    if(!was_open)
        file.open("froced_oscillations_" + imn::double_to_string(beta) + ".dat");

    std::vector<double> u_vals(grid_density, 0);
    std::vector<double> v_vals(grid_density, 0);
    std::vector<double> a_vals(grid_density, 0);

    wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, WaveProblemType::FORCED, file, -1, nullptr, 1, 0, beta, force);


    if(!was_open)
        file.close();
}

void ::imn::resonances_enegry(unsigned grid_density, double dx, double t_min, double t_max, double dt, double t_res,
                       double beta, double om_min, double om_max, double dom, std::ofstream &file, func2dOm force)
{
    auto was_open = file.is_open();
    if(!was_open)
        file.open("resonances_" + imn::double_to_string(beta) + ".dat");

    std::vector<double> u_vals(grid_density);
    std::vector<double> v_vals(grid_density);
    std::vector<double> a_vals(grid_density);
    std::vector<double> u_prev(grid_density);

    for(auto omega = om_min; omega <= om_max; omega += dom){
        auto integral = 0.;

        u_vals.assign(grid_density, 0);
        v_vals.assign(grid_density, 0);
        a_vals.assign(grid_density, 0);

        for(auto t = t_min; t <= t_max; t += dt){

            u_prev = u_vals;

            for(auto i = 0u; i < grid_density; ++i){
                v_vals[i] = v_vals[i] + 0.5 * dt * a_vals[i];
                u_vals[i] = u_vals[i] + dt * v_vals[i];
            }

            for(auto i = 1u; i < grid_density-1; ++i)
                a_vals[i] = ((u_vals[i+1] + u_vals[i-1] - 2 * u_vals[i]) / (dx*dx))
                            - 2 * beta * ((u_vals[i] - u_prev[i]) / dt) + force(i*dx, t, omega);

            for(auto i = 0u; i < grid_density; ++i)
                v_vals[i] = v_vals[i] + 0.5 * dt * a_vals[i];

            if(t >= t_res){
                auto sum1 = 0.;
                auto sum2 = 0.;

                for(auto i = 0; i < grid_density-1; ++i){
                    sum1 += v_vals[i] * v_vals[i];
                    sum2 += ((u_vals[i+1] - u_vals[i]) / dx) * ((u_vals[i+1] - u_vals[i]) / dx);
                }
                sum1 += v_vals[grid_density-1] * v_vals[grid_density-1];
                integral += 0.5 * sum1 * dx + 0.5 * sum2 * dx;

            }

        }

        file << omega << " " << 0.25 * integral * dt << std::endl;

    }

    if(!was_open)
        file.close();
}

void ::imn::wave_solver(std::vector<double> &u_vals, std::vector<double> &v_vals, std::vector<double> &a_vals, double t_min,
                 double t_max, double dt, double dx, WaveProblemType type, std::ofstream &ofile, double t_w, func2d u_sol,
                 double rho, double boundary, double beta, func2d force)
{
    auto iter = 0l;
    auto size = u_vals.size();
    auto dens = 1.;

    auto force_f = force;

    if(!force_f)
        force_f = [](double x, double t){ return 0; };

    std::vector<double> u_prev(size, 0);

    for(auto t = t_min; t <= t_max; t += dt, ++iter){

        if(type == WaveProblemType::DAMPED || type == WaveProblemType::FORCED){
            u_prev = u_vals;
        }

        for(auto i = 0u; i < size; ++i){
            v_vals[i] = v_vals[i] + 0.5 * dt * a_vals[i];
            u_vals[i] = u_vals[i] + dt * v_vals[i];
        }

        if(type == WaveProblemType::BOC_L){
            u_vals[0] = u_vals[1];
            u_vals[size-1] = u_vals[size-2];
        }

        for(auto i = 1u; i < size-1; ++i) {

            if(type == WaveProblemType::BOC_TC)
                i*dx < boundary ? dens = 1 : dens = rho;

            a_vals[i] = (u_vals[i+1] + u_vals[i-1] - 2 * u_vals[i]) / (dens * dx * dx)
                        - 2 * beta * ((u_vals[i] - u_prev[i]) / dt)
                        + force_f(i*dx, t);
        }

        for(auto i = 0u; i < size; ++i)
            v_vals[i] = v_vals[i] + 0.5 * dt * a_vals[i];

        if(type == WaveProblemType::VERLET && t_w >= 0 && fabs(t_w - t) < 1e-9){

            for(auto i = 0u; i < size; ++i){

                ofile << i*dx << " " << u_vals[i];
                if(u_sol)
                    ofile << " " << u_sol(i*dx, t_w+dt);
                ofile << std::endl;

            }
            break;
        }

        else if(type != WaveProblemType::VERLET){

            for(auto i = 0u; i < size; ++i)
               ofile << t << " " << i*dx << " " << u_vals[i] << std::endl;
            ofile << std::endl;

        }

    }

    if(type == WaveProblemType::VERLET && t_w < 0){

        for(auto i = 0u; i < size; ++i){

            ofile << i*dx << " " << u_vals[i];
            if(u_sol)
                ofile << " " << u_sol(i*dx, t_w+dt);
            ofile << std::endl;

        }
    }

}

