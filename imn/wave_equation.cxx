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

    wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, ProblemType::VERLET, file, t_write,
                u_sol, 1, 0, 0);

    if(!was_open)
        file.close();

}


void imn::boundary_of_center(unsigned grid_density, double dx, double t_min, double t_max, double dt, initFunc u,
                        std::ofstream &file, InitType condition_type, double density, double boundary_x)
{
    auto was_open = file.is_open();
    if(!was_open) {
        switch (condition_type){
            case InitType::RIGID:       file.open("boundary_rigid.dat"); break;
            case InitType::LOOSE:       file.open("boundary_loose.dat"); break;
            case InitType::TWIN_CENTER: file.open("boundary_twin.dat");  break;
        }
    }

    std::vector<double> u_vals(grid_density);
    std::vector<double> v_vals(grid_density, 0);
    std::vector<double> a_vals(grid_density, 0);

    for(auto i = 0u; i < grid_density; ++i){
        u_vals[i] = u(i*dx);
    }

    if(condition_type == InitType::RIGID || condition_type == InitType::TWIN_CENTER)
        u_vals[0] = u_vals[grid_density-1] = 0;

    if(condition_type == InitType::LOOSE){
        u_vals[0] = u_vals[1];
        u_vals[grid_density-1] = u_vals[grid_density-2];
    }

    switch (condition_type){
        case InitType::RIGID:
            wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, ProblemType::BOC_R, file, -1,
                        nullptr, density, boundary_x, 0); break;
        case InitType::LOOSE:
            wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, ProblemType::BOC_L, file, -1,
                        nullptr, density, boundary_x, 0); break;
        case InitType::TWIN_CENTER:
            wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, ProblemType::BOC_TC, file, -1,
                        nullptr, density, boundary_x, 0); break;

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

    wave_solver(u_vals, v_vals, a_vals, t_min, t_max, dt, dx, ProblemType::DAMPED, file, -1, nullptr, 1, 0, beta);

    if(!was_open)
        file.close();

}

void ::imn::wave_solver(std::vector<double> &u_vals, std::vector<double> &v_vals, std::vector<double> &a_vals, double t_min,
                 double t_max, double dt, double dx, ProblemType type, std::ofstream &ofile, double t_w, func2d u_sol,
                 double rho, double boundary, double beta)
{
    auto iter = 0l;
    auto size = u_vals.size();
    auto dens = 1.;

    std::vector<double> u_prev(size, 0);

    for(auto t = t_min; t <= t_max; t += dt, ++iter){

        if(type == ProblemType::DAMPED){
            u_prev = u_vals;
        }

        for(auto i = 0u; i < size; ++i){
            v_vals[i] = v_vals[i] + 0.5 * dt * a_vals[i];
            u_vals[i] = u_vals[i] + dt * v_vals[i];
        }

        if(type == ProblemType::BOC_L){
            u_vals[0] = u_vals[1];
            u_vals[size-1] = u_vals[size-2];
        }

        for(auto i = 1u; i < size-1; ++i) {

            if(type == ProblemType::BOC_TC)
                i*dx < boundary ? dens = 1 : dens = rho;

            a_vals[i] = (u_vals[i+1] + u_vals[i-1] - 2 * u_vals[i]) / (dens * dx * dx) - 2 * beta * ((u_vals[i] - u_prev[i]) / dt);
        }

        for(auto i = 0u; i < size; ++i)
            v_vals[i] = v_vals[i] + 0.5 * dt * a_vals[i];

        if(type == ProblemType::VERLET && t_w >= 0 && fabs(t_w - t) < 1e-9){

            for(auto i = 0u; i < size; ++i){

                ofile << i*dx << " " << u_vals[i];
                if(u_sol)
                    ofile << " " << u_sol(i*dx, t_w+dt);
                ofile << std::endl;

            }
            break;
        }

        else if(type != ProblemType::VERLET){

            for(auto i = 0u; i < size; ++i)
               ofile << t << " " << i*dx << " " << u_vals[i] << std::endl;
            ofile << std::endl;

        }

    }

    if(type == ProblemType::VERLET && t_w < 0){

        for(auto i = 0u; i < size; ++i){

            ofile << i*dx << " " << u_vals[i];
            if(u_sol)
                ofile << " " << u_sol(i*dx, t_w+dt);
            ofile << std::endl;

        }
    }

}

