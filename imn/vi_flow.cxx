/**
 * @file   vi_flow.cxx
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   June, 2015
 * @brief  Viscous flow numerical methods
 *
 * Source for implementation of numerical methods created for solving viscous flow of liquid.
 * Everything is based on solving Navier-Stokes equasion.
 */

#include "vi_flow.hpp"

auto vortHelper = [](const imn::OGrid& f, const imn::OGrid& v, unsigned x, unsigned y)
{
    // 0.0625 = 1/16
    return 0.25 * (v(x+1, y) + v(x-1, y) + v(x, y-1) + v(x, y+1)) -
           0.0625 * ( (f(x, y+1) - f(x, y-1)) * (v(x+1, y) - v(x-1, y)) -
           (f(x+1, y) - f(x-1, y)) * (v(x, y+1) - v(x, y-1)) );
};


auto fluxHelper = [](const imn::OGrid& f, const imn::OGrid& v, unsigned x, unsigned y)
{
    return (f(x+1, y) + f(x-1, y) + f(x, y+1) + f(x, y-1) - v(x, y) *
            v.dx() * v.dy()) * 0.25;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imn::poiseuille(OGrid& flux_grid, OGrid& vort_grid, func2d flux, func2d vort, std::ofstream &file, double prec,
                     int min_iter)
{
    // actually it's only wrapper
    imn::navier_stokes(flux_grid, vort_grid, flux, vort, file, prec, min_iter, false);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imn::ns_obstacle(OGrid& flux_grid, OGrid& vort_grid, func2d flux, func2d vort, std::ofstream &file, double prec,
                      int min_iter)
{
    imn::navier_stokes(flux_grid, vort_grid, flux, vort, file, prec, min_iter, true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imn::navier_stokes(OGrid& flux_grid, OGrid& vort_grid, func2d flux, func2d vort, std::ofstream &file, double prec, int min_iter,
                        bool ns_obs)
{
    flux_grid.clear();
    vort_grid.clear();

    if(ns_obs){
        flux_grid.apply_to_single(flux, Grid::Edge::DOWN);
        flux_grid.apply_to_single(flux, Grid::Edge::UP);
        auto ymin = flux_grid.y_min();
        auto ymax = flux_grid.y_max();
        switch(flux_grid.obstacle_type()){

            case OGrid::Obstype::DOWN:
                flux_grid.apply_to_obstackle(
                        [&flux, &flux_grid, &ymin](double x, double y){ return flux(x, ymin); }
                );
                break;
            case OGrid::Obstype::UP:
                flux_grid.apply_to_obstackle(
                        [&flux, &flux_grid, &ymax](double x, double y){ return flux(x, ymax); }
                );
                break;
            case OGrid::Obstype::MID:
                flux_grid.apply_to_obstackle(
                        [&flux, &flux_grid, &ymin, &ymax](double x, double y){ return flux(x, 0.5 * (ymin + ymax)); }
                );
                break;
        }
        flux_grid.apply_to_edges(flux);
        flux_grid.apply_not_to_obstackle(flux);
        vort_grid.apply_to_edges(vort);
        vort_grid.apply_not_to_obstackle(vort);
        auto zero = [](double x, double y) { return 0; };
        vort_grid.apply_to_single(zero, Grid::Edge::UP);
        vort_grid.apply_to_single(zero, Grid::Edge::DOWN);
    }
    else{
        flux_grid.apply_to_edges(flux);
        vort_grid.apply_to_edges(vort);
    }

    // some special values - we check condition in 2/3 x and 0.5 y point
    auto p_x = static_cast<unsigned>(2. * flux_grid.x_size() / 3.);
    auto p_y = flux_grid.y_size() >> 1;
    auto prev_f = 0.;
    auto prev_v = 0.;
    auto iter = 0L;

    auto obs = flux_grid.get_obstacle();

    while(iter < min_iter || fabs(flux_grid(p_x, p_y) - prev_f) > prec || fabs(vort_grid(p_x, p_y) - prev_v) > prec){

        ++iter;
        prev_f = flux_grid(p_x, p_y);
        prev_v = vort_grid(p_x, p_y);


        if(ns_obs)
            vort_grid.auto_von_Neumann(OGrid::Condtype::VIFLOW, &flux_grid);

        //if(iter == 1) std::cout << vort_grid << std::endl;


        // vorticity comes first
        for(auto i = 1u; i < vort_grid.x_size()-1; ++i){
            for(auto j = 1u; j < vort_grid.y_size()-1; ++j){
                if( !obs || !obs(vort_grid.xpos(i), vort_grid.ypos(j)) )
                    vort_grid(i, j) = vortHelper(flux_grid, vort_grid, i, j);
            }
        }

        // flux is next (they cannot be done in same loop)
        for(auto i = 1u; i < flux_grid.x_size()-1; ++i){
            for(auto j = 1u; j < flux_grid.y_size()-1; ++j){
                if( !obs || !obs(flux_grid.xpos(i), flux_grid.ypos(j)) )
                    flux_grid(i, j) = fluxHelper(flux_grid, vort_grid, i, j);
            }
        }

    }
    std::cout << iter << std::endl;

    // checking if file is was open and eventually opening it
    auto was_open = file.is_open();
    if(!was_open)
        ns_obs ? file.open("ns_obstackle.dat") : file.open("poiseuille.dat");

    file << "# x y u(x,y) v(x,y)" << std::endl;

    for(auto i = 0u; i < flux_grid.x_size(); ++i){
        for(auto j = 0u; j < flux_grid.y_size(); ++j){

            file << flux_grid.xpos(i) << " " << flux_grid.ypos(j) << " ";

            if(i == flux_grid.x_size()-1 && j == flux_grid.y_size()-1)
                file << 0 << " " << 0;

            if(i == flux_grid.x_size()-1 && j != flux_grid.y_size()-1)
                file << (flux_grid(i, j+1) - flux_grid(i, j)) / flux_grid.dy() << " " << 0;

            if(j == flux_grid.y_size()-1 && i != flux_grid.x_size()-1)
                file << 0 << " " << -(flux_grid(i+1, j) - flux_grid(i, j)) / flux_grid.dx();

            if(i != flux_grid.x_size()-1 && j != flux_grid.y_size()-1)
                file <<  (flux_grid(i, j+1) - flux_grid(i, j)) / flux_grid.dy() << " "
                     << -(flux_grid(i+1, j) - flux_grid(i, j)) / flux_grid.dx();

            file << std::endl;
        }
        file << std::endl;
    }

    //for(auto j = 0u; j < flux_grid.y_size()-1; ++j){
    //    file << flux_grid.ypos(j) << " " << (flux_grid(150, j+1) - flux_grid(150, j)) / flux_grid.dx() << std::endl;
    //}

    //if file wasn't open closing it
    if(!was_open)
        file.close();
}
