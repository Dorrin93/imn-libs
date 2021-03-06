/**
 * @file   p_flow.cxx
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Potential flow numerical methods.
 *
 * Source of implementation of numerical methods created for solving potential flow of liquids (Laplace equasion).
 */
#include "p_flow.hpp"

auto sqr = [](double x){ return x * x; };

void imn::stream_lines(OGrid &grd, double u0, std::ofstream &file, double eps, bool w_grid)
{
    base_lines(grd, u0, file, eps, w_grid, true);
}

void imn::potential_lines(OGrid& grd, double u0, std::ofstream& file, double eps, bool w_grid)
{
    base_lines(grd, u0, file, eps, w_grid, false);
}


void ::imn::base_lines(imn::OGrid &grd, double u0, std::ofstream &file, double eps, bool w_grid, bool stream)
{
    grd.clear();

    if(stream) {

        grd.apply_to_single([u0](double x, double y) { return u0 * y; }, Grid::Edge::LEFT);
        grd.apply_to_single([u0](double x, double y) { return u0 * y; }, Grid::Edge::RIGHT);
        grd.apply_to_single([&grd](double x, double y) { return grd(0, 0); }, Grid::Edge::DOWN);
        grd.apply_to_single([&grd](double x, double y) { return grd(0, grd.y_size() - 1); }, Grid::Edge::UP);

        switch (grd.obstacle_type()) {

            case OGrid::Obstype::DOWN:
                grd.apply_to_obstackle([&grd](double x, double y) { return grd(0, 0); });
                break;
            case OGrid::Obstype::UP:
                grd.apply_to_obstackle([&grd](double x, double y) { return grd(0, grd.y_size() - 1); });
                break;
            case OGrid::Obstype::MID:
                grd.apply_to_obstackle([&grd](double x, double y) { return grd(0, (grd.y_size() - 1) / 2); });
                break;
        }
    }

    else
        grd.apply_to_all([u0](double x, double y){ return u0 * x; });

    // file operations
    auto was_open = file.is_open();

    if(!was_open)
        if(stream)
            file.open("sl_integral.dat");
        else
            file.open("pl_integral.dat");

    file << "# n a" << std::endl;

    std::ofstream grdfile;
    if(w_grid)
        if(stream)
            grdfile.open("sl_grid.dat");
        else
            grdfile.open("pl_grid.dat");

    // calculations
    auto a = 0.;
    auto a_bf = 1.;

    auto obs = grd.get_obstacle();

    for(auto iter = 1L; fabs(a - a_bf) > eps; ++iter){
        a_bf = a; a = 0;

        if(!stream)
            grd.auto_von_Neumann(OGrid::Condtype::PFLOW, nullptr);

        grd.apply_index_not_to_obstackle(
                [&grd](unsigned x, unsigned y){ return 0.25 * (grd(x-1, y) + grd(x+1, y) + grd(x, y-1) + grd(x, y+1)); }
        );

        for(auto i = 1u; i < grd.x_size()-1; ++i){
            for(auto j = 1u; j < grd.y_size()-1; ++j){
                if(!obs(grd.xpos(i), grd.ypos(j)))
                    a += 0.125 * (sqr(grd(i+1, j) - grd(i-1, j)) + sqr(grd(i, j+1) - grd(i, j-1)));
            }
        }

        file << iter << " " << a << std::endl;

    }

    if(w_grid) {
        grd.apply_to_obstackle([](double x, double y){ return 0; });
        grdfile << grd;
        grdfile.close();
    }

    if(!was_open)
        file.close();

}
