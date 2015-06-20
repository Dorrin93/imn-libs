/**
 * \file   vi_flow.cxx
 * \Author Bartłomiej Meder (bartem93@gmail.com)
 * \date   June, 2015
 * \brief  Viscous flow numerical methods
 *
 * Source for implementation of numerical methods created for solving viscous flow of liquid.
 * Everything is based on solving Navier-Stokes equasion.
 */

#include "vi_flow.hpp"

double imn::vortHelper(double **f, double **v, long x, long y)
{
    // 0.0625 is 1/16
    return 0.25 * (v[x+1][y] + v[x-1][y] + v[x][y-1] + v[x][y+1]) -
           0.0625 * ( (f[x][y+1] - f[x][y-1]) * (v[x+1][y] - v[x-1][y]) -
           (f[x+1][y] - f[x-1][y]) * (v[x][y+1] - v[x][y-1]) );
}


double imn::fluxHelper(double **f, double **v, long x, long y,
        double dx, double dy)
{
    return (f[x+1][y] + f[x-1][y] + f[x][y+1] + f[x][y-1] - v[x][y] * 
            dx * dy) * 0.25;
}


bool imn::inside(double x, double y, double dx, double dy, const std::vector<Wall*> *walls)
{
    bool l, r, u, d;
    l = r = u = d = false;
    // casting to long in order to avoid numerical errors
    long X, Y, Ax, Ay, Bx, By;
    X = lround(x / dx);
    Y = lround(y / dy);
    for(auto& it : *walls){

        Ax = lround(it->A.x / dx);
        Ay = lround(it->A.y / dy);
        Bx = lround(it->B.x / dx);
        By = lround(it->B.y / dy);

        switch(it->type){

            case Wall::UP:
                if(X >= Ax && X <= Bx && Y <= Ay)
                    u = true;
                break;
            case Wall::DOWN:
                if(X >= Ax && X <= Bx && Y >= Ay)
                    d = true;
                break;
            case Wall::LEFT:
                if(Y >= Ay && Y <= By && X >= Ax)
                    l = true;
                break;
            case Wall::RIGHT:
                if(Y >= Ay && Y <= By && X <= Ax)
                    r = true;
                break;
        }
    }
    return l && r && u && d;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imn::poiseuille(double xmin, double xmax, double ymin, double ymax, double dx, double dy, imn::func flux,
                     imn::func vort, std::ofstream &file, double prec, int min_iter)
{
    // actually it's only wrapper
    imn::navier_stokes(xmin, xmax, ymin, ymax, dx, dy, flux, vort, file, prec, min_iter);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imn::ns_obstacle(double xmin, double xmax, double ymin, double ymax, double dx, double dy, func flux, func vort,
                  const std::vector<Wall*> &walls, std::ofstream &file, double prec, int min_iter, bool upper)
{
    imn::navier_stokes(xmin, xmax, ymin, ymax, dx, dy, flux, vort, file, prec, min_iter, true, &walls, upper);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imn::navier_stokes(double xmin, double xmax, double ymin, double ymax, double dx, double dy, func flux, func vort,
                       std::ofstream &file, double prec, int min_iter, bool obstacle, const std::vector<Wall*> *walls,
                        bool top)
{
    // 'gettin grid size
    long x_size = lround((xmax - xmin) / dx);
    long y_size = lround((ymax - ymin) / dy);

    // gird allocation
    double** f_grid = new double*[x_size+1];
    double** v_grid = new double*[x_size+1];
    for(long i = 0; i <= x_size; ++i){
        f_grid[i] = new double[y_size+1];
        v_grid[i] = new double[y_size+1];
    }

    // initial conditions
    for(long i = 0; i <= x_size; ++i){
        for(long j = 0; j <= y_size; ++j){

            if(!i || !j || i == x_size || j == y_size){
                f_grid[i][j] = flux(xmin + i*dx, ymin + j*dy);
                v_grid[i][j] = vort(xmin + i*dy, ymin + j*dy);
            }
            else{
                f_grid[i][j] = 0;
                v_grid[i][j] = 0;
            }

            if(obstacle){
                if(imn::inside(xmin + i*dx, ymin + j*dy, dx, dy, walls))

                    if(top)
                        f_grid[i][j] = flux(xmin + i*dx, ymax);
                    else
                        f_grid[i][j] = flux(xmin + j*dx, ymin);

                else {
                    f_grid[i][j] = flux(xmin + i*dx, ymin + j*dy);
                    v_grid[i][j] = vort(xmin + i*dx, ymin + j*dy);
                }
            }

        }
    }

    // some special values - we check condition in 2/3 x and 0.5 y point
    long p_x = lround((2./3.) * x_size);
    long p_y = y_size >> 1;
    double prev_f = 0;
    double prev_v = 0;
    long iter = 0;

    while(iter < min_iter || fabs(f_grid[p_x][p_y] - prev_f) > prec || fabs(v_grid[p_x][p_y] - prev_v) > prec){

        ++iter;
        prev_f = f_grid[p_x][p_y];
        prev_v = v_grid[p_x][p_y];

        //von Neumann conditions
        if(obstacle){
            for(long i = 0; i <= x_size; ++i){
                v_grid[i][0] = (2 * (f_grid[i][1] - f_grid[i][0])) / (dx * dy);
                v_grid[i][y_size] = (2 * (f_grid[i][y_size-1] - f_grid[i][y_size])) / (dx * dy);
            }

            for(auto& it : *walls){
                //A.x <= B.x && A.y <= B.y
                long left = lround( (it->A.x - xmin) / dx );
                long right = lround( (it->B.x - xmin) / dx);
                long down = lround( ( it->A.y - ymin) / dy);
                long up = lround ( ( it->B.y - ymin) / dy);
                switch(it->type){
                    case Wall::LEFT:
                        for(long j = down; j <= up; ++j)
                           v_grid[left][j] = (2 * (f_grid[left-1][j] - f_grid[left][j])) / (dx*dy);
                        break;
                    case Wall::RIGHT:
                        for(long j = down; j <= up; ++j)
                            v_grid[right][j] = (2 * (f_grid[right+1][j] - f_grid[right][j])) / (dx*dy);
                        break;
                    case Wall::UP:
                        for(long i = left; i <= right; ++i)
                            v_grid[i][up] = (2 * (f_grid[i][up+1] - f_grid[i][up])) / (dx*dy);
                        break;
                    case Wall::DOWN:
                        for(long i = left; i <= right; ++i)
                            v_grid[i][down] = (2 * (f_grid[i][down-1] - f_grid[i][down])) / (dx*dy);
                        break;
                }
            }
        }

        // vorticity comes first
        for(long i = 1; i < x_size; ++i){
            for(long j = 1; j < y_size; ++j){
                if(!obstacle || !imn::inside(xmin + i*dx, ymin + j*dy, dx, dy, walls))
                    v_grid[i][j] = imn::vortHelper(f_grid, v_grid, i, j);
            }
        }
        // flux is next (they cannot be done in same loop)
        for(long i = 1; i < x_size; ++i){
            for(long j = 1; j < y_size; ++j){
                if(!obstacle || !imn::inside(xmin + i*dx, ymin + j*dy, dx, dy, walls))
                    f_grid[i][j] = imn::fluxHelper(f_grid, v_grid, i, j , dx, dy);
            }
        }
    }

    // checking if file is was open and eventually opening it
    bool was_open = file.is_open();
    if(!was_open)
        obstacle ? file.open("ns_obstackle.dat") : file.open("poiseuille.dat");

    file << "# x y u(x,y) v(x,y)" << std::endl;

    for(long i = 0; i <= x_size; ++i){
        for(long j = 0; j <= y_size; ++j){

            file << xmin + i*dx << " " << ymin + j*dy << " ";

            if(obstacle && inside(xmin + i*dx, ymin + j*dy, dx, dy, walls)){
                file << 0 << " " << 0;
            }

            else{
                if(i != x_size && j != y_size)
                    file << (f_grid[i][j+1] - f_grid[i][j]) / dy << " " << -(f_grid[i+1][j] - f_grid[i][j]) / dx;

                else if(i == x_size && j != y_size)
                    file << (f_grid[i][j+1] - f_grid[i][j]) / dy << " " << 0;

                else if(j == y_size && i != x_size)
                    file << 0 << " " << -(f_grid[i+1][j] - f_grid[i][j]) / dx;

                else
                    file << 0 << " " << 0;
            }

            file << std::endl;
        }
        file << std::endl;
    }

    //if file wasn't open closing it
    if(!was_open)
        file.close();


    // grid deallocation
    for(long i = 0; i <= x_size; ++i){
        delete[] f_grid[i];
        delete[] v_grid[i];
    }
    delete[] f_grid;
    delete[] v_grid;
}
