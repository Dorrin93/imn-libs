#ifndef __VI_FLUX__
#define __VI_FLUX__

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

namespace imn{

typedef double(*func)(double, double);

enum WallType{UP, DOWN, LEFT, RIGHT};

struct Point{
    double x;
    double y;
};

struct Wall{
    Point A;
    Point B;
    WallType type;
};

void poiseuille(double xmin, double xmax, double ymin,
        double ymax, double dx, double dy, func flux,
        func vort, std::ofstream &file,
        double prec = 1e-7, int min_iter = 200);

void obstacle(double xmin, double xmax, double ymin,
        double ymax, double dx, double dy, func flux,
        func vort, const std::vector<Wall*> &walls, 
        std::ofstream &file,
        double prec = 1e-7, int min_iter = 500);



double vortHelper(double **f, double **v, int x, int y);

double fluxHelper(double **f, double **v, int x, int y,
        double dx, double dy);

bool inside(double x, double y, const std::vector<Wall*> &walls);

}

#endif

