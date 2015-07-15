/**
 * \file   grid.hpp
 * \Author Bartłomiej Meder (bartem93@gmail.com)
 * \date   July, 2015
 * \brief  Header of class containing (hopefully) pretty optimal 2D computational grid and all its parameters.
 */
#ifndef __IMN_GRID__
#define __IMN_GRID__
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

/**
 * \addtogroup imn
 * \{
 */

//* Engineering Numerical Methods functions
namespace imn{

    using func2d = double(*)(double, double);

    class Grid{
        // basicaly yeah, this whole class is justy fancy packed vector of doubles
        // I've been thinking about using std::array but nah, to much bothering with constexpr
        // this vector should be fixed size anyway
        using vMatrix = std::vector<double>;

    public:

        Grid(int x_min, int x_max, int y_min, int y_max, double dx, double dy, func2d initial = nullptr);

        void apply_to_all(func2d function) noexcept;
        void apply_to_edge(func2d function) noexcept;
        void apply_not_to_edge(func2d function) noexcept;
        void clear() noexcept;
        bool write_to_file(std::ofstream &file, bool zeros = true) const;
        std::string grid_to_string() const;

        double& operator()(int x, int y) noexcept {return _mtx[x * _x_size + y];}
        double operator()(int x, int y) const noexcept {return _mtx[x * _x_size + y];}

        friend std::ostream& operator<<(std::ostream& os, const Grid& obj);

        std::size_t x_size() const noexcept {return _x_size;}
        std::size_t y_size() const noexcept {return _y_size;}
        int x_min()  const noexcept {return _x_min;}
        int x_max() const noexcept {return _x_max;}
        int y_min() const noexcept {return _y_min;}
        int y_max() const noexcept {return _y_max;}
        double dx() const noexcept {return _dx;}
        double dy() const noexcept {return _dy;}

        Grid(const Grid&) = delete;
        Grid& operator=(const Grid&) = delete;

    private:
        vMatrix _mtx;
        int _x_min;
        int _x_max;
        int _y_min;
        int _y_max;
        double _dx;
        double _dy;
        std::size_t _x_size;
        std::size_t _y_size;
    };

    //another declaration (damn you, compiler :< )
    std::ostream& operator<<(std::ostream& os, const Grid& obj);

}

/** \} End of group */


#endif
