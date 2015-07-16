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

    /**
     * Standard 2D funtion.
     * It could be, for example, density of potential.
     */
    using func2d = double(*)(double, double);

    /**
     * @class Grid
     * @brief 2D computational grid.
     *
     * Class represents two-dimensional grid of double values. It delivers interface for changing values inside of
     * grid, so analogy to matrix is rather misguided.
     * Grid sizes are calculated as (n_max - n_min) / dn, where n is x or y. If one needs to have slightly larger grid,
     * there is widen parameter. If true, sizes will be increased by one.
     */
    class Grid{

        // basicaly yeah, this whole class is justy fancy packed vector of doubles
        // I've been thinking about using std::array but nah, to much bothering with constexpr
        // this vector should be fixed size anyway
        using vMatrix = std::vector<double>;

    public:

        /**
         * @param x_min Minimal horizontal grid value.
         * @param x_max Maximal horizontal grid value.
         * @param y_min Minimal vertical grid value.
         * @param y_max Maximal vertical grid value.
         * @param dx Horizontal step.
         * @param dy Vertical step.
         * @param initial Function we want to use for initial values computation. If none given, values will be set to
         * zeros. Default to nullptr.
         * @param widen Indicator, if grid will be widen by one (sometimes it's useful). Default to false.
         */
        Grid(int x_min, int x_max, int y_min, int y_max, double dx, double dy, func2d initial = nullptr, bool widen = false);

        /**
         * @brief Method for aplying function on whole grid.
         *
         * This method takes 2D double funtion, and applies it for the whole grid.
         *
         * @param function Function to apply.
         */
        void apply_to_all(func2d function) noexcept;

        /**
         * @brief Method for aplying function on the grid edges.
         *
         * This method takes 2D double funtion, and applies it only on grid edges.
         *
         * @param function Function to apply.
         */
        void apply_to_edge(func2d function) noexcept;

        /**
         * @brief Method for aplying function on whole grid except edges.
         *
         * This method takes 2D double funtion, and applies it the whole grid, but omits its edges.
         *
         * @param function Function to apply.
         */
        void apply_not_to_edge(func2d function) noexcept;

        /**
         * @brief Method setting whole grid to zero.
         */
        void clear() noexcept;

        /**
         * @brief Method writing grid to file.
         *
         * Method taking otput file stream and writing grid to it, in gnuplot friendly format:
         * x_point y_point grid_value
         * File passed to method must be open. Method is also able to omit values equal to zero.
         *
         * @param file Output file stream we want write to.
         * @param zeros Indicator for omiting zero grid values. Defalut to false.
         * @return False if ofstream closed, true if everything went fine.
         */
        bool write_to_file(std::ofstream &file, bool zeros = true) const;

        /**
         * @brief Grid-to-string conversion.
         *
         * Funtion writes grid to string in 2D matix manner, rather than gnuplot manner, so only grid values are taken.
         * @return String representation of grid.
         */
        std::string grid_to_string() const;

        /**
         * @param x Horizontal index
         * @param y Vertical index
         * @return Matrix value in given point
         */
        double& operator()(unsigned x, unsigned y) noexcept {return _mtx[x * _x_size + y];}
         /**
         * @param x Horizontal index
         * @param y Vertical index
         * @return Matrix value in given point
         */
        double operator()(unsigned x, unsigned y) const noexcept {return _mtx[x * _x_size + y];}

        /**
         * @brief x size getter.
         * @return size of x dimension
         */
        std::size_t x_size() const noexcept {return _x_size;}
        /**
         * @brief y size getter
         * @return size of y dimension
         */
        std::size_t y_size() const noexcept {return _y_size;}
        /**
         * @brief x min getter
         * @return minimal x value
         */
        int x_min()  const noexcept {return _x_min;}
        /**
         * @brief x max getter
         * @return maximal x value
         */
        int x_max() const noexcept {return _x_max;}
        /**
         * @brief y min getter
         * @return minimal y value
         */
        int y_min() const noexcept {return _y_min;}
        /**
         * @brief y max getter
         * @return maximal y value
         */
        int y_max() const noexcept {return _y_max;}
        /**
         * @brief x step getter
         * @return x step value
         */
        double dx() const noexcept {return _dx;}
        /**
         * @brief y step getter
         * @return y step value
         */
        double dy() const noexcept {return _dy;}

        friend std::ostream& operator<<(std::ostream& os, const Grid& obj);
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

    //next operator<< declaration (damn you, compiler :< )
    std::ostream& operator<<(std::ostream& os, const Grid& obj);

}

/** \} End of group */


#endif
