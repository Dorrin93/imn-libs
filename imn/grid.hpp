/**
 * @file   grid.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Header of class containing (hopefully) pretty optimal 2D computational grid and all its parameters.
 */
#ifndef __IMN_GRID__
#define __IMN_GRID__
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include "utils.hpp"

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn{

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
        using gridFunc = std::function<double(unsigned, unsigned)>;

        /**
         * @brief Edge type indicator enum class
         * @var UP
         * @var DOWN
         * @var LEFT
         * @var RIGHT
         */
        enum class Edge{ UP, DOWN, LEFT, RIGHT };

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
        Grid(double x_min, double x_max, double y_min, double y_max, double dx, double dy, func2d initial = nullptr,
             bool widen = false);

        /**
         * @brief Method for aplying function on whole grid.
         *
         * This method takes 2D double funtion, and applies it for the whole grid.
         *
         * @param function Function to apply
         */
        void apply_to_all(const func2d function) noexcept;

        /**
         * @brief Method for applying function on the grid edges.
         *
         * This method takes 2D double funtion, and applies it only on grid edges.
         *
         * @param function Function to apply
         */
        void apply_to_edges(const func2d function) noexcept;

        /**
         * @brief Method for applying function on whole grid except edges.
         *
         * This method takes 2D double funtion, and applies it the whole grid, but omits its edges.
         *
         * @param function Function to apply
         */
        void apply_not_to_edges(const func2d function) noexcept;

        /**
         * @brief Method for applying function on whole grid except edges by index instead of position.
         *
         * Function takes two unsigned "matrix" coordinates and assign function result to those coordinates.
         *
         * @param function Function to apply, defined as <double(unsigned, unsigned)>.
         */
        void apply_index_not_to_edges(const gridFunc function) noexcept;

        /**
         * @brief Method for applying function on the single grid edge.
         * @param function Function to apply
         * @param type Enum edge indication. Could be Edge::left, Edge::right, Edge::up and Edge::down
         */
        void apply_to_single(const func2d function, const Edge type) noexcept;

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
         * @return bool False if ofstream closed, true if everything went fine.
         */
        bool write_to_file(std::ofstream &file, bool zeros = true) const;

        /**
         * @brief Grid-to-string conversion.
         *
         * Funtion writes grid to string in 2D matix manner, rather than gnuplot manner, so only grid values are taken.
         * @return std::string String representation of grid.
         */
        std::string grid_to_string() const;

        /**
         * @brief Value on x-axis from given index
         * @param i Index
         * @return double value on x-axis
         */
        double xpos(const unsigned i) const noexcept { return x_min_ + i * dx_; }

        /**
         * @brief Value on y-axis from given index
         * @pram j Index
         * @return double value on y-axis
         */
        double ypos(const unsigned j) const noexcept { return y_min_ + j * dy_; }

        /**
         * @param x Horizontal index
         * @param y Vertical index
         * @return double Matrix value in given point
         */
        double& operator()(const unsigned x, const unsigned y) noexcept {return mtx_[x * y_size_ + y];}
         /**
         * @param x Horizontal index
         * @param y Vertical index
         * @return double Matrix value in given point
         */
        double operator()(const unsigned x, const unsigned y) const noexcept {return mtx_[x * y_size_ + y];}

        /**
         * @brief Exception safe value return
         * @param x Horizontal index
         * @param y Vertical index
         * @return double Matrix value in given point
         * @throws out_of_range If x * y_size + y > vector size
         */
        double at(const unsigned x, const unsigned y) const { return mtx_.at(x * y_size_ + y); }

        /**
         * @brief Exception safe value return
         * @param x Horizontal index
         * @param y Vertical index
         * @return double Matrix value in given point
         * @throws out_of_range If x * y_size + y > vector size
         */
        double& at(const unsigned x, const unsigned y) { return mtx_.at(x * y_size_ + y); }

        /**
         * @brief x size getter.
         * @return double size of x dimension
         */
        unsigned x_size() const noexcept {return x_size_;}
        /**
         * @brief y size getter
         * @return double size of y dimension
         */
        unsigned y_size() const noexcept {return y_size_;}
        /**
         * @brief x min getter
         * @return int minimal x value
         */
        double x_min() const noexcept {return x_min_;}
        /**
         * @brief x max getter
         * @return int maximal x value
         */
        double x_max() const noexcept {return x_max_;}
        /**
         * @brief y min getter
         * @return int minimal y value
         */
        double y_min() const noexcept {return y_min_;}
        /**
         * @brief y max getter
         * @return int maximal y value
         */
        double y_max() const noexcept {return y_max_;}
        /**
         * @brief x step getter
         * @return double x step value
         */
        double dx() const noexcept {return dx_;}
        /**
         * @brief y step getter
         * @return double y step value
         */
        double dy() const noexcept {return dy_;}

        friend std::ostream& operator<<(std::ostream& os, const Grid& obj);
        Grid(const Grid&) = delete;
        Grid& operator=(const Grid&) = delete;

    protected:

        vMatrix mtx_;
        void increm(unsigned& i, unsigned& j, std::ostream* stm = nullptr) const noexcept;

    private:
        const double x_min_;
        const double x_max_;
        const double y_min_;
        const double y_max_;
        const double dx_;
        const double dy_;
        unsigned x_size_;
        unsigned y_size_;

    };

    std::ostream& operator<<(std::ostream& os, const Grid& obj);

}

/** @} End of group */


#endif
