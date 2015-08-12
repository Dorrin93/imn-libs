/**
 * @file   ogrid.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Header of the class containing 2D computational grid with obstackle.
 */
#ifndef __IMN_OGIRD__
#define __IMN_OGIRD__
#include "grid.hpp"
#include <iostream>
#include <tuple>

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn{

    /**
     * @class OGrid
     * @brief 2D computational grid with obstacle.
     *
     * Class represents two-dimensional grid of double values. It delivers interface for changing values inside of
     * grid, and / or obstacle.
     * Obstacle is defined as function taking two coordinates in double format and returning bool: true if position
     * is inside obstacle, or false if it is not.
     */
    class OGrid: public Grid{

        using ofunc = std::function<bool(double, double)>;

    public:
        /**
         * @brief Obstacle point type indicator.
         * @var LEFT Left edge
         * @var RIGHT Right edge
         * @var UP Upper edge
         * @var DOWN Lower edge
         * @var LU Left upper corner
         * @var LD Left lower corner
         * @var RU Right upper corner
         * @var RD Right lower corner
         * @var IN Inside
         */
        enum class Ptype { LEFT, RIGHT, UP, DOWN, LU, LD, RU, RD, IN };

        /**
         * @brief Obstacle placement indicator.
         * @var DOWN Obstacle by the lower edge of grid
         * @var UP Obstacle by the upper edge of grid
         * @var MID "Floating" obstacle
         */
        enum class Obstype { DOWN, UP, MID };

        /**
         * @brief Von Neumann condition type.
         * @var PFLOW Condition for potential flow function
         * @var VIFLOW Condition for viscous flow function
         */
        enum class Condtype { PFLOW, VIFLOW };

        /**
         * @param x_min Minimal horizontal grid value.
         * @param x_max Maximal horizontal grid value.
         * @param y_min Minimal vertical grid value.
         * @param y_max Maximal vertical grid value.
         * @param dx Horizontal step.
         * @param dy Vertical step.
         * @param widen Indicator, if grid will be widen by one (usualy useful). Default to true.
         */
        OGrid(double x_min, double x_max, double y_min, double y_max, double dx, double dy, bool widen = true);

        /**
         * @brief Obstacle setter.
         *
         * Setting obstacle is the first thing user should do after creating an object.
         * By the default obstacle is set to nullptr, so not setting it probably will cause
         * std::bad_function_call exception.
         * This funtion also automatically sets obstacle type.
         * @param obs Obstacle funtion, type <bool(double, double)>.
         */
        void set_obstackle(const ofunc obs);

        /**
         * @brief Obstacle getter.
         * @return std::function<bool(double, double)> Obstacle function.
         */
        ofunc get_obstacle() const noexcept { return obstacle_p_m; }

        /**
         * @brief Obstacle placement type getter.
         * @return imn::OGrid::Obstype Placement type.
         */
        Obstype obstacle_type() const noexcept { return obs_type_; }

        /**
         * @brief Method for applying 2D function to the obstacle.
         *
         * This method takes 2D double funtion, and applies it on elements, which posision on grid is true for
         * obstacle function.
         *
         * @param function Function to apply
         */
        void apply_to_obstackle(const func2d function) noexcept;

        /**
         * @brief Method for applying 2D function everywhere, but not to obstacle and grid edges.
         *
         * This method takes 2D double funtion, and applies it on elements, which posision on grid is false for
         * obstacle function and edges.
         *
         * @param function Function to apply
         */
        void apply_not_to_obstackle(const func2d function) noexcept;

        /**
         * @brief Method for applying vector index function everywhere, but not to obstacle and grid edges.
         *
         * Method takes function opearating on two vector coordinates, rather than grid coordinates. It may be
         * used for direct grid vector modyfication. Then, method iterates everywhere, but obstacle and edges.
         * gridFunc is defined as std::function<double(unsigned, unsigned)>.
         * @param function Index-based function to apply.
         */
        void apply_index_not_to_obstackle(const gridFunc function) noexcept;

        /**
         * @brief Method autmaticaly applying von Neumann condition to the obstacle.
         *
         * Depending on problem, this method can apply Von Neumann condition autmatically to the obstacle edges.
         * At this moment, there are two types of such possibilities:
         * - Potential flow
         * - Viscous flow
         * It is recommended to use this method for such problem. If problem is different, user should use
         * von_Neumann_cond method separatly for every obstacle edge / corner.
         *
         * @param type Type of problem we want to apply von Neumann conditions.
         * @param other Used only in Viscous flow, for passing flux grid.
         */
        void auto_von_Neumann(Condtype type, const OGrid *other = nullptr);

        /**
         * @brief Method for applying von Neuman condition for one specyfic point / edge type.
         * @param function 2D double funtion we want to apply to edge or point type.
         * @param pointType Type of edge or point we want condition to apply.
         */
        void von_Neumann_cond(const func2d function, const Ptype pointType);

    private:
        using vNfunc = std::function<double(int, int, Ptype, const OGrid*)>;
        using vMatrixObsCoord = std::vector< std::tuple<unsigned, unsigned, Ptype> >;
        using vMatrixCoord = std::vector< std::pair<unsigned, unsigned> >;

        ofunc obstacle = nullptr;
        ofunc obstacle_p_m = nullptr;
        Obstype obs_type_;

        vMatrixObsCoord obs_in_;
        vMatrixCoord obs_out_;

        Ptype point_type(unsigned i, unsigned j) const noexcept;
        double pFlowFunc(unsigned i, unsigned j, Ptype type, const OGrid *other = nullptr) const noexcept;
        double viFlowFunc(unsigned i, unsigned j, Ptype type, const OGrid *other) const noexcept;


    };

}

/** @} End of group */

#endif
