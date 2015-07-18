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

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn{

    class OGrid: public Grid{

        using ofunc = std::function<bool(double, double)>;
        using vNfunc = std::function<double(int, int, Ptype)>;

    public:
        enum class Ptype {left, right, up, down, lu, ld, ru, rd, in};

        enum class VonNeumann {pFlow, viFlow};

        OGrid(int x_min, int x_max, int y_min, int y_max, double dx, double dy, bool widen = true) :
                Grid(x_min, x_max, y_min, y_max, dx, dy, nullptr, widen){}

        void set_obstackle(ofunc&& obstackle) noexcept { clear(); obstacle = std::move(obstackle); }

        void apply_to_obstackle(func2d function) noexcept;

        void apply_not_to_obstackle(func2d function) noexcept;

        void auto_von_Neumann(VonNeumann type);

        void von_Neumann_cond(func2d function, Ptype pointType);

    private:
        ofunc obstacle = nullptr;
        Ptype type(unsigned i, unsigned j) const noexcept;
        double pFlowFunc(unsigned i, unsigned j, Ptype type) const noexcept;
        double viFlowFunc(unsigned i, unsigned j, Ptype type) const noexcept;

    };

}

/** @} End of group */

#endif
