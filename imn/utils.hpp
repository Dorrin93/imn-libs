/**
 * @file   utils.hpp
 * @author Bartłomiej Meder (bartem93@gmail.com)
 * @date   July, 2015
 * @brief  Some useful functions.
 */
#ifndef __IMN_UTILS__
#define __IMN_UTILS__
#include <string>
#include <functional>

/**
 * @addtogroup imn
 * @{
 */

//* Engineering Numerical Methods functions
namespace imn{

    /**
     * @brief Two-dimensional function
     *
     * It could be function of location and time, or x-location and y-location etc.
     */
    using func2d = std::function<double(double, double)>;

    /**
     * @brief Initial funtion
     *
     * Function of location when time = 0
     */
    using initFunc = std::function<double(double)>;

    /**
     * @brief Casting double to string
     *
     * Function takes double and returns its string representation, without trailing zeros.
     * If last digit after dot is zero, it will stay.
     * Examples:
     * 0.004500 -> 0.0045
     * 1.000000 -> 1.0
     *
     * @param x Value we want to cast
     * @return std::string String representation of x
     */
    std::string double_to_string(double x);

}

/** @} End of group */

#endif
