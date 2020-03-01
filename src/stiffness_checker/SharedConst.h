#pragma once
#include <string>
#include <iostream>

#ifndef NDEBUG
// https://stackoverflow.com/questions/3767869/adding-message-to-assert
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

namespace conmech {

/**
 * @brief numerical zero
 * See: https://github.com/libigl/libigl/blob/master/include/igl/EPS.h
 */
const static double DOUBLE_EPS    = 1.0e-14;
const static double DOUBLE_EPS_SQ = 1.0e-28;
const static double EPSILON = 1.0e-14;

/**
 * @brief unit system used inside conmech. Input unit systems will be converted to these before computation.
 * 
 */
const static std::string LENGTH_UNIT = "meter";
const static std::string ROT_ANGLE_UNIT = "rad";
const static std::string FORCE_UNIT = "kN";
const static std::string MOMENT_UNIT = "kN-m";

// deduced from above
const static std::string PRESSURE_UNIT = "kN/m^2"; // modulus

} // ns conmech