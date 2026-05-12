#pragma once

// Compile-time grid resolution. Changing RES requires recompilation because
// all solver arrays are statically allocated with this size.
// All other parameters are read at runtime from bennu_input.txt.
constexpr int    RES = 200;
constexpr double PI  = 3.14159265358979323846;
