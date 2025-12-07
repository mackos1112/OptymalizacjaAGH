#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);

// Problem rzeczywisty - zbiorniki wody
matrix ff_tanks(matrix, matrix = NAN, matrix = NAN);
matrix df_tanks(double, matrix, matrix = NAN, matrix = NAN);