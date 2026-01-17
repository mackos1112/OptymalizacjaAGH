#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix x1, matrix x2, matrix ud1);
matrix ff_ball(matrix, matrix = NAN, matrix = NAN);
matrix ff3T(matrix x1, matrix x2, matrix ud1);
matrix ff3R_sym(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff3R(matrix x1, matrix x2, matrix ud1);
// Problem rzeczywisty - zbiorniki wody
matrix ff_tanks(matrix, matrix = NAN, matrix = NAN);
matrix df_tanks(double, matrix, matrix = NAN, matrix = NAN);