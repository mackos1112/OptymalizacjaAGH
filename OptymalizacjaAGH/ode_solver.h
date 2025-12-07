//Ten plik nie powinien byc edytowany

#pragma once

#include"matrix.h"
#include"user_funs.h"

// Solve ODE with signature: diff(t, Y, ud1, ud2)
matrix* solve_ode(matrix (*diff)(double, matrix, matrix, matrix), double t0, double dt, double tend, matrix Y0, matrix ud1 = matrix(), matrix ud2 = matrix()); // throw (string);
