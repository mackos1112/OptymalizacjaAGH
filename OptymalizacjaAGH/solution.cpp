//Ten plik nie powinien byæ edytowany

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

void solution::clear_calls()
{
	f_calls = 0;
	g_calls = 0;
	H_calls = 0;
}

solution::solution(double L)
{
	x = L;
	g = NAN;
	H = NAN;
	y = NAN;
	ud = NAN;
	flag = -1;
}

solution::solution(const matrix& A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
	ud = NAN;
	flag = -1;
}

solution::solution(int n, double* A)
{
	try
	{
		x = matrix(n, A);
		g = NAN;
		H = NAN;
		y = NAN;
		ud = NAN;
		flag = -1;
	}
	catch (string ex_info)
	{
		throw ("solution::solution(int,double*):\n" + ex_info);
	}
}

solution::solution(const solution& A)
{
	x = A.x;
	g = A.g;
	H = A.H;
	y = A.y;
	if (!isnan(A.ud(0, 0)))
		ud = A.ud;
	flag = A.flag;
}

solution& solution::operator=(const solution& A)
{
	if (&A == this)
		return *this;
	x = A.x;
	g = A.g;
	H = A.H;
	y = A.y;
	if (!isnan(A.ud(0, 0)))
		ud = A.ud;
	flag = A.flag;
	return *this;
}

matrix solution::fit_fun(matrix(*ff)(matrix, matrix, matrix), matrix ud1, matrix ud2)
{
	try
	{
		++f_calls;
		y = ff(x, ud1, ud2);
		return y;
	}
	catch (string ex_info)
	{
		throw ("matrix solution::fit_fun(...):\n" + ex_info);
	}
}

matrix solution::grad(matrix(*gf)(matrix, matrix, matrix), matrix ud1, matrix ud2)
{
	try
	{
		++g_calls;
		g = gf(x, ud1, ud2);
		return g;
	}
	catch (string ex_info)
	{
		throw ("matrix solution::grad(...):\n" + ex_info);
	}
}

matrix solution::hess(matrix(*Hf)(matrix, matrix, matrix), matrix ud1, matrix ud2)
{
	try
	{
		++H_calls;
		H = Hf(x, ud1, ud2);
		return H;
	}
	catch (string ex_info)
	{
		throw ("matrix solution::hess(...):\n" + ex_info);
	}
}

ostream& operator<<(ostream& S, const solution& A)
{
	S << "x = " << A.x << endl;
	S << "y = " << A.y << endl;
	S << "f_calls = " << solution::f_calls << endl;
	if (solution::g_calls > 0)
		S << "g_calls = " << solution::g_calls << endl;
	if (solution::H_calls > 0)
		S << "H_calls = " << solution::H_calls << endl;
	S << "Exit flag: " << A.flag << endl;
	return S;
}

int get_dim(const solution& A)
{
	try
	{
		return get_len(A.x);
	}
	catch (string ex_info)
	{
		throw ("int get_dim(const solution&):\n" + ex_info);
	}
}