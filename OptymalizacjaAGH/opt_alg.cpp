#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo³añ funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix),
	double x0, double d, double alpha, int Nmax,
	matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0.0, 0.0 };

		// pomocnicza funkcja f(x)
		auto f = [&](double x) -> double
		{
			matrix X(1, 1);
			X(0, 0) = x;
			matrix Y = ff(X, ud1, ud2);
			return Y(0, 0);
		};

		int f_calls = 0;

		int i = 0;
		double x_0 = x0;          // x^(0)
		double x1 = x0 + d;      // x^(1)

		double fx0 = f(x_0);  ++f_calls;
		double fx1 = f(x1);   ++f_calls;

		// linie 3–4
		if (fx1 == fx0)
		{
			p[0] = x_0;
			p[1] = x1;
			return p;
		}

		// linie 6–12
		if (fx1 > fx0)
		{
			d = -d;
			x1 = x0 + d;
			fx1 = f(x1);  ++f_calls;

			if (fx1 >= fx0)
			{
				p[0] = x1;
				p[1] = x0 - d;
				return p;
			}
		}

		// mamy: x^(0) = x_0, x^(1) = x1
		double x_im1 = x_0;   // x^(i-1)
		double x_i = x1;    // x^(i)
		double x_ip1;         // x^(i+1)

		double f_im1 = fx0;
		double f_i = fx1;
		double f_ip1;

		// linie 13–19  (pêtla repeat ... until)
		while (true)
		{
			if (f_calls > Nmax)
			{
				delete[] p;
				throw std::string("expansion: przekroczono maksymalna liczbe wywolan funkcji");
			}

			i = i + 1;
			x_ip1 = x0 + std::pow(alpha, i) * d;   // x^(i+1)
			f_ip1 = f(x_ip1);  ++f_calls;

			if (f_i <= f_ip1)
				break;

			// przesuwamy indeksy: (i-1) <- i, i <- i+1
			x_im1 = x_i;   f_im1 = f_i;
			x_i = x_ip1; f_i = f_ip1;
		}

		// linie 20–23 – zwracamy przedzia³ zawieraj¹cy minimum
		if (d > 0.0)
		{
			p[0] = x_im1;   // x^(i-1)
			p[1] = x_ip1;   // x^(i+1)
		}
		else
		{
			p[0] = x_ip1;   // x^(i+1)
			p[1] = x_im1;   // x^(i-1)
		}

		return p;
	}
	catch (std::string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}


solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
