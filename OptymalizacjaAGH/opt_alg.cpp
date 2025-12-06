#include"opt_alg.h"

#define DEBUG false

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejściowe:
	// ff - wskaźnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zakłądana dokładność rozwiązania
	// Nmax - maksymalna liczba wywołań funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosując rozkład jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwiązanie do przedziału [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartość funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwiązanie z zadaną dokładnością
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywołań funkcji celu
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
		double* p = new double[3]{0.0, 0.0, 0.0};

		// pomocnicza funkcja celu f(x)
		auto f = [&](double x) -> double {
			matrix X(1, 1);
			X(0, 0) = x;
			matrix Y = ff(X, ud1, ud2);
			return Y(0, 0);
		};

		int f_calls = 0;                          // licznik wywołań funkcji celu

		// punkty startowe
		double x0_curr = x0;                      // x^(0)
		double x1_curr = x0 + d;                  // x^(1)

		double f0 = f(x0_curr); ++f_calls;        // f(x^(0))
		double f1 = f(x1_curr); ++f_calls;        // f(x^(1))

		if( DEBUG )
		std::cout << "[DEBUG] i=0, x0=" << x0_curr << " f(x0)=" << f0
			<< ", x1=" << x1_curr << " f(x1)=" << f1 << std::endl;

		// przypadek: wartości równe
		if (std::abs(f1 - f0) < 1e-12) {
			p[0] = x0_curr;
			p[1] = x1_curr;

			if (DEBUG)
			std::cout << "[DEBUG] f(x1) == f(x0) → zwracam [" << p[0] << ", " << p[1] << "]" << std::endl;
			
			p[2] = f_calls;
			return p;
		}

		// przypadek: funkcja rośnie → zmiana kierunku
		if (f1 > f0) {
			d = -d;
			x1_curr = x0 + d;
			f1 = f(x1_curr); ++f_calls;

			if (DEBUG)
			std::cout << "[DEBUG] f(x1) > f(x0), zmiana kierunku: x1=" << x1_curr << " f(x1)=" << f1 << std::endl;

			if (f1 >= f0) {
				p[0] = x1_curr;
				p[1] = x0 + d;   // poprawka zgodna z pseudokodem

				if (DEBUG)
				std::cout << "[DEBUG] f(x1) >= f(x0) → zwracam [" << p[0] << ", " << p[1] << "]" << std::endl;
				
				p[2] = f_calls;
				return p;
			}
		}

		// pętla ekspansji
		int i = 0;
		double x_next, f_next;
		while (true) {
			if (f_calls > Nmax) {
				delete[] p;
				throw std::string("expansion: przekroczono maksymalną liczbę wywołań funkcji celu");
			}

			i++;
			x_next = x0 + std::pow(alpha, i) * d; // x^(i+1) = x^(0) + α^i * d
			f_next = f(x_next); ++f_calls;

			if (DEBUG)
			std::cout << "[DEBUG] i=" << i << ", x_next=" << x_next << " f(x_next)=" << f_next << std::endl;

			if (f1 <= f_next) {
				if (DEBUG)
				std::cout << "[DEBUG] warunek stopu: f(x(i)) <= f(x(i+1))" << std::endl;
				
				break;
			}

			// przesuwamy okno: (i-1) ← i
			x0_curr = x1_curr; f0 = f1;
			x1_curr = x_next; f1 = f_next;
		}

		// zwracamy przedział zawierający minimum
		if (d > 0.0) {
			p[0] = x0_curr;   // x^(i-1)
			p[1] = x_next;    // x^(i+1)
		}
		else {
			p[0] = x_next;    // x^(i+1)
			p[1] = x0_curr;   // x^(i-1)
		}

		if (DEBUG)
		std::cout << "[DEBUG] Zwracam przedział [" << p[0] << ", " << p[1] << "]" << std::endl;

		p[2] = f_calls;
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
		int f_calls = 0;

		// helper f(x)
		auto f = [&](double x) {
			matrix X(1,1);
			X(0,0) = x;
			f_calls++;
			return ff(X, ud1, ud2)(0,0);
			};

		// build fibonacci seq until φ_k > (b-a)/ε
		std::vector<long long> F = {1,1};
		while (F.back() < (b - a) / epsilon)
			F.push_back(F[F.size()-1] + F[F.size()-2]);
		int k = F.size() - 1;

		double a_i = a;
		double b_i = b;

		double c_i = b_i - (double)F[k-1] / F[k] * (b_i - a_i);
		double d_i = a_i + b_i - c_i;

		for (int i = 0; i < k-2; i++)
		{
			if (f(c_i) < f(d_i)) {
				b_i = d_i;
			} else {
				a_i = c_i;
			}

			c_i = b_i - (double)F[k-i-2] / F[k-i-1] * (b_i - a_i);
			d_i = a_i + b_i - c_i;
		}

		Xopt.x = c_i;
		Xopt.y = f(c_i);
		Xopt.f_calls = f_calls;


		// zakładamy lokalne minimum, bo fib() nie sprawdza całej funkcji
		Xopt.flag = 0;       // 0 = lokalne, 1 = globalne

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lagrange(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		// Implementation based on provided pseudocode

		int i = 0;
		double ai = a;
		double bi = b;
		double ci = (a + b) / 2.0;  // Initial c(0) as midpoint
		int f_calls = 0;

		// Helper function for evaluating the original function
		auto f = [&](double x) -> double {
			matrix X(1, 1);
			X(0, 0) = x;
			matrix Y = ff(X, ud1, ud2);
			f_calls++;
			return Y(0, 0);
			};

		while (f_calls < Nmax) {
			// Step 4
			double l = f(ai) * ((bi * bi) - (ci * ci))
				+ f(bi) * ((ci * ci) - (ai * ai))
				+ f(ci) * ((ai * ai) - (bi * bi));

			// Step 5
			double m = f(ai) * (bi - ci)
				+ f(bi) * (ci - ai)
				+ f(ci) * (ai - bi);

			// Step 6: error if m <= 0
			if (m <= 0) {
				throw std::string("lag: m <= 0, interpolation error");
			}

			// Step 9
			double di = 0.5 * l / m;

			// Step 10: if ai < di < ci
			if (ai < di && di < ci) {
				// Step 11
				if (f(di) < f(ci)) {
					// Step 12-14
					double ai_next = ai;
					double ci_next = di;
					double bi_next = ci;
					ai = ai_next;
					bi = bi_next;
					ci = ci_next;
				}
				else {
					// Step 16-18
					double ai_next = di;
					double ci_next = ci;
					double bi_next = bi;
					ai = ai_next;
					bi = bi_next;
					ci = ci_next;
				}
			}
			else {
				// Step 20: exit if interval too small
				break;
			}

			// Stopping criterion (not specified; could use |ci - ai| < epsilon)
			if (fabs(ci - ai) < epsilon) {
				break;
			}

			i++;
		}

		// Save best found so far
		matrix xopt_mat(1, 1);
		xopt_mat(0, 0) = ci;
		Xopt.x = xopt_mat;
		Xopt.y = f(ci);
		Xopt.flag = 1;

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
