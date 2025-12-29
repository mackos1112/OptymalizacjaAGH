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
				std::cout << "[DEBUG] f(x1) == f(x0) -> zwracam [" << p[0] << ", " << p[1] << "]" << std::endl;

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
					std::cout << "[DEBUG] f(x1) >= f(x0) -> zwracam [" << p[0] << ", " << p[1] << "]" << std::endl;

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

		// helper f(x) with caching
		auto f = [&](double x) {
			matrix X(1,1);
			X(0,0) = x;
			f_calls++;
			return ff(X, ud1, ud2)(0,0);
			};

		// build fibonacci sequence until φ_k > (b-a)/epsilon
		std::vector<long long> F = {1,1};
		while (F.back() < (b - a) / epsilon)
			F.push_back(F[F.size()-1] + F[F.size()-2]);
		int k = F.size() - 1;

		double a_i = a;
		double b_i = b;

		double c_i = b_i - (double)F[k-1] / F[k] * (b_i - a_i);
		double d_i = a_i + b_i - c_i;

		double fc = f(c_i);
		double fd = f(d_i);

		for (int i = 0; i < k-2; i++)
		{
			if (fc < fd) {
				b_i = d_i;
				d_i = c_i;
				fd = fc;
				c_i = b_i - (double)F[k-i-2] / F[k-i-1] * (b_i - a_i);
				fc = f(c_i);
			} else {
				a_i = c_i;
				c_i = d_i;
				fc = fd;
				d_i = a_i + b_i - c_i;
				fd = f(d_i);
			}
		}

		Xopt.x = c_i;
		Xopt.y = fc;
		Xopt.f_calls = f_calls;
		Xopt.flag = 0; // local minimum

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

		double ai = a;
		double bi = b;
		double ci = (a + b) / 2.0;  // initial midpoint
		int f_calls = 0;

		// helper f(x) with counter
		auto f = [&](double x) -> double {
			matrix X(1,1);
			X(0,0) = x;
			f_calls++;
			return ff(X, ud1, ud2)(0,0);
			};

		// cache initial values
		double fa = f(ai);
		double fb = f(bi);
		double fc = f(ci);

		while (f_calls < Nmax)
		{
			double l = fa * ((bi*bi) - (ci*ci))
				+ fb * ((ci*ci) - (ai*ai))
				+ fc * ((ai*ai) - (bi*bi));

			double m = fa * (bi - ci)
				+ fb * (ci - ai)
				+ fc * (ai - bi);

			if (m <= 1e-12) break;  // avoid division by zero / negative

			double di = 0.5 * l / m;
			double fdi = f(di);

			if (ai < di && di < ci)
			{
				if (fdi < fc)
				{
					bi = ci; fb = fc;
					ci = di; fc = fdi;
				}
				else
				{
					ai = di; fa = fdi;
				}
			}
			else
			{
				// if di outside interval, reduce interval
				if (fabs(ci - ai) < epsilon) break;
				// optionally shrink interval slightly
				ai = ai; // keep same for now
				break;
			}

			if (fabs(ci - ai) < epsilon) break;
		}

		Xopt.x = ci;
		Xopt.y = fc;
		Xopt.f_calls = f_calls;
		Xopt.flag = 0; // local minimum

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

		// Nelder-Mead simplex method (sympleks Neldera-Meada) implementation
		int n = 2; // dimension of the problem

		// initialize simplex points p0 = x0, pi = x0 + s*ei for i=1..n
		std::vector<matrix> P(n + 1);
		P[0] = x0;
		for (int i = 1; i <= n; ++i)
		{
			P[i] = x0;
			// add s to the (i-1)-th coordinate
			P[i](i - 1, 0) = P[i](i - 1, 0) + s;
		}

		int f_calls = 0;
		auto f = [&](const matrix& X) -> double {
			ud2(0, 0) = f_calls;
			matrix Y = ff(X, ud1, ud2);
			++f_calls;
			return Y(0, 0);
			};

		// compute initial function values
		std::vector<double> fv(n + 1);
		for (int i = 0; i <= n; ++i) fv[i] = f(P[i]);

		while (true)
		{
			// determine pmin and pmax
			int idx_min = 0, idx_max = 0;
			for (int i = 1; i <= n; ++i)
			{
				if (fv[i] < fv[idx_min]) idx_min = i;
				if (fv[i] > fv[idx_max]) idx_max = i;
			}

			// if all values equal, break (can't progress)
			if (idx_min == idx_max) break;

			// compute centroid p = (sum_{i != max} pi) / n
			matrix p(n, 1, 0.0);
			for (int i = 0; i <= n; ++i)
			{
				if (i == idx_max) continue;
				for (int j = 0; j < n; ++j)
					p(j, 0) += P[i](j, 0);
			}
			for (int j = 0; j < n; ++j) p(j, 0) /= static_cast<double>(n);

			// reflection: podb = p + alpha*(p - pmax)
			matrix temp = p - P[idx_max];
			for (int j = 0; j < n; ++j) temp(j, 0) *= alpha;
			matrix podb = p + temp;
			double f_podb = f(podb);

			if (f_podb < fv[idx_min])
			{
				// expansion: pe = p + gamma*(podb - p)
				matrix temp2 = podb - p;
				for (int j = 0; j < n; ++j) temp2(j, 0) *= gamma;
				matrix pe = p + temp2;
				double f_pe = f(pe);

				if (f_pe < f_podb)
				{
					P[idx_max] = pe;
					fv[idx_max] = f_pe;
				}
				else
				{
					P[idx_max] = podb;
					fv[idx_max] = f_podb;
				}
			}
			else
			{
				if (fv[idx_min] <= f_podb && f_podb < fv[idx_max])
				{
					// accept reflection
					P[idx_max] = podb;
					fv[idx_max] = f_podb;
				}
				else
				{
					// contraction: pz = p + beta*(pmax - p)
					matrix temp3 = P[idx_max] - p;
					for (int j = 0; j < n; ++j) temp3(j, 0) *= beta;
					matrix pz = p + temp3;
					double f_pz = f(pz);

					if (f_pz >= fv[idx_max])
					{
						// shrink towards pmin: pi = pmin + delta*(pi - pmin) for i != min
						for (int i = 0; i <= n; ++i)
						{
							if (i == idx_min) continue;
							matrix diff = P[i] - P[idx_min];
							for (int j = 0; j < n; ++j) diff(j, 0) *= delta;
							P[i] = P[idx_min] + diff;
							fv[i] = f(P[i]);
						}
					}
					else
					{
						P[idx_max] = pz;
						fv[idx_max] = f_pz;
					}
				}
			}

			if (f_calls > Nmax)
				throw std::string("sym_NM: przekroczono maksymalna liczbe wywolan funkcji celu");

			// check termination: max_i ||pmin - pi||2 < epsilon
			double max_dist = 0.0;
			for (int i = 0; i <= n; ++i)
			{
				matrix diff = P[i] - P[idx_min];
				double d = norm(diff);
				if (d > max_dist) max_dist = d;
			}
			if (max_dist < epsilon) break;
		}

		// return best found point
		int best = 0;
		for (int i = 1; i <= n; ++i) if (fv[i] < fv[best]) best = i;

		Xopt.x = P[best];
		Xopt.y = ff(Xopt.x, ud1, ud2);
		Xopt.f_calls = f_calls;
		Xopt.flag = 0; // local minimum

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

		// Implementation based on provided pseudocode.
		// Parameters:
		// N       - number of variables
		// lb, ub  - lower/upper bounds (matrices of length N)
		// mi      - population size (mu)
		// lambda  - number of offspring
		// sigma0  - initial sigma (either scalar or vector of length N)
		// epsilon - target function value threshold
		// Nmax    - max allowed function evaluations
		//
		// P(i) contains pairs (solution, sigma matrix)

		// helper to extract scalar from 1x1 matrix safely
		auto mat2d_safe = [&](const matrix& M)->double {
			int* sz = get_size(M);
			if (sz[0] != 1 || sz[1] != 1)
				throw string("EA: expected scalar matrix");
			return M();
		};

		struct Individual {
			solution s;
			matrix sigma; // Nx1 vector
		};

		// compute alpha, beta (as in pseudocode: α = N^-0.5, β = (2N)^-0.25)
		double alpha = std::pow((double)N, -0.5);
		double beta = std::pow(2.0 * (double)N, -0.25);

		// random generators
		std::random_device rd;
		std::default_random_engine gen(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
		std::uniform_real_distribution<double> unif01(0.0, 1.0);
		std::normal_distribution<double> normal01(0.0, 1.0);

		// initialize population P(0)
		std::vector<Individual> P;
		P.reserve(mi);

		for (int i = 0; i < mi; ++i) {
			Individual ind;
			// random x in [lb, ub]
			matrix r = rand_mat(N, 1);
			ind.s.x = matrix(N, 1);
			for (int k = 0; k < N; ++k) {
				double rv = r(k, 0);
				ind.s.x(k, 0) = (ub(k) - lb(k)) * rv + lb(k);
			}
			// sigma initialization: if sigma0 is scalar use that, else copy vector
			int* ssz = get_size(sigma0);
			if (ssz[0] == 1 && ssz[1] == 1) {
				ind.sigma = matrix(N, 1);
				for (int k = 0; k < N; ++k) ind.sigma(k, 0) = sigma0();
			}
			else {
				// try to copy regardless of orientation (assume Nx1 or 1xN)
				ind.sigma = matrix(N, 1);
				for (int k = 0; k < N; ++k) {
					// if sigma0 has matching dimension, use it, otherwise fallback to first element
					int* ssz2 = get_size(sigma0);
					if (ssz2[0] == N && ssz2[1] == 1)
						ind.sigma(k, 0) = sigma0(k, 0);
					else if (ssz2[0] == 1 && ssz2[1] == N)
						ind.sigma(k, 0) = sigma0(0, k);
					else
						ind.sigma(k, 0) = sigma0(0, 0);
				}
			}

			// evaluate fitness
			ind.s.fit_fun(ff, ud1, ud2);

			P.push_back(ind);
			if (solution::f_calls > Nmax) {
				Xopt = P[0].s;
				Xopt.flag = 0;
				return Xopt;
			}
		}

		int iter = 0;
		while (true) {
			// compute φ_j = 1 / f(x_j)
			// to avoid division by zero or negative values, shift by min_f if needed
			double min_f = std::numeric_limits<double>::infinity();
			for (int j = 0; j < mi; ++j) {
				double fv = m2d(P[j].s.y);
				if (fv < min_f) min_f = fv;
			}
			double shift = 0.0;
			if (min_f <= 0.0) shift = -min_f + 1e-12;

			std::vector<double> phi(mi);
			double Phi = 0.0;
			for (int j = 0; j < mi; ++j) {
				double fv = m2d(P[j].s.y);
				phi[j] = 1.0 / (fv + shift + 1e-12);
				Phi += phi[j];
			}

			// build cumulative distribution q (roulette)
			std::vector<double> q(mi);
			double cum = 0.0;
			for (int j = 0; j < mi; ++j) {
				cum += phi[j] / Phi;
				q[j] = cum;
			}
			// ensure last q is exactly 1
			q[mi-1] = 1.0;

			// draw global a (as in pseudocode line 14)
			double a_global = normal01(gen);

			// generate lambda offspring T
			std::vector<Individual> T;
			T.reserve(lambda);
			for (int j = 0; j < lambda; ++j) {
				// select parent A
				double r = unif01(gen);
				int k = 0;
				while (k < mi && r > q[k]) ++k;
				if (k >= mi) k = mi - 1;
				Individual A = P[k];

				// select parent B
				r = unif01(gen);
				k = 0;
				while (k < mi && r > q[k]) ++k;
				if (k >= mi) k = mi - 1;
				Individual B = P[k];

				// crossover scalar r_c
				double r_c = unif01(gen);
				Individual child;
				child.s = A.s;
				child.s.x = matrix(N, 1);
				child.sigma = matrix(N, 1);
				// combine x and sigma
				for (int d = 0; d < N; ++d) {
					double Ax = A.s.x(d, 0);
					double Bx = B.s.x(d, 0);
					child.s.x(d, 0) = r_c * Ax + (1.0 - r_c) * Bx;

					double As = A.sigma(d, 0);
					double Bs = B.sigma(d, 0);
					child.sigma(d, 0) = r_c * As + (1.0 - r_c) * Bs;
				}

				// mutation: update sigma with global a and per-dimension b
				for (int d = 0; d < N; ++d) {
					double b_sigma = normal01(gen);
					child.sigma(d, 0) = child.sigma(d, 0) * std::exp(alpha * a_global + beta * b_sigma);
					// ensure sigma stays positive and not too small
					if (child.sigma(d, 0) < 1e-16) child.sigma(d, 0) = 1e-16;
				}

				// mutation: mutate x using per-dimension normal
				for (int d = 0; d < N; ++d) {
					double b_x = normal01(gen);
					child.s.x(d, 0) = child.s.x(d, 0) + b_x * child.sigma(d, 0);
					// enforce bounds
					if (child.s.x(d, 0) < lb(d)) child.s.x(d, 0) = lb(d);
					if (child.s.x(d, 0) > ub(d)) child.s.x(d, 0) = ub(d);
				}

				// evaluate child
				child.s.fit_fun(ff, ud1, ud2);

				T.push_back(child);

				if (solution::f_calls > Nmax) {
					// return best found so far with flag = 0 (error / exceeded)
					// find best in current union P U T
					std::vector<Individual> unionPT = P;
					for (auto &it : T) unionPT.push_back(it);
					// find best (min y)
					int bestIdx = 0;
					double bestVal = m2d(unionPT[0].s.y);
					for (size_t ii = 1; ii < unionPT.size(); ++ii) {
						double vv = m2d(unionPT[ii].s.y);
						if (vv < bestVal) { bestVal = vv; bestIdx = ii; }
					}
					Xopt = unionPT[bestIdx].s;
					Xopt.flag = 0;
					return Xopt;
				}
			} // end offspring generation

			// P(i+1) = select mu best from P U T
			std::vector<Individual> unionPT = P;
			for (auto &it : T) unionPT.push_back(it);

			// sort union by objective (ascending)
			std::sort(unionPT.begin(), unionPT.end(), [&](const Individual& A, const Individual& B) {
				return m2d(A.s.y) < m2d(B.s.y);
			});

			// take first mi individuals as new population
			P.clear();
			for (int i = 0; i < mi; ++i) P.push_back(unionPT[i]);

			// find best x* in new population
			int bestIdx = 0;
			double bestVal = m2d(P[0].s.y);
			for (int i = 1; i < mi; ++i) {
				double vv = m2d(P[i].s.y);
				if (vv < bestVal) { bestVal = vv; bestIdx = i; }
			}
			Xopt = P[bestIdx].s;

			// check stopping criteria
			if (m2d(Xopt.y) < epsilon) {
				Xopt.flag = 1;
				return Xopt;
			}
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}

			iter++;
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
