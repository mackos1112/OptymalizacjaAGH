//Ten plik nie powinien byæ edytowany

#include"matrix.h"

matrix::matrix(double L)
{
	n = m = 1;
	M = new double*[1];
	M[0] = new double[1];
	M[0][0] = L;
}

matrix::matrix(int nv, int mv, double L)
{
	if (nv <= 0 || mv <= 0)
		throw string("matrix::matrix(int,int,double):\nwymiary macierzy musza byc dodatnie");
	n = nv;
	m = mv;
	M = new double*[n];
	for (int i = 0; i < n; ++i)
	{
		M[i] = new double[m];
		for (int j = 0; j < m; ++j)
			M[i][j] = L;
	}
}

matrix::matrix(int nv, double* A)
{
	if (nv <= 0)
		throw string("matrix::matrix(int,double*):\ndlugosc wektora musi byc dodatnia");
	n = nv;
	m = 1;
	M = new double*[n];
	for (int i = 0; i < n; ++i)
	{
		M[i] = new double[1];
		M[i][0] = A[i];
	}
}

matrix::matrix(int nv, int mv, double** A)
{
	if (nv <= 0 || mv <= 0)
		throw string("matrix::matrix(int,int,double**):\nwymiary macierzy musza byc dodatnie");
	n = nv;
	m = mv;
	M = new double*[n];
	for (int i = 0; i < n; ++i)
	{
		M[i] = new double[m];
		for (int j = 0; j < m; ++j)
			M[i][j] = A[i][j];
	}
}

matrix::matrix(const matrix& A)
{
	n = A.n;
	m = A.m;
	M = new double*[n];
	for (int i = 0; i < n; ++i)
	{
		M[i] = new double[m];
		for (int j = 0; j < m; ++j)
			M[i][j] = A.M[i][j];
	}
}

matrix::~matrix()
{
	for (int i = 0; i < n; ++i)
		delete[] M[i];
	delete[] M;
}

matrix& matrix::operator=(const matrix& A)
{
	if (&A == this)
		return *this;
	for (int i = 0; i < n; ++i)
		delete[] M[i];
	delete[] M;
	n = A.n;
	m = A.m;
	M = new double*[n];
	for (int i = 0; i < n; ++i)
	{
		M[i] = new double[m];
		for (int j = 0; j < m; ++j)
			M[i][j] = A.M[i][j];
	}
	return *this;
}

matrix matrix::operator[](int nv) const
{
	if (nv >= m || nv < 0)
		throw string("matrix matrix::operator[](int nv) const:\nnumer kolmuny jest poza zakresem");
	matrix A(n, 1);
	for (int i = 0; i < n; ++i)
		A.M[i][0] = M[i][nv];
	return A;
}

double& matrix::operator()(int nv, int mv)
{
	if (nv >= n || mv >= m || nv < 0 || mv < 0)
		throw string("double& matrix::operator()(int,int):\nindeks jest poza zakresem");
	return M[nv][mv];
}

double matrix::operator()(int nv, int mv) const
{
	if (nv >= n || mv >= m || nv < 0 || mv < 0)
		throw string("double matrix::operator()(int,int) const:\nindeks jest poza zakresem");
	return M[nv][mv];
}

void matrix::set_col(const matrix& c, int mv)
{
	if (mv >= m || mv < 0)
		throw string("void matrix::set_col(const matrix&,int):\nnumer kolmuny jest poza zakresem");
	if (n != c.n)
		throw string("void matrix::set_col(const matrix&,int):\nliczba wierszy macierzy musi byc rowna dlugosci wektora");
	if (c.m != 1)
		throw string("void matrix::set_col(const matrix&,int):\nwstawiana kolumna musi miec postac wektora pionowego");
	for (int i = 0; i < n; ++i)
		M[i][mv] = c(i);
}

void matrix::set_row(const matrix& c, int nv)
{
	if (nv >= n || nv < 0)
		throw string("void matrix::set_row(const matrix&,int):\nnumer wiersza jest poza zakresem");
	if (m != c.m)
		throw string("void matrix::set_row(const matrix&,int):\nliczba kolumn macierzy musi byc rowna dlugosci wektora");
	if (c.n != 1)
		throw string("void matrix::set_row(const matrix&,int):\nwstawiany wiersz musi miec postac wektora poziomego");
	for (int i = 0; i < m; ++i)
		M[nv][i] = c(0, i);
}

void matrix::add_col(double L)
{
	matrix A(n, m + 1);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			A(i, j) = M[i][j];
		A(i, m) = L;
	}
	*this = A;
}

void matrix::add_row(double L)
{
	matrix A(n + 1, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			A(i, j) = M[i][j];
	for (int j = 0; j < m; ++j)
		A(n, j) = L;
	*this = A;
}

void matrix::add_col(const matrix& c)
{
	try
	{
		matrix A(n, m + 1);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				A(i, j) = M[i][j];
		}
		A.set_col(c, m);
		*this = A;
	}
	catch (string ex_info)
	{
		throw ("void matrix::add_col(const matrix&):\n" + ex_info);
	}
}

void matrix::add_row(const matrix& c)
{
	try
	{
		matrix A(n + 1, m);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				A(i, j) = M[i][j];
		A.set_row(c, n);
		*this = A;
	}
	catch (string ex_info)
	{
		throw ("void matrix::add_row(const matrix&):\n" + ex_info);
	}
}

matrix operator+(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] == 1 && nA[1] == 1)
	{
		matrix C(B);
		for (int i = 0; i < nB[0]; ++i)
			for (int j = 0; j < nB[1]; ++j)
				C(i, j) += A();
		return C;
	}
	else if (nB[0] == 1 && nB[1] == 1)
	{
		matrix C(A);
		for (int i = 0; i < nA[0]; ++i)
			for (int j = 0; j < nA[1]; ++j)
				C(i, j) += B();
		return C;
	}
	else if (nA[0] == nB[0] && nA[1] == nB[1])
	{
		matrix C(A);
		for (int i = 0; i < nA[0]; ++i)
			for (int j = 0; j < nA[1]; ++j)
				C(i, j) += B(i, j);
		return C;
	}
	else
		throw string("matrix operator+(const matrix&, const matrix&):\nwymiary macierzy nie sa zgodne");
}

matrix operator-(const matrix& A, const matrix& B)
{
	try
	{
		return A + (-B);
	}
	catch (string ex_info)
	{
		throw ("matrix operator-(const matrix&, const matrix&):\n" + ex_info);
	}
}

matrix operator*(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] == 1 && nA[1] == 1)
	{
		matrix C(B);
		for (int i = 0; i < nB[0]; ++i)
			for (int j = 0; j < nB[1]; ++j)
				C(i, j) *= A();
		return C;
	}
	else if (nB[0] == 1 && nB[1] == 1)
	{
		matrix C(A);
		for (int i = 0; i < nA[0]; ++i)
			for (int j = 0; j < nA[1]; ++j)
				C(i, j) *= B();
		return C;
	}
	else if (nA[1] == nB[0])
	{
		matrix C(nA[0], nB[1]);
		for (int i = 0; i < nA[0]; ++i)
			for (int j = 0; j < nB[1]; ++j)
				for (int k = 0; k < nA[1]; ++k)
					C(i, j) += A(i, k) * B(k, j);
		return C;
	}
	else
		throw string("matrix operator*(const matrix&, const matrix&):\nwymiary macierzy nie sa zgodne");
}

matrix operator/(const matrix& A, const matrix& B)
{
	try
	{
		return A * inv(B);
	}
	catch (string ex_info)
	{
		throw ("matrix operator/(const matrix&, const matrix&):\n" + ex_info);
	}
}

matrix operator-(const matrix& A)
{
	int* n = get_size(A);
	matrix B(A);
	for (int i = 0; i < n[0]; ++i)
		for (int j = 0; j < n[1]; ++j)
			B(i, j) = -A(i, j);
	return B;
}

bool operator<(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] != 1 || nA[1] != 1 || nB[0] != 1 || nB[1] != 1)
		throw string("bool operator<(const matrix&, const matrix&):\noperator relacji jest zdefiniwany tylko dla macierzy 1x1");
	return A() < B();
}

bool operator>(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] != 1 || nA[1] != 1 || nB[0] != 1 || nB[1] != 1)
		throw string("bool operator>(const matrix&, const matrix&):\noperator relacji jest zdefiniwany tylko dla macierzy 1x1");
	return A() > B();
}

bool operator<=(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] != 1 || nA[1] != 1 || nB[0] != 1 || nB[1] != 1)
		throw string("bool operator<=(const matrix&, const matrix&):\noperator relacji jest zdefiniwany tylko dla macierzy 1x1");
	return A() <= B();
}

bool operator>=(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] != 1 || nA[1] != 1 || nB[0] != 1 || nB[1] != 1)
		throw string("bool operator>=(const matrix&, const matrix&):\noperator relacji jest zdefiniwany tylko dla macierzy 1x1");
	return A() >= B();
}

bool operator==(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] != 1 || nA[1] != 1 || nB[0] != 1 || nB[1] != 1)
		throw string("bool operator==(const matrix&, const matrix&):\noperator relacji jest zdefiniwany tylko dla macierzy 1x1");
	return A() == B();
}

bool operator!=(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] != 1 || nA[1] != 1 || nB[0] != 1 || nB[1] != 1)
		throw string("bool operator!=(const matrix&, const matrix&):\noperator relacji jest zdefiniwany tylko dla macierzy 1x1");
	return A() != B();
}

matrix ident_mat(int nv)
{
	try
	{
		matrix A(nv, nv);
		for (int i = 0; i < nv; ++i)
			A(i, i) = 1;
		return A;
	}
	catch (string ex_info)
	{
		throw ("matrix ident_mat(int):\n" + ex_info);
	}
}

matrix rand_mat(int nv, int mv)
{
	try
	{
		matrix A(nv, mv);
		random_device R;
		for (int i = 0; i < nv; ++i)
			for (int j = 0; j < mv; ++j)
				A(i, j) = 1.0 * R() / R.max();
		return A;
	}
	catch (string ex_info)
	{
		throw ("matrix rand_mat(int,int):\n" + ex_info);
	}
}

matrix randn_mat(int nv, int mv)
{
	try
	{
		matrix A(nv, mv);
		random_device rd;
		default_random_engine gen;
		gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
		normal_distribution<double> distr(0.0, 1.0);
		for (int i = 0; i < nv; ++i)
			for (int j = 0; j < mv; ++j)
				A(i, j) = distr(gen);
		return A;
	}
	catch (string ex_info)
	{
		throw ("matrix randn_mat(int,int):\n" + ex_info);
	}
}

double m2d(const matrix& A)
{
	int* nA = get_size(A);
	if (nA[0] != 1 || nA[1] != 1)
		throw string("double m2d(const matrix&):\nzamiana macierzy na liczbe mozliwa jest tylko dla skalarow");
	return A(0, 0);
}

double det(const matrix& A)
{
	int* nA = get_size(A);
	if (nA[0] != nA[1])
		throw string("double det(const matrix&):\nmacierz musi byc kwadratowa");
	double D = 0;
	if (nA[0] == 1)
		D = A();
	else
	{
		for (int k = 0; k < nA[0]; ++k)
		{
			matrix T(nA[0] - 1, nA[0] - 1);
			for (int i = 0; i < nA[0] - 1; ++i)
				for (int j = 0; j < nA[1] - 1; ++j)
					T(i, j) = A(i + 1, j >= k ? j + 1 : j);
			D = D + A(0, k) * pow(-1.0, k) * det(T);
		}
	}
	return D;
}

matrix inv(const matrix& A)
{
	try
	{
		double D = det(A);
		if (D == 0)
			throw string("matrix inv(const matrix&):\nwyznacznik macierzy wynosi 0");
		int* nA = get_size(A);
		matrix I(nA[0], nA[0]);
		if (nA[0] == 1)
			I() = 1 / A();
		else
		{
			for (int k = 0; k < nA[0]; ++k)
				for (int l = 0; l < nA[1]; ++l)
				{
					matrix T(nA[0] - 1, nA[0] - 1);
					for (int i = 0; i < nA[0] - 1; ++i)
						for (int j = 0; j < nA[1] - 1; ++j)
							T(i, j) = A(i >= k ? i + 1 : i, j >= l ? j + 1 : j);
					I(k, l) = pow(-1.0, k + l) * det(T);
				}
			I = 1 / D * trans(I);
		}
		return I;
	}
	catch (string ex_info)
	{
		throw ("matrix inv(const matrix&):\n" + ex_info);
	}
}

matrix trans(const matrix& A)
{
	int* nA = get_size(A);
	matrix B(nA[1], nA[0]);
	for (int i = 0; i < nA[1]; ++i)
		for (int j = 0; j < nA[0]; ++j)
			B(i, j) = A(j, i);
	return B;
}

matrix pow(const matrix& A, int n)
{
	if (n < 0)
		throw string("matrix pow(const matrix&,int):\nwykladnik potegi nie moze byc ujemny");
	int* nA = get_size(A);
	if (nA[0] != nA[1])
		throw string("matrix pow(const matrix&,int):\npotegowanie jest mozliwe tylko dla macierzy kwadratowych");
	matrix B = ident_mat(nA[0]);
	for (int i = 1; i <= n; ++i)
		B = B * A;
	return B;
}

double norm(const matrix& A)
{
	int* nA = get_size(A);
	if (nA[1] != 1)
		throw string("double norm(const matrix&):\nnorma jest zdefiniowana tylko dla wektorow pionowych");
	double N = 0;
	for (int i = 0; i < nA[0]; ++i)
		N += pow(A(i), 2);
	return sqrt(N);
}

matrix hcat(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[0] != nB[0])
		throw string("matrix hcat(const matrix&,const matrix&):\nliczba wierszy macierzy musi byc taka sama");
	matrix C(nA[0], nA[1] + nB[1]);
	for (int i = 0; i < nA[0]; ++i)
	{
		for (int j = 0; j < nA[1]; ++j)
			C(i, j) = A(i, j);
		for (int j = 0; j < nB[1]; ++j)
			C(i, j + nA[1]) = B(i, j);
	}
	return C;
}

matrix vcat(const matrix& A, const matrix& B)
{
	int* nA = get_size(A);
	int* nB = get_size(B);
	if (nA[1] != nB[1])
		throw string("matrix vcat(const matrix&,const matrix&):\nliczba kolumn macierzy musi byc taka sama");
	matrix C(nA[0] + nB[0], nA[1]);
	for (int i = 0; i < nA[0]; ++i)
		for (int j = 0; j < nA[1]; ++j)
			C(i, j) = A(i, j);
	for (int i = 0; i < nB[0]; ++i)
		for (int j = 0; j < nB[1]; ++j)
			C(i + nA[0], j) = B(i, j);
	return C;
}

matrix get_col(const matrix& A, int mv)
{
	int* n = get_size(A);
	if (mv >= n[1] || mv < 0)
		throw string("matrix get_col(const matrix&,int):\nnumer kolmuny jest poza zakresem");
	matrix B(n[0], 1);
	for (int i = 0; i < n[0]; ++i)
		B(i, 0) = A(i, mv);
	return B;
}

matrix get_row(const matrix& A, int nv)
{
	int* n = get_size(A);
	if (nv >= n[0] || nv < 0)
		throw string("matrix get_row(const matrix&,int):\nnumer wiersza jest poza zakresem");
	matrix B(1, n[1]);
	for (int j = 0; j < n[1]; ++j)
		B(0, j) = A(nv, j);
	return B;
}

ostream& operator<<(ostream& OS, const matrix& A)
{
	int* nA = get_size(A);
	ostringstream OSS;
	string S;
	string::size_type p;
	for (int i = 0; ; ++i)
	{
		for (int j = 0; j < nA[1]; ++j)
		{
			OSS << A(i, j);
			S = OSS.str();
			OSS.str("");
			p = S.find('.');
			if (p != string::npos)
				S[p] = SEP_SYMBOL;
			OS << S << "; ";
		}
		if (i == nA[0] - 1)
			return OS;
		OS << endl;
	}
}

istream& operator>>(istream& IS, matrix& A)
{
	istringstream ISS;
	string S;
	string::size_type p;
	int* nA = get_size(A);
	for (int i = 0; i < nA[0]; ++i)
		for (int j = 0; j < nA[1]; ++j)
		{
			getline(IS, S, ';');
			p = S.find(SEP_SYMBOL);
			if (p != string::npos)
				S[p] = '.';
			ISS.str(S);
			ISS >> A(i, j);
			if (ISS.fail())
				throw string("istream& operator>>(istream&,matrix&):\nblad podczas odczytu macierzy");
			ISS.clear();
			if (IS.eof())
				throw string("istream& operator>>(istream&,matrix&):\nzbyt malo liczb");
		}
	return IS;
}

int* get_size(const matrix& A)
{
	int* s = new int[2]{ A.n, A.m };
	return s;
}

int get_len(const matrix& A)
{
	if (A.m != 1)
		throw string("int get_len(const matrix&):\ndlugosc jest zwracana tylko dla wektorow pionowych");
	return A.n;
}
