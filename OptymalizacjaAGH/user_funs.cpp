#include"user_funs.h"
#include <cmath>

#define MATH_PI 3.1415926

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartosc funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wspolrzedne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartosc funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki poczatkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment sily dzialajacy na wahadlo oraz czas dzialania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwiazujemy rownanie rozniczkowe
	int n = get_len(Y[0]);									// dlugosc rozwiazania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahadla
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartosc funkcji celu (ud1 to załozone maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamieci rozwiazanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z polozenia to predkosc
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z predkosci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;   // y zawiera wartosc funkcji celu

	// pobieramy wartosc x jako double
	double xv = m2d(x);

	// obliczamy wartosc funkcji:
	// f(x) = -cos(0.1*x) * exp((-0.1*x - 2*pi)^2) + 0.002*(0.1*x)^2
	double fx = -cos(0.1 * xv) * exp(-pow(0.1 * xv - 2 * 3.14, 2))
		+ 0.002 * pow(0.1 * xv, 2);

	y = matrix(1, 1);   // wynik jako macierz 1x1
	y(0, 0) = fx;

	return y;
}

// Funkcja rozniczkowa dla modelu zbiornikow
// Y = [VA, VB, TA, TB] - objetosci i temperatury w zbiornikach A i B
// ud2 = [DA, DB] - pola przekroju otworow (w cm², konwertowane na m²)
matrix df_tanks(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(4, 1);

	// Parametry stałe
	double PA = 2.0;          // m² - pole podstawy zbiornika A
	double PB = 1.0;          // m² - pole podstawy zbiornika B
	double FBin = 0.01;       // m³/s - szybkość dopływu wody zewnętrznej (10 l/s = 0.01 m³/s)
	double TBin = 20.0;       // °C - temperatura wody zewnętrznej
	double a = 0.98;          // współczynnik lepkości
	double b = 0.63;          // współczynnik zwężenia strumienia
	double g = 9.81;          // m/s² - przyspieszenie ziemskie

	// Konwersja pól przekroju z cm² na m²
	// Sprawdzenie czy ud2 ma odpowiedni rozmiar
	int* size_ud2 = get_size(ud2);
	double DA, DB;
	if (size_ud2[0] >= 2 && size_ud2[1] >= 1)
	{
		DA = ud2(0) / 10000.0;  // cm² -> m²
		DB = ud2(1) / 10000.0;  // cm² -> m²
	}
	else
	{
		// Wartości domyślne jeśli ud2 jest puste
		DA = 50.0 / 10000.0;  // domyślnie 50 cm²
		DB = 36.5665 / 10000.0;  // domyślnie 36.5665 cm²
	}
	delete[] size_ud2;

	// Stan układu
	double VA = Y(0);  // objętość wody w zbiorniku A
	double VB = Y(1);  // objętość wody w zbiorniku B
	double TA = Y(2);  // temperatura wody w zbiorniku A
	double TB = Y(3);  // temperatura wody w zbiorniku B

	// Wypływ z zbiornika A
	double dVA_out = -a * b * DA * sqrt(2.0 * g * VA / PA);
	if (VA <= 0) dVA_out = 0;  // zabezpieczenie przed ujemnymi objętościami

	// Wypływ z zbiornika B
	double dVB_out = -a * b * DB * sqrt(2.0 * g * VB / PB);
	if (VB <= 0) dVB_out = 0;

	// Zmiana objętości w zbiorniku A
	dY(0) = dVA_out;

	// Zmiana objętości w zbiorniku B
	dY(1) = -dVA_out + FBin + dVB_out;  // dopływ z A (-dVA_out jest dodatnie) + dopływ zewnętrzny + wypływ z B (dVB_out jest ujemne)

	// Temperatura w zbiorniku A się nie zmienia (tylko wypływ)
	dY(2) = 0.0;

	// Zmiana temperatury w zbiorniku B przez mieszanie
	if (VB > 1e-10)  // unikamy dzielenia przez zero
	{
		// Objętość wpływająca z A na sekundę
		double Vin_A = -dVA_out;  // dodatnia wartość dopływu z A
		if (Vin_A < 0) Vin_A = 0;

		// Całkowita objętość wpływająca
		double Vin_total = Vin_A + FBin;

		// Średnia ważona temperatura wpływającej wody
		double Tin_avg = 0.0;
		if (Vin_total > 1e-10)
		{
			Tin_avg = (Vin_A * TA + FBin * TBin) / Vin_total;
		}
		else
		{
			Tin_avg = TB;  // jeśli brak dopływu, temperatura się nie zmienia
		}

		// Zmiana temperatury zgodnie ze wzorem: dT/dt = (Vin/V) * (Tin - T)
		dY(3) = (Vin_total / VB) * (Tin_avg - TB);
	}
	else
	{
		dY(3) = 0.0;
	}

	return dY;
}

// Funkcja celu dla optymalizacji pola przekroju DA
// x - pole przekroju DA w cm²
// Celem jest znalezienie DA takiego, że max(TB) = 50°C
matrix ff_tanks(matrix x, matrix ud1, matrix ud2)
{
	matrix y;

	// Parametry
	double PA = 2.0;          // m²
	double PB = 1.0;          // m²
	double VA0 = 5.0;         // m³ - początkowa objętość w A
	double VB0 = 1.0;         // m³ - początkowa objętość w B
	double TA0 = 95.0;        // °C - początkowa temperatura w A
	double TB0 = 20.0;        // °C - początkowa temperatura w B
	double DB = 36.5665;      // cm² - pole przekroju otworu w B
	double t0 = 0.0;          // s
	double dt = 1.0;          // s
	double tend = 2000.0;     // s
	double T_target = 50.0;   // °C - docelowa maksymalna temperatura

	// DA z wektora optymalizacji (w cm²)
	double DA = m2d(x);

	// Warunki początkowe: [VA, VB, TA, TB]
	matrix Y0(4, 1);
	Y0(0) = VA0;
	Y0(1) = VB0;
	Y0(2) = TA0;
	Y0(3) = TB0;

	// Parametry dla funkcji różniczkowej: [DA, DB] w cm²
	matrix params(2, 1);
	params(0) = DA;
	params(1) = DB;

	// Rozwiązanie równania różniczkowego
	matrix* Y = solve_ode(df_tanks, t0, dt, tend, Y0, ud1, params);

	// Znajdź maksymalną temperaturę w zbiorniku B
	// Po transpozycji w solve_ode: Y[1] ma format (kroki × zmienne)
	// Zmienne: 0=VA, 1=VB, 2=TA, 3=TB
	int n = get_len(Y[0]);
	double T_max = Y[1](0, 3);  // TB w pierwszym kroku (krok 0, zmienna TB=3)
	for (int i = 1; i < n; ++i)
	{
		if (T_max < Y[1](i, 3))  // krok i, zmienna TB=3
			T_max = Y[1](i, 3);
	}

	// Funkcja celu: różnica między maksymalną temperaturą a docelową
	y = abs(T_max - T_target);

	// Zwolnienie pamięci
	Y[0].~matrix();
	Y[1].~matrix();

	return y;
}

matrix ff2T(matrix x1, matrix x2, matrix ud1)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartosc funkcji celu
	matrix pom = MATH_PI * pow(pow(m2d(x1) / MATH_PI, 2) + pow(m2d(x2) / MATH_PI, 2), 0.5);
	y = sin(m2d(pom)) / m2d(pom);

	return y;
}

matrix ff3T(matrix x1, matrix x2, matrix ud1)
{
	matrix y;
	y = pow(m2d(x1), 2) + pow(m2d(x2), 2) - cos(2.5 * MATH_PI * m2d(x1)) - cos(2.5 * MATH_PI * m2d(x2)) + 2;
	return y;
}