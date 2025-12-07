/*********************************************
Kod stanowi uzupenienie materiałow do cwiczen
w ramach przedmiotu metody optymalizacji.
Kod udostepniony na licencji CC BY-SA 3.0
Autor: dr inz. Lukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Gorniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include<ctime>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dokladnosc
	int Nmax = 10000;										// maksymalna liczba wywolan funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dokladne rozwiazanie optymalne
	solution opt;											// rozwiazanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywolanie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznikow

	//Wahadlo
	Nmax = 1000;											// dokladnosc
	epsilon = 1e-2;											// maksymalna liczba wywolan funkcji celu
	lb = 0, ub = 5;											// dolne oraz gorne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahadla
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywolanie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznikow

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki poczatkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment sily dzialajacy na wahadlo oraz czas dzialania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwiazujemy rownanie rozniczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumien do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumień
	Y[0].~matrix();											// usuwamy z pamieci rozwiazanie RR
	Y[1].~matrix();
}


void lab1()
{
	matrix ud1, ud2;   // macierze pomocnicze

	std::srand(std::time(0));
	//double x0 = 62.0;   // wyjsciowa dziedzina [-100,100]
	double a = -100.0;
	double b = 100.0;
	double d = 1.0;    // rozmiar kroku wyjsciowego
	//double alpha = 1.1; // wspolczynnik ekspansji
	double* alphat = new double[3]{1.1, 1.5, 2.0};
	int Nmax = 1000;   // maksymalna liczba wywolan espansji
	double epsilon = 1e-6; // dokladnosc dla metod lokalnych
	double gamma = 0.618; // wspolczynnik zlagrania dla metody Lagrange'a

	ofstream Sout("ekspansja_lab1.csv");// definiujemy strumien do pliku .csv

	for (int e = 0; e < 3; e++ )
		for (int i = 0; i < 100; i++)
		{
			double x0 = (double)( (double)(rand() % 20000)/100 - 100); //losowa liczba pomiedzy -100 a 100

			double* interval = expansion(ff1T, x0, d, alphat[e], Nmax, ud1, ud2);
			Sout << x0 << ";"<< interval[0] << ";" << interval[1] << ";"<< interval[2] << ";";

			solution sol1 = fib(ff1T, interval[0], interval[1], epsilon, ud1, ud2);
			Sout << m2d(sol1.x) << ";" << m2d(sol1.y) << ";" << sol1.f_calls << ";" << sol1.flag << ";";
			solution sol2 = lag(ff1T, interval[0], interval[1], epsilon, gamma, Nmax, ud1, ud2);
			Sout << m2d(sol2.x) << ";" << m2d(sol2.y) << ";" << sol2.f_calls << ";" << sol2.flag << "\n";
		}
	Sout.close();
	solution sol1 = fib(ff1T, -100, 100, epsilon, ud1, ud2);
	cout << m2d(sol1.x) << ";" << m2d(sol1.y) << ";" << sol1.f_calls << ";" << sol1.flag << ";";
	solution sol2 = lag(ff1T, -100, 100, epsilon, gamma, Nmax, ud1, ud2);
	cout << m2d(sol2.x) << ";" << m2d(sol2.y) << ";" << sol2.f_calls << ";" << sol2.flag << "\n";

	// ========== Problem rzeczywisty - zbiorniki wody ==========
	cout << "\n=== Problem rzeczywisty - optymalizacja DA ===\n" << endl;
	solution::clear_calls();

	// Parametry optymalizacji
	double DA_min = 1.0;      // cm2 - dolne ograniczenie
	double DA_max = 100.0;    // cm2 - gorne ograniczenie
	double DA_start = 50.0;   // cm2 - punkt startowy (dla testu)
	double d_exp = 1.0;       // cm2 - rozmiar kroku ekspansji
	double alpha_exp = 1.5;   // wspolczynnik ekspansji
	int Nmax_exp = 1000;      // maksymalna liczba wywolan ekspansji
	double epsilon_opt = 1e-3; // dokladnosc optymalizacji
	double gamma_lag = 0.618; // wspolczynnik dla metody Lagrange'a
	int Nmax_opt = 1000;      // maksymalna liczba wywolan dla Lagrange'a

	// Test poprawnosci implementacji dla DA = 50 cm2
	cout << "Test poprawnosci dla DA = 50 cm2:" << endl;
	matrix DA_test(1, 1);
	DA_test(0) = 50.0;
	matrix y_test = ff_tanks(DA_test, ud1, ud2);

	// Symulacja dla DA = 50 cm2
	matrix Y0_test(4, 1);
	Y0_test(0) = 5.0;   // VA0
	Y0_test(1) = 1.0;   // VB0
	Y0_test(2) = 95.0;  // TA0
	Y0_test(3) = 20.0;  // TB0
	matrix params_test(2, 1);
	params_test(0) = 50.0;      // DA w cm2
	params_test(1) = 36.5665;   // DB w cm2
	matrix* Y_test = solve_ode(df_tanks, 0.0, 1.0, 2000.0, Y0_test, ud1, params_test);
	// Po transpozycji w solve_ode: Y[1] ma format (kroki × zmienne)
	// Zmienne: 0=VA, 1=VB, 2=TA, 3=TB
	int n_test = get_len(Y_test[0]);
	double T_max_test = Y_test[1](0, 3);  // krok 0, zmienna TB=3
	for (int i = 1; i < n_test; ++i)
	{
		if (T_max_test < Y_test[1](i, 3))  // krok i, zmienna TB=3
			T_max_test = Y_test[1](i, 3);
	}
	cout << "  Maksymalna temperatura TB = " << T_max_test << " degC (oczekiwana ~62.5degC)" << endl;
	Y_test[0].~matrix();
	Y_test[1].~matrix();

	// Optymalizacja - ekspansja przedzialu
	cout << "\nEkspansja przedzialu poszukiwan:" << endl;
	double* interval_tanks = expansion(ff_tanks, DA_start, d_exp, alpha_exp, Nmax_exp, ud1, ud2);
	cout << "  Przedzial znaleziony: [" << interval_tanks[0] << ", " << interval_tanks[1] << "] cm2" << endl;
	cout << "  Liczba wywolan funkcji celu: " << (int)interval_tanks[2] << endl;

	// Optymalizacja metoda Fibonacciego
	cout << "\nOptymalizacja metoda Fibonacciego:" << endl;
	solution::clear_calls();
	solution sol_fib = fib(ff_tanks, interval_tanks[0], interval_tanks[1], epsilon_opt, ud1, ud2);
	cout << "  DA_opt = " << m2d(sol_fib.x) << " cm2" << endl;
	cout << "  f(DA_opt) = " << m2d(sol_fib.y) << " (roznica od 50degC)" << endl;
	cout << "  Liczba wywolan funkcji celu: " << sol_fib.f_calls << endl;
	cout << "  Flaga: " << sol_fib.flag << endl;

	// Optymalizacja metoda Lagrange'a
	cout << "\nOptymalizacja metoda Lagrange'a:" << endl;
	solution::clear_calls();
	solution sol_lag = lag(ff_tanks, interval_tanks[0], interval_tanks[1], epsilon_opt, gamma_lag, Nmax_opt, ud1, ud2);
	cout << "  DA_opt = " << m2d(sol_lag.x) << " cm2" << endl;
	cout << "  f(DA_opt) = " << m2d(sol_lag.y) << " (roznica od 50degC)" << endl;
	cout << "  Liczba wywolan funkcji celu: " << sol_lag.f_calls << endl;
	cout << "  Flaga: " << sol_lag.flag << endl;

	// Zapis wynikow optymalizacji do pliku
	ofstream Sout_tanks("optymalizacja_zbiorniki.csv");
	Sout_tanks << "Metoda;DA_opt[cm2];f(DA_opt);f_calls;flag\n";
	Sout_tanks << "Fibonacci;" << m2d(sol_fib.x) << ";" << m2d(sol_fib.y) << ";" << sol_fib.f_calls << ";" << sol_fib.flag << "\n";
	Sout_tanks << "Lagrange;" << m2d(sol_lag.x) << ";" << m2d(sol_lag.y) << ";" << sol_lag.f_calls << ";" << sol_lag.flag << "\n";
	Sout_tanks.close();

	// Symulacja dla znalezionego optymalnego DA (uzywamy wyniku z Fibonacciego)
	double DA_opt = m2d(sol_fib.x);
	cout << "\nSymulacja dla optymalnego DA = " << DA_opt << " cm2:" << endl;

	matrix Y0_sim(4, 1);
	Y0_sim(0) = 5.0;   // VA0
	Y0_sim(1) = 1.0;   // VB0
	Y0_sim(2) = 95.0;  // TA0
	Y0_sim(3) = 20.0;  // TB0
	matrix params_sim(2, 1);
	params_sim(0) = DA_opt;    // DA w cm2
	params_sim(1) = 36.5665;   // DB w cm2
	matrix* Y_sim = solve_ode(df_tanks, 0.0, 1.0, 2000.0, Y0_sim, ud1, params_sim);

	// Zapis symulacji do pliku CSV
	// Format: t, VA, VB, TA, TB
	// Po transpozycji w solve_ode: Y[1] ma format (kroki × zmienne)
	// Zmienne: 0=VA, 1=VB, 2=TA, 3=TB
	ofstream Sout_sim("symulacja_zbiorniki.csv");
	Sout_sim << "t[s];VA[m3];VB[m3];TA[degC];TB[degC]\n";
	int n_sim = get_len(Y_sim[0]);
	for (int i = 0; i < n_sim; ++i)
	{
		Sout_sim << Y_sim[0](i) << ";" 
			<< Y_sim[1](i, 0) << ";"   // krok i, zmienna VA=0
			<< Y_sim[1](i, 1) << ";"   // krok i, zmienna VB=1
			<< Y_sim[1](i, 2) << ";"   // krok i, zmienna TA=2
			<< Y_sim[1](i, 3) << "\n"; // krok i, zmienna TB=3
	}
	Sout_sim.close();

	// Znajdz maksymalna temperature
	double T_max_sim = Y_sim[1](0, 3);  // krok 0, zmienna TB=3
	for (int i = 1; i < n_sim; ++i)
	{
		if (T_max_sim < Y_sim[1](i, 3))  // krok i, zmienna TB=3
			T_max_sim = Y_sim[1](i, 3);
	}
	cout << "  Maksymalna temperatura TB = " << T_max_sim << " degC" << endl;
	cout << "  Wyniki zapisane do: symulacja_zbiorniki.csv" << endl;

	Y_sim[0].~matrix();
	Y_sim[1].~matrix();
	delete[] interval_tanks;
	delete[] alphat;
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
