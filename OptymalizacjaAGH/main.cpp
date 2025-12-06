/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include"ctime"

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
	double epsilon = 1e-2;									// dok³adnoœæ
	int Nmax = 10000;										// maksymalna liczba wywo³añ funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz górne ograniczenie
		a(2, 1);											// dok³adne rozwi¹zanie optymalne
	solution opt;											// rozwi¹zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Wahadlo
	Nmax = 1000;											// dok³adnoœæ
	epsilon = 1e-2;											// maksymalna liczba wywo³añ funkcji celu
	lb = 0, ub = 5;											// dolne oraz górne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad³a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumieñ do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumieñ
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
}


void lab1()
{
	matrix ud1, ud2;   // macierze pomocnicze

	std::srand(std::time(0));
	//double x0 = 62.0;   // wyjsciowa dziedzina [-100,100]
	
	double d = 1.0;    // rozmiar kroku wyjsciowego
	//double alpha = 1.1; // wspolczynnik ekspansji
	double* alphat = new double[3]{1.1, 1.5, 2.0};
	int Nmax = 1000;   // maksymalna liczba wywolan espansji

	ofstream Sout("ekspansja_lab1.csv");// definiujemy strumieñ do pliku .csv

	for (int e = 0; e < 3; e++ )
	for (int i = 0; i < 100; i++)
	{
		double x0 = (double)( (double)(rand() % 20000)/100 - 100); //losowa liczba pomiedzy -100 a 100

		double* interval = expansion(ff1T, x0, d, alphat[e], Nmax, ud1, ud2);
		Sout << x0 << ";"<< interval[0] << ";" << interval[1] << ";"<< interval[2] << ";";

		solution sol1 = fib(ff1T, interval[0], interval[1], 1e-6, ud1, ud2);
		Sout << sol1.x << ";" << sol1.y << ";" << sol1.f_calls << ";";
		//double* interval = lagrange(ff1T, x0, d, alphat[e], Nmax, ud1, ud2);
		//Sout << x0 << ";" << interval[0] << ";" << interval[1] << ";" << interval[2] << "\n";
	}
	Sout.close();
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
