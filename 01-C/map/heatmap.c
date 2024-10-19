#include <stdio.h>

// Definieren Sie ein 3x3-Array Namens map, das Werte vom Typ double enthält
static double map[3][3] = {0}; // Initialisierung aller Werte mit 0


// Die Funktion set_temperature soll an Position [x, y] den Wert dir in das Array map eintragen
// Überprüfen Sie x und y, um mögliche Arrayüberläufe zu verhindern
void set_temperature (int x, int y, double temperature)
{
	if (x >= 0 && x < 3 && y >= 0 && y < 3)
	{
		map[x][y] = temperature;
	}
	else
	{
		printf("Error: Arrayueberlauf bei x=%d, y=%d. t=%f wird nicht eingefügt\n", x, y, temperature); // %d für integer und %f für doubles
	}
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
	for (int i = 2; i >= 0; i--) // Subtraktiv, weil Y-Ache nach oben positiver wird
	{
		printf(" ");
		for (int j = 0; j < 3; j++) // Additiv, weil X-Ache nach rechts positiver wird
		{
			printf(" %6.2f ", map[j][i]); // 6 Stellen insg. und 2 Nachkommastellen
		}
		printf("\n");
	}
	printf("\n");
}

// Die Funktion average_value soll an Position [x, y] den Durchschnitt der 8 umgebenen
// Temperaturen in das Array map eintragen. 
// Für Werte außerhalb des Arrays soll der Wert 0 angenommen werden.
// Verwenden Sie hierfür auch die Funktion set_temperature.
void set_average (int x, int y)
{
	double sum = 0;
	int count = 0;
	
	for (int i = x - 1; i <= x + 1; i++)
	{
		for (int j = y - 1; j <= y + 1; j++)
		{
			if (i >= 0 && i < 3 && j >= 0 && j < 3)
			{
				sum += map[i][j];
				count++;
			}
		}
	}
	
	set_temperature(x, y, sum / count);
}

// In dieser Funktion darf nichts verändert werden!
int main (void)
{
	set_temperature(0, 1, 40);
	set_temperature(1, 0, 160);
	set_temperature(1, 4, 75);
	set_temperature(1, 2, 80);
	set_temperature(2, 1, 120);

	show_map();

	set_temperature(0, 0, 20.5);
	set_temperature(0, 2, 14.8);
	set_temperature(0, 2, 22.7);
	set_temperature(2, 0, 100.2);
	set_temperature(2, 2, 20.6);
	set_temperature(2, 2, 200.2);
	set_temperature(1, 3, 200.06);
	set_temperature(1, 1, 50.5);

	show_map();
  
  set_average(0,0);
  set_average(2,0);
  set_average(1,2);
  
  show_map();

	return 0;
}
