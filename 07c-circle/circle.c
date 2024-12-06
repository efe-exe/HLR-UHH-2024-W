#include <stdio.h> // Standard I/O Funktionen
#include <stdlib.h> // Standard Bibliothek für Speicherverwaltung, Zufallszahlen, etc.
#include <mpi.h> // MPI Bibliothek für parallele Programmierung
#include <time.h> // Zeitfunktionen für Zufallszahlengenerierung
#include <string.h> // Für memset und memcpy

// Funktion zur Initialisierung des lokalen Arrays
int* init_local_array(int local_size, int rank) {
    srand(time(NULL) + rank); // Unterschiedlicher Seed pro Prozess
    int* local_buf = (int*)malloc(sizeof(int) * local_size); // Speicher für lokales Array allokieren
    if (!local_buf) { // Überprüfen, ob die Speicherallokation erfolgreich war
        perror("Memory allocation failed"); // Fehlermeldung ausgeben
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // MPI-Programm abbrechen
    }

    for (int i = 0; i < local_size; i++) { // Array mit Zufallszahlen füllen
        local_buf[i] = rand() % 12; // Zufallswerte zwischen 0 und 11
    }

    return local_buf; // Lokales Array zurückgeben
}

// Funktion zum synchronisierten Drucken der Arrays
void printSynced(int rank, int nprocs, int *numbers, int max_local_size) {
    MPI_Status status; // MPI Statusvariable

    if (rank == 0) { // Wenn der Prozess der erste ist
        for (int i = 0; i < max_local_size; i++) { // Array drucken
            int num = numbers[i];
            if (num >= 0) {
                printf("%02d ", num); // Zahl drucken
            } else {
                printf("   "); // Leerraum drucken
            }
        }
        printf(" | "); // Trennzeichen drucken
        int l = 1;
        MPI_Ssend(&l, 1, MPI_INT, 1, 0, MPI_COMM_WORLD); // Nachricht an den nächsten Prozess senden
        return;
    }

    int l = 0;
    MPI_Recv(&l, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status); // Nachricht vom vorherigen Prozess empfangen

    for (int i = 0; i < max_local_size; i++) { // Array drucken
        int num = numbers[i];
        if (num >= 0) {
            printf("%02d ", num); // Zahl drucken
        } else {
            printf("   "); // Leerraum drucken
        }
    }

    if (rank != nprocs - 1) { // Wenn der Prozess nicht der letzte ist
        printf(" | "); // Trennzeichen drucken
        MPI_Ssend(&l, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD); // Nachricht an den nächsten Prozess senden
    } else {
        printf("\n"); // Neue Zeile drucken
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv); // MPI-Umgebung initialisieren

    int rank, nprocs, N;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rang des Prozesses ermitteln
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); // Anzahl der Prozesse ermitteln

    if (argc < 2) { // Überprüfen, ob die Array-Größe als Argument übergeben wurde
        if (rank == 0) { // Wenn der Prozess der erste ist
            fprintf(stderr, "Usage: %s <array_size>\n", argv[0]); // Fehlermeldung ausgeben
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // MPI-Programm abbrechen
    }

    N = atoi(argv[1]); // Array-Größe aus dem Argument lesen
    if (N < nprocs) { // Überprüfen, ob die Anzahl der Prozesse die Array-Größe überschreitet
        if (rank == 0) { // Wenn der Prozess der erste ist
            fprintf(stderr, "Error: Number of processes cannot exceed array size\n"); // Fehlermeldung ausgeben
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // MPI-Programm abbrechen
    }

    // Lokale Array-Größe berechnen
    int local_size = N / nprocs + (rank < N % nprocs ? 1 : 0);
    int* local_buf = init_local_array(local_size, rank); // Lokales Array initialisieren

    // Maximale lokale Größe ermitteln
    int max_local_size;
    MPI_Allreduce(&local_size, &max_local_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); // Maximale lokale Größe berechnen

    // Sendepuffer vorbereiten und bei Bedarf mit Nullen auffüllen
    int* send_buf = (int*)calloc(max_local_size, sizeof(int)); // Speicher für Sendepuffer allokieren
    memcpy(send_buf, local_buf, sizeof(int) * local_size); // Lokales Array in den Sendepuffer kopieren

    // Ringkommunikation einrichten
    int left = (rank == 0) ? nprocs - 1 : rank - 1; // Linker Nachbarprozess
    int right = (rank + 1) % nprocs; // Rechter Nachbarprozess
    int* recv_buf = (int*)malloc(sizeof(int) * max_local_size); // Speicher für Empfangspuffer allokieren
    int iterations = 0; // Iterationszähler
    int initial_first_element = 0; // Variable zum Speichern des initialen ersten Elements

    // Initiale Arrays sammeln und drucken
    int* all_arrays = NULL;
    if (rank == 0) { // Wenn der Prozess der erste ist
        all_arrays = (int*)malloc(sizeof(int) * max_local_size * nprocs); // Speicher für alle Arrays allokieren
    }

    MPI_Gather(send_buf, max_local_size, MPI_INT, all_arrays, max_local_size, MPI_INT, 0, MPI_COMM_WORLD); // Arrays sammeln

    if (rank == 0) { // Wenn der Prozess der erste ist
        printf("Iteration %d: ", iterations); // Iterationsnummer drucken
        for (int i = 0; i < nprocs * max_local_size; i++) { // Alle Arrays drucken
            printf("%02d ", all_arrays[i]); // Zahlen mit führenden Nullen drucken
            if ((i + 1) % max_local_size == 0 && i != nprocs * max_local_size - 1) {
                printf("| "); // Trennzeichen drucken
            }
        }
        printf("\n");
        // Initiales erstes Element speichern
        initial_first_element = all_arrays[0];
    }

    MPI_Barrier(MPI_COMM_WORLD); // Sicherstellen, dass alle Prozesse das initiale erste Element haben

    // Terminationsflag initialisieren
    int terminate_flag = 0;

    // Schleife mit Terminationsbedingung basierend auf dem initialen ersten Element
    do {
        // Arrays senden und empfangen wegen Deadlock vermeiden
        if (rank % 2 == 0) {
            MPI_Ssend(send_buf, max_local_size, MPI_INT, right, 0, MPI_COMM_WORLD);
            MPI_Recv(recv_buf, max_local_size, MPI_INT, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(recv_buf, max_local_size, MPI_INT, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Ssend(send_buf, max_local_size, MPI_INT, right, 0, MPI_COMM_WORLD);
        }

        // Sendepuffer für die nächste Iteration aktualisieren
        memcpy(send_buf, recv_buf, sizeof(int) * max_local_size);
        iterations++;

        // Arrays nach jeder Iteration sammeln und drucken
        MPI_Gather(send_buf, max_local_size, MPI_INT, all_arrays, max_local_size, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0) { // Wenn der Prozess der erste ist
            printf("Iteration %d: ", iterations); // Iterationsnummer drucken
            for (int i = 0; i < nprocs * max_local_size; i++) { // Alle Arrays drucken
                printf("%02d ", all_arrays[i]);
                if ((i + 1) % max_local_size == 0 && i != nprocs * max_local_size - 1) {
                    printf("| "); // Trennzeichen drucken
                }
            }
            printf("\n");

            // Terminationsbedingung nur nach der ersten Iteration überprüfen
            if (iterations > 0) {
                int last_rank_first = all_arrays[(nprocs - 1) * max_local_size]; // Erstes Element des letzten Prozesses
                if (last_rank_first == initial_first_element) { // Überprüfen, ob das erste Element des letzten Prozesses dem initialen ersten Element entspricht
                    printf("Termination condition met: last rank's first element matches initial first element.\n");
                    terminate_flag = 1; // Terminationsflag setzen
                }
            }
        }

        // Terminationsflag an alle Prozesse senden
        MPI_Bcast(&terminate_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (terminate_flag) { // Wenn das Terminationsflag gesetzt ist
            break; // Schleife beenden
        }

    } while (iterations < nprocs); // Schleife fortsetzen, bis die maximale Anzahl an Iterationen erreicht ist

    if (rank == 0) { // Wenn der Prozess der erste ist
        printf("Abbruch nach %d Iterationen\n", iterations); // Anzahl der Iterationen drucken
        printf("Abbruch wegen %d\n", initial_first_element); // Grund für den Abbruch drucken
    }

    // Allokierten Speicher freigeben
    free(local_buf);
    free(send_buf);
    free(recv_buf);
    if (rank == 0) { // Wenn der Prozess der erste ist
        free(all_arrays);
    }

    MPI_Finalize(); // MPI-Umgebung beenden
    return 0; // Programm beenden
}