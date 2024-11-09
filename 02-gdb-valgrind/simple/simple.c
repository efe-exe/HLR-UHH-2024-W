/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> //für memcpy

int *mistakes1(void) {
  int buf[] = {1, 1, 2, 3, 4, 5};
  int* out = malloc(sizeof(buf));                     //reserviert speicher auf heap und gibt pointer zu diesem speicher, sizeof(buf) größe des arrays in bytes, 24bytes
  memcpy(out, buf, sizeof(buf));                // (destination, source, length) ist die signatur, kopiert lengthbyte von source zu destination. beides pointer, array ist ein pointer
  return out;                                        //heap nicht an den stack gebunden, daher "returbar"
}

int *mistakes2(void) {
  int *buf = malloc(sizeof(char) * 4);
  buf[1] = 2;                                        // andere zahlen für buf geben 0
  return buf;
}

int *mistakes3(void) {
  /* In dieser Funktion darf kein Speicher direkt d.h. explizit allokiert werden. */ //pointer auf mistakes2 gecastet, macht keinen sinn für Speicher aufrufen weil statisch
  int *buf = mistakes2();
  buf[0] = 3;
  return buf;
}

int *mistakes4(void) {
  int *buf = malloc(sizeof(char) * 4);
  buf[0] = 4;                                  // free buf musste weg // andere zahlen für buf geben 0 aus
  return buf;
}

int *mistakes5(void) {
  int *buf = malloc(4 * 5);
  buf[4] = 5;                      //tippfehler 44
  return buf;
}

int main(void) {
  /* Diese Zeile darf NICHT verändert werden! */
  int *p[5] = {&mistakes1()[1], &mistakes2()[1], mistakes3(), mistakes4(), mistakes5()+4};

  printf("1: %d\n", *p[0]);
  printf("2: %d\n", *p[1]);
  printf("3: %d\n", *p[2]);
  printf("4: %d\n", *p[3]);
  printf("5: %d\n", *p[4]);

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  free(p[0]-1);
  free(p[1]-1);
  free(p[2]);
  free(p[3]);
  free(p[4]-4);                          // 1 1 und 4 aus zeile 44 zurückrechnen


  return 0;
}
