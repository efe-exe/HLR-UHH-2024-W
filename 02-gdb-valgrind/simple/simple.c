/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

int *mistakes1(void) {
  int buf[] = {1, 1, 2, 3, 4, 5};
  return buf;
}

int *mistakes2(void) {
  int *buf = malloc(sizeof(char) * 4);
  buf[2] = 2;
  return buf;
}

int *mistakes3(void) {
  /* In dieser Funktion darf kein Speicher direkt d.h. explizit allokiert werden. */
  int mistakes2_ = 0;
  int *buf = (int *)&mistakes2;
  buf[0] = 3;
  return buf;
}

int *mistakes4(void) {
  int *buf = malloc(sizeof(char) * 4);
  buf[4] = 4;
  free(buf);
  return buf;
}

int *mistakes5(void) {
  int *buf = malloc(4 * 5);
  buf[44] = 5;
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
  free(p[1]); /* welcher Pointer war das doch gleich?, TODO: Fixme... ;-) */

  return 0;
}
