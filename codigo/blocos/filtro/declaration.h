#include "filtro.h"

#include <math.h>

Filtro_t filtro;

unsigned int tipo;
unsigned int ordem;
float freqCorte;
float freqAmostragem;
float zeta;

unsigned int estaInicializado;
