// Inclui o arquivo de cabeçalho relativo ao local do modelo
#include "srf_pll.h"

// Instância Global da PLL
// 'static' garante que ela seja visível apenas neste bloco
static SRF_PLL_t pll;

// Flag de inicialização
static int estaInicializado = 0;