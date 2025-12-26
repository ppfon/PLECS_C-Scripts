/**
 * @file filtro.h
 * @brief Implementação de filtro IIR universal (biquad)
 * Suporta LPF, HPF, BPF e Notch de 1ª e 2ª ordens.
 * @author Pedro Paulo Fontolan de Faria
 */

#ifndef FILTRO_H
#define FILTRO_H

#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

// Definição dos tipos de filtro
typedef enum {
    FILTRO_PB,   // passa-baixa
    FILTRO_PA,   // passa-alta
    FILTRO_PF,   // passa-faixa - apenas 2ª ordem
    FILTRO_RF  // rejeita-faixa - apenas 2ª ordem
} FiltroTipo_t;

typedef struct {
    float b0, b1, b2; // Coeficientes das entradas
    float a1, a2;     // Coeficientes das saídas
} Coeficientes_t;

typedef struct {
    float x1, x2; // Entradas passadas
    float y1, y2; // Saídas passadas
} Estados_t;

typedef struct {
    Coeficientes_t coef;
    Estados_t est;
} Filtro_t;

/**
 * @brief Calcula coeficientes para qualquer filtro IIR até 2ª ordem.
 * @param filtro        Ponteiro para o filtro
 * @param tipo      FILTRO_PB, FILTRO_PA, FILTRO_PF, FILTRO_RF
 * @param ordem     1 ou 2 (BPF e Notch ignoram isso e forçam 2)
 * @param freqCorte        Frequência de corte ou central (Hz)
 * @param freqAmostragem        Frequência de amostragem (Hz)
 * @param zeta      Fator de amortecimento. 
 * Use 0.7071 para Butterworth LPF/HPF.
 * Para BPF/Notch, quanto menor o zeta, mais estreita a banda (maior o Q).
 */
static inline void Filtro_inicia(Filtro_t *filtro, FiltroTipo_t tipo, unsigned int ordem, float freqCorte, float freqAmostragem, float zeta)
{
    // Limpa estados e coeficientes
    filtro->est.x1 = 0; filtro->est.x2 = 0;
    filtro->est.y1 = 0; filtro->est.y2 = 0;
    filtro->coef.b0 = 0; filtro->coef.b1 = 0; filtro->coef.b2 = 0;
    filtro->coef.a1 = 0; filtro->coef.a2 = 0;

    // Pre-warping comum a todos (Tustin)
    float K = tanf(M_PI * freqCorte / freqAmostragem);
    float K2 = K * K;
    float div = 1.0f; // Divisor de normalização

    if (ordem == 1 && (tipo == FILTRO_PB || tipo == FILTRO_PA)) {
        // --- 1ª ORDEM ---
        div = K + 1.0f;
        
        if (tipo == FILTRO_PB) {
            filtro->coef.b0 = K / div;
            filtro->coef.b1 = filtro->coef.b0;
        } else { // HPF
            filtro->coef.b0 = 1.0f / div;
            filtro->coef.b1 = -filtro->coef.b0;
        }
        filtro->coef.a1 = (K - 1.0f) / div;
        filtro->coef.a2 = 0.0f;     // Não usa termo z^-2
        filtro->coef.b2 = 0.0f;

    } else {
        // --- 2ª ORDEM (genérica) ---
        // Denominador comum: s^2 + 2*zeta*w*s + w^2 -> mapeado em Z
        div = K2 + (2.0f * zeta * K) + 1.0f;

        // Cálculo dos pólos (a1, a2) - igual para todos dessa família
        filtro->coef.a1 = (2.0f * (K2 - 1.0f)) / div;
        filtro->coef.a2 = (K2 - (2.0f * zeta * K) + 1.0f) / div;

        switch (tipo) {
            case FILTRO_PB:
                // H(s) = w^2 / D(s)
                filtro->coef.b0 = K2 / div;
                filtro->coef.b1 = 2.0f * filtro->coef.b0;
                filtro->coef.b2 = filtro->coef.b0;
                break;

            case FILTRO_PA:
                // H(s) = s^2 / D(s)
                filtro->coef.b0 = 1.0f / div;
                filtro->coef.b1 = -2.0f * filtro->coef.b0;
                filtro->coef.b2 = filtro->coef.b0;
                break;

            case FILTRO_PF:
                // H(s) = 2*zeta*w*s / D(s)  (ganho unitário na ressonância)
                filtro->coef.b0 = (2.0f * zeta * K) / div;
                filtro->coef.b1 = 0.0f;
                filtro->coef.b2 = -filtro->coef.b0;
                break;

            case FILTRO_RF:
                // H(s) = (s^2 + w^2) / D(s)
                filtro->coef.b0 = (K2 + 1.0f) / div;
                filtro->coef.b1 = (2.0f * (K2 - 1.0f)) / div; 
                filtro->coef.b2 = filtro->coef.b0;
                break;
        }
    }
}

/**
 * @brief Executa o Biquad (forma direta I transposta)
 */
static inline float Filtro_calcula(Filtro_t *filtro, float entrada)
{
    float saida = (filtro->coef.b0 * entrada) + 
              (filtro->coef.b1 * filtro->est.x1) + 
              (filtro->coef.b2 * filtro->est.x2) - 
              (filtro->coef.a1 * filtro->est.y1) - 
              (filtro->coef.a2 * filtro->est.y2);

    // Armazena o estado atual 
    filtro->est.x2 = filtro->est.x1;
    filtro->est.x1 = entrada;
    filtro->est.y2 = filtro->est.y1;
    filtro->est.y1 = saida;

    return saida;
}

#endif