/**
 * @file srf_pll.h
 * @brief Implementação de PLL em referencial síncrono (SRF-PLL).
 *
 * Estrutura modular que utiliza as bibliotecas de Transformadas e PID
 * para sincronização com a rede trifásica.
 *
 * @author Pedro Paulo
 */

#ifndef SRF_PLL_H
#define SRF_PLL_H

#include "/home/unknow/Desktop/PLECS/codigo/include/transformadas.h"
#include "/home/unknow/Desktop/PLECS/codigo/include/PID.h"
#include <math.h>

#define DOIS_PI 6.28318530718f

/**
 * @brief Objeto principal da SRF-PLL.
 */
typedef struct {
    // Sub-módulos
    PID_t filtroMalha;           /**< Controlador PI do filtro de loop */
    
    // Estados internos
    float theta;                /**< Ângulo estimado atual (rad) */
    float omegaFF;             /**< Frequência feedforward (rad/s) */
    float omegaEstimado;       /**< Saída do filtro de malha + feedforward */
    
    // Detector de fase
    DqZero_t vDQ;              /**< Tensões no referencial síncrono */
    float sinTheta;            /**< Seno do ângulo atual (otimização) */
    float cosTheta;            /**< Cosseno do ângulo atual (otimização) */

} SRF_PLL_t;

/**
 * @brief Inicializa a PLL.
 * @param pll Ponteiro para a estrutura.
 * @param kp Ganho proporcional do controlador.
 * @param ki Ganho integral do controlador.
 * @param freq_rede Frequência nominal da rede.
 * @param periodo Período de amostragem.
 */
static inline void SRF_PLL_inicia(SRF_PLL_t *pll, float kp, float ki, float freqRede, float periodo)
{
    // Inicializa o PID (Filtro de malha).
    // SRF-PLL geralmente não usa derivativo (Kd=0), nem filtro (N=1).
    // Saturação simétrica para evitar desvios grandes de frequência (ex: +/- 10Hz em rad/s)
    float maxFreqDev = DOIS_PI * 10.0f; 
    PID_inicia(&pll->filtroMalha, kp, ki, 0.0f, -maxFreqDev, maxFreqDev, periodo, PID_DISC_BILINEAR, 1.0f);

    pll->omegaFF = DOIS_PI * freqRede;
    pll->theta = 0.0f;
    pll->omegaEstimado = pll->omegaFF;
    pll->sinTheta = 0.0f;
    pll->cosTheta = 1.0f;
}

/**
 * @brief Executa um passo da PLL.
 * Deve ser chamado na frequência de amostragem.
 * * @note O código que chamar esta função deve realizar a Transformada de Park externamente 
 * utilizando o seno/cosseno fornecidos por esta PLL no passo anterior,
 * e passar o vetor resultante v_dq0 aqui.
 * * @param pll Ponteiro para a estrutura.
 * @param v_dq0 Tensão no referencial síncrono (já transformada externamente).
 */
static inline void SRF_PLL_executa(SRF_PLL_t *pll, DqZero_t vDQ)
{
    // 1. Atualiza o estado interno com a entrada externa
    // Assume-se que vDQ foi calculado usando pll->sinTheta e pll->cosTheta
    pll->vDQ = vDQ;

    // 2. Cálculo da magnitude para normalização
    // (Pode-se usar apenas v_dq.d se a PLL já estiver travada, mas sqrt é mais seguro na partida)
    float magnitude = sqrtf(pll->vDQ.d * pll->vDQ.d + pll->vDQ.q * pll->vDQ.q);
    
    float erroFase;

    // 3. Cálculo do Erro de Fase Normalizado
    // O objetivo da PLL é alinhar o eixo D com o vetor tensão (Vq -> 0).
    // Proteção contra divisão por zero (se a rede cair ou na inicialização)
    if (magnitude > 1.0f) { // Limiar de segurança (ex: 1V)
        // O erro agora é adimensional (radianos), desacoplando a amplitude da tensão da dinâmica da malha
        erroFase = pll->vDQ.q / magnitude; 
    } else {
        erroFase = 0.0f;
    }

    // 4. Filtro de malha (controlador PI)
    PID_calcula(&pll->filtroMalha, erroFase);

    // 5. VCO (Oscilador Controlado por Tensão)
    // Frequência = Nominal + correção do PI
    pll->omegaEstimado = pll->omegaFF + pll->filtroMalha.saida;

    // 6. Integrador (frequência -> ângulo)
    // Método: Euler Progressivo (Forward Euler) / Acumulador com atraso.
    // Equação: theta[k+1] = theta[k] + omegaEstimado[k] * periodo.
    //
    // Por que NÃO usar Transformada Bilinear (Tustin) ou Euler Regressivo aqui?
    // 1. Loop Algébrico: Métodos implícitos (como Tustin) requerem o ângulo atual 
    //    para calcular o erro de fase atual, mas o erro é necessário para calcular 
    //    o ângulo, criando um travamento matemático (dependência circular).
    // 2. Custo vs benefício: A entrada deste integrador (frequência) é praticamente 
    //    constante (sinal CC). Para sinais CC, a precisão do Euler é matematicamente 
    //    equivalente à do Tustin, mas com menor custo computacional.
    pll->theta += pll->omegaEstimado * pll->filtroMalha.periodo;

    // 7. Normalização do ângulo (Wrap around 2*pi)
    if (pll->theta > DOIS_PI) {
        pll->theta -= DOIS_PI;
    } else if (pll->theta < 0.0f) {
        pll->theta += DOIS_PI;
    }

    // 8. Pré-calcula sin/Ccs para o próximo ciclo (otimização)
    // Estes valores serão usados externamente para calcular o próximo vDQ
    pll->sinTheta = sinf(pll->theta);
    pll->cosTheta = cosf(pll->theta);
}

#endif // SRF_PLL_H