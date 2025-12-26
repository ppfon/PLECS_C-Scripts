/**
 * @file transformadas.h
 * @brief Transformadas de Coordenadas 
 *
 * Este arquivo implementa as transformadas de Clarke (ABC <-> AlphaBeta) e 
 * Park (AlphaBeta <-> DQ) mais ou menos como elas são realmente implementadas 
 * no DSP.
 *
 * @note As funções utilizam ponto flutuante de precisão simples (float) para
 * aproveitar a FPU nativa de 32-bits do DSP.
 * @note As funções trigonométricas (sin/cos) devem ser calculadas externamente
 * para garantir eficiência de ciclos de clock. Se for calculado dentro da função,
 * pode ser computado várias vezes ao mesmo tempo desnecessariamente. 
 * A transformada de Clarke, por exemplo, pode ser usada tanto na regulação de corrente usando 
 * um controlador proporcional-ressonante quanto em uma SRF-PLL. Para otimizar, a função trigonométrica
 * é calculada apenas uma vez a cada interrupção.
 *
 * @author Pedro Paulo
 * @date Dezembro 2025
 */

#ifndef TRANSFORMADAS_H
#define TRANSFORMADAS_H

#include <stdbool.h>
#include <stdint.h>

/* Constantes matemáticas */ 

/** @brief Raiz quadrada de 3. */
#define CONST_RAIZ_3      1.7320508076f 

/** @brief Raiz quadrada de 2. */
#define CONST_RAIZ_2      1.4142135624f 

/** @brief Inverso da raiz quadrada de 2 (1/sqrt(2)). Usado para evitar divisões. */
#define CONST_INV_RAIZ_2  0.7071067812f 

/** * @brief Constante de Clarke para invariância em amplitude.
 * Garante que a amplitude do vetor alpha/beta seja igual à amplitude no referencial natural.
 */
#define CONST_INV_AMP     0.66666666666f 

/** * @brief Constante de Clarke para invariância em potência.
 * Garante que a potência no domínio estacionário seja igual à do domínio natural.
 */
#define CONST_INV_POT     0.81649658092f 

/* Estruturas de Dados */

/** * @brief Estrutura para o referencial natural trifásico.
 * Representa as grandezas instantâneas das fases A, B e C.
 */
typedef struct {
    float a; /**< Valor instantâneo da Fase A */
    float b; /**< Valor instantâneo da Fase B */
    float c; /**< Valor instantâneo da Fase C */
} Abc_t;

/** * @brief Estrutura para o referencial estacionário (Clarke).
 * Representa o vetor espacial em um plano ortogonal fixo.
 */
typedef struct {
    float alpha; /**< Componente no eixo real (alinhado com a fase A) */
    float beta;  /**< Componente no eixo imaginário (90 graus da fase A) */
    float gamma; /**< Componente de sequência zero (homopolar) */
} AlphaBetaGamma_t;

/** * @brief Estrutura para o referencial síncrono (Park).
 * Representa o vetor espacial em um referencial girante (dq0).
 */
typedef struct {
    float d;    /**< Componente direta (alinhada com o fluxo/tensão) */
    float q;    /**< Componente em quadratura (torque/potência reativa) */
    float zero; /**< Componente de sequência zero (inalterada pela rotação) */
} DqZero_t;


/* Funções de Transformadas */

/**
 * @brief Executa a Transformada de Clarke (ABC -> AlphaBetaGamma).
 *
 * Converte componentes trifásicas naturais para componentes ortogonais estacionárias.
 *
 * Matriz de Transformação (Geral):
 * \f[
 * \begin{bmatrix} \alpha \\ \beta \\ \gamma \end{bmatrix} = K \cdot 
 * \begin{bmatrix} 
 * 1 & -1/2 & -1/2 \\ 
 * 0 & \sqrt{3}/2 & -\sqrt{3}/2 \\ 
 * 1/2 & 1/2 & 1/2 
 * \end{bmatrix} 
 * \cdot \begin{bmatrix} a \\ b \\ c \end{bmatrix}
 * \f]
 *
 * Onde K depende da invariância escolhida.
 *
 * @param natural Ponteiro para a estrutura de entrada (ABC).
 * @param estacionario Ponteiro para a estrutura de saída (AlphaBeta).
 * @param invariancia Flag de controle:
 * - 0: Invariância em amplitude (K = 2/3). Vetor resultante tem mesma amplitude da fase.
 * - 1: Invariância em potência (K = sqrt(2/3)). Transformação ortonormal.
 */
static inline void nat_para_estacionario(Abc_t *natural, AlphaBetaGamma_t *estacionario, bool invariancia)
{
    // Soma comum usada no cálculo do Gamma
    float soma_abc = natural->a + natural->b + natural->c;

    // Transformada invariante em amplitude
    // Objetivo do Gamma: Resultar na média aritmética = 1/3 * (a+b+c)
    if (invariancia == 0) 
    {
        estacionario->alpha = CONST_INV_AMP * (natural->a - 0.5f * (natural->b + natural->c)); 
        estacionario->beta  = CONST_INV_AMP * (CONST_RAIZ_3 * 0.5f * (natural->b - natural->c)); 
        
        // K da matriz é 2/3. A linha da sequência zero na matriz original é [0.5 0.5 0.5].
        // Portanto: (2/3) * 0.5 * (a+b+c) = 1/3 * (a+b+c).
        estacionario->gamma = CONST_INV_AMP * 0.5f * soma_abc;
    }
    // Transformada invariante em potência 
    // Objetivo do Gamma: Resultar em 1/sqrt(3) * (a+b+c) para manter a norma do vetor.
    else 
    {
        estacionario->alpha = CONST_INV_POT * (natural->a - 0.5f * (natural->b + natural->c)); 
        estacionario->beta  = CONST_INV_POT * (CONST_RAIZ_3 * 0.5f * (natural->b - natural->c));
        
        // K da matriz é sqrt(2/3). A linha da sequência zero na matriz ortonormal é [1/sqrt(2) ...].
        // Portanto: sqrt(2/3) * 1/sqrt(2) = sqrt(1/3) = 1/sqrt(3).
        estacionario->gamma = CONST_INV_POT * CONST_INV_RAIZ_2 * soma_abc;
    }
}

/**
 * @brief Executa a Transformada inversa de Clarke (AlphaBetaGamma -> ABC).
 *
 * Converte componentes estacionárias de volta para o sistema trifásico natural.
 *
 * @warning A matriz inversa muda drasticamente dependendo da invariância escolhida.
 * Na invariância de amplitude, a inversa não utiliza o ganho 2/3.
 * Na invariância de potência, a inversa é a transposta da direta.
 *
 * @param estacionario Ponteiro para a estrutura de entrada (AlphaBetaGamma_t).
 * @param natural Ponteiro para a estrutura de saída (Abc_t).
 * @param invariancia Flag de controle (deve ser a mesma usada na transformada direta): 0 para amplitude e 1 para potência.
 */
static inline void estacionario_para_nat(AlphaBetaGamma_t *estacionario, Abc_t *natural, bool invariancia)
{ 
    // Invariante em Amplitude
    // Matriz Inversa (Ganho Unitário para projeção Alpha):
    // A = alpha + gamma
    // B = -0.5*alpha + (sqrt(3)/2)*beta + gamma
    // C = -0.5*alpha - (sqrt(3)/2)*beta + gamma
    if (invariancia == 0) 
    {
        natural->a = estacionario->alpha + estacionario->gamma;
        
        natural->b = -0.5f * estacionario->alpha + (CONST_RAIZ_3 * 0.5f) * estacionario->beta + estacionario->gamma;
        
        natural->c = -0.5f * estacionario->alpha - (CONST_RAIZ_3 * 0.5f) * estacionario->beta + estacionario->gamma;
    }       
    // Invariante em Potência (Transposta da Direta * K)
    else 
    {
        // O termo gamma precisa ser escalado por 1/sqrt(2) dentro da matriz transposta
        // para respeitar a ortogonalidade.
        float gamma_norm = estacionario->gamma * CONST_INV_RAIZ_2;

        natural->a = CONST_INV_POT * (estacionario->alpha + gamma_norm);
        
        natural->b = CONST_INV_POT * (-0.5f * estacionario->alpha + (CONST_RAIZ_3 * 0.5f) * estacionario->beta + gamma_norm);
        
        natural->c = CONST_INV_POT * (-0.5f * estacionario->alpha - (CONST_RAIZ_3 * 0.5f) * estacionario->beta + gamma_norm);
    }
}

/**
 * @brief Executa a Transformada de Park (AlphaBeta -> DQ).
 *
 * Rotaciona o referencial estacionário por um ângulo theta para alinhá-lo
 * com o referencial síncrono.
 *
 * Equações:
 * \f$ d = \alpha \cos(\theta) + \beta \sin(\theta) \f$
 * \f$ q = -\alpha \sin(\theta) + \beta \cos(\theta) \f$
 *
 * @param estacionario Ponteiro para entrada (AlphaBetaGamma_t).
 * @param sincrono Ponteiro para saída (DqZero_t).
 * @param sin_theta Seno do ângulo de rotação (calculado externamente/PLL).
 * @param cos_theta Cosseno do ângulo de rotação (calculado externamente/PLL).
 */
static inline void estacionario_para_sincrono(AlphaBetaGamma_t *estacionario, DqZero_t *sincrono, float sin_theta, float cos_theta)
{
    // Eixo direto (projeção no vetor girante)
    sincrono->d = (estacionario->alpha * cos_theta) + (estacionario->beta * sin_theta);
    
    // Eixo quadratura (projeção ortogonal ao vetor girante - adiantado em 90 graus)
    sincrono->q = (-estacionario->alpha * sin_theta) + (estacionario->beta * cos_theta);
    
    // Eixo Zero (Homopolar) é o eixo Z, ortogonal ao plano de rotação. Não sofre alteração.
    sincrono->zero = estacionario->gamma;
}

/**
 * @brief Executa a Transformada inversa de Park (DQ -> AlphaBeta).
 *
 * Rotaciona o referencial síncrono de volta para o estacionário (ângulo -theta).
 *
 * Equações:
 * \f$ \alpha = d \cos(\theta) - q \sin(\theta) \f$
 * \f$ \beta = d \sin(\theta) + q \cos(\theta) \f$
 *
 * @param sincrono Ponteiro para entrada (DqZero_t).
 * @param estacionario Ponteiro para saída (AlphaBetaGamma_t).
 * @param sin_theta Seno do ângulo de rotação original.
 * @param cos_theta Cosseno do ângulo de rotação original.
 */
static inline void sincrono_para_estacionario(DqZero_t *sincrono, AlphaBetaGamma_t *estacionario, float sin_theta, float cos_theta)
{
    estacionario->alpha = (sincrono->d * cos_theta) - (sincrono->q * sin_theta);
    
    estacionario->beta  = (sincrono->d * sin_theta) + (sincrono->q * cos_theta);
    
    // O componente zero passa direto 
    estacionario->gamma = sincrono->zero;
}

#endif // TRANSFORMADAS_H