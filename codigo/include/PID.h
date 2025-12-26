/**
 * @file pid.h
 * @brief Controlador PID Digital Otimizado com Filtro Derivativo.
 *
 * Implementação de um controlador Proporcional-Integral-Derivativo (PID)
 *
 * Principais Características:
 * - Estrutura paralela: PID(s) = Kp + Ki/s + Kd*s/(tau*s+1).
 * - Termo derivativo com filtro passa-baixa de primeira ordem (essencial para evitar ruído).
 * - Pré-cálculo de coeficientes para evitar divisões no loop de controle.
 * - Suporte a múltiplos métodos de discretização (Tustin/Bilinear e Euler regressivo/progressivo).
 * - Anti-windup através de clamping condicional.
 *
 * @author Pedro Paulo
 * @date Dezembro 2025
 */

#ifndef PID_H
#define PID_H

#include <stdint.h>

/**
 * @brief Métodos de discretização (mapeamento s -> z).
 * Determina como as equações diferenciais contínuas são convertidas em equações de diferenças.
 */
typedef enum {
    /** * @brief Método de Euler Progressivo (Forward Euler).
     * Aproximação: s ~= (z-1)/T
     * @warning Pode ser instável se o período de amostragem for grande em relação à constante de tempo.
     */
    PID_DISC_EULER_PROGRESSIVO = -1,

    /** * @brief Transformada Bilinear (Tustin).
     * Aproximação: s ~= 2/T [ (z-1)/(z+1) ] 
     * @note Método padrão da indústria. Preserva a estabilidade dentro do círculo unitário.
     * Mapeia o eixo imaginário jw no círculo unitário sem distorção de estabilidade.
     */
    PID_DISC_BILINEAR = 0,

    /** * @brief Método de Euler Regressivo (Backward Euler).
     * Aproximação: s ~= (z-1)(zT) $
     * @note Método incondicionalmente estável para funções de transferência estáveis.
     * Tende a introduzir mais amortecimento (atraso de fase) que o Tustin.
     */
    PID_DISC_EULER_REGRESSIVO = 1
} PID_Metodo_e;

/**
 * @brief Estrutura de Ganhos do Controlador.
 */
typedef struct {
    float kp; /**< Ganho Proporcional */
    float ki; /**< Ganho Integral */
    float kd; /**< Ganho Derivativo */
    
    /** * @brief Constante de tempo do filtro derivativo (Tau).
     * Define o polo do filtro passa-baixa aplicado ao termo D.
     * D(s) = (kd * s)/(tau_{filtro} * s + 1)
     * @note Tipicamente tauFiltroDerivativo ~= (Td/10).
     */
    float tauFiltroDerivativo; 
} GanhosPID_t;

/**
 * @brief Objeto principal do Controlador PID.
 * Contém os estados, coeficientes pré-calculados e configurações.
 */
typedef struct {
    GanhosPID_t ganhos; /**< Cópia dos ganhos configurados */
    
    /* -------------------------------------------------------------------------
       Coeficientes Pré-calculados (Otimização)
       ------------------------------------------------------------------------- */
    
    /** * @brief Coeficiente 'b0' do termo Integral.
     * Multiplicador do erro atual erro[k].
     */
    float bi0;

    /** * @brief Coeficiente 'b1' do termo Integral.
     * Multiplicador do erro anterior erro[k-1].
     */
    float bi1;

    /** * @brief Coeficiente 'a1' do filtro Derivativo.
     * Multiplicador da saída derivativa anterior D[k-1]. Representa o polo do filtro.
     */
    float ad1;

    /** * @brief Coeficiente 'b0' do filtro Derivativo.
     * Multiplicador da diferença de erro (erro[k] - erro[k-1]). Representa o ganho do zero.
     */
    float bd0;

    /* -------------------------------------------------------------------------
       Memórias de Estado
       ------------------------------------------------------------------------- */
    float erro;             /**< Erro atual e[k] */
    float erroAnterior;     /**< Erro passo anterior e[k-1] */
    float saida;            /**< Saída total do PID u[k] */
    float integral;         /**< Estado do acumulador integral I[k] */
    float derivativo;       /**< Estado do termo derivativo filtrado D[k] */
    
    /* -------------------------------------------------------------------------
       Configurações
       ------------------------------------------------------------------------- */
    float satMin;           /**< Saturação mínima da saída */
    float satMax;           /**< Saturação máxima da saída */
    float periodo;          /**< Período de amostragem em segundos */
    PID_Metodo_e metodo;    /**< Método de discretização selecionado */

} PID_t;

/**
 * @brief Inicializa o PID, configura o filtro derivativo e pré-calcula coeficientes.
 * * Realiza o mapeamento dos ganhos contínuos (domínio s) para coeficientes discretos (domínio z)
 * baseados no método escolhido. Isso remove operações de divisão da função periódica.
 *
 * A relação do filtro derivativo segue a regra prática da indústria:
 * tauFiltroDerivativo = alpha T_d \f$, onde tipicamente \f$ \alpha = 0.1 \f$ (ou seja, N=10).
 * * @param pid Ponteiro para a estrutura do PID.
 * @param kp Ganho Proporcional.
 * @param ki Ganho Integral.
 * @param kd Ganho Derivativo.
 * @param satMin Limite inferior da saída.
 * @param satMax Limite superior da saída.
 * @param periodo Período de amostragem em segundos.
 * @param metodo Método de discretização (Euler ou Bilinear).
 * @param N_filtro Fator de divisão do filtro derivativo (Tf = Td / N). Recomendado: 10.0f.
 */
static inline void PID_inicia(PID_t *pid, float kp, float ki, float kd, 
                              float satMin, float satMax, float periodo, 
                              PID_Metodo_e metodo, float N_filtro)
{
    pid->ganhos.kp = kp;
    pid->ganhos.ki = ki;
    pid->ganhos.kd = kd;
    
    // Calcula constante de tempo do filtro derivativo: Tau = Kd / (Kp * N)
    // Se o controlador for puramente I ou D (Kp=0), assume Tau direto pelo Kd.
    if (kp > 1e-6f) {
        pid->ganhos.tauFiltroDerivativo = (kd / kp) / N_filtro;
    } else {
        pid->ganhos.tauFiltroDerivativo = kd / N_filtro;
    }

    pid->satMin = satMin;
    pid->satMax = satMax;
    pid->periodo = periodo;
    pid->metodo = metodo;

    // -----------------------------------------------------------------
    // Pré-cálculo dos coeficientes discretos (mapeamento para o domínio z)
    // -----------------------------------------------------------------
    switch (metodo)
    {
        case PID_DISC_EULER_PROGRESSIVO: 
            // Aproximação: s ~ (z-1)/Ts
            // Integral: I[k] = I[k-1] + Ki*Ts * e[k-1] 
            // (Euler Progressivo usa o erro anterior para integrar "retangular à esquerda")
            pid->bi0 = 0.0f; 
            pid->bi1 = ki * pid->periodo;
            
            // Derivativo Filtrado: D(s) = Kds / (Tfs + 1)
            // Equação de Diferenças: D[k] = (1 - Ts/Tf)*D[k-1] + (Kd/Tf)*(e[k] - e[k-1])
            if (pid->ganhos.tauFiltroDerivativo > 1e-9f) {
                pid->ad1 = 1.0f - (pid->periodo / pid->ganhos.tauFiltroDerivativo);
                pid->bd0 = kd / pid->ganhos.tauFiltroDerivativo;
            } else {
                pid->ad1 = 0.0f; pid->bd0 = 0.0f;
            }
            break;

        case PID_DISC_EULER_REGRESSIVO: 
            // Aproximação: s ~ (z-1)/(z*Ts)
            // Integral: I[k] = I[k-1] + Ki*Ts * e[k] 
            // (Euler Regressivo usa o erro atual "retangular à direita")
            pid->bi0 = ki * pid->periodo;
            pid->bi1 = 0.0f;

            // Derivativo Filtrado:
            // Equação: D[k] = (Tf / (Tf + Ts)) * D[k-1] + (Kd / (Tf + Ts)) * (e[k] - e[k-1])
            {
                float den = pid->ganhos.tauFiltroDerivativo + pid->periodo;
                pid->ad1 = pid->ganhos.tauFiltroDerivativo / den;
                pid->bd0 = kd / den;
            }
            break;

        case PID_DISC_BILINEAR: 
        default:
            // Aproximação: s ~ (2/Ts) * (z-1)/(z+1)
            // Integral (Trapezoidal): I[k] = I[k-1] + (Ki*Ts/2)*(e[k] + e[k-1])
            pid->bi0 = ki * pid->periodo * 0.5f;
            pid->bi1 = pid->bi0;

            // Derivativo Filtrado:
            // Equação: D[k] = ( (2*Tf - Ts)/(2*Tf + Ts) ) * D[k-1] + ( 2*Kd/(2*Tf + Ts) ) * (e[k] - e[k-1])
            {
                float den = 2.0f * pid->ganhos.tauFiltroDerivativo + pid->periodo;
                pid->ad1 = (2.0f * pid->ganhos.tauFiltroDerivativo - pid->periodo) / den;
                pid->bd0 = (2.0f * kd) / den;
            }
            break;
    }

    // Reinicia os estados
    pid->saida = 0.0f;
    pid->integral = 0.0f;
    pid->erroAnterior = 0.0f;
    pid->derivativo = 0.0f;
    pid->erro = 0.0f;
}

/**
 * @brief Executa o cálculo do PID para um passo de tempo.
 * * Esta função deve ser chamada na interrupção de controle (ISR).
 * Implementa a lei de controle: u[k] = P[k] + I[k] + D[k].
 *
 * Estratégia Anti-Windup: clamping (congelamento).
 * A integração é interrompida se a saída saturar E o erro tentar aumentar a saturação.
 * @param pid Ponteiro para a estrutura do PID.
 * @param erro O valor do erro medido no instante atual.
 */
static inline void PID_calcula(PID_t *pid, float erro)
{
    pid->erro = erro;

    /* -------------------------------------------------------------------------
       1. Termo Proporcional: P[k] = Kp * e[k]
       ------------------------------------------------------------------------- */
    float P = pid->ganhos.kp * pid->erro;

    /* -------------------------------------------------------------------------
       2. Termo Integral: I[k] = I[k-1] + bi0*e[k] + bi1*e[k-1]
       ------------------------------------------------------------------------- */
    // Calcula tentativa de integração. Só será efetivada se não houver saturação (clamping).
    float I_tentativa = pid->integral + pid->bi0 * pid->erro + pid->bi1 * pid->erroAnterior;

    /* -------------------------------------------------------------------------
       3. Termo Derivativo filtrado: D[k] = ad1*D[k-1] + bd0*(e[k] - e[k-1])
       ------------------------------------------------------------------------- */
    // Implementa filtro passa-alta na banda do sinal e passa-baixa para ruído (band-limited differentiator).
    float delta_erro = pid->erro - pid->erroAnterior;
    float D = (pid->ad1 * pid->derivativo) + (pid->bd0 * delta_erro);
    
    // Atualiza estado do derivativo imediatamente
    pid->derivativo = D;

    /* -------------------------------------------------------------------------
       4. Cálculo da Saída Total
       ------------------------------------------------------------------------- */
    float saidaCalculada = P + I_tentativa + D;

    /* -------------------------------------------------------------------------
       5. Saturação e Anti-Windup (Método Clamping)
       ------------------------------------------------------------------------- */
    if (saidaCalculada > pid->satMax)
    {
        pid->saida = pid->satMax;
        
        // Anti-Windup: Se saturou no topo, só permite atualizar a integral 
        // se o erro for negativo (ajudando a sair da saturação).
        if (erro < 0.0f) { 
            pid->integral = I_tentativa;
        }
    }
    else if (saidaCalculada < pid->satMin)
    {
        pid->saida = pid->satMin;
        
        // Anti-Windup: Se saturou no fundo, só permite atualizar a integral
        // se o erro for positivo.
        if (erro > 0.0f) {
            pid->integral = I_tentativa;
        }
    }
    else
    {
        // Região linear: Atualiza saída e integral normalmente.
        pid->saida = saidaCalculada;
        pid->integral = I_tentativa;
    }

    // Atualiza memória do erro para o próximo ciclo
    pid->erroAnterior = pid->erro;
}

#endif // PID_H