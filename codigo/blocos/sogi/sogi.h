#ifndef SOGI_H
#define SOGI_H

typedef struct {
    float entrada;
    float omegaEstimado;
    float saidaEmFase;      // v_alpha
    float saidaEmQuadratura; // v_beta
} Conexoes_t;

typedef struct {
    float entrada;         // Input atual (u[k])
    float entradaAnterior; // Input anterior (u[k-1])
    float saida;           // Output atual (y[k])
    float saidaAnterior;   // Output anterior (y[k-1])
    float periodo;         // Ts
    float satMin;
    float satMax;
} Integrador_t;

typedef struct {
    Conexoes_t conexoes;
    Integrador_t integradorEmFase;
    Integrador_t integradorEmQuadratura;
    float ganho; // k (Amortecimento, tipicamente 1.414)
} SOGI_t;

// --- Função do Integrador (Tustin/Trapezoidal) ---
static inline void Integrador_calcula(Integrador_t *integrador, float entrada)
{
    // y[n] = y[n-1] + (Ts/2)*(u[n] + u[n-1])
    integrador->saida = integrador->saidaAnterior + 
                        (0.5f * integrador->periodo) * (entrada + integrador->entradaAnterior);

    // Saturação
    if (integrador->saida > integrador->satMax) {
        integrador->saida = integrador->satMax;
    } else if (integrador->saida < integrador->satMin) {
        integrador->saida = integrador->satMin;
    }
    // Else implícito: mantém o valor calculado

    // Atualiza estados para o próximo ciclo
    integrador->saidaAnterior = integrador->saida;
    integrador->entradaAnterior = entrada;
}

// Inicialização (Mantive a sua, apenas removendo redundâncias)
static inline void Integrador_inicia(Integrador_t *integrador, float periodo, float satMin, float satMax)
{
    integrador->entrada = 0.0f;
    integrador->entradaAnterior = 0.0f;
    integrador->saida = 0.0f;
    integrador->saidaAnterior = 0.0f;
    integrador->periodo = periodo;
    integrador->satMin = satMin;
    integrador->satMax = satMax;
}

static inline void SOGI_inicia(SOGI_t *sogi, float ganho, float periodo, float satMin, float satMax)
{
    sogi->ganho = ganho;
    sogi->conexoes.omegaEstimado = 0.0f;
    sogi->conexoes.entrada = 0.0f;
    sogi->conexoes.saidaEmFase = 0.0f;
    sogi->conexoes.saidaEmQuadratura = 0.0f;

    Integrador_inicia(&sogi->integradorEmFase, periodo, satMin, satMax);
    Integrador_inicia(&sogi->integradorEmQuadratura, periodo, satMin, satMax);
}

// --- Função SOGI Corrigida ---
static inline void SOGI_calcula(SOGI_t *sogi, float entrada, float omegaEstimado)
{
    float erro;
    sogi->conexoes.omegaEstimado = omegaEstimado;
    float w = omegaEstimado; // Atalho para clareza

    sogi->conexoes.entrada = entrada;

    // 1. Cálculo do Erro usando a realimentação ANTERIOR (evita loop algébrico)
    // Erro = V_in - V_alpha_ant
    erro = entrada - sogi->integradorEmFase.saidaAnterior;

    // 2. Cálculo das entradas dos integradores (Derivadas)
    // Abordagem Simétrica: dot_alpha = w * (k*erro - v_beta)
    sogi->integradorEmFase.entrada = w * ((sogi->ganho * erro) - sogi->integradorEmQuadratura.saidaAnterior);

    // dot_beta = w * v_alpha
    sogi->integradorEmQuadratura.entrada = sogi->integradorEmFase.saidaAnterior;

    // 3. EXECUTA OS INTEGRADORES 
    Integrador_calcula(&sogi->integradorEmFase, sogi->integradorEmFase.entrada);
    Integrador_calcula(&sogi->integradorEmQuadratura, sogi->integradorEmQuadratura.entrada);

    // 4. Atualiza as saídas da estrutura principal
    sogi->conexoes.saidaEmFase = sogi->integradorEmFase.saida;
    sogi->conexoes.saidaEmQuadratura = sogi->integradorEmQuadratura.saida;
}

#endif //SOGI_H