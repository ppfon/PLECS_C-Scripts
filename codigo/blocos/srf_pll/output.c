// 1. Inicialização
// Fazemos aqui dentro para garantir que os parâmetros estejam corretos
if (estaInicializado == 0) 
{
    // Lê as variáveis do ambiente do PLECS  
    float kpPLL = ParamRealData(0,0);  
    float kiPLL = ParamRealData(1,0); 
    float freqRede = ParamRealData(2,0);
    float periodo = ParamRealData(3,0); 

    SRF_PLL_inicia(&pll, kpPLL, kiPLL, freqRede, periodo);
    estaInicializado = 1;
}

// 2. Leitura das entradas
DqZero_t tensaoRede;
tensaoRede.d = InputSignal(0, 0); // Porta 1, canal 1
tensaoRede.q = InputSignal(0, 1); // Porta 1, canal 2
tensaoRede.zero = 0; // Não utilizado

// 3. Execução da Lógica de Controle
SRF_PLL_executa(&pll, tensaoRede);

// 4. Escrita das saídas (Portas do PLECS)
OutputSignal(0, 0) = pll.vDQ.d;
OutputSignal(0, 1) = pll.omegaEstimado / DOIS_PI; // Em Hz
OutputSignal(0, 2) = pll.theta;