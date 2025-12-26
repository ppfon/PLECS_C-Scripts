entrada = InputSignal(0, 0);
omegaEstimado = InputSignal(0, 1);

SOGI_calcula(&sogi, entrada, omegaEstimado);


OutputSignal(0, 0) = sogi.conexoes.saidaEmFase; // Saída em fase
OutputSignal(0, 1) = sogi.conexoes.saidaEmQuadratura; // Saída em quadratura