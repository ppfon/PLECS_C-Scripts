unsigned int tipo = (unsigned int)ParamRealData(0, 0);
ordem = ParamRealData(1, 0);
freqCorte = ParamRealData(2, 0);
freqAmostragem = ParamRealData(3, 0);
zeta = ParamRealData(4, 0);

FiltroTipo_t tipoEnum = (FiltroTipo_t)tipo;
Filtro_inicia(&filtro, tipoEnum, ordem, freqCorte, freqAmostragem, zeta);  