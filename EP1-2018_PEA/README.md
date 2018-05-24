# MAP3121 e PEA3301
## Exercício Programa 1

Daniel Nery Silva de Oliveira - 9349051

23/05/2018

-------------------

## Introdução

O objetivo desse exercício programa foi utilizar o método de Newton para achar a raiz de sistemas de equações, atravé da resolução de sistemas lineares com decomposição LU.

## Arquivos

O código está dividido em arquivos `.c` e `.h`, nas pastas `src` e `inc` respectivamente. Aqui um breve resumo do que é encontrado em cada arquivo, mais informações podem ser encontradas nos .h:
* `log.c`:

  Logging para o terminal, retirado de https://github.com/rxi/log.c, usado para printar boa parte do código

* `erro.c`:

  Função de erro para ser chamada pelos outros arquivos, envolve o logger acima, adicionando mais informações.

* `lu.c`:

  Arquivo que contém a implementação da decomposição LU, como pedido no enunciado e a resolução de sistem Ax = b, com A sendo a matriz decomposta, também inclui uma versão paralela da decomposição, utilizando a biblioteca OpenMP para casos em que a matriz é muito grande.

* `matrix.c` e `vector.c`:

  Funções auxiliares para criação e manipulação de vetores e matrizes

* `newton.c`:

  Função que acha as raizes de um sistema de equações pelo método de Newton, utilizando sucessivas decomposições LU e soluções de sistemas lineares.

* `rede.c`:

  Funções para leitura, cálculo de F e J e finalização das redes propostas.

* `testes.c`:

    Funções para o cálculo de F e J dos testes iniciais propostos

Além dos arquivos fonte, há na raiz o enunciado do EP, que explica todos os métodos numéricos utilizados, e um Makefile, que deve ser utilizado na compilação do exercício, além deste arquivo.

## Testes Iniciais

Todos os testes iniciais foram executados com sucesso e suas saídas podem ser conferidas no arquivo `saida_testes.txt`, o teste 1 necessitou de 2 iterações para convergir, independente do valor incial, dado que são retas.

## Cálculo das Redes

Os cálculos referentes às redes propostas foram realizados com sucesso, abaixo os resultados pedidos:

### Rede 1
| Barra |  Modulo (pu)  |  Angulo (g)  |  Tensão Complexa  (V) |
| ----- | ------------- | ------------ | ------------------ |
|     0 |     1.000000  |     +0.0000  |        132790.562  |
|     1 |     0.939956  |     -5.2930  |        124817.271  |
|     2 |     1.000000  |     -2.3990  |        132790.562  |
|     3 |     0.911989  |     -8.9880  |        121103.490  |
|     4 |     0.945681  |     -5.2480  |        125577.508  |


| Barra Inicial | Barra Final  |  Potencia Ativa (kW)  |  Perda Ativa (kW)  |
| ------- | ------- | --------------------- | ------------------ |
|       0 |      1  |           +57547.523  |         +1627.615  |
|       0 |      4  |           +75462.604  |         +2003.022  |
|       1 |      2  |           -45684.555  |         +1105.219  |
|       2 |      3  |           +35898.252  |         +1386.122  |
|       2 |      4  |           +27311.974  |          +597.499  |
|       3 |      4  |           -23708.505  |          +448.990  |


|                                           |    Valor (kW)  |
| ----------------------------------------  | ---------------|
| Potência ativa total gerada               |   +378010.127  |
| Potência ativa total de carga (absorvida) |   +370841.661  |
| Perda ativa total                         |     +7168.466  |


### Rede 2
| Barra |  Modulo (pu)  |  Angulo (g)  |  Tensão Complexa  (V) |
| ----- | ------------- | ------------ | ------------------ |
|     2 |     0.969477  |     -2.3034  |           123.140  |
|    11 |     0.970124  |     -2.3430  |           123.222  |
|    25 |     0.973131  |     -1.9674  |           123.604  |
|    28 |     0.972731  |     -1.9701  |           123.553  |
|    30 |     0.969597  |     -2.3468  |           123.155  |
|    42 |     0.991272  |     -0.7801  |          7897.894  |
|    43 |     0.991236  |     -0.7801  |          7897.605  |
|    47 |     0.991905  |     -0.7744  |          7902.941  |
|    48 |     0.991666  |     -0.7766  |          7901.034  |
|    49 |     0.991480  |     -0.7784  |          7899.554  |


| Barra Inicial | Barra Final  |  Potencia Ativa (kW)  |  Perda Ativa (kW)  |
| ------- | ------- | --------------------- | ------------------ |
|       3 |      4  |            +1106.332  |            +0.253  |
|       6 |      5  |             -516.354  |            +0.056  |
|      12 |     28  |          +116901.478  |       +114883.571  |
|      13 |      9  |            -2072.832  |       +114937.309  |
|      17 |      1  |               -7.314  |            +0.041  |
|      18 |      2  |               +7.399  |            +0.041  |
|      19 |     20  |            +1040.765  |            +0.228  |
|      24 |     52  |            +6680.793  |            -0.000  |
|      60 |     62  |          +116922.663  |       +114851.424  |
|      75 |      2  |          +116913.444  |       +114849.410  |


|                                           |    Valor (kW)  |
| ----------------------------------------  | ---------------|
| Potência ativa total gerada               |     +6680.793  |
| Potência ativa total de carga (absorvida) |  -2865946.746  |
| Perda ativa total                         |  +2872627.540  |


### Rede 3
| Barra |  Modulo (pu)  |  Angulo (g)  |  Tensão Complexa  (V) |
| ----- | ------------- | ------------ | ------------------ |
|     0 |     1.000000  |     +0.0000  |          7967.434  |
|     1 |     0.950002  |     -1.6971  |          7569.075  |
|    47 |     0.949670  |     -1.6986  |          7566.435  |
|   633 |     0.952944  |     -1.4924  |          7592.516  |
|  1414 |     0.959761  |     -1.1210  |          7646.833  |
|  1429 |     0.958100  |     -1.1904  |          7633.600  |
|  1528 |     0.955043  |     -1.3299  |          7609.242  |
|  1607 |     0.997308  |     -0.1625  |          7945.985  |
|  1609 |     0.999228  |     -0.0489  |          7961.285  |
|  1636 |     0.997308  |     -0.1625  |          7945.985  |


| Barra Inicial | Barra Final  |  Potencia Ativa (kW)  |  Perda Ativa (kW)  |
| ------- | ------- | --------------------- | ------------------ |
|       0 |   1185  |            +1853.159  |            +0.156  |
|       1 |      2  |             +315.603  |            +0.018  |
|       1 |     92  |             -315.603  |            +0.018  |
|      47 |      6  |             -294.432  |            +0.029  |
|      47 |     31  |             +294.432  |            +0.017  |
|     633 |    632  |             +673.094  |            +0.069  |
|     633 |    634  |             -673.094  |            +0.076  |
|    1414 |   1415  |            +1544.286  |            +0.138  |
|    1607 |    286  |               +0.000  |            +0.000  |
|    1621 |   1622  |            +1850.432  |            +0.148  |


|                                           |    Valor (kW)  |
| ----------------------------------------  | ---------------|
| Potência ativa total gerada               |     +1853.159  |
| Potência ativa total de carga (absorvida) |     +1781.390  |
| Perda ativa total                         |       +71.769  |


### Rede 4
| Barra |  Modulo (pu)  |  Angulo (g)  |  Tensão Complexa  (V) |
| ----- | ------------- | ------------ | ------------------ |
|     3 |     0.951245  |     -1.6736  |          7578.981  |
|   990 |     0.994312  |     -0.2686  |          7922.114  |
|  1310 |     0.952854  |     -1.1879  |          7591.798  |
|  1466 |     0.956551  |     -1.2977  |          7621.260  |
|  3947 |     0.929398  |     -2.8625  |           203.903  |
|  4105 |     0.936458  |     -1.4129  |           205.452  |
|  4188 |     0.905968  |     -2.5287  |           198.763  |
|  5820 |     0.928335  |     -2.4081  |           235.829  |
|  5830 |     0.947113  |     -1.4556  |           240.599  |
|  5840 |     0.947269  |     -1.4574  |           240.638  |


| Barra Inicial | Barra Final  |  Potencia Ativa (kW)  |  Perda Ativa (kW)  |
| ------- | ------- | --------------------- | ------------------ |
|       0 |   1185  |            +1810.442  |            +0.149  |
|     710 |    543  |              +51.533  |            +0.002  |
|     776 |   1748  |           +13473.327  |        +13086.688  |
|    1543 |   1542  |             +660.532  |            +0.039  |
|    1600 |   1387  |            +1517.924  |            +0.149  |
|    1631 |   1630  |              +25.681  |            +0.000  |
|    1748 |    776  |             -386.638  |        +13086.688  |
|    2867 |   2868  |               +5.605  |            +0.010  |
|    2878 |   2877  |               -4.107  |            +0.014  |
|    3640 |   3947  |               +6.548  |            +0.001  |


|                                           |    Valor (kW)  |
| ----------------------------------------  | ---------------|
| Potência ativa total gerada               |     +1810.442  |
| Potência ativa total de carga (absorvida) |  -1089232.069  |
| Perda ativa total                         |  +1091042.511  |
