/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <complex.h>

#include "rede.h"
#include "matrix.h"
#include "vector.h"
#include "log.h"
#include "error.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#define len(array) ((&array)[1] - array)

static int caso = 0;

// G B
static matrix_t* Y[2] = { MAT_NULL, MAT_NULL };

// [thetaPQ ... thetaPV ... thetaSwing]
static vector_t* thetas = VEC_NULL;

// [VPQ ... VPV ... VSwing]
static vector_t* Vs = VEC_NULL;

// [PespPQ ... PespPV ... PespSwing]
static vector_t* Pesp  = VEC_NULL;

// Tensao nominal
static vector_t* V_nom = VEC_NULL;

// Vetor que mapeia posições originais das barras para
// a ordem [PQ PV Swing]
static vector_t* map = VEC_NULL;

// Numero total de barras de cada tipo
// 0 -> PQ   1 -> PV    2 -> SW
static int n_barras[3] = { 0 };

// Barras e trechos pedidos
// Para as tabelas
static const int barras[4][10] = {
    { 0,   1,    2,    3,    4,   -1,   -1,   -1,   -1,   -1 },
    { 2,  11,   25,   28,   30,   42,   43,   47,   48,   49 },
    { 0,   1,   47,  633, 1414, 1429, 1528, 1607, 1609, 1636 },
    { 3, 990, 1310, 1466, 3947, 4105, 4188, 5820, 5830, 5840 },
};

static const int trechos[4][10][2] = {
    {{ 0,    1 }, {   0,   4 }, {   1,    2 }, {    2,    3 }, {    2,    4 }, {    3,    4 }, {   -1,  -1 }, {   -1,   -1 }, {   -1,   -1 }, {   -1,   -1 }},
    {{ 3,    4 }, {   6,   5 }, {  12,   28 }, {   13,    9 }, {   17,    1 }, {   18,    2 }, {   19,  20 }, {   24,   52 }, {   60,   62 }, {   75,    2 }},
    {{ 0, 1185 }, {   1,   2 }, {   1,   92 }, {   47,    6 }, {   47,   31 }, {  633,  632 }, {  633, 634 }, { 1414, 1415 }, { 1607,  286 }, { 1621, 1622 }},
    {{ 0, 1185 }, { 710, 543 }, { 776, 1748 }, { 1543, 1542 }, { 1600, 1387 }, { 1631, 1630 }, { 1748, 776 }, { 2867, 2868 }, { 2878, 2877 }, { 3640, 3947 }}
};

vector_t* prepara_rede(const int n) {
    FILE* fp;

    char* dados_barras = "";
    char* dados_ynodal = "";

    switch (n) {
        case 1:
            dados_barras = "Redes/1_Stevenson/1_Stevenson_DadosBarras.txt";
            dados_ynodal = "Redes/1_Stevenson/1_Stevenson_Ynodal.txt";
            break;

        case 2:
            dados_barras = "Redes/2_Reticulada/2_Reticulada_DadosBarras.txt";
            dados_ynodal = "Redes/2_Reticulada/2_Reticulada_Ynodal.txt";
            break;

        case 3:
            dados_barras = "Redes/3_Distribuicao_Primaria/3_Distribuicao_Primaria_DadosBarras.txt";
            dados_ynodal = "Redes/3_Distribuicao_Primaria/3_Distribuicao_Primaria_Ynodal.txt";
            break;

        case 4:
            dados_barras = "Redes/4_Distribuicao_Pri_Sec/4_Distribuicao_Primaria_Secundaria_DadosBarras.txt";
            dados_ynodal = "Redes/4_Distribuicao_Pri_Sec/4_Distribuicao_Primaria_Secundaria_Ynodal.txt";
            break;

        default:
            error(ERR_OOB, "prepara_rede");
    }

    caso = n - 1;

    if ((fp = fopen(dados_barras, "r")) == NULL)
        error(ERR_OPEN_FILE, "prepara_rede");

    int total_barras;
    if (fscanf(fp, "%d", &total_barras) != 1)
        error(ERR_READ_FILE, "prepara_rede");
    int pos_file = ftell(fp);

    thetas = vector_create(total_barras);
    Vs     = vector_create(total_barras);
    V_nom  = vector_create(total_barras);
    Pesp   = vector_create(total_barras);
    map    = vector_create(total_barras);

    int ind, tipo;
    double c3, c4, c5;

    // Ve quantidade de cada barra
    for (int i = 0; i < total_barras; i++) {
        if (fscanf(fp, "%d %d %lf %lf %lf", &ind, &tipo, &c3, &c4, &c5) != 5)
            error(ERR_READ_FILE, "prepara_rede");

        n_barras[tipo]++;
        vector_set(map, i, ind);
    }

    fseek(fp, pos_file, SEEK_SET);

    // [thetaPQ thetaPV VPQ]
    vector_t* x0 = vector_create(2*n_barras[0] + n_barras[1]);

    // Posicoes nos vetores
    int pos_pq = 0;
    int pos_pv = n_barras[0];
    int pos_sw = n_barras[0] + n_barras[1];
    int pos_x = n_barras[0] + n_barras[1];

    for (int i = 0; i < total_barras; i++) {
        if (fscanf(fp, "%d %d %lf %lf %lf", &ind, &tipo, &c3, &c4, &c5) != 5)
            error(ERR_READ_FILE, "prepara_rede");

        switch (tipo) {
            case 0: // PQ
                vector_set(map, i, pos_pq);
                vector_set(V_nom, pos_pq, c3); // Vnominal
                vector_set(Vs, pos_pq++, c3); // V
                vector_set(x0,  pos_x++, c3); // V
                break;

            case 1: // PV
                vector_set(map, i, pos_pv);
                vector_set(V_nom, pos_pv, c3); // Vnominal
                vector_set(Vs, pos_pv, c5); // V
                vector_set(Pesp, pos_pv++, c4); // Pesp
                break;

            case 2: // Swing
                vector_set(map, i, pos_sw);
                vector_set(V_nom, pos_sw, c3); // Vnominal
                vector_set(Vs, pos_sw, c4); // V
                vector_set(thetas, pos_sw++, c5*M_PI/180.0); // thetas em rad
                break;
        }
    }

    fclose(fp);

    if ((fp = fopen(dados_ynodal, "r")) == NULL)
        error(ERR_OPEN_FILE, "prepara_rede");

    Y[0] = matrix_create(total_barras, total_barras);
    Y[1] = matrix_create(total_barras, total_barras);

    int lines;
    if (fscanf(fp, "%d", &lines) != 1)
        error(ERR_READ_FILE, "prepara_rede");

    int j, k;
    for (int i = 0; i < lines; i++) {
        if (fscanf(fp, "%d %d %lf %lf", &j, &k, &c3, &c4) != 4)
            error(ERR_READ_FILE, "prepara_rede");

        matrix_set(Y[0], vector_get(map, j), vector_get(map, k), c3);
        matrix_set(Y[1], vector_get(map, j), vector_get(map, k), c4);
    }

    fclose(fp);

    return x0;
}

vector_t* F_rede(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "F_rede");

    vector_t* F = vector_create(2*n_barras[0] + n_barras[1]);

    // Atualiza Vs e thetas
    for (int i = 0; i < n_barras[0] + n_barras[1]; i++) {
        vector_set(thetas, i, vector_get(x, i));
        if (i < n_barras[0])
            vector_set(Vs, i, vector_get(x, n_barras[0] + n_barras[1] + i));
    }

    double sum = 0;
    double theta_kj = 0;

    // fps
    for (int i = 0; i < n_barras[0] + n_barras[1]; i++) {
        sum = 0;
        for (int j = 0; j < Y[0]->l; j++) {
            if (i == j) continue;

            theta_kj = vector_get(thetas, j) - vector_get(thetas, i);
            sum += vector_get(Vs, j) *
                (matrix_get(Y[0], i, j) * cos(theta_kj) -
                 matrix_get(Y[1], i, j) * sin(theta_kj));
        }

        vector_set(F, i, pow(vector_get(Vs, i), 2)*matrix_get(Y[0], i, i) + vector_get(Vs, i) * sum - vector_get(Pesp, i));
    }

    // fqs
    for (int i = 0; i < n_barras[0]; i++) {
        sum = 0;
        for (int j = 0; j < Y[0]->l; j++) {
            if (i == j) continue;

            theta_kj = vector_get(thetas, j) - vector_get(thetas, i);
            sum += vector_get(Vs, j) *
                (matrix_get(Y[0], i, j) * sin(theta_kj) +
                matrix_get(Y[1], i, j) * cos(theta_kj));
        }

        vector_set(F, n_barras[0] + n_barras[1] + i, -pow(vector_get(Vs, i), 2)*matrix_get(Y[1], i, i)-vector_get(Vs, i) * sum);
    }

    // print_vector(F);

    return F;
}

matrix_t* J_rede(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "J_rede");

    // Numero de barras PQ e PV
    int n = n_barras[0] + n_barras[1];

    // Atualiza Vs e thetas
    for (int i = 0; i < n; i++) {
        vector_set(thetas, i, vector_get(x, i));
        if (i < n_barras[0])
            vector_set(Vs, i, vector_get(x, n + i));
    }

    matrix_t* J = matrix_create(2*n_barras[0] + n_barras[1], 2*n_barras[0] + n_barras[1]);

    double sum = 0;
    double theta_kj = 0;


    // dfps
    for (int i = 0; i < n; i++) {
        // /dThetas
        for (int j = 0; j < n; j++) {
            // Propria barra
            if (i == j) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(thetas, k) - vector_get(thetas, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * sin(theta_kj) +
                        matrix_get(Y[1], i, k) * cos(theta_kj));
                }

                matrix_set(J, i, j, sum * vector_get(Vs, i));
            } else {
                theta_kj = vector_get(thetas, j) - vector_get(thetas, i);
                matrix_set(J, i, j,
                    -vector_get(Vs, i) * vector_get(Vs, j) *
                    (matrix_get(Y[0], i, j) * sin(theta_kj) +
                     matrix_get(Y[1], i, j) * cos(theta_kj)));
            }
        }

        // /dVs
        for (int j = n; j < J->c; j++) {
            // Propria barra
            if (j - n == i) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(thetas, k) - vector_get(thetas, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * cos(theta_kj) -
                         matrix_get(Y[1], i, k) * sin(theta_kj));
                }

                matrix_set(J, i, j, sum + 2 * vector_get(Vs, i) * matrix_get(Y[0], i, i));
            } else {
                theta_kj = vector_get(thetas, j - n) - vector_get(thetas, i);
                matrix_set(J, i, j, vector_get(Vs, i) *
                    (matrix_get(Y[0], i, j - n) * cos(theta_kj) -
                     matrix_get(Y[1], i, j - n) * sin(theta_kj)));
            }
        }
    }

    // dfqs
    for (int i = 0; i < n_barras[0]; i++) {
        // /dThetas
        for (int j = 0; j < n; j++) {
            // Propria barra
            if (i == j) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(thetas, k) - vector_get(thetas, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * cos(theta_kj) -
                        matrix_get(Y[1], i, k) * sin(theta_kj));
                }

                matrix_set(J, i + n, j, sum * vector_get(Vs, i));
            } else {
                theta_kj = vector_get(thetas, j) - vector_get(thetas, i);
                matrix_set(J, i + n, j,
                    -vector_get(Vs, i) * vector_get(Vs, j) *
                    (matrix_get(Y[0], i, j) * cos(theta_kj) -
                     matrix_get(Y[1], i, j) * sin(theta_kj)));
            }
        }

        // /dVs
        for (int j = n; j < J->c; j++) {
            // Propria barra
            if (j - n == i) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(thetas, k) - vector_get(thetas, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * sin(theta_kj) +
                         matrix_get(Y[1], i, k) * cos(theta_kj));
                }

                matrix_set(J, i + n, j, -sum - 2 * vector_get(Vs, i) * matrix_get(Y[1], i, i));
            } else {
                theta_kj = vector_get(thetas, j - (n)) - vector_get(thetas, i);
                matrix_set(J, i + n, j, -vector_get(Vs, i) *
                    (matrix_get(Y[0], i, j - (n)) * sin(theta_kj) +
                     matrix_get(Y[1], i, j - (n)) * cos(theta_kj)));
            }
        }
    }

    // print_matrix(J);

    return J;
}

void finaliza_rede(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "finaliza_rede");

    log_debug("Finalizando Rede");

    log_debug("Atualizando Vs e Thetas");
    // Atualiza Vs e thetas
    for (int i = 0; i < n_barras[0] + n_barras[1]; i++) {
        vector_set(thetas, i, vector_get(x, i));
        if (i < n_barras[0])
            vector_set(Vs, i, vector_get(x, n_barras[0] + n_barras[1] + i));
    }
    log_debug("Concluido");

    log_debug("Calculando dados da tabela 1");
    // Imprime tabela 1
    printf(" ___________________________________________________________ \n");
    printf("|       |       Tensão Complexa        |  Modulo da Tensão  |\n");
    printf("| Barra |  Modulo (pu)  |  Angulo (g)  |    Complexa  (V)   |\n");

    int ind;
    double graus = 0;
    for (int i = 0; i < len(barras[caso]); i++) {
        if (barras[caso][i] == -1)
            break;
        ind = vector_get(map, barras[caso][i]);
        graus = vector_get(thetas, ind)*180.0/M_PI;
        printf("| %5d |  %11.6f  |  %+10.4f  |  %16.3f  |\n", barras[caso][i], vector_get(Vs, ind)/vector_get(V_nom, ind), graus, vector_get(Vs, ind));
    }
    printf(" ");

    // Borda inferior
    // for (int i = 0; i < 59; i++)
    //     printf("\u203E");
    printf(" \n");
    log_debug("Concluido");

    // Calculo das potencias nos trechos pedidos para a tabela 2
    log_debug("Calculando dados da tabela 2");

    // Imprime tabela 2
    printf(" ________________________________________________________________ \n");
    printf("|       Trecho      |                       |                    |\n");
    printf("| Inicial |  Final  |  Potencia Ativa (kW)  |  Perda Ativa (kW)  |\n");

    int ind1, ind2;
    for (int i = 0; i < len(trechos[caso]); i++) {
        if (trechos[caso][i][0] == -1)
            break;

        ind1 = vector_get(map, trechos[caso][i][0]);
        ind2 = vector_get(map, trechos[caso][i][1]);

        double complex Vi = vector_get(Vs, ind1)*cos(vector_get(thetas, ind1)) + vector_get(Vs, ind1)*sin(vector_get(thetas, ind1))*I;
        double complex Vj = vector_get(Vs, ind2)*cos(vector_get(thetas, ind2)) + vector_get(Vs, ind2)*sin(vector_get(thetas, ind2))*I;

        double complex Yij = matrix_get(Y[0], ind1, ind2) + matrix_get(Y[1], ind1, ind2)*I;
        double complex Iij = (Vi - Vj) * (-Yij);
        double complex Sij = 0.001 * Vi * conj(Iij) * 3.0;
        double perda = pow(cabs(Vi - Vj), 2.0) * (-creal(Yij)) * 3.0 / 1000.0; // k

        printf("| %7d |  %5d  |  %+19.3f  |  %+16.3f  |\n", trechos[caso][i][0], trechos[caso][i][1], creal(Sij), perda);
    }
    printf(" ");

    // Borda inferior
    // for (int i = 0; i < 64; i++)
    //     printf("\u203E");
    printf(" \n");
    log_debug("Concluido");

    log_debug("Calculando dados da tabela 3");
    // Calculo Geral
    double P = 0.0;
    double soma_ps = 0.0;
    double soma_perda = 0;
    double theta_ji = 0.0;
    for (int i = 0; i < map->size; i++) {
        P = 0.0;

        for (int j = 0; j < map->size; j++) {
            theta_ji = vector_get(thetas, j) - vector_get(thetas, i);
            P += vector_get(Vs, i) * vector_get(Vs, j) * (matrix_get(Y[0], i, j)*cos(theta_ji) - matrix_get(Y[1], i, j)*sin(theta_ji));

            if (j > i) {
                double complex Vi = vector_get(Vs, i)*cos(vector_get(thetas, i)) + vector_get(Vs, i)*sin(vector_get(thetas, i))*I;
                double complex Vj = vector_get(Vs, j)*cos(vector_get(thetas, j)) + vector_get(Vs, j)*sin(vector_get(thetas, j))*I;

                soma_perda += 3 * pow(cabs(Vi - Vj), 2.0) * (-matrix_get(Y[0], i, j)) / 1000.0;
            }
        }

        soma_ps += 3 * P / 1000.0;
    }

    // Imprime tabela 3
    printf(" ____________________________________________________________ \n");
    printf("|                                           |    Valor (kW)  |\n");
    printf("| Potência ativa total gerada               |  %+12.3f  |\n", soma_ps);
    printf("| Potência ativa total de carga (absorvida) |  %+12.3f  |\n", soma_ps - soma_perda);
    printf("| Perda ativa total                         |  %+12.3f  |\n ", soma_perda);
    // Borda inferior
    // for (int i = 0; i < 60; i++)
    //     printf("\u203E");
    printf(" \n");
    log_debug("Concluido");
}
