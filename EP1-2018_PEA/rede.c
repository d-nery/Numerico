#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "rede.h"
#include "matrix.h"
#include "vector.h"
#include "log.h"
#include "error.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

static matrix_t* Y[2] = { MAT_NULL, MAT_NULL };

// [thetaPQ thetaPV thetaSwing]
static vector_t* theta = VEC_NULL;

// [VPQ VPV VSwing]
static vector_t* Vs = VEC_NULL;

static vector_t* tipos = VEC_NULL;
static vector_t* Pesp  = VEC_NULL;
static vector_t* map = VEC_NULL;

static int n_barras[3] = { 0 };
static double V_ref;

vector_t* prepara_rede(const int n) {
    FILE* fp;

    char* name;

    switch (n) {
        case 1:
            name = "Redes/1_Stevenson/1_Stevenson_DadosBarras.txt";
            break;

        case 2:
            name = "Redes/2_Reticulada/2_Reticulada_DadosBarras.txt";
            break;

        case 3:
            name = "Redes/3_Distribuicao_Primaria/3_Distribuicao_Primaria_DadosBarras.txt";
            break;

        case 4:
            name = "Redes/4_Distribuicao_Pri_Sec/4_Distribuicao_Primaria_Secundaria_DadosBarras.txt";
            break;

        default:
            error(ERR_OOB, "prepara_rede");
    }

    if ((fp = fopen(name, "r")) == NULL)
        error(ERR_OPEN_FILE, "prepara_rede");

    int size;
    fscanf(fp, "%d", &size);
    int pos_file = ftell(fp);

    tipos = vector_create(size);
    theta = vector_create(size);
    Vs    = vector_create(size);
    Pesp  = vector_create(size);

    map = vector_create(size);

    int ind, tipo;
    double c3, c4, c5;

    // Ve quantidade de cada barra
    for (int i = 0; i < tipos->size; i++) {
        if (fscanf(fp, "%d %d %lf %lf %lf", &ind, &tipo, &c3, &c4, &c5) != 5)
            error(ERR_READ_FILE, "prepara_rede");

        n_barras[tipo]++;
        vector_set(tipos, ind, tipo);
        vector_set(map, i, ind);
    }

    fseek(fp, pos_file, SEEK_SET);

    // [thetaPQ thetaPV VPQ]
    vector_t* x = vector_create(2*n_barras[0] + n_barras[1]);

    // Posicoes nos vetores
    int pos_pq = 0;
    int pos_pv = n_barras[0];
    int pos_sw = n_barras[0] + n_barras[1];
    int pos_x = n_barras[0] + n_barras[1];

    for (int i = 0; i < tipos->size; i++) {
        if (fscanf(fp, "%d %d %lf %lf %lf", &ind, &tipo, &c3, &c4, &c5) != 5)
            error(ERR_READ_FILE, "prepara_rede");

        switch (tipo) {
            case 0: // PQ
                vector_set(map, i, pos_pq);
                vector_set(Vs, pos_pq++, c3); // V
                vector_set(x,  pos_x++, c3); // V
                break;

            case 1: // PV
                vector_set(map, i, pos_pv);
                vector_set(Vs, pos_pv, c5); // V
                vector_set(Pesp, pos_pv++, c4); // Pesp
                break;

            case 2: // Swing
                vector_set(map, i, pos_sw);
                V_ref = c4;
                vector_set(Vs, pos_sw, c4); // V
                vector_set(theta, pos_sw++, c5); // theta
                break;
        }
    }

    fclose(fp);

    switch (n) {
        case 1:
            name = "Redes/1_Stevenson/1_Stevenson_Ynodal.txt", "r";
            break;

        case 2:
            name = "Redes/2_Reticulada/2_Reticulada_Ynodal.txt", "r";
            break;

        case 3:
            name = "Redes/3_Distribuicao_Primaria/3_Distribuicao_Primaria_Ynodal.txt", "r";
            break;

        case 4:
            name = "Redes/4_Distribuicao_Pri_Sec/4_Distribuicao_Primaria_Secundaria_Ynodal.txt", "r";
            break;

        default:
            error(ERR_OOB, "prepara_rede");
    }

    if ((fp = fopen(name, "r")) == NULL)
        error(ERR_OPEN_FILE, "prepara_rede");

    Y[0] = matrix_create(tipos->size, tipos->size);
    Y[1] = matrix_create(tipos->size, tipos->size);

    int lines;
    fscanf(fp, "%d", &lines);

    int j, k;
    for (int i = 0; i < lines; i++) {
        if (fscanf(fp, "%d %d %lf %lf", &j, &k, &c3, &c4) != 4)
            error(ERR_READ_FILE, "prepara_rede");

        matrix_set(Y[0], vector_get(map, j), vector_get(map, k), c3);
        matrix_set(Y[1], vector_get(map, j), vector_get(map, k), c4);
    }

    fclose(fp);

    return x;
}

vector_t* F_rede(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "F_rede");

    vector_t* F = vector_create(2*n_barras[0] + n_barras[1]);

    // Atualiza Vs e thetas
    for (int i = 0; i < n_barras[0] + n_barras[1]; i++) {
        vector_set(theta, i, vector_get(x, i));
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

            theta_kj = vector_get(theta, j) - vector_get(theta, i);
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

            theta_kj = vector_get(theta, j) - vector_get(theta, i);
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

    // Atualiza Vs e thetas
    for (int i = 0; i < n_barras[0] + n_barras[1]; i++) {
        vector_set(theta, i, vector_get(x, i));
        if (i < n_barras[0])
            vector_set(Vs, i, vector_get(x, n_barras[0] + n_barras[1] + i));
    }

    matrix_t* J = matrix_create(2*n_barras[0] + n_barras[1], 2*n_barras[0] + n_barras[1]);

    double sum = 0;
    double theta_kj = 0;

    // dfps
    for (int i = 0; i < n_barras[0] + n_barras[1]; i++) {
        // /dThetas
        for (int j = 0; j < n_barras[0] + n_barras[1]; j++) {
            // Propria barra
            if (i == j) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(theta, k) - vector_get(theta, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * sin(theta_kj) +
                        matrix_get(Y[1], i, k) * cos(theta_kj));
                }

                matrix_set(J, i, j, sum * vector_get(Vs, i));
            } else {
                theta_kj = vector_get(theta, j) - vector_get(theta, i);
                matrix_set(J, i, j,
                    -vector_get(Vs, i) * vector_get(Vs, j) *
                    (matrix_get(Y[0], i, j) * sin(theta_kj) +
                     matrix_get(Y[1], i, j) * cos(theta_kj)));
            }
        }

        // /dVs
        for (int j = n_barras[0] + n_barras[1]; j < J->c; j++) {
            // Propria barra
            if (j - (n_barras[0] + n_barras[1]) == i) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(theta, k) - vector_get(theta, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * cos(theta_kj) -
                         matrix_get(Y[1], i, k) * sin(theta_kj));
                }

                matrix_set(J, i, j, sum + 2 * vector_get(Vs, i) * matrix_get(Y[0], i, i));
            } else {
                theta_kj = vector_get(theta, j - (n_barras[0] + n_barras[1])) - vector_get(theta, i);
                matrix_set(J, i, j, vector_get(Vs, i) *
                    (matrix_get(Y[0], i, j - (n_barras[0] + n_barras[1])) * cos(theta_kj) -
                     matrix_get(Y[1], i, j - (n_barras[0] + n_barras[1])) * sin(theta_kj)));
            }
        }
    }

    // dfqs
    for (int i = 0; i < n_barras[0]; i++) {
        // /dThetas
        for (int j = 0; j < n_barras[0] + n_barras[1]; j++) {
            // Propria barra
            if (i == j) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(theta, k) - vector_get(theta, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * cos(theta_kj) -
                        matrix_get(Y[1], i, k) * sin(theta_kj));
                }

                matrix_set(J, i + n_barras[0] + n_barras[1], j, sum * vector_get(Vs, i));
            } else {
                theta_kj = vector_get(theta, j) - vector_get(theta, i);
                matrix_set(J, i + n_barras[0] + n_barras[1], j,
                    -vector_get(Vs, i) * vector_get(Vs, j) *
                    (matrix_get(Y[0], i, j) * cos(theta_kj) -
                     matrix_get(Y[1], i, j) * sin(theta_kj)));
            }
        }

        // /dVs
        for (int j = n_barras[0] + n_barras[1]; j < J->c; j++) {
            // Propria barra
            if (j - (n_barras[0] + n_barras[1]) == i) {
                sum = 0;
                for (int k = 0; k < Y[0]->l; k++) {
                    if (k == i) continue;

                    theta_kj = vector_get(theta, k) - vector_get(theta, i);
                    sum += vector_get(Vs, k) *
                        (matrix_get(Y[0], i, k) * sin(theta_kj) +
                         matrix_get(Y[1], i, k) * cos(theta_kj));
                }

                matrix_set(J, i + n_barras[0] + n_barras[1], j, -sum - 2 * vector_get(Vs, i) * matrix_get(Y[1], i, i));
            } else {
                theta_kj = vector_get(theta, j - (n_barras[0] + n_barras[1])) - vector_get(theta, i);
                matrix_set(J, i + n_barras[0] + n_barras[1], j, -vector_get(Vs, i) *
                    (matrix_get(Y[0], i, j - (n_barras[0] + n_barras[1])) * sin(theta_kj) +
                     matrix_get(Y[1], i, j - (n_barras[0] + n_barras[1])) * cos(theta_kj)));
            }
        }
    }

    // print_matrix(J);

    return J;
}

void finaliza_rede(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "finaliza_rede");

    // Atualiza Vs e thetas
    for (int i = 0; i < n_barras[0] + n_barras[1]; i++) {
        vector_set(theta, i, vector_get(x, i));
        if (i < n_barras[0])
            vector_set(Vs, i, vector_get(x, n_barras[0] + n_barras[1] + i));
    }

    // Radiano -> Grau
    for (int i = 0; i < theta->size; i++)
        vector_set(theta, i, vector_get(theta, i)*180/M_PI);

    // Imprime tabela 1
    printf(" ____________________________________________________ \n");
    printf("|       |       Tensão Complexa        |             |\n");

    unsigned char* t =" Barra |  Modulo (pu)  |  Angulo (g)  |   |V| (V)   ";
    printf("|");
    for (int i = 0; i < strlen(t); i++)
        printf("\u0332%c", t[i]);          // Com underlines
    printf("|\n");

    int ind;
    for (int i = 0; i < theta->size; i++) {
        ind = vector_get(map, i);
        printf("| %5d |    %8.6f   |   %s%6.4f    |  %10.3f |\n", i, vector_get(Vs, ind)/V_ref, vector_get(theta, ind) < 0 ? "" : " ", vector_get(theta, ind), vector_get(Vs, ind));
    }
    printf(" ");
    for (int i = 0; i < 52; i++)
        printf("\u203E");
    printf(" \n");
}
