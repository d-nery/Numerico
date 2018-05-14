/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 2
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#include "chua.h"

double g(double v) {
    return
        (v) <= -Emax ? (Gc*(v) + Emax*(Gc - Gb) + E*(Gb - Ga)) :
        (v) <= -E    ? (Gb*(v) + (Gb - Ga)*E)                  :
        (v) < E      ? (Ga*(v))                                :
        (v) < Emax   ? (Gb*(v) + (Ga - Gb)*E)                  :
                       (Gc*(v) + Emax*(Gb - Gc) + E*(Ga - Gb));
}
