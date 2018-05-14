/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#include <iostream>
#include <iomanip>

#include "Matrix.hpp"

using namespace MAP3121;
using namespace std;

int main(int argc, char** argv) {
    try {
        Matrix A(5, 10);

        for (int i = 0; i < A.m; i++) {
            for (int j = 0; j < A.n; j++) {
                std::cout << std::setw(5) << A.data.at(i).at(j) << " ";
            }
            std::cout << std::endl;
        }
    } catch (exception& e) {
        cout << e.what() << endl;
    }

    return 0;
}