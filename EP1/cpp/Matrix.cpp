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
#include <stdexcept>

#include "Matrix.hpp"

namespace MAP3121 {
    Matrix::Matrix(int m, int n) : m(m), n(n) {
        if (m <= 0 || n <= 0)
            throw std::out_of_range("[ERR]: " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + " negative or 0 index");

        data.resize(m);
        for (int i = 0; i < m; i++)
            data.at(i).resize(n);
    }
}

