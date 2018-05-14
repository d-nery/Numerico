/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <vector>

namespace MAP3121 {
    class Matrix {
        private:
        
        public:
            int m, n;
            std::vector< std::vector<double> > data;
            Matrix();
            Matrix(int m, int n);
    };
}

#endif