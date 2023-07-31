//#include <iostream>
#include "matrix.h"

int main() {

    SquareMatrix<6, Residue<7>> s ({{1, 1, 1, 1, 1, 1},
                                        {1, 3, 2, 6, 4, 5},
                                        {1, 2, 4, 1, 2, 4},
                                        {1, 6, 1, 6, 1, 6},
                                        {1, 4, 2, 1, 4, 2},
                                        {1, 5, 4, 6, 2, 3}});

    Matrix<6, 1, Residue<7>> a ({{-1},
                                     {2},
                                     {0},
                                     {1},
                                     {0},
                                     {0}});


    SquareMatrix<6, Residue<7>> s1({{1, 1, 1, 1, 1, 1},
                                        {1, 5, 4, 6, 2, 3},
                                        {1, 4, 2, 1, 4, 2},
                                        {1, 6, 1, 6, 1, 6},
                                        {1, 2, 4, 1, 2, 4},
                                        {1, 3, 2, 6, 4, 5}});

    Matrix<6, 1, Residue<7>> a1 ({{1},{1},{1},{0},{0},{0}});

    a = s * a;
    a1 = s * a1;

    std::vector<std::vector<Residue<7>>> ae(6, std::vector<Residue<7>>(1));
    for(int i = 0; i < 6; i++)
        ae[i][0] = a[i][0] * a1[i][0];

    Matrix<6, 1, Residue<7>> a2(ae);

    std::cout << s1 * a2 * static_cast<Residue<7>>(6);

    return 0;
}