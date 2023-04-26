#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;
class Matrix {
protected:
    int numberRows;
    int numberColumns;
    vector<vector<double>> curMatrix;
public:
    Matrix() {
        this->numberRows = 0;
        this->numberColumns = 0;
    }
    explicit Matrix(vector<double>& v) {
        this->numberRows = (int)v.size();
        this->numberColumns = 1;
        this->curMatrix.resize(this->numberRows, vector<double>(this->numberColumns));
        for (int i = 0; i < this->numberRows; ++i) {
            for (int j = 0; j < this->numberColumns; ++j) {
                this->curMatrix[i][j] = v[i];
            }
        }
    }
    Matrix(int n, vector<double>& v) {
        this->numberRows = (int)v.size();
        this->numberColumns = n + 1;
        this->curMatrix.resize(this->numberRows, vector<double>(this->numberColumns));
        for (int i = 0; i < this->numberRows; ++i) {
            for (int j = 0; j < this->numberColumns; ++j) {
                this->curMatrix[i][j] = pow(v[i], j);
            }
        }
    }
    Matrix(int n, int m) {
        this->numberRows = n;
        this->numberColumns = m;
        this->curMatrix.resize(n, vector<double>(m));
    }
    int getIndForSwap(int cJ) {
        int ind = cJ;
        double mx = abs(this->curMatrix[ind][cJ]);
        for (int i = cJ + 1; i < this->numberRows; ++i) {
            if (mx < abs(this->curMatrix[i][cJ])) {
                mx = abs(this->curMatrix[i][cJ]);
                ind = i;
            }
        }
        return ind;
    }
    [[nodiscard]] vector<vector<double>>& getCurMatrix() {
        return curMatrix;
    }
    [[nodiscard]] int getNumberRows() const {
        return this->numberRows;
    }
    [[nodiscard]] int getNumberColumns() const {
        return this->numberColumns;
    }
    Matrix& transposeMatrix() {
        auto *newMatrix = new Matrix(this->numberColumns, this->numberRows);
        for (int i = 0; i < this->numberRows; ++i) {
            for (int j = 0; j < this->numberColumns; ++j) {
                newMatrix->curMatrix[j][i] = this->curMatrix[i][j];
            }
        }
        return *newMatrix;
    }
    Matrix& operator*(Matrix& matrix) {
        auto *newMatrix = new Matrix(this->numberRows, matrix.numberColumns);
        for (int i = 0; i < this->numberRows; ++i) {
            for (int j = 0; j < matrix.numberColumns; ++j) {
                for (int z = 0; z < this->numberColumns; ++z) {
                    newMatrix->curMatrix[i][j] += this->curMatrix[i][z] * matrix.curMatrix[z][j];
                    if (abs(newMatrix->curMatrix[i][j]) < 1e-10) {
                        newMatrix->curMatrix[i][j] = 0;
                    }
                }
            }
        }
        return *newMatrix;
    }
    friend istream& operator>>(istream& in, Matrix& matrix) {
        for (int i = 0; i < matrix.numberRows; ++i) {
            for (int j = 0; j < matrix.numberColumns; ++j) {
                in >> matrix.curMatrix[i][j];
            }
        }
        return in;
    }
    friend ostream& operator<<(ostream& out, Matrix& matrix) {
        for (int i = 0; i < matrix.numberRows; ++i) {
            for (int j = 0; j < matrix.numberColumns; ++j) {
                if (abs(matrix.curMatrix[i][j]) < 1e-10) {
                    matrix.curMatrix[i][j] = 0;
                }
                out << fixed << setprecision(4) << matrix.curMatrix[i][j] << ' ';
            }
            out << '\n';
        }
        return out;
    }
};
class IdenticalMatrix : public Matrix {
public:
    IdenticalMatrix() : Matrix() {}
    explicit IdenticalMatrix(int n) : Matrix(n, n) {
        for (int i = 0; i < this->numberRows; ++i) {
            this->curMatrix[i][i] = 1;
        }
    }
};
class EliminationMatrix : public  IdenticalMatrix {
public:
    EliminationMatrix() : IdenticalMatrix() {}
    explicit EliminationMatrix(int n) : IdenticalMatrix(n) {}
    void elimination(Matrix* A, int rI, int rJ) {
        this->curMatrix[rI][rJ] = (A->getCurMatrix()[rJ][rJ] != 0) ? -A->getCurMatrix()[rI][rJ] / A->getCurMatrix()[rJ][rJ] : 0;
    }
};
class PermutationMatrix : public  IdenticalMatrix {
public:
    PermutationMatrix() : IdenticalMatrix() {}
    explicit PermutationMatrix(int n) : IdenticalMatrix(n) {}
    void permutation(int rI, int rJ) {
        this->curMatrix[rI][rI] = 0;
        this->curMatrix[rI][rJ] = 1;
        this->curMatrix[rJ][rJ] = 0;
        this->curMatrix[rJ][rI] = 1;
    }
};
class AugmentedMatrix {
private:
    Matrix* A;
    Matrix* B;
public:
    AugmentedMatrix() = default;
    AugmentedMatrix(Matrix *a, Matrix *b) : A(a), B(b) {}
    void directWay() {
        int numRows = A->getNumberRows();
        for (int i = 0; i < numRows; ++i) {
            PermutationMatrix P(numRows);
            int swapInd = A->getIndForSwap(i);
            if (swapInd != i) {
                P.permutation(i, swapInd);
                *A = P * *A;
                *B = P * *B;
            }
            for (int j = i + 1; j < numRows; ++j) {
                EliminationMatrix E(numRows);
                E.elimination(A, j, i);
                if (E.getCurMatrix()[j][i] != 0) {
                    *A = E * *A;
                    *B = E * *B;
                }
            }
        }
    }
    void wayBack() {
        int numRows = A->getNumberRows();
        for (int i = numRows - 1; i >= 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                EliminationMatrix E(numRows);
                E.elimination(A, j, i);
                if (E.getCurMatrix()[j][i] != 0) {
                    *A = E * *A;
                    *B = E * *B;
                }
            }
        }
    }
    void diagonalNormalization() {
        int numRows = A->getNumberRows();
        for (int i = 0; i < numRows; ++i) {
            double coefficient = A->getCurMatrix()[i][i];
            for (int j = 0; j < numRows; ++j) {
                B->getCurMatrix()[i][j] /= coefficient;
            }
            A->getCurMatrix()[i][i] = 1;
        }
    }
    [[nodiscard]] Matrix *getB() const {
        return B;
    }
};
int main() {
    int m, n;
    cin >> m;
    vector<vector<double>> v(2, vector<double>(m));
    for (int i = 0; i < m; ++i) {
        cin >> v[0][i] >> v[1][i];
    }
    cin >> n;
    Matrix A(n, v[0]), ATA = A.transposeMatrix() * A;
    Matrix B(v[1]), ATB = A.transposeMatrix() * B;
    cout << "A:\n" << A;
    cout << "A_T*A:\n" << ATA;
    IdenticalMatrix I(ATA.getNumberRows());
    AugmentedMatrix augmentedMatrix(&ATA, &I);
    augmentedMatrix.directWay();
    augmentedMatrix.wayBack();
    augmentedMatrix.diagonalNormalization();
    Matrix invATA = *augmentedMatrix.getB();
    cout << "(A_T*A)^-1:\n" << invATA;
    cout << "A_T*b:\n" << ATB;
    cout << "x~:\n" << invATA * ATB;
}