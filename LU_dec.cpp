#include <iostream>

double** new_matrix(int n, int m, int byte = 0){
    double** A = new double*[n];
    A[0] = new double[n*m];
    for(int i = 1; i < n; ++i){
        A[i] = A[i-1] + m;
    }
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            if(byte!=0){
                std::cin >> A[i][j];
            }
            else{
                A[i][j] = byte;
            }
        }
    }
    return A;
}

void show(double** A, int x, int y){
    std::cout << "\n";
    for(int i = 0; i < x; ++i){
        for(int j = 0; j < y; ++j){
            std::cout << A[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void LU_dec(double** A, double** L, double** U, int n){
    for(int k = 0; k < n; ++k){
        U[k][k] = A[k][k];
        for(int i = k; i < n; ++i){
            L[i][k] = A[i][k]/U[k][k];
            U[k][i] = A[k][i];
        }
        for(int i = k; i < n; ++i){
            for(int j = k; j < n; ++j){
                A[i][j] = A[i][j]-L[i][k]*U[k][j];
            }
        }
    }
}

int main(){
    int n, m;
    std::cin >> n >> m;
    double** A = new_matrix(n,m,1);
    double** L = new_matrix(n,m);
    double** U = new_matrix(n,m);
    LU_dec(A,L,U,n);
    std::cout << "L: " << "\n";
    show(L, n, n);
    std::cout << "U: " << "\n";
    show(U, n, n);
    return 0;
}
