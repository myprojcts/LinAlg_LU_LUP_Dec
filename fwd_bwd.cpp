#include <iostream>
#include <cstring>

void del_matrix(double** A, int row){
    for(int i = 1; i < row; ++i){
        delete[] (A[i - 1] + row) ;
    }
    delete[] A[0];
    A = nullptr;
}

double* lup_solve(double** L, double** U, int* p, double* b,
                  int n){
    double* x = new double [n];
    double* y = new double [n];
    std::memset(y, 0, sizeof(double)*n);
    std::memset(x, 0, sizeof(double)*n);
    for(int i = 0; i < n; ++i){
        double sum = 0;
        for(int j = 0; j < i; ++j){
            sum += L[i][j] * y[j];
        }
        y[i] = b[p[i]] - sum;
    }
    for(int i = 0; i < n; ++i){
        std::cout << y[i] << " ";
    }
    std::cout << "\n";
    for(int i = n-1; i >= 0; --i){
        double sum = 0;
        for(int j = n-1; j > i; --j){
            sum += U[i][j]*x[j];
        }
        x[i] = (y[i] - sum)/U[i][i];
    }
    delete [] y;
    y = nullptr;
    return x;
}

int main(){
    int n = 3;
    double* b = new double[n];
    int* p = new int[n];
    double** L = new double* [n];
    double** U = new double* [n];
    L[0] = new double[n * n]; // выделение двумерного массива без сегментации памяти
    U[0] = new double[n * n];
    for(int i = 1; i < n; ++i){
        L[i] = L[i - 1] + n ;
        U[i] = U[i - 1] + n;
    }
    std::cout << "b:" << "\n";
    for(int i = 0; i < n; ++i){
        std::cin >> b[i];
    }
    std::cout << "p:" << "\n";
    for(int i = 0; i < n; ++i){
        std::cin >> p[i];
    }
    std::cout << "L:" << "\n";
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            std::cin >> L[i][j];
        }
    }
    std::cout << "U:" << "\n";
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            std::cin >> U[i][j];
        }
    }
    double* ans = lup_solve(L,U,p,b,n);
    for(int i = 0; i < n; ++i){
        std::cout << ans[i] << " ";
    }
    delete [] ans;
    ans = nullptr;
    delete [] p;
    p = nullptr;
    delete [] b;
    b = nullptr;
    del_matrix(L, n);
    del_matrix(U, n);
    return 0;
}
