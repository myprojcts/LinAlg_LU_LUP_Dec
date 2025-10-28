#include <iostream>
#include <cmath>
#include <cstring>

void del_matrix(double** A){
    delete[] A[0];
    delete[] A;
    A = nullptr;
}

int rang_U(double** U, int n){
    double eps = 1e-10;
    int rang = 0;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if(std::abs(U[i][j]) > eps){
                rang++;
                break;
            }
        }
    }
    return rang;
}

double** matrix_realloc(double** A, int* n, int* m){
    int max = std::max(*n,*m);
    double** C = new double*[max];
    C[0] = new double[max * max];
    for(int i = 1; i != max; ++i){
        C[i] = C[i-1] + max;
    }
    for(int i = 0; i < max; ++i){
        for(int j = 0; j < max; ++j){
            C[i][j] = 0;
            if(i < *n && j < *m){
                C[i][j] = A[i][j];
            }
        }
    }
    delete[] A[0];
    delete[] A;
    *n = max;
    *m = max;
    return C;
}

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

void LUP_dec(double** A, double** L, double** U, int* pi, int n, int m){
    double eps = 1e-12;
    for(int i = 0; i < n; ++i){
        pi[i] = i;
    }
    int k1 = 0;
    for(int k = 0; k < n; ++k){
        U[k][k] = A[k][k];
        // вычисление максимального элемента по столбцу
        double max = 0.0;
        for(int i = k; i < n; ++i){
            if((std::abs(A[i][k])) > max){
                max = (std::abs(A[i][k]));
                k1 = i;
            }
        }
        // обработка случая нулевого столбца eps >= машинной точности
        if(max < eps){
            throw "error: singular matrix";
        }
        std::swap(pi[k], pi[k1]);
        for(int i = 0; i < n; ++i){
            std::swap(A[k][i], A[k1][i]);
            std::swap(L[k][i], L[k1][i]);
            for(int p = k; p < n; ++p){
                L[p][k] = A[p][k]/U[k][k];
                U[k][p] = A[k][p];
            }
        }
        for(int i = k; i < n; ++i){
            for(int j = k; j < n; ++j){
                A[i][j] = A[i][j] - L[i][k]*U[k][j];
            }
        }
        std::cout << "A: " << "\n";
        show(A, n, n);
    }
}

double* LUP_solve(double** L, double** U, int* p, double* b, int n){
    double* x = new double [n];
    double* y = new double [n];
    memset(y, 0, sizeof(double)*n);
    memset(x, 0, sizeof(double)*n);
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
    int n = 0;
    int m = 0;
    std::cout << "input n & m" << "\n";
    std::cin >> n >> m;
    std::cout << "input matrix" << "\n";
    double** A = new_matrix(n,m,1);
    double** L = nullptr;
    double** U = nullptr;
    if(n!=m){
        A = matrix_realloc(A,&n,&m);
    }
    L = new_matrix(n,m);
    U = new_matrix(n,m);
    int* pi = new int[n];
    std::cout << "input b size of " << n << "\n";
    double* b = new double[n];
    for(int i = 0; i < n; ++i){
        std::cin >> b[i];
    }
    try{
        LUP_dec(A,L,U,pi,n,m);
    }catch(const char* e){
        std::cout<<e<<"\n";
        return 0;
    }

    std::cout << "L: " << "\n";
    show(L, n, n);
    std::cout << "U: " << "\n";
    show(U, n, n);
    std::cout << "rang = " << rang_U(U,n) << "\n";
    for(int i = 0; i < n; ++i){
        std::cout << pi[i] << " ";
    }
    double* solve = LUP_solve(L,U,pi,b,n);
    for(int i = 0; i < n; ++i){
        std::cout << "x" << i+1 << " = "<< solve[i] << " ";
    }
    del_matrix(A);
    del_matrix(L);
    del_matrix(U);
    delete[] pi;
    delete[] solve;
    delete[] b;
    return 0;
}

