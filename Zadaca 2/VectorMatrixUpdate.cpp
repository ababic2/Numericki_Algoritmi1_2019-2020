#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <climits>
#include <iomanip>

class Vector {
    std::vector<double>vektor;
    public:
        explicit Vector(int n) {
            if(n <= 0) throw std::range_error("Bad dimension");
            vektor.resize(n);
            std::fill(vektor.begin(), vektor.end(), 0);
        }
        Vector(std::initializer_list<double> l) {
            if(l.size() == 0) throw std::range_error("Bad dimension");
            vektor.resize(l.size());
            std::copy(l.begin(), l.end(), vektor.begin());
        }
        int NElems() const {
            return vektor.size();
        }
        double operator[] (int i) const {
            return vektor[i];
        }
        double &operator[] (int i) {
            return vektor[i];
        }
        double operator() (int i) const {
            if(i < 1 || i > vektor.size()) throw std::range_error("Invalid index");
            return vektor.at(i - 1);
        }
        double &operator() (int i) {
            if(i < 1 || i > vektor.size()) throw std::range_error("Invalid index");
            return vektor.at(i - 1);
        }
        double Norm() const {
            double sum(0);
            for(int i = 0; i < vektor.size(); i++) sum += vektor.at(i) * vektor.at(i);
            return sqrt(sum);
        }
        friend double VectorNorm(const Vector &v) {
            return v.Norm();                                                    //sakrivena tiii metoda amina
        }
        double GetEpsilon() const {
            return 10 * Norm() * std::numeric_limits<double>::epsilon();
        }
        void Print(char separator = '\n', double eps = -1) const {
            double prag;
            if(eps < 0) prag = GetEpsilon();
            else prag = eps;
            for(int i = 0; i < vektor.size(); i++)
            {
                if(i != vektor.size() - 1) {
                    if(fabs(vektor.at(i)) < prag) std::cout << '0' << separator;
                    else std::cout << vektor.at(i) << separator;
                }
                else if(separator == '\n') std::cout << vektor.at(i) << separator;
                else std::cout << vektor.at(i);
            }
        }
        friend void PrintVector(const Vector &v, char separator = '\n', double eps = -1) {
            v.Print(separator, eps);
        }
        friend Vector operator +( const Vector &v1, const Vector & v2);
        Vector &operator +=(const Vector &v);
        friend Vector operator -( const Vector &v1, const Vector &v2);
        Vector &operator -=(const Vector &v);
        friend Vector operator *(double s, const Vector &v);
        friend Vector operator *(const Vector &v, double s);
        Vector &operator *=(double s);
        friend double operator *(const Vector &v1, const Vector &v2);
        friend Vector operator /(const Vector &v, double s);
        Vector &operator /=(double s);

        void Chop(double eps = -1) {
            if(eps < 0) eps = GetEpsilon();
            for(int i = 0; i < NElems(); i++)
                if(std::fabs(vektor[i]) < eps) vektor[i] = 0;
        }
        bool EqualTo(const Vector &v, double eps = -1) const;
};


inline bool Vector::EqualTo(const Vector &v, double eps) const {
    if(NElems() != v.NElems()) return false;
    if(eps < 0) eps = GetEpsilon();

    for(int i = 0; i < NElems(); i++)
        if(vektor[i] - v[i] > eps) return false;
    return true;
}

double operator *(const Vector &v1, const Vector &v2) {
     if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
     double skalarni_proizvod = 0;
     for(int i = 0; i < v1.NElems(); i++) skalarni_proizvod +=  v1[i] * v2[i];
     return skalarni_proizvod;
}

inline Vector& Vector::operator /=(double s) {
    if(s < std::numeric_limits<double>::epsilon() ) throw std::domain_error("Division by zero");
    for(int i = 0; i < NElems(); i++) vektor[i] /= s;
    return *this;
}

inline Vector& Vector::operator *=(double s) {
    for(int i = 0; i < NElems(); i++) vektor[i] *= s;
    return *this;
}

inline Vector& Vector::operator +=(const Vector &v) {
    if(NElems() != v.NElems()) throw std::domain_error("Incompatible formats");
    for(int i = 0; i < v.NElems(); i++) vektor[i] += v[i];
    return *this;
}

inline Vector& Vector::operator -=(const Vector &v) {
    if(NElems() != v.NElems()) throw std::domain_error("Incompatible formats");
    for(int i = 0; i < v.NElems(); i++) vektor[i] -= v[i];
    return *this;
}

Vector operator /(const Vector &v, double s) {
    if(s < std::numeric_limits<double>::epsilon()) throw std::domain_error("Division by zero");
    Vector v_n(v.NElems());
    for(int i = 0; i < v.NElems(); i++) v_n[i] = v[i] / s;
    return v_n;
}

Vector operator *(double s, const Vector &v) {
    Vector v_n(v.NElems());
    for(int i = 0; i < v.NElems(); i++) v_n[i] = v[i] * s;
    return v_n;
}

Vector operator *(const Vector &v, double s) {
    Vector v_n(v.NElems());
    for(int i = 0; i < v.NElems(); i++) v_n[i] = v[i] * s;
    return v_n;
}

Vector operator -( const Vector &v1, const Vector & v2) {
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i = 0; i < v1.NElems(); i++) v3[i] = v1[i] - v2[i];
    return v3;
}

Vector operator +( const Vector &v1, const Vector & v2) {
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector v3(v1.NElems());
    for(int i = 0; i < v1.NElems(); i++) v3[i] = v1[i] + v2[i];
    return v3;
}

class Matrix {
    std::vector<std::vector<double>> matrica;
    Matrix VektorKolona(Vector v) {      // dodana moja metoda, vidi teba li biti private
            Matrix mat(v.NElems(), 1);
            for(int i = 0; i < v.NElems(); i++)
                mat[i][0] = v[i];
            return mat;
        }


    public:
        Matrix(int m, int n) {
            if( m <= 0 || n <= 0) throw std::range_error("Bad dimension");
            matrica.resize(m);
            for(int i = 0; i < m; i++) {
                matrica[i].resize(n);
                std::fill(matrica[i].begin(), matrica[i].end(), 0);
            }
        }

        Matrix(const Vector &v) {
            matrica.resize(1);
            matrica[0].resize(v.NElems());
            for(int i = 0; i < v.NElems(); i++)
                matrica[0][i] = v[i];

        }

        Matrix(std::initializer_list< std::vector<double> >l) {
            if(l.size() == 0) throw std::range_error("Bad dimension");
            matrica.resize(l.size());
            int i = 0;
            int velicina = 0;
            for(const auto &x : l)  {
                if(i == 0) velicina = x.size();
                if(velicina != x.size() || velicina == 0)
                    throw std::logic_error("Bad matrix");
                matrica[i++] = x;
            }
        }
        int NRows() const {
            return matrica.size();
        }
        int NCols() const {
            return matrica[0].size();
        }

        double *operator[] (int i) {

            return &matrica[i][0];
        }
        const double *operator[] (int i) const {            //vraca pokazivac na red matrice Aminaaaaaa!!
             return &matrica[i][0];
        }
        double &operator() (int i, int j) {
            if(i < 1 || j < 1 || i > matrica.size() || j > matrica[0].size())
                throw std::range_error("Invalid index");
            return matrica[i - 1][j - 1];
        }
        double operator ()(int i, int j) const {
            if(i < 1 || j < 1 || i > matrica.size() || j > matrica[0].size())
                throw std::range_error("Invalid index");
            return matrica[i - 1][j - 1];
        }
        double Norm() const;
        friend double MatrixNorm(const Matrix &m);
        double GetEpsilon() const {
            return 10 * Norm() * std::numeric_limits<double>::epsilon();
        }
        void Print(int width = 2, double eps = -1) const;
        friend void PrintMatrix(const Matrix &m, int width = 2, double eps = -1) {
            m.Print(width, eps);
        }
        friend Matrix operator +(const Matrix &m1, const Matrix &m2);
        Matrix &operator +=(const Matrix &m);
        friend Matrix operator -(const Matrix &m1, const Matrix &m2);
        Matrix &operator -=(const Matrix &m);
        friend Matrix operator *(double s, const Matrix &m);
        friend Matrix operator *(const Matrix &m, double s);
        Matrix &operator *=(double s);
        friend Matrix operator *(const Matrix &m1, const Matrix &m2);
        Matrix &operator *=(const Matrix &m);
        friend Vector operator *(const Matrix &m, const Vector &v);
        friend Matrix Transpose(const Matrix &m);
        void Transpose();

        //new added

        void Chop(double eps = -1) {
            if(eps < 0) eps = GetEpsilon();
            for(int i = 0; i < NRows(); i++)
                for(int j = 0; j < NCols(); j++)
                if(std::fabs(matrica[i][j]) < eps) matrica[i][j] = 0;
        }
        bool EqualTo(const Matrix &m, double eps = -1) const;
        friend Matrix LeftDiv(Matrix m1, Matrix m2);

        friend Vector LeftDiv(Matrix m, Vector v);

        friend Matrix operator/(const Matrix& m, double s);
        Matrix &operator/=(double s);

        friend Matrix operator/(Matrix m1, Matrix m2);
        Matrix &operator/=(Matrix m);

        double Det() const;
        friend double Det(Matrix m);

        void Invert();
        friend Matrix Inverse(Matrix m);

        friend Matrix RREF(Matrix m);
        void ReduceToRREF();

        int Rank() const;
        friend int Rank(Matrix m);
};

int Rank(Matrix m) {
    return m.Rank();
}

int Matrix::Rank() const {
    Matrix X(NRows(), NCols());
    X += *this;
    int k = -1;
    int l =  -1;
    int n = X.NCols();
    int m = X.NRows();
    std::vector<bool> w;
    w.resize(n);
    int p;

    for(int j = 0; j < n; j++) {
        w.at(j) = false;
    }

        while( k < m && l < n)
        {
            l++;
            k++;
            double v = 0;

            while( v < X.GetEpsilon() && l < n) {
                p = k;
                for(int i = k; i < m; i++) {
                    if(std::fabs(X.matrica[i][l]) > v) {
                        v = std::fabs(X.matrica[i][l]);
                        p = i;
                    }
                }
                if(v < X.GetEpsilon()) l++;
            }

            if(l < n) {
                w[l] = true;
                if(p != k) { // razmijeni k ti i pti red matrice
                    X.matrica[k].swap(X.matrica[p]);
                }

                double mi = X.matrica[k][l];
                for(int j = l; j < n; j++)
                    X.matrica[k][j] /= mi;
                for(int i = 0; i < m; i++) {
                    if(i != k) {
                        mi = X.matrica[i][l];
                        for(int j = l; j < n; j++)
                            X.matrica[i][j] -= mi * X.matrica[k][j];
                    }
                }
            }
        }
    return k;
}

Matrix RREF(Matrix m) {
    Matrix result(m.NRows(), m.NCols());
    result += m;
    result.ReduceToRREF();
    return result;
}

void Matrix::ReduceToRREF() {
    int k = -1;
    int l =  -1;
    int n = NCols();
    int m = NRows();
    std::vector<bool> w;
    w.resize(n);
    int p;

    for(int j = 0; j < n; j++) {
        w.at(j) = false;
    }

        while( k < m && l < n)
        {
            l++;
            k++;
            double v = 0;

            while( v < GetEpsilon() && l < n) {
                p = k;
                for(int i = k; i < m; i++) {
                    if(std::fabs(matrica[i][l]) > v) {
                        v = std::fabs(matrica[i][l]);
                        p = i;
                    }
                }
                if(v < GetEpsilon()) l++;
            }

            if(l < n) {
                w[l] = true;
                if(p != k) { // razmijeni k ti i pti red matrice
                    matrica[k].swap(matrica[p]);
                }

                double mi = matrica[k][l];
                for(int j = l; j < n; j++)
                    matrica[k][j] /= mi;
                for(int i = 0; i < m; i++) {
                    if(i != k) {
                        mi = matrica[i][l];
                        for(int j = l; j < n; j++)
                            matrica[i][j] -= mi * matrica[k][j];
                    }
                }
            }
        }
}

Matrix operator/(Matrix m1, Matrix m2) {
    if(m2.NCols() != m2.NRows()) throw std::domain_error("Divisor matrix is not square");
    else if(m1.NCols() != m2.NCols()) throw std::domain_error("Incompatible formats");
    else if(std::fabs(m2.Det()) < std::numeric_limits<double>::epsilon() ) throw std::domain_error("Divisor matrix is singular");

    Matrix X (m1.NRows(), m1.NCols());

    int n = m2.NRows();
    int m = m1.NRows();

    for(int k = 0; k < n; k++) {

        //pivotizacija
        int p = k;
        for(int i = k + 1; i < n; i++)
            if(std::fabs(m2[k][i]) > std::fabs(m2[k][p])) p = i;

        if(p != k) {
            for (int i = 0; i < m2.NRows(); i++) {
                double t = m2[i][p];
                m2[i][p] = m2[i][k];
                m2[i][k] = t;
            }
        }

        // svodjenje na trougaonu matrica m1 i paralelno m2
        double mi;
        for(int i = k + 1; i < n; i++) {
            mi = m2[k][i] / m2[k][k];
            for(int j = k + 1; j < n; j++)
                m2[j][i] -= mi * m2[j][k];

            for(int j = 0; j < m; j++)
                m1[j][i] -= mi * m1[j][k];
        }
    }

     // supstitucija unazad
    double s;
    for(int k = 0; k < m; k++) {
        for(int i = n - 1; i >= 0; i--) {
            s = m1[k][i];
            for(int j = i + 1; j < n; j++)
                s -= m2[j][i] * X[k][j];
            X[k][i] = s / m2[i][i];
        }
    }
    return X;
}

Matrix & Matrix::operator/=(Matrix m) {
    *this = *this / m;
    return *this;
}

Matrix Inverse(Matrix m) {
   Matrix kopija(m.NRows(), m.NCols());
   kopija += m;
   kopija.Invert();
   return kopija;
}

void Matrix::Invert() {

    if(NCols() != NRows()) throw std::domain_error("Matrix is not square");
    if(std::fabs(Det()) < std::numeric_limits<double>::epsilon()) throw std::domain_error("Matrix is singular");

    int n = NRows();
    double mi;
    for(int k = 0; k < n; k++) {
        mi = matrica[k][k];

        matrica[k][k] = 1;
        for(int j = 0; j < n; j++)
            matrica[k][j] /= mi;
        for(int i = 0; i < n; i++) {
            if(i != k) {
                mi = matrica[i][k];
                matrica[i][k] = 0;
                for(int j = 0; j < n; j++)
                    matrica[i][j] -= mi * matrica[k][j];
            }
        }
    }
}

double Det(Matrix m) {
    return m.Det();
}

double Matrix::Det() const{
    if(NCols() != NRows()) throw std::domain_error("Matrix is not square");
    double d =1;
    Matrix X(NRows(), NCols());
    X = *this;
    double p;
    for(int k = 0; k < NRows(); k++) {
        p = k;
        for(int i = k + 1; i < NRows(); i++)
            if(fabs(X.matrica[i][k]) > fabs(X.matrica[p][k])) p = i;

        if(fabs(X.matrica[p][k]) < X.GetEpsilon()) return 0;
        if(std::fabs(p - k) > std::numeric_limits<double>::epsilon()) {
            d = (-1) * d;
            X.matrica[p].swap(X.matrica[k]);
        }
        d *= X.matrica[k][k];

        //sada nuliraj dole ove i predji dalje

        for(int i = k + 1; i < NRows(); i++) {
            double mi = X.matrica[i][k] / X.matrica[k][k];
            for(int j = k + 1; j < NRows(); j++)
                X.matrica[i][j] -= mi * X.matrica[k][j];
        }
    }
    return d;
}

Matrix operator/(const Matrix& m, double s) {
    Matrix result(m.NRows(), m.NCols());
    if(s < std::numeric_limits<double>::epsilon() ) throw std::domain_error("Division by zero");
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        result[i][j] = m[i][j] / s;
    return result;
}

Matrix & Matrix::operator/=(double s) {
    if(s < std::numeric_limits<double>::epsilon() ) throw std::domain_error("Division by zero");
     for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
        matrica[i][j] /= s;
    return *this;
}

Vector LeftDiv(Matrix m, Vector v) {
    Matrix m2 (v.NElems(),1);
    m2 = m2.VektorKolona(v);
    m2 = LeftDiv(m,m2);
    Vector result(v.NElems());
    for(int i = 0; i < v.NElems(); i++)
        result[i] = m2[i][0];
    return result;
}


Matrix LeftDiv(Matrix m1, Matrix m2) {
    if(m1.NCols() != m1.NRows()) throw std::domain_error("Divisor matrix is not square");
    else if(m1.NCols() != m2.NRows()) throw std::domain_error("Incompatible formats");

    Matrix X (m2.NRows(), m2.NCols());
    int n = m1.NRows();
    int m = m2.NCols();
    int p; // za razmjenu kolona u slucaju pivotizacije
    for(int k = 0; k < n; k++) {

        //pivotizacija
        p = k;
        for(int i = k + 1; i < n; i++)
            if(std::fabs(m1[i][k]) > std::fabs(m1[p][k])) p = i;
        if(std::fabs(m1[p][k]) < m1.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");

        if(p != k) {// razmjeni redove matrice m1 i m2 // i guess
            m1.matrica[p].swap(m1.matrica[k]);
            m2.matrica[p].swap(m2.matrica[k]);
        }

        // svodjenje na trougaonu matrica m1 i paralelno m2
        double mi;
        for(int i = k + 1; i < n; i++) {
            mi = m1[i][k] / m1[k][k];
            for(int j = k + 1; j < n; j++)
                m1[i][j] -= mi * m1[k][j];

            for(int j = 0; j < m; j++)
                m2[i][j] -= mi * m2[k][j];
        }

    }

     // supstitucija unazad
    double s;
    for(int k = 0; k < m; k++) {
        for(int i = n - 1; i >= 0; i--) {
            s = m2[i][k];
            for(int j = i + 1; j < n; j++)
                s -= m1[i][j] * X[j][k];
            X[i][k] = s / m1[i][i];
        }
    }
    return X;
}


bool Matrix:: EqualTo(const Matrix &m, double eps) const {
    if(NCols() != m.NCols() || NRows() != m.NRows()) return false;
    if(eps < 0) eps = GetEpsilon();
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
            if(matrica[i][j] - m[i][j] > eps) return false;
    return true;
}

void Matrix::Transpose() {

   if(NCols() == NRows()) {
        for(int i = 0; i < NRows(); i++)
            for(int j = i + 1; j < NCols(); j++)
               std::swap(matrica[i][j], matrica[j][i]);
   }
    else {
        Matrix pomocna(NCols(), NRows());
        for(int i = 0; i < NRows(); i++)
            for(int j = 0; j < NCols(); j++)
                pomocna[j][i] = matrica[i][j];
        *this = pomocna;
    }
}

Matrix Transpose(const Matrix &m) {
    Matrix trans_matrix(m.NCols(), m.NRows());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
            trans_matrix[j][i] = m[i][j];
    return trans_matrix;
}

Vector operator *(const Matrix &m, const Vector &v) {
    if(m.NCols() != v.NElems()) throw std::domain_error("Incompatible formats");
    Vector vec (m.NRows());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < v.NElems(); j++)
            vec[i] += m[i][j] * v[j];
    return vec;
}

inline Matrix& Matrix::operator *=(const Matrix &m) {
    *this = *this * m;
    return *this;
}

Matrix operator *(const Matrix &m1, const Matrix &m2) {
    if(m1.NCols() != m2.NRows()) throw std::domain_error("Incompatible formats");
    Matrix proizvod_matrica (m1.NRows(), m2.NCols());
    for(int i = 0; i < m1.NRows(); i++)
        for(int j = 0; j < m2.NCols(); j++)
            for(int k = 0; k < m1.NCols(); k++)
                proizvod_matrica[i][j] += m1[i][k] * m2[k][j];
    return proizvod_matrica;
}

void Matrix::Print(int width, double eps) const {
    double prag = eps;
    if(eps < 0) prag = GetEpsilon();
    for(int i = 0; i < NRows(); i++) {
        for(int j = 0; j < NCols(); j++) {
            if(fabs(matrica[i][j]) < prag) std::cout << std::setw(width) << '0';
            else if(matrica[i][j] < 0) std::cout << std::setw(width + 1) << matrica[i][j];
            else std::cout << std::setw(width) << matrica[i][j];
        }
    std::cout << std::endl;
    }
}

double Matrix::Norm() const {
    double sumOfSquares(0);
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
        sumOfSquares += matrica[i][j] * matrica[i][j];
    //return sqrt(sumOfSquares);
    sumOfSquares = sqrtl(sumOfSquares);
    return sumOfSquares;
}

double MatrixNorm(const Matrix &m) {
    return m.Norm();
}

Matrix& Matrix::operator *=(double s) {
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
        matrica[i][j] *= s;
    return *this;
}

Matrix& Matrix::operator +=(const Matrix &m) {
    if((NRows() != m.NRows() ) || (NCols() != m.NCols())) throw std::domain_error("Incompatible formats");
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        matrica[i][j] += m[i][j];
    return *this;
}

Matrix& Matrix::operator -=(const Matrix &m) {
    if((NRows() != m.NRows() ) || (NCols() != m.NCols())) throw std::domain_error("Incompatible formats");
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        matrica[i][j] -= m[i][j];
    return *this;
}

Matrix operator *(double s, const Matrix &m) {
    Matrix m3(m.NRows(), m.NCols());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        m3[i][j] = m[i][j] * s;
    return m3;
}

Matrix operator *(const Matrix &m, double s) {
    Matrix m3(m.NRows(), m.NCols());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        m3[i][j] = m[i][j] * s;
    return m3;
}

Matrix operator +(const Matrix &m1, const Matrix &m2) {
    if((m1.NRows() != m2.NRows() ) || (m1.NCols() != m2.NCols())) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(), m1.NCols());
    for(int i = 0; i < m1.NRows(); i++)
        for(int j = 0; j < m1.NCols(); j++)
        m3[i][j] = m1[i][j] + m2[i][j];
    return m3;
}

Matrix operator -(const Matrix &m1, const Matrix &m2) {
    if((m1.NRows() != m2.NRows() ) || (m1.NCols() != m2.NCols())) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(), m1.NCols());
    for(int i = 0; i < m1.NRows(); i++)
        for(int j = 0; j < m1.NCols(); j++)
        m3[i][j] = m1[i][j] - m2[i][j];
    return m3;
}


class LUDecomposer {
    Matrix matrix;
    Vector w; //cuva izmjene redova
    void Solve_vectors(std::vector<double> &b, std::vector<double> &x);
    public:
        LUDecomposer(Matrix m) : matrix(m), w(m.NCols()){
            if(matrix.NCols() != matrix.NRows()) throw std::domain_error("Matrix is not square");

            int n = matrix.NRows();


            for(int j = 0; j < n; j++) {

                for(int i = 0; i <= j; i++) {

                    double s = matrix[i][j];
                    for(int k = 0; k < i; k++)
                        s -= matrix[i][k] * matrix[k][j];

                    matrix[i][j] = s;
                }

                int p = j;
                for(int i = j + 1; i < n; i++) {
                    double s = matrix[i][j];
                    for(int k = 0; k <=  j - 1; k++)
                        s -= matrix[i][k] * matrix[k][j];

                    matrix[i][j] = s;
                    if(std::fabs(s) > std::fabs(matrix[p][j]))
                        p = i;
                }
                if(std::fabs(matrix[p][j]) < matrix.GetEpsilon())
                    throw std::domain_error("Matrix is singular");
                if(p != j) {
                    for(int r = 0; r < n; r++) {
                        double temp = matrix[p - 1][r];
                        matrix[p - 1][r] = matrix[j][r];
                        matrix[j][r] = temp;
                    }
                }

                w[j] = p;
                double mi = matrix[j][j];
                for(int i = j + 1; i < n; i++) {
                    matrix[i][j] /= mi;
                }
            }

        }

        void Solve(const Vector &b, Vector &x) const;
        Vector Solve(Vector b)const;
        void Solve(Matrix &b, Matrix &x) const;
        Matrix Solve(Matrix b) const;

        Matrix GetCompactLU() const {
            return matrix;
        }
        Matrix GetL() const;

        Matrix GetU() const;

        Vector GetPermuation() const {
            Vector kopija(w.NElems());
            for(int i = 0; i < w.NElems(); i++)
                kopija[i] = w[i];
            return kopija;
        }
};

Matrix LUDecomposer::Solve(Matrix b) const {
    if(b.NRows() != matrix.NRows() || b.NCols() != matrix.NCols())
        throw std::domain_error("Incompatible formats");
    Matrix x (b.NRows(), b.NCols());
    Solve(b, x);
    return x;
}

inline void LUDecomposer::Solve(Matrix &b, Matrix &x) const {
    if(b.NRows() != x.NRows() || b.NCols() != x.NCols())
        throw std::domain_error("Incompatible formats");

    for(int m = 0; m < b.NCols(); m++) {

        for(int i = 0; i < b.NRows(); i++) {
            int p = w[i];
            long double s = b[p][m];
            b[p][m] = b[i][m];

            for(int j = 0; j < i; j++)
                s -= matrix[i][j] * b[j][m];
            b[i][m] = s;
        }

        //supstitucija unazad
        for(int i = b.NRows() - 1; i >= 0; i--) {
            double s = b[i][m];
            for(int j = i + 1; j < b.NRows(); j++)
                s -= matrix[i][j] * x[j][m];
            x[i][m] = s /  matrix[i][i];
        }
    }
}

Vector LUDecomposer::Solve(Vector b) const {
   if(b.NElems() != matrix.NRows()) throw std::domain_error("Incompatible formats");

   for(int i = 0; i < b.NElems(); i++) {
        int p = w[i];
        long double s = b[p];
        b[p] = b[i];

        for(int j = 0; j < i; j++)
            s -= matrix[i][j] * b[j];
        b[i] = s;
    }

    //supstitucija unazad
    for(int i = b.NElems() - 1; i >= 0; i--) {
        double s = b[i];
        for(int j = i + 1; j < b.NElems(); j++)
            s -= matrix[i][j] * b[j];
        b[i] = s /  matrix[i][i];
    }
    return b;
}

void LUDecomposer::Solve(const Vector &b, Vector &x) const {
    if(b.NElems() != matrix.NRows()) throw std::domain_error("Incompatible formats");
    else if(b.NElems() != x.NElems()) throw std::domain_error("Incompatible formats");
    x = std::move(Solve(b));
}

inline Matrix LUDecomposer::GetU() const {
    Matrix U (matrix.NRows(), matrix.NCols());
    for(int i = 0; i < U.NRows(); i++)
        for(int j = 0; j < U.NRows(); j++)
            if(j >= i)
                U[i][j] = matrix[i][j];
    return U;
}

inline Matrix LUDecomposer::GetL() const {
    Matrix L (matrix.NRows(), matrix.NCols());
    for(int i = 0; i < L.NRows(); i++)
        for(int j = 0; j <= i; j++) {
            if(i == j) L[i][j] = 1;
            else L[i][j] = matrix[i][j];
        }
    return L;
}

class QRDecomposer {

    Matrix matrix; // R matrica se dobije gornji R elementi
    Matrix V;   //cuva vektore
    Vector diagonal;
    int rows;
    int n;

    public:

        QRDecomposer(Matrix m);
        void Solve(const Vector &b, Vector &x) const;
        Vector Solve(Vector b) const;
        void Solve(Matrix &b, Matrix &x) const;
        Matrix Solve(Matrix b) const;
        Vector MulQWith(Vector v) const;
        Matrix MulQWith(Matrix m) const;
        Vector MulQTWith(Vector v) const;
        Matrix MulQTWith(Matrix m) const;
        Matrix GetQ() const;
        Matrix GetR() const;
};

Matrix QRDecomposer::MulQTWith(Matrix m) const {
    if(matrix.NCols() != m.NRows())
        throw std::domain_error("Incompatible formats");

     for(int j = 0; j < m.NRows(); j++) {
        for(int k = 0; k < m.NRows(); k++) {
        double s = 0;
        for(int i = k; i < rows; i++)
            s += matrix[i][k] * m[i][j];

        for(int i = k; i < rows; i++)
            m[i][j] -= s * matrix[i][k];
        }
    }
    return m;
}

Matrix QRDecomposer::MulQWith(Matrix m) const {
      if(matrix.NCols() != m.NRows())
        throw std::domain_error("Incompatible formats");

     for(int j = 0; j < m.NRows(); j++) {
        for(int k = n - 1; k >= 0; k--) {
        double s = 0;
        for(int i = k; i < rows; i++)
            s += matrix[i][k] * m[i][j];

        for(int i = k; i < rows; i++)
            m[i][j] -= s * matrix[i][k];
        }
    }
    return m;
}

Matrix QRDecomposer::GetQ()const {
    Matrix Q (matrix.NRows(), matrix.NRows());

    for(int j = 0; j < rows; j++) {
        for(int i = 0; i < rows; i++)
            Q[i][j] = 0;
        Q[j][j] = 1;

        for(int k = n - 1; k >= 0; k--) {
        double s = 0;
        for(int i = k; i < rows; i++)
            s += matrix[i][k] * Q[i][j];

        for(int i = k; i < rows; i++)
            Q[i][j] -= s * matrix[i][k];
        }
    }
    return Q;
}

Matrix QRDecomposer::GetR()const {
    Matrix R(matrix.NRows(), matrix.NCols());
    for(int i = 0; i < matrix.NRows(); i++)
        for(int j = 0; j < matrix.NCols(); j++) {
            if(i == j) R[i][j] = diagonal[i];
            else if(i > j) R[i][j] = 0;
            else R[i][j] = matrix[i][j];
        }
    return R;
}

Vector QRDecomposer::MulQWith(Vector v) const {
    if(matrix.NCols() != matrix.NRows())
        throw std::domain_error("Matrix is not square");

    else if(v.NElems() != matrix.NRows())
        throw std::domain_error("Incompatible formats");

    for(int k = v.NElems() - 1; k >= 0; k--) {
        double s = 0;
        for(int i = k; i < rows; i++)
            s += matrix[i][k] * v[i];

        for(int i = k; i < rows; i++)
            v[i] -= s * matrix[i][k];
    }
    return v;
}

Vector QRDecomposer::MulQTWith(Vector v) const {
    if(matrix.NCols() != matrix.NRows())
        throw std::domain_error("Matrix is not square");

    else if(v.NElems() != matrix.NRows())
        throw std::domain_error("Incompatible formats");

    for(int k = 0; k < v.NElems(); k++) {
        double s = 0;
        for(int i = k; i < rows; i++)
            s += matrix[i][k] * v[i];

        for(int i = k; i < rows; i++)
            v[i] -= s * matrix[i][k];
    }
    return v;
}

Matrix QRDecomposer::Solve(Matrix b) const {

    if(matrix.NCols() != matrix.NRows())
        throw std::domain_error("Matrix is not square");
     if(b.NRows() != matrix.NRows() || b.NCols() != matrix.NCols())
        throw std::domain_error("Incompatible formats");

    Matrix x(b.NRows(), b.NCols());

    for(int i = 0; i < b.NCols(); i++) {
       Vector v(b.NRows());
       for(int j = 0; j < b.NRows(); j++) {
        v[j] = b[j][i];
       }

        v = Solve(v); // rijesimo kolonu
        for(int j = 0; j < b.NRows(); j++) {
            x[j][i] = v[j]; }
    }
    return x;
}


void QRDecomposer::Solve(Matrix &b, Matrix &x) const {

    if(matrix.NCols() != matrix.NRows())
        throw std::domain_error("Matrix is not square");
    if(b.NRows() != matrix.NRows() || b.NCols() != matrix.NCols())
        throw std::domain_error("Incompatible formats");
    if(b.NRows() != x.NRows() || b.NCols() != x.NCols())
        throw std::domain_error("Incompatible formats");
    x = Solve(b);
}

Vector QRDecomposer::Solve(Vector b) const {
    if(matrix.NCols() != matrix.NRows())
        throw std::domain_error("Matrix is not square");

    else if(b.NElems() != matrix.NRows())
        throw std::domain_error("Incompatible formats");

    for(int k = 0; k < b.NElems(); k++) {
        double s = 0;
        for(int i = k; i < rows; i++)
            s += matrix[i][k] * b[i];

        for(int i = k; i < rows; i++)
            b[i] -= s * matrix[i][k];

    }

    //suspstitucija unazad
    for(int i = b.NElems() - 1; i >= 0; i--) {
        double s = b[i];
        for(int j = i + 1; j < b.NElems(); j++)
            s -= matrix[i][j] * b[j];

        b[i] = s / diagonal[i];
    }
    return b;
}

void QRDecomposer::Solve(const Vector &b, Vector &x) const {
    if(b.NElems() != matrix.NRows()) throw std::domain_error("Incompatible formats");
    else if(b.NElems() != x.NElems()) throw std::domain_error("Incompatible formats");
    x = std::move(Solve(b));
}

QRDecomposer::QRDecomposer(Matrix m) : matrix(m), V(m.NRows(),m.NCols()), diagonal(matrix.NRows()), rows(m.NRows()), n(m.NCols()) {

    if(rows < n) throw std::domain_error("Invalid matrix format");
    //n jee cols
    for(int k = 0; k < n; k++) {
        long double s = 0;
        for(int i = k; i < rows; i++)
            s += matrix[i][k] * matrix[i][k];

        s = std::sqrt(s);
        long double mi = std::sqrt(s * (s + std::fabs(matrix[k][k])));

        if(std::fabs(mi) < matrix.GetEpsilon())
            throw std::domain_error("Matrix is singular");

        if(matrix[k][k] < 0) s = -s;
        matrix[k][k] = (matrix[k][k] + s) / mi;

        for(int i = k + 1; i < n; i++)
            matrix[i][k] = matrix[i][k] / mi;

        diagonal[k] = -s;

        for(int j = k + 1; j < n; j++) {
            s = 0;
            for(int i = k; i < rows; i++)
                s += matrix[i][k] * matrix[i][j];

            for(int i = k; i < rows; i++)
                matrix[i][j] -= s * matrix[i][k];
        }
    }
}



int main ()
{
    {
        std::cout << "Vector EqualTo: ";

        Vector v1{1.1, 2.05, 3.7, 1.5};
        if (!v1.EqualTo({1, 2, 3.75, 1.5}, 0.11)) {
            std::cout << "NOT OKAY EQUALTO" << std::endl;
            return false;
        }

        if (v1.EqualTo({1, 2, 3.75, 1.5})) {
            std::cout << "NOT OKAY EQUALTO" << std::endl;
            return false;
        }

        std::cout << "OK" << std::endl;

           std::cout << "Vector Chop: ";

        Vector v{0.1, 2, 0.3, 0.4, 0.2, -1, -0.1};
        v.Chop(0.25);

        if (!v.EqualTo({0, 2, 0.3, 0.4, 0, -1, 0})) {
            std::cout << "NOT OKAY CHOP" << std::endl;
            return false;
        }

        {
            Vector v1{2, 5, 7, 2, 3, 0.1};
            v1.Chop(5);
            if (!v1.EqualTo({0, 5, 7, 0, 0, 0})) {
                std::cout << "NOT OKAY CHOP" << std::endl;
            }

            Vector v2{2, 5, 7, 2, 3, 0.1};
            v2.Chop(5);
            if (v2.EqualTo({0, 0, 7, 0, 0, 0})) {
                std::cout << "NOT OKAY CHOP" << std::endl;
            }

            std::cout << "OK" << std::endl;
        }

    }

    {
        std::cout<<"\n\n<<<<<<<PRVA RUNDA TESTOVA ZA LU>>>>>>>\n\n";
        Matrix A{{0,2,3},{9,8,7},{13,21,64}};
        LUDecomposer solver(A);

        try{
            LUDecomposer(Matrix{{1,2}});
            std::cout<<"Ne hvata nekvadratne"<<std::endl;
        }catch(std::domain_error e){
            std::cout << "LU konstruktor hvata: "<<e.what()<<std::endl;
        }

        try{
            LUDecomposer(Matrix{{1,1},{2,2}});
            std::cout<<"Ne hvata singularne"<<std::endl;
        }catch(std::domain_error e){
            std::cout << "LU konstruktor hvata: "<<e.what()<<std::endl;
        }

        try{
            LUDecomposer(Matrix{{1,1},{2,3}}).Solve(Vector({1,2,3,4,5,6}));
            std::cout<<"Ne hvata los format"<<std::endl;
        }catch(std::domain_error e){
            std::cout << "LU Solve hvata: "<<e.what()<<std::endl;
        }

        try{
            LUDecomposer(Matrix{{1,1},{2,3}}).Solve(Matrix({{1,2,3},{4,5,6},{7,8,9}}));
            std::cout<<"Ne hvata los format"<<std::endl;
        }catch(std::domain_error e){
            std::cout << "LU Solve hvata: "<<e.what()<<std::endl;
        }

    }

    //TEST LU2
    {
        std::cout<<"\n\n<<<<<<<DRUGA RUNDA TESTOVA ZA LU>>>>>>>\n\n";
        Matrix A{{2340,215452,514513},{1359,153458,5147},{11323,5421,143564}};
        LUDecomposer solver(A);

        if((solver.GetL()+solver.GetU() - Matrix({{1,0,0},{0,1,0},{0,0,1}})).EqualTo(solver.GetCompactLU()))
        std::cout<<"CompactLU test USPJESAN\n";
        else std::cout<<"CompactLU test NIJE USPJESAN\n";

        Vector w = solver.GetPermuation();
        if(std::abs(w[0]-1)<5*std::numeric_limits<double>::epsilon() || std::abs(A[0][0]) > A.GetEpsilon()) std::cout<<"GetPermutation test USPJESAN\n";
        else std::cout<<"GetPermutation test NIJE USPJESAN\n";

        {
            Matrix A{{0,7,2},{4,6,5},{2,1,7}};
            Matrix x{{4,2,3},{2,3,1},{8,11,9}};
            Matrix rez=A*x;
            LUDecomposer lu(A);
            lu.Solve(rez,rez);

            if(rez.EqualTo(lu.GetL()*lu.GetU())) std::cout<<"LU test GetL() i GetU() USPJESAN\n";
             else std::cout<<"LU test NIJE USPJESAN\n";
        }

        {
            Matrix A{{0,7,2},{4,6,5},{2,1,7}};
            Matrix x{{4,2,3},{2,3,1},{8,11,9}};
            Matrix rez=A*x;
            LUDecomposer lu(A);
            lu.Solve(rez);
            if(x.EqualTo(rez)) std::cout<<"Test1 Solve Matrix uspjesan\n";
            else std::cout<<"Test Solve Matrix neuspjesan ";
        }

        {
            Matrix A{{0,7,2},{4,6,5},{2,1,7}};
            Matrix x{{4,1,5},{2,3,1},{8,7,9}};
            Matrix rez=A*x;
            LUDecomposer lu(A);
            lu.Solve(rez,rez);
            if(x.EqualTo(rez)) std::cout<<"Test Solve Matrix uspjesan\n";
            else std::cout<<"Test Solve Matrix neuspjesan ";
        }

        {
            Matrix A{{0,5,2},{4,6,2},{3,1,7}};
            Vector x{1,2,4};
            Vector rez=A*x;
            LUDecomposer lu(A);
            if(x.EqualTo(lu.Solve(rez))) std::cout<<"Test SOlve vector uspjesan\n"<<std::endl;
            else std::cout << "Test SOlve vector neuspjesan";
        }
    }

    //QR IZUZECI
    std::cout << "\n<<<< QR IZUZECI >>>>\n" << std::endl;
    try {
        Matrix Test{ {1,2,3,4,5}, {1,2,3,4,5}, {1,2,3,4,5}, {1,2,3,4,5} };  // test konstruktora format
        QRDecomposer q(Test);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {
        Matrix Test{ {1,1,1}, {0,1,0}, {1,0,1} };                       // test konstruktora singular
        QRDecomposer q(Test);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                               // test solve(vektor, vektor) 1. slucaj
        Matrix Test{ {12,-51}, {6,167}, {7,8} };
        QRDecomposer q(Test);
        Vector b{1,2,3,4};
        Vector x{1,2,3,4};
        q.Solve(b, x);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                           // test solve(vektor, vektor) 2. slucaj
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9} };
        QRDecomposer q(Test);
        Vector b{1,2,3};
        Vector x{1,2,3,4};
        q.Solve(b, x);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                         // test solve(vektor) 1.slucaj
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9}, {14,58,-33} };
        QRDecomposer q(Test);
        Vector b{1,2,3,4};
        Vector x = q.Solve(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                         // test solve(vektor) 2.slucaj
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9} };
        QRDecomposer q(Test);
        Vector b{1,2,3,4};
        Vector x = q.Solve(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }


    try {                                                         // test solve(matrix, matrix) 1.slucaj
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9}, {14,58,-33} };
        QRDecomposer q(Test);
        Matrix b { {-3,5,7}, {43,68,-21}, {13,36,23} };
        Matrix x { {6,3,5}, {-8,21,45}, {31,4,-12} };
        q.Solve(b, x);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                         // test solve(matrix, matrix) 2.slucaj
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9} };
        QRDecomposer q(Test);
        Matrix b { {-3,5,7}, {43,68,-21} };
        Matrix x { {6,3,5}, {-8,21,45}, {31,4,-12} };
        q.Solve(b, x);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                         // test MulQWith(vektor)
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9}, {14,58,-33} };
        QRDecomposer q(Test);
        Vector b {-3,5,7,21};
        Vector x = q.MulQWith(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                         // test MulQWith 2.slucaj
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9} };
        QRDecomposer q(Test);
        Vector b {-3,5,7,21};
        Vector x = q.MulQWith(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }


    try {                                                         // test MulQWith(matrix)
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9}, {14,58,-33} };
        QRDecomposer q(Test);
        Matrix b { {-3,5,7}, {43,68,-21} };
        Matrix x = q.MulQWith(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                         // test MulQTWith(vektor)
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9}, {14,58,-33} };
        QRDecomposer q(Test);
        Vector b {-3,5,7,21};
        Vector x = q.MulQTWith(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }

    try {                                                         // test MulQTWith 2.slucaj
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9} };
        QRDecomposer q(Test);
        Vector b {-3,5,7,21};
        Vector x = q.MulQTWith(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }


    try {                                                         // test MulQTWith(matrix)
        Matrix Test{ {12,-51,4}, {6,167,-68}, {7,8,9}, {14,58,-33} };
        QRDecomposer q(Test);
        Matrix b { {-3,5,7}, {43,68,-21} };
        Matrix x = q.MulQTWith(b);
        std::cout << " NOT OK\n";
    } catch (std::domain_error e) { std::cout << e.what() << " OK\n"; }


    std::cout << "\n<<<QR TESTOVI >>>\n" << std::endl;

    //Q i R test
    Matrix A{{2.34,3.45,4.56},{5.67,6.78,7.89},{8.91,9.12,1.23}};
    QRDecomposer qr(A);
    if ((qr.GetQ()*qr.GetR()).EqualTo(A)) std::cout<<"TEST 1 USPJESAN"<<std::endl;
    else std::cout<<"TEST 1 NIJE USPJESAN"<<std::endl;

    //Matrix solve test1
    Matrix A1{{2,2,3}, {4,5,5}, {7,7,8}};
    Matrix B1{{1,3,3}, {2,2,4}, {5,6,6}};
    QRDecomposer qr1(A1);
    Matrix x1=qr1.Solve (B1);
    if (x1.EqualTo({{2,1,-1},{-0.6,-2.2,-0.2},{-0.6,1.8,1.8}})) std::cout<<"TEST 2.Matrix solve USPJESAN"<<std::endl; //ova matrica sa kojom se poredi se dobije kao A1\B1 u scilabu
    else std::cout<<"TEST 2 NIJE USPJESAN"<<std::endl;

    //Matrix solve test2
    Matrix x11(3,3);
    qr1.Solve (B1,x11);
    if (x11.EqualTo({{2,1,-1},{-0.6,-2.2,-0.2},{-0.6,1.8,1.8}})) std::cout<<"TEST 3.Matrix Solve USPJESAN"<<std::endl;
    else std::cout<<"TEST 3 NIJE USPJESAN"<<std::endl;

    //Vector solve test1
    Vector v{1,5,1};
    Vector x2=qr1.Solve(v);
    if (x2.EqualTo({-5,4,1})) std::cout<<"TEST 4.Vector Solve USPJESAN"<<std::endl;
    else std::cout<<"TEST 4 NIJE USPJESAN"<<std::endl;

    //Vector solve test2
    Vector x22(3);
    qr1.Solve(v,x22);
    if (x22.EqualTo({-5,4,1})) std::cout<<"TEST 5.VectorSolve USPJESAN"<<std::endl;
    else std::cout<<"TEST 5 NIJE USPJESAN"<<std::endl;

    //MulQWith vector test
    if ((qr1.GetQ()*v).EqualTo(qr1.MulQWith(v))) std::cout<<"TEST 6.MulQWith vector USPJESAN"<<std::endl;
    else std::cout<<"TEST 6 NIJE USPJESAN"<<std::endl;

    //MulQTWith vector test
    if ((Transpose(qr1.GetQ())*v).EqualTo(qr1.MulQTWith(v))) std::cout<<"TEST 7.MulQTWith vector USPJESAN"<<std::endl;
    else std::cout<<"TEST 7 NIJE USPJESAN"<<std::endl;

    //MulQWith matrix test
    Matrix m{{1,1,1},{2,2,2},{3,3,3}};
    if ((qr1.GetQ()*m).EqualTo(qr1.MulQWith(m))) std::cout<<"TEST 8.MulQWith matrix USPJESAN"<<std::endl;
    else std::cout<<"TEST 8 NIJE USPJESAN"<<std::endl;

    //MulQTWith matrix test
    if ((Transpose(qr1.GetQ())*m).EqualTo(qr1.MulQTWith(m))) std::cout<<"TEST 9.MulQTWith matrix USPJESAN"<<std::endl;

    else std::cout<<"TEST 9 NIJE USPJESAN"<<std::endl;


    std::cout << "\n<<<MATRIX TESTOVI >>>\n" << std::endl;


   //TEST ZA LEFT DIVISION i operator=/, operator /
    {
        {
            //test za  /=  i / sa skalarom
            Matrix matrica{{2,4,6},{1,1,1}};
            matrica = matrica / 2;
            Matrix result{{1,2,3},{0.5,0.5,0.5}};
            if(matrica.EqualTo(result)) std::cout <<  "Test za /=(scalar) uspješan"<<std::endl;

            if(matrica.EqualTo(result)) std::cout <<  "Test za / (scalar) uspješan"<<std::endl;
            else std::cout << "nije uspjesan test";

            try {
                matrica /= 0;
            } catch(std::domain_error e) { std::cout << "Test za /=(scalar) uspješan -izuzetak-"<<e.what()<<std::endl;}

            try {
                matrica = matrica / 0;
            } catch(std::domain_error e) { std::cout << "Test za / (scalar) uspješan -izuzetak-"<<e.what()<<std::endl; }
        }

        { //singular divisior
            Matrix m2{{1,2,3},{3,4,5},{6,7,8}};
            Matrix m1{{1,1,1},{3,6,9},{2,1,2}};
            try {
                m1/=m2;
            } catch(std::domain_error e) { std::cout << "Test za operator /=(Matrix) uspjesan-izuzetak-"<<e.what()<<std::endl;}

            try {
                m1=m1 / m2;
            } catch(std::domain_error e) { std::cout << "Test za operator / (Matrix) uspjesan-izuzetak-"<<e.what()<<std::endl;}
        }

        {
            //div matrix is not square
            Matrix m2{{1,2,3},{3,4,5}};
            Matrix m1{{1,1,1},{3,6,9},{2,1,2}};
            try {
                m1/=m2;
            } catch(std::domain_error e) { std::cout << "Test za operator /=(Matrix) uspjesan-izuzetak-"<<e.what()<<std::endl;}
            try {
                m1 = m1 / m2;
            } catch(std::domain_error e) { std::cout << "Test za operator / (Matrix) uspjesan-izuzetak-"<<e.what()<<std::endl;}
        }

        {
            //div matrix is not square
            Matrix m2{{1,2,3,4},{3,4,5,6},{7,8,9,10},{2,4,6,8}};
            Matrix m1{{1,9,11},{3,6,9},{2,1,2}};
            try {
                m1/=m2;
            } catch(std::domain_error e) { std::cout << "Test za operator /=(Matrix) uspjesan-izuzetak-"<<e.what()<<std::endl;}
            try {
                m1 = m1 / m2;
            } catch(std::domain_error e) { std::cout << "Test za operator / (Matrix) uspjesan-izuzetak-"<<e.what()<<std::endl;}
        }

        {
            //izuzetak za not square
            Matrix m1{{1,2,3},{1,2,3}};
            Vector v{1,2,3};
            try {
                Matrix result = LeftDiv(m1,v);
            } catch(std::domain_error e) {
                std::cout << "Test za LeftDiv() uspješan - izuzetak-"<<e.what()<<std::endl;
            }

        }

        {
            Matrix m1{{1,2,3},{1,2,3},{1,2,3}};
            Vector v{1,2,3,4};
            try {
                Matrix result = LeftDiv(m1,v);
            } catch(std::domain_error e) {
                std::cout << "Test za LeftDiv() uspješan - izuzetak-"<<e.what()<<std::endl;
            }
        }

        {
            Matrix m1{{1,2,3},{3,2,1},{2,1,3}};
            Vector v{2,1,2};
            Vector result = LeftDiv(m1,v);
            Vector compare_to{0.0833333,0.0833333,0.583333};
            if(result.EqualTo(compare_to,1)) std::cout << "Test za LeftDiv() uspješan\n";
            else std::cout << "Test za LeftDiv() neuspješan";
        }

        {
            //izuzetak za not square
            Matrix m1{{1,2,3},{1,2,3}};
            Matrix m2 {{1,2},{2,3}};
            try {
                Matrix result = LeftDiv(m1,m2);
            } catch(std::domain_error e) {
                std::cout << "Test za LeftDiv() uspješan - izuzetak-"<<e.what()<<std::endl;
            }
        }

         {
             //izuzetak za is singular
             Matrix m1{{1,-2},{-3,6}};
             Matrix m2 {{1,2},{2,3}};
              try {
                Matrix result = LeftDiv(m1,m2);
            } catch(std::domain_error e) {
                std::cout << "Test za LeftDiv() uspješan - izuzetak-"<<e.what()<<std::endl;
            }
         }

          {
             Matrix m1{{1,2,3},{6,5,4},{7,6,7}};
             Matrix m2 {{1,0,0},{0,1,0},{0,0,1}};
             Matrix result = LeftDiv(m1,m2);

             if(result.EqualTo(m1) ) std::cout << "Test za LeftDiv uspjesan\n";
                else std::cout << "Test za LeftDiv neuspjesan";
         }

          {
              // izuzetak Inc.formats
             Matrix m1{{1,2,3},{6,5,4},{7,6,7}};
             Matrix m2 {{1,0},{1,1}};
             try {
                 Matrix result = LeftDiv(m1,m2);
             } catch(std::domain_error e) {
                  std::cout << "Test za LeftDiv uspjesan-" << e.what() << std::endl;
             }
         }
    }


    //TEST ZA DETERMINANTU

    {
        {
            Matrix n_square{{1,2,3},{4,5,6}};
            try {
                n_square.Det();
            }catch(std::domain_error e) {
                std::cout << "Test za Det() uspjesan - izuzetak: " << e.what() << std::endl;
            }
        }

        {
            Matrix square{{1,2,3},{1.1,1.2,3},{4,5,6}};
            if(square.Det() - 5.1 < 0.00000001) std::cout << "Test za Det uspjesan\n";
            else std::cout << "Test za det neuspjesan";
        }

    }

    //TEST RREF

    {
        {
            Matrix matrix{{1,2,3},{4,5,6},{7,8,9}};
            matrix.ReduceToRREF();
            Matrix m{{1,0,-1},{0,1,2},{0,0,0}};
             if(matrix.EqualTo(m)) std::cout << "Test ReduceToRREF() USPJESAN\n";
                else std::cout << "TEST ZA RREF() NIJE USPJESAN";
        }

        {
            Matrix matrix{{1,2,3},{4,5,6},{7,8,9}};
            Matrix m{{1,0,-1},{0,1,2},{0,0,0}};
             if(RREF(matrix).EqualTo(m)) std::cout << "Test friend RREF() USPJESAN\n";
                else std::cout << "TEST ZA RREF() NIJE USPJESAN";
        }
    }


    //TEST CHOP
    {
        double num =std::numeric_limits<double>::epsilon();

        Matrix matrix{{num,2,3},{4,-num,6},{11,-123,3.4}};
        Matrix ChopM{{0,2,3},{4,0,6},{11,-123,3.4}};

        if(matrix.EqualTo(ChopM)) std::cout << "Test Chop() USPJESAN\n";
            else std::cout << "TEST ZA Chop() NIJE USPJESAN";
    }

    //TEST EQUAL
    {
        Matrix matrix{{1,2,3},{4,5,6},{11,123,3.4}};
        Matrix matrix2{{1,2,3},{4,5,6},{11,123,3.4}};
        if(matrix.EqualTo(matrix2)) std::cout << "Test  EqualTo() USPJESAN\n";
            else std::cout << "TEST ZA Equal() NIJE USPJESAN";
    }

    //TEST RANK
    {
        Matrix matrix{{1,2,3},{5,1,8},{7,9,1}};
        if(Rank(matrix) == 3) std::cout << "Test  Rank(Matrix m) USPJESAN\n";
            else std::cout << "TEST ZA Rank(Matrix m) NIJE USPJESAN";

        if(matrix.Rank() == 3) std::cout << "Test  Rank() USPJESAN\n";
            else std::cout << "TEST ZA Rank() NIJE USPJESAN";
    }


 //TEST ZA INVERZNU MATRICU (funcija INVERT())- provjera rezultata + testiranje izuzetaka
    {
        {
            Matrix matrix {{3,0,2},{0,0,1},{2,-2,1}};
            Matrix inverse {{0.2,-0.2, 0.2},{0.2,0.3,-0.3},{0,1,0}};
            matrix.Invert();
            if(inverse.EqualTo(matrix)) std::cout << "Test  Invert() USPJESAN";
            else std::cout << "TEST ZA INV.MATRICU NIJE USPJESAN";
        }

         {
            Matrix matrix {{3,0,2},{0,0,1}};
            try {
                matrix.Invert();
            } catch(std::domain_error e) {
                std::cout << "\nTest za Invert() uspjesan - izuzetak : " << e.what();
            }
        }

        {
            Matrix matrix {{1,1,1},{1,1,1},{1,1,1}};
            try {
                matrix.Invert();
            } catch(std::domain_error e) {
                std::cout << "\nTest za Invert() uspjesan - izuzetak : " << e.what();
            }

        }
    }

    {
        {
            Matrix matrix {{3,0,2},{0,0,1},{2,-2,1}};
            Matrix inverse {{0.2,-0.2, 0.2},{0.2,0.3,-0.3},{0,1,0}};

            if(inverse.EqualTo(Inverse(matrix))) std::cout << "\nTEST ZA Inverse() USPJESAN";
            else std::cout << "TEST ZA INV.MATRICU NIJE USPJESAN";
        }

         {
            Matrix matrix {{3,0,2},{0,0,1}};
            try {
                Inverse(matrix);
            } catch(std::domain_error e) {
                std::cout << "\nTest za Inverse() uspjesan - izuzetak : " << e.what();
            }
        }

        {
            Matrix matrix {{1,1,1},{1,1,1},{1,1,1}};
            try {
                Inverse(matrix);
            } catch(std::domain_error e) {
                std::cout << "\nTest za Inverse() uspjesan - izuzetak : " << e.what();
            }

        }
    }
 	return 0;
}
