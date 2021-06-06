//NA 2019/2020: Zadaća 1, Zadatak 1
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
};

double operator *(const Vector &v1, const Vector &v2) {
     if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
     double skalarni_proizvod = 0;
     for(int i = 0; i < v1.NElems(); i++) skalarni_proizvod +=  v1[i] * v2[i];
     return skalarni_proizvod;
}

inline Vector& Vector::operator /=(double s) {
    if(s == 0 ) throw std::domain_error("Division by zero");
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
    if(s == 0) throw std::domain_error("Division by zero");
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
            if(l.size() == 1) throw std::range_error("Bad dimension");
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
        const double *operator[] (int i) const {
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
        void Print(int width = 10, double eps = -1) const;
        friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1) {
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
};

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
    return sqrt(sumOfSquares);
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

int main ()
{
  // TEST ZA KONSTRUKTOR VEKTORA
  {
      std::cout << "TEST ZA KOSNTRUKTOR VEKTORA:\n";
      Vector probni_v(5);
      std::cout << "Saljemo kao parametar 5 i sada velicinu vektora provjeravamo pozivom metode NElems(): ";
      std::cout << probni_v.NElems() << std::endl;
      std::cout << "Sada kao parametar konstruktoru šaljemo respektivno 0 i negativan broj."<< std::endl;
      try {
          Vector null_v(0);
      } catch(std::range_error e) {
          std::cout << "Null dimension - Izuzetak: " << e.what() << std::endl; }

        try {
          Vector neg_v(-5);
      } catch(std::range_error e) {
          std::cout << "Negative dimension - Izuzetak: " << e.what() << std::endl; }
  }

  //TEST ZA SEKVENCIJSKI KONSTRUKTOR

  {
      std::cout << "\nTEST ZA SEKVENCIJSKI KONSTRUKTOR:\n";
      Vector probni_v { 1.1, 1.2, 1.3, 2.2, 3.9 };
      std::cout << "Zadali ste listu velicine: " << probni_v.NElems() << std::endl;
      try {
         Vector null_v{};
      } catch(std::range_error e) { std::cout << "Zadata prazna inic.lista: " << e.what() << std::endl; }
  }

  //TEST ZA OPERATOR[]

  {
      std::cout << "\nTEST ZA operator[]:\n";
      Vector probni_v { 1.1, 1.2, 1.3, 2.2, 3.9 };
      std::cout << "Zadat je vektor v { 1.1, 1.2, 1.3, 2.2, 3.9 }\n";
      std::cout << "1-element je " << probni_v[0] << std::endl
                << "5-ti element je " << probni_v[4] << std::endl;

      std::cout << "Dosad smo ispitivali inspektore, probajmo promijenit vrijednost preko operatora():\n";
      std::cout << "Pomnozimo elemente vektora sa 2\n";
      for(int i = 0; i < 5; i++) probni_v[i] *= 2;
      for(int i = 0; i < 5; i++) std::cout << probni_v[i] << " ";
      std::cout << std::endl;
  }

  //TEST ZA OPERATOR() i NELEMS()

  {
      Vector probni_v { 1.1, 1.2, 1.3, 2.2, 3.9 };
      std::cout << "\nTEST ZA OPERATOR():\n";
      std::cout << "Zadat je vektor v { 1.1, 1.2, 1.3, 2.2, 3.9 }\n";
      std::cout << "1-element je " << probni_v(1) << std::endl
                 << "5-ti element je " << probni_v(5) << std::endl;

      try {
          std::cout << "6-ti element je " << probni_v(6) << std::endl;
      } catch(std::range_error e) {
          std::cout << "Izvan opsega: " << e.what() << std::endl;
      }

      try {
          std::cout << "(-5)-ti element je " << probni_v(6) << std::endl;
      } catch(std::range_error e) {
          std::cout << "Negativna vrijedost: " << e.what() << std::endl;
      }

      std::cout << "Dosad smo ispitivali inspektore, probajmo promijenit vrijednost preko operatora():\n";
      std::cout << "Pomnozimo elemente vektora sa 2\n";
      for(int i = 1; i < 6; i++) probni_v(i) *= 2;
      for(int i = 1; i < 6; i++) std::cout << probni_v(i) << " "; std::cout << std::endl;

      std::cout << "\nTEST ZA NElems(): \n";
      std::cout << "Broj elemenara vektora je: " << probni_v.NElems();
      std::cout << std::endl;
  }

  //TEST ZA Norm(), VctorNorm(), GetEpsilon(), Print(), PrintVector()

  {
      Vector probni_v { 1.1, 1.2, 1.3, 2.2, 3.9 };
      std::cout << "\nTEST ZA Norm() :\n";
      std::cout << "Euklidska norma vektora: " << probni_v.Norm() << std::endl;

      std::cout << "\nTEST ZA VectorNorm() :\n";
      std::cout << "Euklidska norma vektora: " << VectorNorm(probni_v) << std::endl;

      std::cout << "\nTEST ZA GetEpsilon():\n";
      std::cout << "Tolerancija za greske iznosi: " << probni_v.GetEpsilon() << std::endl;

      std::cout << "\nTEST ZA Print(): \n";
      probni_v.Print();

      std::cout << "\nTEST ZA PrintVector():\n";
      PrintVector(probni_v);
  }


  //TEST ZA OPERATORE +, +=, -, -=,/, /=, *, *=

  {
      Vector vec1 { 1.1, 1.2, 1.3, 2.2, 3.9 };
      Vector vec2 { 1, 1, 1, 1 ,1 };
      Vector manji { 1, 2.2 };

      std::cout << "\nTEST ZA operator+ :\n";
      {
          std::cout << "Zbir vec1 + vec2 je (bez promjene vec1) \n ";
          (vec1 + vec2).Print();
          std::cout << "Sabiranje vektora nejednakih dimenzija\n";
          try{
              (vec1 + manji).Print();
          } catch(std::domain_error e) { std::cout <<"IZUZETAK: "<< e.what() << std::endl;}

      }

      std::cout << "\nTEST ZA operator+= :\n";
      {
          vec1 += vec2;
          std::cout << "Dodavanjem vec2 na vec1, vec1 je:\n";
          vec1.Print();
          std::cout << "Dodajemo na vec1 vektor manjeg formata\n";
          try{
              vec1 += manji;
          } catch(std::domain_error e) { std::cout <<"IZUZETAK: "<< e.what() << std::endl;}
      }

      std::cout << "\nTEST ZA operator- :\n";
      {
          std::cout << "Razlika vec1 - vec2 je \n";
          (vec1 - vec2).Print();
          std::cout << "Sabiranje vektora nejednakih dimenzija\n";
          try{
              (vec1 - manji).Print();
          } catch(std::domain_error e) { std::cout <<"IZUZETAK: "<< e.what() << std::endl;}
      }

      std::cout << "\nTEST ZA operator-= :\n";
      {
          vec1 -= vec2;
          std::cout << "Oduzimanjem vec2 od vec1, vec1 je:\n";
          vec1.Print();
          std::cout << "Oduzimamo od vec1 vektor manjeg formata\n";
          try{
              vec1 -= manji;
          } catch(std::domain_error e) { std::cout <<"IZUZETAK: "<< e.what() << std::endl;}
      }

      std::cout << "\nTEST ZA operator* :\n";
      {
          std::cout << "Rezultat 2 * vec1 :\n";
          (2 * vec1).Print();

          std::cout <<  "Rezultat vec1 * 2:\n";
          (vec1 * 2).Print();

          std::cout << "\nTEST ZA operator* :\n";
          std::cout << "Skalarni proizvod vec1 i vec2:\n";
          std::cout << vec1 * vec2 << std::endl;

          std::cout << "Skalarni proizvod vec1 i vec manjeg formata :\n";
          try { vec1 * manji; }
            catch(std::domain_error e) { std::cout << e.what() << std::endl; }

            std::cout << "\nTEST ZA operator*= :\n";
            std::cout << "Mnozimo vec1 sa 2 i vec1 je sada:\n";
            vec1 *= 2;
            vec1.Print();
      }

      std::cout << "\nTEST ZA operator/ :\n";
      {
          std::cout << "Rezultat vec1 / 2 je: \n";
          (vec1 / 2).Print();
          std::cout << "Dijeljenjem sa nulom ocekujemo izuzetak: ";
          try {
            (vec1 / 0);
          }
          catch(std::domain_error e) { std::cout << "IZUZETAK: "<< e.what() << std::endl; }
      }
      std::cout << "\nTEST ZA operator/= :\n";
      {
          std::cout << "Dijeljenjem vec1 sa dva, vec1 je sada: \n";
          vec1 /= 2;
          vec1.Print();
          std::cout << "Pokusaj dijeljenja sa nulom baca izuzetak.\n";
          try{
              vec1 /= 0;
            } catch(std::domain_error e){ std::cout << "IZUZETAK: " << e.what()<< std::endl; }
      }
  }


  //TESTIRAMO KONSTRUKTORE MATRICE

  {
      std::cout << "\nTESTIRAMO KONSTRUKTORE MATRICE:\n";

      try{
        Matrix matrica(5, 2);
        std::cout << "Postavljanje dimenzija matrice uspjelo.\n";
      } catch(...){}
      std::cout << "\nZadavanjem neke od dimenzija da je nula ili negativna:\n";
      try{
          Matrix matrica(0, 2);
      } catch(std::range_error e) { std::cout << "IZUZETAK: "<< e.what() << std::endl;}

      try{
          Matrix matrica(2, 0);
      } catch(std::range_error e) { std::cout << "IZUZETAK: "<< e.what() << std::endl;}

      try{
          Matrix matrica(2,-1);
      } catch(std::range_error e) { std::cout << "IZUZETAK: "<< e.what() << std::endl;}

      std::cout << "\nTESTIRAMO KONSTRUKTOR MATRICE KOJI PRIMA VEKTOR:\n";
      {
          Vector vektor{1.1, 1.2, 1.3};


            Matrix matrica(vektor);
            std::cout << "Pretvorba vektora u matricu uspjela.\n";

      }

      std::cout << "\nTESTIRAMO SEKVENCIJSKI KONSTRUKTOR:\n";
      {
          try{
            Matrix matrica{{1,2,3}, {3,4,5}};
            std::cout << "Pretvorba vektora u matricu uspjela.\n";
          } catch(...){}

          try{
              Matrix matrica{{1,2}, {3,4,5}};
          } catch(std::logic_error e) { std::cout << "IZUZETAK: " << e.what() << std::endl;}

          try{
              Matrix matrica{{}};
          } catch(std::range_error e) { std::cout << "IZUZETAK: "<< e.what() << std::endl;}
      }
  }

  //TESTIRAMO METODE i OPERATORE

  {
      Matrix matrica {{1,2,3}, {4,5,6}, {7,8,9}};

      std::cout << "\nTESTIRAMO METODU PRINT():\n";
      matrica.Print();

      std::cout << "\nTESTIRAMO METODE NRows() I NCols():\n";
      std::cout << "Matrica ima "<< matrica.NRows() <<" reda i " << matrica.NCols() << " kolona\n";

      std::cout << "\nTESTIRAMO METODU Norm():\n";
      std::cout << "Frobeniusova norma matrice je: " << matrica.Norm() << std::endl;

      std::cout << "\nTESTIRAMO METODU MatrixNorm():\n";
      std::cout << "Frobeniusova norma matrice je: " << MatrixNorm(matrica) << std::endl;

      std::cout << "\nTESTIRAMO METODU GetEpsilon():\n";
      std::cout << "Prag je " << matrica.GetEpsilon() << std::endl;

      Matrix matrica_2{{1,2,1}, {4,1,6}, {5,8,9}};
      Matrix manja {{1,2}, {3,4}};

      std::cout << "\nTESTIRAMO OPERATOR +\n";

      std::cout << "Zbir sljedecih matrica :\n";
      matrica.Print();
      std::cout << std::endl;
      PrintMatrix(matrica_2);

      std::cout << "Zbir je: \n";
      (matrica + matrica_2).Print();

      std::cout << "Ukoliko pokusamo sabrati sljedece matrice razlicitog formata:\n";

      matrica_2.Print(); std::cout << std::endl;
      PrintMatrix(manja);
      try{
            matrica_2 + manja;
        } catch(std::domain_error e) { std::cout << "IZUZETAK: " << e.what() << std::endl; }


      std::cout << "\nTESTIRAMO OPERATOR -\n";

      std::cout << "Razlika sljedecih matrica :\n";
      matrica.Print();
      std::cout << std::endl;
      PrintMatrix(matrica_2);

      std::cout << "Razlika je: \n";
      (matrica - matrica_2).Print();

      std::cout << "Ukoliko pokusamo oduzeti sljedece matrice razlicitog formata:\n";

      matrica_2.Print(); std::cout << std::endl;
      PrintMatrix(manja);
      try{
            matrica_2 - manja;
        } catch(std::domain_error e) { std::cout << "IZUZETAK: " << e.what() << std::endl; }


      std::cout << "\nTESTIRAMO OPERATOR +=\n";

      std::cout << "Dodavanjem matrice :\n";
      matrica.Print();
      std::cout << std::endl << "Na matricu \n";
      PrintMatrix(matrica_2);

      std::cout << "Druga matrica ce biti: \n";
      matrica_2 += matrica;
      (matrica_2).Print();

      std::cout << "Ukoliko pokusamo na matricu :\n";

      matrica_2.Print();
      std::cout << std::endl << "Dodati matricu manjeg formata\n";
      PrintMatrix(manja);
      try{
            matrica_2 += manja;
        } catch(std::domain_error e) { std::cout << "IZUZETAK: " << e.what() << std::endl; }


      std::cout << "\nTESTIRAMO OPERATOR -=\n";

      std::cout << "Oduzimanjem matrice :\n";
      matrica.Print();
      std::cout << std::endl << "Od matrice \n";
      PrintMatrix(matrica_2);

      std::cout << "Druga matrica ce biti: \n";
      matrica_2 -= matrica;
      (matrica_2).Print();

      std::cout << "Ukoliko pokusamo od matrice :\n";

      matrica_2.Print();
      std::cout << std::endl << "Oduzeti matricu manjeg formata\n";
      PrintMatrix(manja);
      try{
            matrica_2 -= manja;
        } catch(std::domain_error e) { std::cout << "IZUZETAK: " << e.what() << std::endl; }


      std::cout << "\nTESTIRAMO OPERATOR *=\n";
      std::cout << "Mnozenjem ispod matrice sa 2:\n";
      matrica.Print();
      matrica *= 2;
      std::cout << std::endl << "Ona ce biti: \n";
      matrica.Print();

      std::cout << "\nTESTIRAMO OPERATORE *\n";
      manja.Print();
      std::cout << "Mnozenjem date matrice sa 2 (operator ciji je prvi parametar matrica, a drugi skalar):\n";
      (manja * 2).Print();
      std::cout << "Mnozenjem date matrice sa 2 (operator ciji je prvi parametar skalar, a drugi matrica):\n";
      (2 * manja).Print();

      std::cout << "\nTESTIRAMO MNOZENJE MATRICA\n";
      std::cout << "Rezultat mnozenja matrica (bez promjene ijedne od njih):\n";
      matrica.Print();
      std::cout << std::endl;
      matrica_2.Print();
      std::cout << "Result: \n";
      (matrica * matrica_2).Print();

      std::cout << std::endl << "Mnozenjem prve matrice drugom, prva ce biti:\n";

      matrica *= matrica_2;
      matrica.Print();
      std::cout << "ako sada nju zelimo pomnoziti operatorom *= sa : \n";
      manja.Print();
      try{
          matrica *= manja;
      } catch(std::domain_error e) { std::cout << "IZUZETAK: " << e.what() << std::endl;}

      std::cout << " a ako je zelimo pomnoziti operatorom * sa : \n";
      manja.Print();
      try{
          matrica * manja;
      } catch(std::domain_error e) { std::cout << "IZUZETAK: " << e.what() << std::endl;}


      std::cout << "\nTESTIRAMO TRANSPONOVANJE MATRICA: \n";
      std::cout << "Rezultat transponovanja matrice:\n";
      matrica.Print();
      std::cout << "ce biti\n";
      (Transpose(matrica)).Print();

      std::cout << "\nTransponovana matrica: \n";
      matrica_2.Print();
      matrica_2.Transpose();
      std::cout << std::endl << "ce sada postati: \n";
      matrica_2.Print();

      std::cout << "\nTransponovana matrica: \n";
      {
          Matrix pravougaona{{11,12,13},{24,45,67}};
          pravougaona.Print();
          pravougaona.Transpose();
          std::cout << std::endl << "ce sada postati: \n";
          pravougaona.Print();
      }

      std::cout << "\nTESTIRAMO MNOZENJE MATRICE I VEKTORA: \n";

      {
          Vector vektor{1,2,3};
          std::cout << "Mnozimo sljedecu matricu i vektor:\n";
          matrica_2.Print();
          std::cout << std::endl;
          vektor.Print();
          std::cout << std::endl << "Result je: \n";
          (matrica_2 * vektor).Print();

          std::cout << "Ako pokusamo pomnoziti sljedecu matricu i vektor\n";
          manja.Print();
          std::cout << std::endl;
          vektor.Print();
          std::cout << std::endl << "Result je: \n";
          try{
            (manja* vektor).Print();
          }catch(...) { std::cout << "IZUZETAK !"; }
      }

      std::cout << "\nTESTIRAMO OPERATOR []: \n";
      {
          matrica.Print();

          std::cout << "Dijagonalni elementi su: " << matrica[0][0] <<" " <<matrica[1][1] << " " << matrica[2][2] << std::endl;


          std::cout << "\nTESTIRAMO OPERATOR (): \n";
          matrica_2.Print();

          std::cout << "Dijagonalni elementi su: " << matrica_2(1,1) <<" " <<matrica_2(2,2) << " " << matrica_2(3,3) << std::endl;

          std::cout << "Medjutim ako probamo pristupiti sa invalidnim indeksom (negativan ili izvan opsega): \n";
          manja.Print();
          std::cout << "\nElement (-1,3) je \n";
          try {
              matrica(-1,3);
          } catch (std::range_error e){ std::cout << "IZUZETAK: " << e.what() << std::endl;}

          std::cout << "\nElement (4,4) je \n";
          try {
              matrica(4,4);
          } catch (std::range_error e){ std::cout << "IZUZETAK: " << e.what() << std::endl;}

          std::cout << "Isprobali smo inspektore, sada probbajmo promijeniti neke elemente matrice:\n";
          manja.Print();
          std::cout << std::endl;
          manja[0][0] *= 2;
          manja(1,2) *= 5;
          manja.Print();
      }
  }
    return 0;
}
