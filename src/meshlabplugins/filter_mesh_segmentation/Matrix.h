#include <vector>
#include <algorithm>

class Matrix;
class SparseMatrix;

void save(std::vector<double> towrite, char *filename);
Matrix Heaviside(Matrix rhs);
Matrix onehot(std::vector<int> labels);
Matrix onehot(std::vector<double> labels);
double squared_distance(std::vector<double> a, std::vector<double> b);
std::vector<int> kmeans(Matrix u, int k);
std::vector<int> kmeans(std::vector<double> u, int k);
Matrix operator*(const double& scalar, Matrix matrix);
SparseMatrix operator*(const double& scalar, SparseMatrix matrix);
SparseMatrix sparse_exp(SparseMatrix rhs);
SparseMatrix sparse_abs(SparseMatrix rhs);
SparseMatrix sparse_log(SparseMatrix rhs);
SparseMatrix sparse_pow(SparseMatrix rhs, const double& p);
void print(std::vector<double> x);
void print(std::vector<int> x);

Matrix sign(Matrix rhs);

double norm(std::vector<double> rhs);
std::vector<int> arange(int a, int b);
std::vector<int> arange(int b);
std::vector<int> nonzero(std::vector<double> rhs);

std::vector<double> splice(std::vector<double> w1, std::vector<double> w2, std::vector<bool> J);
std::vector<double> sum(Matrix rhs, int axis);
std::vector<int> argmax(Matrix rhs, int axis);
double max(const std::vector<double>& v);
int max(const std::vector<int>& v);

std::vector<double> sqrt(std::vector<double> rhs);
std::vector<double> sign(std::vector<double> rhs);
std::vector<int> Heaviside(std::vector<double> rhs);

Matrix conjgrad(Matrix& A, Matrix& b, double tol = 1e-10, int T=100000);
Matrix conjgrad(SparseMatrix& A, Matrix& b, double tol = 1e-10, int T=100000);
std::vector<double> conjgrad(SparseMatrix& A, std::vector<double>& b, double tol = 1e-10, int T=100000);

std::vector<bool> operator!(std::vector<bool> a);
std::vector<bool> operator>(std::vector<double> a, double b);
std::vector<bool> operator<(std::vector<double> a, double b);
std::vector<double> operator-(const std::vector<double>& v, double scalar);
std::vector<double> operator/(const std::vector<double>& v, double scalar);
std::vector<double> operator*(const std::vector<double>& v, double scalar);
std::vector<double> operator*(double scalar, std::vector<double> v);
std::vector<double> operator*(const std::vector<bool>& v, double scalar);
std::vector<double> operator*(double scalar, std::vector<bool> v);
std::vector<double> operator/(const std::vector<double>& v, const std::vector<double>& w);
std::vector<double> operator+(const std::vector<double>& v, const std::vector<double>& w);
std::vector<double> operator-(const std::vector<double>& v, const std::vector<double>& w);


class Matrix {
   friend class SparseMatrix;
   private:
      std::vector<double> data;
      std::vector<double*> mat;
      int rows, cols;

   protected:
   
      Matrix(int _rows, int _cols, const double& min, const double& max);
      Matrix(int _rows, int _cols, const double& initial);
      Matrix(int _rows);

   public:
      static Matrix eye(int _rows) {return Matrix(_rows);};
      static Matrix ones(int _rows, int _cols) {return Matrix(_rows,_cols,1);};
      static Matrix zeros(int _rows, int _cols) {return Matrix(_rows,_cols,0);};
      static Matrix rand(int _rows, int _cols) {return Matrix(_rows,_cols,0,1);};
      Matrix(int _rows, int _cols);
      Matrix();

      Matrix(const std::vector<double>& rhs);
      Matrix(const Matrix& rhs);
      virtual ~Matrix();

      //Operator overloading for standard mathematical matrix operations
      Matrix& operator=(const Matrix& rhs);

      //Matrix mathematical operations
      Matrix operator+(const Matrix& rhs);
      Matrix& operator+=(const Matrix& rhs);
      Matrix operator-(const Matrix& rhs);
      Matrix& operator-=(const Matrix& rhs);
      Matrix operator*(const Matrix& rhs);
      Matrix operator*(const SparseMatrix& rhs);
      Matrix& operator*=(const Matrix& rhs);
      Matrix transpose();

      //Matrix/scalar operations
      Matrix operator+(const double& rhs);
      Matrix operator-(const double& rhs);
      Matrix operator*(const double& rhs);
      Matrix operator/(const double& rhs);

      //Matrix vector operations
      Matrix operator+(const std::vector<double>& rhs);
      Matrix& operator+=(const std::vector<double>& rhs);
      Matrix operator-(const std::vector<double>& rhs);
      Matrix& operator-=(const std::vector<double>& rhs);

      //Matrix/vector<double> operations
      std::vector<double> operator*(const std::vector<double>& rhs);
      std::vector<double> diag();
      std::vector<double> tovector();

      Matrix row_divide(const std::vector<double>& rhs);
      Matrix row_multiply(const std::vector<double>& rhs);
      Matrix col_multiply(const std::vector<double>& rhs);
      Matrix multiply(Matrix& rhs);
      Matrix divide(Matrix& rhs);
      Matrix sign();
      Matrix Heaviside();

      //Access individual elements
      double& operator()(const int& i, const int& j);
      const double& operator()(const int& i, const int& j) const;

      //Extract submatrix
      Matrix operator()(const std::vector<int>& I, const std::vector<int>& J);

      //Save to file
      void save(char *filename);


      //Append column
      void append_col(const Matrix& w);

      //Access the row and column sizes
      int get_rows() const;
      int get_cols() const;
      std::vector<int> shape() const;

      //Print matrix
      void print();

      //Operations
      std::vector<int> argmin(int axis);
      std::vector<int> argmax(int axis);
      double min();
      std::vector<double> min(int axis);
      double max();
      std::vector<double> max(int axis);
      double norm();
      std::vector<double> norm(int axis);
      double sum();
      std::vector<double> sum(int axis);
      double mean();
      std::vector<double> mean(int axis);


};

class SparseMatrix {
   friend class Matrix;
   private:
      struct node {
         double val;
         int col;

         bool operator<( const node& rhs ) const { return col < rhs.col; }
         node operator+( const node& rhs ) const { return node{val+rhs.val,col}; }
      };
      int rows, cols;      
      std::vector<std::vector<node>> mat;

   protected:

      SparseMatrix(int _rows);
      SparseMatrix(const std::vector<double>& diags);

   public:

      static SparseMatrix speye(int _rows) {return SparseMatrix(_rows);};
      static SparseMatrix spdiags(const std::vector<double>& diags) {return SparseMatrix(diags);};
      SparseMatrix(int _rows, int _cols);
      SparseMatrix();
      
      SparseMatrix(const SparseMatrix& rhs);
      virtual ~SparseMatrix();

      SparseMatrix& operator=(const SparseMatrix& rhs);
      SparseMatrix operator+(const SparseMatrix& rhs);
      SparseMatrix& operator+=(const SparseMatrix& rhs);
      SparseMatrix operator-(const SparseMatrix& rhs);
      SparseMatrix& operator-=(const SparseMatrix& rhs);
      SparseMatrix operator*(const SparseMatrix& rhs);
      SparseMatrix& operator*=(const SparseMatrix& rhs);
      SparseMatrix transpose() const;

      //Multiplication with dense matrix on right
      Matrix operator*(const Matrix& rhs);

      //Convert to dense matrix
      Matrix todense();

      //Matrix/scalar operations
      SparseMatrix operator+(const double& rhs);
      SparseMatrix operator-(const double& rhs);
      SparseMatrix operator*(const double& rhs);
      SparseMatrix operator/(const double& rhs);
      
      //Matrix/vector<double> operations
      std::vector<double> operator*(const std::vector<double>& rhs);
      std::vector<double> diag();

      //Operations on individual elements
      SparseMatrix pow(const double& p);
      SparseMatrix exp();
      SparseMatrix log();
      SparseMatrix abs();
      SparseMatrix row_divide(const std::vector<double>& rhs);
      SparseMatrix col_divide(const std::vector<double>& rhs);
      std::vector<int> argmin(int axis);
      std::vector<int> argmax(int axis);
     
      //Extract submatrix
      SparseMatrix operator()(const std::vector<int>& I, const std::vector<int>& J);
      SparseMatrix operator()(const std::vector<bool>& I, const std::vector<bool>& J);

      //Add elements to matrix
      void insert(int i, int j, const double& _val);

      //Access elements
      double operator()(const int& i, const int& j);

      int num_nonzero();
      double sum();
      double norm();
      std::vector<double> sum(int axis);

      //Eigenvalue operations
      double LargestEigenvalue(double tol=1E-10);
      //double Largest_NonConst_Eigenvalue(double tol=1E-10);
      double SmallestEigenvalue(double tol=1E-10);
      //double Smallest_NonConst_Eigenvalue(double tol=1E-10);
      //std::vector<double> Smallest_NonConst_Eigenvector(double tol=1E-10);
      //Matrix Smallest_NonConst_Eigenvectors(int k, double tol=1E-10);
      Matrix Largest_Eigenvector(Matrix u, double tol=1E-10);

      //Clean up
      void cleanup();
     
      //Print
      void print();

      //Access the row and column sizes
      int get_rows() const;
      int get_cols() const;
};


