#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <cmath>


std::vector<double> splice(std::vector<double> w1, std::vector<double> w2, std::vector<bool> J){
   std::vector<double> w(J.size());
   int w1_ind=0, w2_ind=0;
   for(int i=0; i<J.size(); i++){
      if(J[i])
         w[i] = w1[w1_ind++];
      else
         w[i] = w2[w2_ind++];
   }
   return w;
}

std::vector<double> operator*(const std::vector<bool>& v, double scalar){
   std::vector<double> result(v.size());
   for(int i=0; i<v.size(); i++)
      result[i] = ((double)v[i])*scalar;
   return result;
}
std::vector<double> operator*(double scalar, std::vector<bool> v){ return v*scalar;}

std::vector<double> operator*(const std::vector<double>& v, double scalar){
   std::vector<double> result(v.size());
   for(int i=0; i<v.size(); i++)
      result[i] = v[i]*scalar;
   return result;
}
std::vector<double> operator*(double scalar, std::vector<double> v){ return v*scalar;}

std::vector<bool> operator!(std::vector<bool> a){
   std::vector<bool> result(a.size());
   for(int i=0; i<a.size(); i++)
      result[i] = !a[i];
   return result;
}

Matrix Heaviside(Matrix rhs){ return rhs.Heaviside(); }

Matrix sign(Matrix rhs) {return rhs.sign();}
//Saves matrix to a file 
void save(std::vector<double> towrite, char *filename){

   double length = (double)towrite.size();
   towrite.insert(towrite.begin(),1.0);
   towrite.insert(towrite.begin(),length);

   std::ofstream fout(filename, std::ios::out | std::ios::binary);
   fout.write((char*)&towrite[0], towrite.size() * sizeof(double));
   fout.close();
}
void Matrix::save(char *filename){

   std::vector<double> towrite = data;
   towrite.insert(towrite.begin(),cols);
   towrite.insert(towrite.begin(),rows);

   std::ofstream fout(filename, std::ios::out | std::ios::binary);
   fout.write((char*)&towrite[0], towrite.size() * sizeof(double));
   fout.close();
}
Matrix onehot(std::vector<double> labels){
   
   int cols = max(labels)+1;
   int rows = labels.size();
   Matrix result = Matrix::zeros(rows,cols);
   for(int i=0; i<rows; i++)
      result(i,(int)(labels[i]))=1.0;
   
   return result;
}
Matrix onehot(std::vector<int> labels){
   
   int cols = max(labels)+1;
   int rows = labels.size();
   Matrix result = Matrix::zeros(rows,cols);
   for(int i=0; i<rows; i++)
      result(i,labels[i])=1.0;
   
   return result;
}

double squared_distance(std::vector<double> a, std::vector<double> b){
   double result = 0.0;
   for(int i=0; i<a.size(); i++)
      result += (a[i]-b[i])*(a[i]-b[i]);
   return result;
}
std::vector<int> kmeans(std::vector<double> u, int k){ return kmeans(Matrix(u), k); }
std::vector<int> kmeans(Matrix u, int k){
   
   int n = u.get_rows();
   int m = u.get_cols();
   std::vector<int> labels(n,0);

   //Convert to vector of vectors
   std::vector<std::vector<double>> P(n);
   for(int i=0; i<n; i++)
      for(int j=0; j<m; j++)
         P[i].push_back(u(i,j));

   //Random intial cluster centers
   std::vector<std::vector<double>> C(k);
   for(int i=0; i<k; i++){
      int j = rand()%n;
      C[i] = P[j];
   }

   int num_changed = 1;
   while(num_changed > 0){

      num_changed = 0;
   
      //Assignment based on distance
      for(int i=0; i<n; i++){
         int prev_label = labels[i];
         labels[i] = 0;
         double min_dist = squared_distance(P[i],C[labels[i]]);
         for(int j=1; j<k; j++){
            double dist = squared_distance(P[i],C[j]);
            if(dist < min_dist){
               labels[i] = j;
               min_dist = dist;
            }   
         }
         if(prev_label != labels[i])
            num_changed++;
      }
      printf("num_changed = %d\n",num_changed);

      //Recompute centroids
      for(int j=0; j<k; j++){
         std::fill(C[j].begin(), C[j].end(), 0.0);
         int num =0;
         for(int i=0; i<n; i++){ 
            if(labels[i] == j){
               C[j] = C[j] + P[i];
               num++;
            }
         }
         printf("%d,",num);
         C[j] = C[j]/(double)num;
      }
      printf("\n");
   }

   return labels;
}
void print(std::vector<int> x){
   printf("\n[");
   for(int i=0; i<x.size()-1; i++)
      printf("%d,\n",x[i]);
   printf("%d]\n",x[x.size()-1]);
}
void print(std::vector<double> x){
   printf("\n[");
   for(int i=0; i<x.size()-1; i++)
      printf("%f,\n",x[i]);
   printf("%f]\n",x[x.size()-1]);
}
std::vector<int> arange(int a, int b){
   std::vector<int> result;
   if(a < b){
      result.resize(b-a);
      for(int i=a; i<b; i++)
         result[i-a] = i;
   }
   return result;
}

std::vector<int> arange(int b){return arange(0,b);}

std::vector<bool> operator>(std::vector<double> a, double b){
   std::vector<bool> result(a.size(),0);
   for(int i=0; i<a.size(); i++)
      if(a[i] > b)
         result[i]=1;
   return result;
}
std::vector<bool> operator<(std::vector<double> a, double b){
   std::vector<bool> result(a.size(),0);
   for(int i=0; i<a.size(); i++)
      if(a[i] < b)
         result[i]=1;
   return result;
}
std::vector<int> nonzero(std::vector<double> rhs){
   std::vector<int> result;
   for(int i=0; i<rhs.size(); i++)
      if(rhs[i] != 0.0)
         result.push_back(i);
   return result;
}
std::vector<double> sign(std::vector<double> rhs){
   std::vector<double> result(rhs.size(),-1.0);
   for(int i=0; i < rhs.size(); i++)
      if(rhs[i]>0)
         result[i] = 1.0;
   return result;
}
std::vector<int> Heaviside(std::vector<double> rhs){
   std::vector<int> result(rhs.size(),0);
   for(int i=0; i < rhs.size(); i++)
      if(rhs[i]>0)
         result[i] = 1;
   return result;
}
std::vector<double> sqrt(std::vector<double> rhs){
   std::vector<double> result(rhs.size());
   for(int i=0; i < rhs.size(); i++)
      result[i] = sqrt(rhs[i]);
   return result;
}

std::vector<double> sum(Matrix rhs, int axis){return rhs.sum(axis);}
std::vector<int> argmax(Matrix rhs, int axis){return rhs.argmax(axis);}
double norm(std::vector<double> rhs){
   double result = 0;
   for(int i=0; i<rhs.size(); i++)
      result += rhs[i]*rhs[i];
   return sqrt(result);
}
std::vector<double> conjgrad(SparseMatrix& A, std::vector<double>& b, double tol, int T){

   Matrix bmat(b);
   Matrix x = conjgrad(A, bmat, tol, T);
   return x.tovector();
}
Matrix conjgrad(SparseMatrix& A, Matrix& b, double tol, int T){

   Matrix x = Matrix::zeros(b.get_rows(),b.get_cols());
   Matrix r = b - A*x;
   Matrix p = r;
   std::vector<double> rsold = sum(r.multiply(r),0);

   int i;
   for(i=0; i<T; i++){
      Matrix Ap = A*p;
      std::vector<double> alpha = rsold/sum(Ap.multiply(p),0);
      x += p.col_multiply(alpha);
      r -= Ap.col_multiply(alpha);

      std::vector<double> rsnew = sum(r.multiply(r),0);
      double err = norm(rsnew);
      printf("i:%d,err:%f\n",i,err);
      if(err < tol)
         break;
      
      p = r + p.col_multiply(rsnew/rsold);
      rsold = rsnew;
   }
   std::cout << i << std::endl;
   return x;
}
Matrix conjgrad(Matrix& A, Matrix& b, double tol, int T){

   Matrix x = Matrix::zeros(b.get_rows(),b.get_cols());
   Matrix r = b - A*x;
   Matrix p = r;
   std::vector<double> rsold = sum(r.multiply(r),0);

   int i;
   for(i=0; i<T; i++){
      Matrix Ap = A*p;
      std::vector<double> alpha = rsold/sum(Ap.multiply(p),0);
      x += p.col_multiply(alpha);
      r -= Ap.col_multiply(alpha);

      std::vector<double> rsnew = sum(r.multiply(r),0);
      if(norm(rsnew) < tol)
         break;
      
      p = r + p.col_multiply(rsnew/rsold);
      rsold = rsnew;
   }
   std::cout << i << std::endl;
   return x;
}

Matrix operator*(const double& scalar, Matrix matrix) {return matrix * scalar;}
std::vector<double> operator-(const std::vector<double>& v, const std::vector<double>& w) {

   std::vector<double> result(v.size());
   for(int i=0; i<v.size(); i++)
      result[i] = v[i]-w[i];

   return result;
}
std::vector<double> operator+(const std::vector<double>& v, const std::vector<double>& w) {

   std::vector<double> result(v.size());
   for(int i=0; i<v.size(); i++)
      result[i] = v[i]+w[i];

   return result;
}
std::vector<double> operator/(const std::vector<double>& v, const std::vector<double>& w) {

   std::vector<double> result(v.size());
   for(int i=0; i<v.size(); i++)
      result[i] = v[i]/w[i];

   return result;
}
std::vector<double> operator-(const std::vector<double>& v, double scalar) {

   std::vector<double> result(v.size());
   for(int i=0; i<v.size(); i++)
      result[i] = v[i]-scalar;

   return result;
}
std::vector<double> operator/(const std::vector<double>& v, double scalar) {

   std::vector<double> result(v.size());
   for(int i=0; i<v.size(); i++)
      result[i] = v[i]/scalar;

   return result;
}
int max(const std::vector<int>& v){
   int result = v[0];
   for(int i=0; i<v.size(); i++)
      result = std::max(v[i],result);
   return result;
}

double max(const std::vector<double>& v){
   double result = v[0];
   for(int i=0; i<v.size(); i++)
      result = std::max(v[i],result);
   return result;
}
SparseMatrix operator*(const double& scalar, SparseMatrix matrix) {return matrix * scalar;}

SparseMatrix sparse_exp(SparseMatrix rhs) {return rhs.exp();}

SparseMatrix sparse_abs(SparseMatrix rhs) {return rhs.abs();}

SparseMatrix sparse_log(SparseMatrix rhs) {return rhs.log();}

SparseMatrix sparse_pow(SparseMatrix rhs, const double& p) {return rhs.pow(p);}


SparseMatrix SparseMatrix::abs() {
   SparseMatrix result(rows, cols);
   cleanup();
   for(int i=0; i<rows; i++) 
      for(int j=0; j<mat[i].size(); j++) 
         result.mat[i].push_back(node{std::abs(mat[i][j].val),mat[i][j].col});

   return result;
}

SparseMatrix SparseMatrix::exp() {
   SparseMatrix result(rows, cols);
   cleanup();
   for(int i=0; i<rows; i++) 
      for(int j=0; j<mat[i].size(); j++) 
         result.mat[i].push_back(node{std::exp(mat[i][j].val),mat[i][j].col});

   return result;
}


SparseMatrix SparseMatrix::log() {
   SparseMatrix result(rows, cols);
   cleanup();
   for(int i=0; i<rows; i++) 
      for(int j=0; j<mat[i].size(); j++) 
         result.mat[i].push_back(node{std::log(mat[i][j].val),mat[i][j].col});

   return result;
}
SparseMatrix SparseMatrix::col_divide(const std::vector<double>& rhs){
   SparseMatrix result(rows, cols);
   cleanup();
   for(int i=0; i<rows; i++) 
      for(int j=0; j<mat[i].size(); j++) 
         result.mat[i].push_back(node{mat[i][j].val/rhs[mat[i][j].col],mat[i][j].col});
   return result;
}
SparseMatrix SparseMatrix::row_divide(const std::vector<double>& rhs){
   SparseMatrix result(rows, cols);
   cleanup();
   for(int i=0; i<rows; i++) 
      for(int j=0; j<mat[i].size(); j++) 
         result.mat[i].push_back(node{mat[i][j].val/rhs[i],mat[i][j].col});
   return result;
}

SparseMatrix SparseMatrix::pow(const double& p) {
   SparseMatrix result(rows, cols);
   cleanup();
   for(int i=0; i<rows; i++) 
      for(int j=0; j<mat[i].size(); j++) 
         result.mat[i].push_back(node{std::pow(mat[i][j].val,p),mat[i][j].col});
   return result;
}


//Convert to dense matrix
Matrix SparseMatrix::todense() {

   Matrix result(rows,cols,0);

   for(int i=0; i<rows; i++) 
      for(int j=0; j<mat[i].size(); j++) 
         result.mat[mat[i][j].col][i] += mat[i][j].val;

   return result;

}

Matrix Matrix::operator()(const std::vector<int>& I, const std::vector<int>& J){
   Matrix result(I.size(),J.size(),0.0);
   for(int j=0; j<J.size(); j++)
      for(int i=0; i<I.size(); i++)
         result.mat[j][i] = mat[J[j]][I[i]];
   return result;
}

//Multiply by a sparse matrix on right
Matrix Matrix::operator*(const SparseMatrix& rhs) {
   Matrix result(rows,rhs.get_cols(),0);

   if(rhs.get_rows() != cols){
      throw std::length_error("Matrices have incompatible sizes for multiplication.");
   }else{
      SparseMatrix rhs_tran = rhs.transpose();
      for(int j=0; j<rhs_tran.get_rows(); j++){
         for(int i=0; i<rows; i++) 
            for(int k=0; k<rhs_tran.mat[j].size(); k++) 
               result.mat[j][i] += mat[rhs_tran.mat[j][k].col][i] * rhs_tran.mat[j][k].val;
      }
   }

   return result;
}
//Extract submatrix
SparseMatrix SparseMatrix::operator()(const std::vector<bool>& I, const std::vector<bool>& J){

   std::vector<int> I1, J1;
   for(int i=0; i<std::min((int)I.size(),rows); i++)
      if(I[i])
         I1.push_back(i);

   for(int j=0; j<std::min((int)J.size(),cols); j++)
      if(J[j])
         J1.push_back(j);

   return (*this)(I1,J1);
   /*int new_rows = std::accumulate(I.begin(), I.end(), 0);
   int new_cols = std::accumulate(J.begin(), J.end(), 0);
   SparseMatrix result(new_rows, new_cols);
   SparseMatrix temp(new_rows, cols);

   cleanup();
   //Subsample rows first
   int r=0;
   for(int i=0; i<rows; i++){
      if(I[i]){
         temp.mat[r] = mat[i];
         r++;
      }
   }

   temp = temp.transpose();

   //Now sample columns
   r=0;
   for(int j=0; j<J.size(); j++){
      if(J[j]){
         result.mat[r] = temp.mat[j];
         r++;
      }
   }

   return result.transpose();*/
}
//Extract submatrix
SparseMatrix SparseMatrix::operator()(const std::vector<int>& I, const std::vector<int>& J){

   SparseMatrix result(I.size(),J.size());
   SparseMatrix temp(I.size(),cols);

   cleanup();
   //Subsample rows first
   for(int i=0; i<I.size(); i++)
      temp.mat[i] = mat[I[i]];

   temp = temp.transpose();

   //Now sample columns
   for(int j=0; j<J.size(); j++)
      result.mat[j] = temp.mat[J[j]];

   return result.transpose();
}
//Multiply by a dense matrix on right
Matrix SparseMatrix::operator*(const Matrix& rhs) {
   Matrix result(rows,rhs.get_cols(),0);

   if(rhs.get_rows() != cols){
      throw std::length_error("Matrices have incompatible sizes for multiplication.");
   }else{
      for(int j=0; j<rhs.get_cols(); j++){
         for(int i=0; i<rows; i++) 
            for(int k=0; k<mat[i].size(); k++) 
               result.mat[j][i] += mat[i][k].val * rhs.mat[j][mat[i][k].col];
      }
   }

   return result;
}

//Multiply a matrix with a vector<double>
std::vector<double> SparseMatrix::operator*(const std::vector<double>& rhs) {
   std::vector<double> result(rows,0.0);

   if(rhs.size() != cols){
      throw std::length_error("Matrix and vector<double> have incompatible sizes for multiplication.");
   }else{
      for(int i=0; i<rows; i++) 
         for(int j=0; j<mat[i].size(); j++) 
            result[i] += mat[i][j].val * rhs[mat[i][j].col];
   }

   return result;
}

//Obtain a vector<double> of the diagonal elements
std::vector<double> SparseMatrix::diag() {

   std::vector<double> result(rows);
   if(rows != cols){
      throw std::length_error("Requested diagonal of non-square matrix.");
   }else{
      for(int i=0; i<rows; i++) 
         result[i] = (*this)(i,i);
   }
   return result;
}
//Left multiplication of this matrix and another
SparseMatrix SparseMatrix::operator*(const SparseMatrix& rhs) {

   SparseMatrix result(rows, rhs.get_cols());

   if(rhs.get_rows() != cols){
      throw std::length_error("Matrices have incompatible sizes for multiplication.");
   }else{
      cleanup();
      SparseMatrix tran = rhs.transpose();
      for(int i=0; i<rows; i++){
         if(!mat[i].empty()){
            for(int j=0; j<rhs.get_cols(); j++){
               if(!tran.mat[j].empty()){
                  std::vector<node> dot_vec = mat[i];
                  dot_vec.insert(dot_vec.end(),tran.mat[j].begin(),tran.mat[j].end());
                  std::sort(dot_vec.begin(),dot_vec.end());
                  double dot = 0;
                  for(int k=0; k<dot_vec.size()-1; k++)
                     if(dot_vec[k].col==dot_vec[k+1].col)
                        dot += dot_vec[k].val*dot_vec[k+1].val;
                  result.mat[i].push_back(node{dot,j});
               }
            }
         }
      }
   }
   result.cleanup();
   return result;
}

//Cumulative left multiplication of this matrix and another
SparseMatrix& SparseMatrix::operator*=(const SparseMatrix& rhs) {
   SparseMatrix result = (*this) * rhs;
   (*this) = result;
   return *this;
}

//Matrix/scalar addition
SparseMatrix SparseMatrix::operator+(const double& rhs) {
   SparseMatrix result(rows, cols);

   for(int i=0; i<rows; i++){
      result.mat[i] = mat[i];
      for(int j=0; j<mat[i].size(); j++)
         result.mat[i][j].val += rhs;
   }

   result.cleanup();
   return result;
}

//Matrix/scalar subtraction
SparseMatrix SparseMatrix::operator-(const double& rhs) {
   SparseMatrix result(rows, cols);

   for(int i=0; i<rows; i++){
      result.mat[i] = mat[i];
      for(int j=0; j<mat[i].size(); j++)
         result.mat[i][j].val -= rhs;
   }

   result.cleanup();
   return result;
}

//Matrix/scalar multiplication
SparseMatrix SparseMatrix::operator*(const double& rhs) {
   SparseMatrix result(rows, cols);

   for(int i=0; i<rows; i++){
      result.mat[i] = mat[i];
      for(int j=0; j<mat[i].size(); j++)
         result.mat[i][j].val *= rhs;
   }

   result.cleanup();
   return result;
}

//Matrix/scalar division
SparseMatrix SparseMatrix::operator/(const double& rhs) {
   SparseMatrix result(rows, cols);

   for(int i=0; i<rows; i++){
      result.mat[i] = mat[i];
      for(int j=0; j<mat[i].size(); j++)
         result.mat[i][j].val /= rhs;
   }

   result.cleanup();
   return result;
}
SparseMatrix& SparseMatrix::operator=(const SparseMatrix& rhs) {

   if (&rhs == this)
      return *this;

   rows = rhs.get_rows();
   cols = rhs.get_cols();
   mat = rhs.mat;

   return *this;
}

SparseMatrix SparseMatrix::operator+(const SparseMatrix& rhs) {

   SparseMatrix result(rows, cols);
   
   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for addition.");
   }else{
      for(int i=0; i<rows; i++){
         result.mat[i] = mat[i];
         result.mat[i].insert(result.mat[i].end(), rhs.mat[i].begin(), rhs.mat[i].end());
      }
   }
   result.cleanup();
   return result;
}

SparseMatrix& SparseMatrix::operator+=(const SparseMatrix& rhs) {

   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for addition.");
   }else{
      for(int i=0; i<rows; i++)
         mat[i].insert(mat[i].end(), rhs.mat[i].begin(), rhs.mat[i].end());
   }
   cleanup();
   return *this;
}

SparseMatrix SparseMatrix::operator-(const SparseMatrix& rhs) {

   SparseMatrix result(rows, cols);

   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for subtraction.");
   }else{
      for(int i=0; i<rows; i++){
         result.mat[i] = mat[i];
         for(int j=0; j<rhs.mat[i].size(); j++)
            result.mat[i].push_back(node{-rhs.mat[i][j].val,rhs.mat[i][j].col});
      }
   }
   result.cleanup();
   return result;
}

SparseMatrix& SparseMatrix::operator-=(const SparseMatrix& rhs) {

   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for subtraction.");
   }else{
      for(int i=0; i<rows; i++){
         for(int j=0; j<rhs.mat[i].size(); j++)
            mat[i].push_back(node{-rhs.mat[i][j].val,rhs.mat[i][j].col});
      }
   }
   cleanup();
   return *this;
}

SparseMatrix SparseMatrix::transpose() const {
   SparseMatrix result(cols, rows);
   
   for(int i=0; i<rows; i++)
      for(int j=0; j<mat[i].size(); j++)
         result.mat[mat[i][j].col].push_back(node{mat[i][j].val,i});

   result.cleanup();
   return result;
}

//Access the individual elements
double SparseMatrix::operator()(const int& i, const int& j) {

   double result = 0;

   for(int k=0; k<mat[i].size(); k++){
      if(mat[i][k].col == j)
         result += mat[i][k].val;
   }

   return result;
}

double SparseMatrix::sum() {
   double result = 0;
   for(int i=0; i<rows; i++)
      result += std::accumulate(mat[i].begin(),mat[i].end(), node{0.0,0}).val;
   return result;
}
double SparseMatrix::norm() {
   double result = 0;
   for(int i=0; i<rows; i++)
      for(int j=0; j<mat[i].size(); j++)
         result += mat[i][j].val*mat[i][j].val;
   return sqrt(result);
}

/*Matrix SparseMatrix::Smallest_NonConst_Eigenvectors(int k, double tol){

   double l = LargestEigenvalue(tol);
   SparseMatrix M = *this - l*speye(rows);

   //Matrix v = Matrix::rand(cols,1);
   //Matrix o = Matrix::ones(1,cols)/std::sqrt(cols);
   //Matrix ot = Matrix::ones(cols,1)/std::sqrt(cols);
   //v = v - (o*v)(0,0)*ot;  //Projection step

   Matrix u = Matrix::ones(cols,1)/std::sqrt(cols);

   //Matrix u = M.Largest_NonConst_Eigenvector(v,tol);
   //u = sign(u);
   //u = u - u.mean();
   //u = u/u.norm();
   for(int i=1; i<k; i++){
      v = Matrix::rand(cols,1);
      v = v - (o*v)(0,0)*ot;  //Projection step
      v = v - u*(u.transpose()*v);
      Matrix w = M.Largest_NonConst_Eigenvector(v,tol);
      w = sign(w);
      w = w - w.mean();
      w = w/w.norm();
      u.append_col(w);
   }
   
   return u;
}*/
//Power method to find smallest eigenvector
/*std::vector<double> SparseMatrix::Smallest_NonConst_Eigenvector(double tol) {

   double l = LargestEigenvalue(tol);
   SparseMatrix M = *this - l*speye(rows);

   //Matrix v = Matrix::rand(cols,1);
   //Matrix o = Matrix::ones(1,cols)/std::sqrt(cols);
   Matrix u = Matrix::ones(cols,1)/std::sqrt(cols);
   //v = v - (o*v)(0,0)*ot;  //Projection step
   return  M.Largest_Eigenvector(u,tol).tovector();
}*/
//Power method to find largest eigenvector
Matrix SparseMatrix::Largest_Eigenvector(Matrix u, double tol) {
   double l=0.0;
   double res = tol + 1.0;
   int i = 0;
   Matrix v = Matrix::rand(rows,1);
   while(res > tol && i < 1E5){
      v = v - u*(u.transpose()*v);
      v = v/v.norm();
      Matrix w = (*this)*v;
      l = (v.transpose()*w)(0,0);
      res = (w - l*v).norm();
      printf("i:%d,res=%f\n",i,res);
      v = w/w.norm();
      i++;
   }
   printf("%d iterations for tol=%f\n",i,tol);
   return v;
}
//Power method to find largest eigenvalue
/*double SparseMatrix::Smallest_NonConst_Eigenvalue(double tol) {

   double l = LargestEigenvalue(tol);
   SparseMatrix M = *this - l*speye(rows);
   double s = M.Largest_NonConst_Eigenvalue(tol);

   return s + l;
}
//Power method to find largest eigenvalue
double SparseMatrix::Largest_NonConst_Eigenvalue(double tol) {
   double l=0.0;
   Matrix v = Matrix::rand(cols,1);
   Matrix o = Matrix::ones(1,cols)/std::sqrt(cols);
   Matrix ot = Matrix::ones(cols,1)/std::sqrt(cols);
   v = v - (o*v)(0,0)*ot;  //Projection step
   double res = tol + 1.0;
   int i = 0;
   while(res > tol){
      Matrix w = (*this)*v;
      l = (v.transpose()*w)(0,0);
      res = (w - l*v).norm();
      v = w/w.norm();
      i++;
   }
   printf("%d iterations for tol=%f\n",i,tol);
   printf("l=%f\n",l);
   return l;
}*/

//Power method to find largest eigenvalue
double SparseMatrix::SmallestEigenvalue(double tol) {

   double l = LargestEigenvalue(tol);
   SparseMatrix M = *this - l*speye(rows);
   double s = M.LargestEigenvalue(tol);

   return s + l;
}

//Power method to find largest eigenvalue
double SparseMatrix::LargestEigenvalue(double tol) {
   double l=0.0;
   Matrix v = Matrix::rand(cols,1);
   double res = tol + 1.0;
   int i = 0;
   while(res > tol){
      Matrix w = (*this)*v;
      l = (v.transpose()*w)(0,0);
      res = (w - l*v).norm();
      v = w/w.norm();
      i++;
   }
   printf("%d iterations for tol=%f\n",i,tol);
   printf("l=%f\n",l);
   return l;
}

std::vector<double> SparseMatrix::sum(int axis) {

   std::vector<double> result;
   SparseMatrix tran(rows, cols);
   switch(axis)
   {
      case 0:
         tran = transpose();
         result.resize(cols,0);
         for(int i=0; i<cols; i++)
            result[i] = std::accumulate(tran.mat[i].begin(),tran.mat[i].end(), node{0.0,0}).val;
         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++)
            result[i] = std::accumulate(mat[i].begin(),mat[i].end(), node{0.0,0}).val;
         break;
      default:
         throw std::invalid_argument("Invalid axis for summation.");
   }
   
   return result; 

}

SparseMatrix::SparseMatrix() {
   rows = 0;
   cols = 0;
}
SparseMatrix::SparseMatrix(const std::vector<double>& diags) {
   rows = diags.size();
   cols = diags.size();
   mat.resize(rows);
   for(int i=0; i<rows; i++)
      if(diags[i] != 0)
         mat[i].push_back(node{diags[i],i});
}
SparseMatrix::SparseMatrix(int _rows) {
   rows = _rows;
   cols = _rows;
   mat.resize(rows);
   for(int i=0; i<rows; i++)
      mat[i].push_back(node{1.0,i});
}
SparseMatrix::SparseMatrix(int _rows, int _cols) {
   rows = _rows;
   cols = _cols;
   mat.resize(rows);
}

//Copy constructor
SparseMatrix::SparseMatrix(const SparseMatrix& rhs) {
   mat = rhs.mat;
   rows = rhs.get_rows();
   cols = rhs.get_cols();
}

int SparseMatrix::num_nonzero() {
   
   int num = 0;
   for(int i=0; i<rows; i++)
      num += mat[i].size();
   return num;
}

//Destructor
SparseMatrix::~SparseMatrix() {}

//Get the number of rows of the matrix
int SparseMatrix::get_rows() const {
   return rows;
}

//Insert entry to sparse matrix
void SparseMatrix::insert(int i, int j, const double& _val){

   if(i >= rows || i < 0 || j >= cols || j < 0){
      throw std::invalid_argument("Invalid index to sparse matrix.");
   }else{
      node a = {_val,j};
      mat[i].push_back(a);
   }
}

//Cleanup sparse matrix (combines multiple entries and removes zeros)
void SparseMatrix::cleanup(){
   for(int i=0; i<rows; i++){
      std::sort(mat[i].begin(),mat[i].end());
      
      std::vector<node> new_row;
      int col = -1;
      node a = {0.0,col};
      for(int j=0; j < mat[i].size(); j++){
         if(col < mat[i][j].col){
            col = mat[i][j].col;
            a.val = mat[i][j].val;
            a.col = col;
            new_row.push_back(a);
         }else{
            new_row.back().val+=mat[i][j].val;
         }
      }
      //Now remove any zero entries
      std::vector<node> newer_row;
      for(int j=0; j< new_row.size(); j++)
         if(new_row[j].val != 0)
            newer_row.push_back(new_row[j]);
      mat[i] = newer_row;
   }
}

//Print a sparse matrix
void SparseMatrix::print(){
   std::cout << std::endl;
   for(int i=0; i<rows; i++){
      for(int j=0; j<mat[i].size(); j++)
         std::cout << "(" << i << "," << mat[i][j].col << "): " << mat[i][j].val << std::endl;
      if(mat[i].size() > 0)
         std::cout << std::endl;
   }
}

//Get the number of columns of the matrix
int SparseMatrix::get_cols() const {
   return cols;
}



Matrix::Matrix(int _rows, int _cols, const double& min, const double& max) {
   std::mt19937_64 rng;
   uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
   std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
   rng.seed(ss);
   std::uniform_real_distribution<double> unif(0, 1);

   rows = _rows;
   cols = _cols;
   data.resize(rows*cols);

   mat.resize(cols);
   mat[0] = &data[0];
   for(int i=1; i<cols; i++)
      mat[i] = mat[i-1] + rows;

   for(int i=0; i<rows*cols; i++)
      data[i] = min + (max-min)*unif(rng);
}

void Matrix::append_col(const Matrix& w){
   cols++;
   data.resize(rows*cols);
   mat.resize(cols);
   mat[0] = &data[0];
   for(int i=1; i<cols; i++)
      mat[i] = mat[i-1] + rows;

   for(int i=0; i<rows; i++)
      mat[cols-1][i] = w(i,0);
}
Matrix::Matrix(int _rows, int _cols, const double& initial) {
   rows = _rows;
   cols = _cols;
   data.resize(rows*cols,initial);

   mat.resize(cols);
   mat[0] = &data[0];
   for(int i=1; i<cols; i++)
      mat[i] = mat[i-1] + rows;
}

Matrix::Matrix() {
   rows = 0;
   cols = 0;
}

Matrix::Matrix(int _rows) {
   rows = _rows;
   cols = _rows;
   data.resize(rows*cols,0);

   mat.resize(cols);
   mat[0] = &data[0];
   for(int i=1; i<cols; i++)
      mat[i] = mat[i-1] + rows;

   for(int i=0; i<rows; i++)
      mat[i][i] = 1;
}

Matrix::Matrix(int _rows, int _cols) {
   rows = _rows;
   cols = _cols;
   data.resize(rows*cols);

   mat.resize(cols);
   mat[0] = &data[0];
   for(int i=1; i<cols; i++)
      mat[i] = mat[i-1] + rows;
}
//Vector copy constructor
Matrix::Matrix(const std::vector<double>& rhs) {
   data = rhs;
   rows = rhs.size();
   cols = 1;
   mat.resize(cols);
   mat[0] = &data[0];
}
//Copy constructor
Matrix::Matrix(const Matrix& rhs) {
   data = rhs.data;
   rows = rhs.get_rows();
   cols = rhs.get_cols();

   mat.resize(cols);
   mat[0] = &data[0];
   for(int i=1; i<cols; i++)
      mat[i] = mat[i-1] + rows;
}

//Destructor
Matrix::~Matrix() {}

double Matrix::norm() {
   double result = 0;
   for(int i=0; i<cols*rows; i++)
      result += data[i]*data[i];

   return sqrt(result);
}
Matrix Matrix::Heaviside(){
   Matrix result(rows, cols, 1.0);
   for(int i=0; i<rows; i++) 
      for(int j=0; j<cols; j++)
         if(mat[j][i] < 0)
            result.mat[j][i] = 0.0;
   return result;
}
Matrix Matrix::sign(){
   Matrix result(rows, cols, 1.0);
   for(int i=0; i<rows; i++) 
      for(int j=0; j<cols; j++)
         if(mat[j][i] < 0)
            result.mat[j][i] = -1.0;
   return result;
}
Matrix Matrix::divide(Matrix& rhs){
   Matrix result(rows, cols);
   for(int i=0; i<rows; i++) 
      for(int j=0; j<cols; j++) 
         result.mat[j][i] = mat[j][i]/rhs.mat[j][i];
   return result;
}
Matrix Matrix::multiply(Matrix& rhs){
   Matrix result(rows, cols);
   for(int i=0; i<rows; i++) 
      for(int j=0; j<cols; j++) 
         result.mat[j][i] = mat[j][i]*rhs.mat[j][i];
   return result;
}
Matrix Matrix::col_multiply(const std::vector<double>& rhs){
   Matrix result(rows, cols);
   for(int i=0; i<rows; i++) 
      for(int j=0; j<cols; j++) 
         result.mat[j][i] = mat[j][i]*rhs[j];
   return result;
}
Matrix Matrix::row_multiply(const std::vector<double>& rhs){
   Matrix result(rows, cols);
   for(int i=0; i<rows; i++) 
      for(int j=0; j<cols; j++) 
         result.mat[j][i] = mat[j][i]*rhs[i];
   return result;
}
Matrix Matrix::row_divide(const std::vector<double>& rhs){
   Matrix result(rows, cols);
   for(int i=0; i<rows; i++) 
      for(int j=0; j<cols; j++) 
         result.mat[j][i] = mat[j][i]/rhs[i];
   return result;
}
std::vector<double> Matrix::norm(int axis) {

   std::vector<double> result;
   switch(axis)
   {
      case 0:
         result.resize(cols,0);
         for(int i=0; i<cols; i++){
            for(int j=0; j<rows; j++)
               result[i] += mat[i][j]*mat[i][j];
            result[i] = sqrt(result[i]);
         }
         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++)
               result[i] += mat[j][i]*mat[j][i];
            result[i] = sqrt(result[i]);
         }
         break;
      default:
         throw std::invalid_argument("Invalid axis for summation.");
   }
   
   return result; 

}
double Matrix::max() {
   double result = data[0];
   for(int i=0; i<rows*cols; i++)
      result = std::max(result,data[i]);
   return result;
}
double Matrix::min() {
   double result = data[0];
   for(int i=0; i<rows*cols; i++)
      result = std::min(result,data[i]);
   return result;
}
std::vector<double> Matrix::max(int axis) {

   std::vector<double> result;
   switch(axis)
   {
      case 0:
         result.resize(cols,0);
         for(int i=0; i<cols; i++){
            result[i] = mat[i][0];
            for(int j=0; j<rows; j++)
               result[i] = std::max(result[i],mat[i][j]);
         }

         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++){
            result[i] = mat[0][i];
            for(int j=0; j<cols; j++)
               result[i] = std::max(result[i],mat[j][i]);
         }

         break;
      default:
         throw std::invalid_argument("Invalid axis for max.");
   }
   
   return result; 

}
std::vector<double> Matrix::min(int axis) {

   std::vector<double> result;
   switch(axis)
   {
      case 0:
         result.resize(cols,0);
         for(int i=0; i<cols; i++){
            result[i] = mat[i][0];
            for(int j=0; j<rows; j++)
               result[i] = std::min(result[i],mat[i][j]);
         }

         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++){
            result[i] = mat[0][i];
            for(int j=0; j<cols; j++)
               result[i] = std::min(result[i],mat[j][i]);
         }

         break;
      default:
         throw std::invalid_argument("Invalid axis for min.");
   }
   
   return result; 

}
std::vector<int> SparseMatrix::argmin(int axis) {

   std::vector<int> result;
   switch(axis)
   {
      case 0:
         throw std::invalid_argument("Invalid axis for argmin.");
         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++){
            result[i] = 0;
            for(int j=0; j<mat[i].size(); j++){
               if(mat[i][j] < mat[i][result[i]])
                  result[i] = j;
            }
            if(!mat[i].empty())
               result[i] = mat[i][result[i]].col;
         }

         break;
      default:
         throw std::invalid_argument("Invalid axis for argmin.");
   }
   
   return result; 

}
std::vector<int> Matrix::argmin(int axis) {

   std::vector<int> result;
   switch(axis)
   {
      case 0:
         result.resize(cols,0);
         for(int i=0; i<cols; i++){
            result[i] = 0;
            for(int j=0; j<rows; j++){
               if(mat[i][j] < mat[i][result[i]])
                  result[i] = j;
            }
         }

         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++){
            result[i] = 0;
            for(int j=0; j<cols; j++){
               if(mat[j][i] < mat[result[i]][i])
                  result[i] = j;
            }
         }

         break;
      default:
         throw std::invalid_argument("Invalid axis for min.");
   }
   
   return result; 

}
std::vector<int> Matrix::argmax(int axis) {

   std::vector<int> result;
   switch(axis)
   {
      case 0:
         result.resize(cols,0);
         for(int i=0; i<cols; i++){
            result[i] = 0;
            for(int j=0; j<rows; j++){
               if(mat[i][j] > mat[i][result[i]])
                  result[i] = j;
            }
         }

         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++){
            result[i] = 0;
            for(int j=0; j<cols; j++){
               if(mat[j][i] > mat[result[i]][i])
                  result[i] = j;
            }
         }

         break;
      default:
         throw std::invalid_argument("Invalid axis for max.");
   }
   
   return result; 

}
double Matrix::sum() {
   double result = std::accumulate(data.begin(),data.end(), 0.0);
   return result;
}

std::vector<double> Matrix::sum(int axis) {

   std::vector<double> result;
   switch(axis)
   {
      case 0:
         result.resize(cols,0);
         for(int i=0; i<cols; i++)
            for(int j=0; j<rows; j++)
               result[i] += mat[i][j];

         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++)
            for(int j=0; j<cols; j++)
               result[i] += mat[j][i];

         break;
      default:
         throw std::invalid_argument("Invalid axis for summation.");
   }
   
   return result; 

}

double Matrix::mean() {
   double result = std::accumulate(data.begin(),data.end(), 0.0);
   return result/double(rows*cols);
}


std::vector<double> Matrix::mean(int axis) {

   std::vector<double> result;
   switch(axis)
   {
      case 0:
         result.resize(cols,0);
         for(int i=0; i<cols; i++){
            for(int j=0; j<rows; j++)
               result[i] += mat[i][j];
            result[i]/=double(rows);
         }

         break;
      case 1:
         result.resize(rows,0);
         for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++)
               result[i] += mat[j][i];
            result[i]/=double(cols);
         }

         break;
      default:
         throw std::invalid_argument("Invalid axis for summation.");
   }
   
   return result; 

}

void Matrix::print() {

   //std::cout << std::setprecision(2);
   std::cout << "[";
   int c = 0;
   for (int i=0; i < rows; i++) {
      if(i==0)
         std::cout << "[";
      else
         std::cout << " [";
      for (int j=0; j < cols-1; j++) {
         std::cout << mat[j][i] << ", ";
      }
      std::cout << mat[cols-1][i] << "]";
      if(i == rows-1)
         std::cout << "]" << std::endl;
      else
         std::cout << "," << std::endl;
   }
}

//Assignment operator
Matrix& Matrix::operator=(const Matrix& rhs) {

   if (&rhs == this)
      return *this;

   rows = rhs.get_rows();
   cols = rhs.get_cols();
   data = rhs.data;

   mat.resize(cols);
   mat[0] = &data[0];
   for(int i=1; i<cols; i++)
      mat[i] = mat[i-1] + rows;

   return *this;
}

//Addition of two matrices
Matrix Matrix::operator+(const Matrix& rhs) {

   Matrix result(rows, cols);
   
   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for addition.");
   }else{
      for(int i=0; i<rows*cols; i++) 
         result.data[i] = data[i] + rhs.data[i];
   }
   return result;
}

//Cumulative addition of this matrix and another
Matrix& Matrix::operator+=(const Matrix& rhs) {

   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for addition.");
   }else{
      for(int i=0; i<rows*cols; i++) 
         data[i] += rhs.data[i];
   }
   return *this;
}
//Addition of this matrix and a vector (adds to rows)
Matrix Matrix::operator+(const std::vector<double>& rhs) {

   Matrix result(rows, cols);

   if(rhs.size() != cols){
      throw std::length_error("Matrix and vector have incompatible sizes for addition.");
   }else{
      for(int i=0; i<rows; i++) 
         for(int j=0; j<cols; j++) 
            result.mat[j][i] = mat[j][i] + rhs[j];
   }
   return result;
}

//Cumulative addition of this matrix and a vector (adds to rows)
Matrix& Matrix::operator+=(const std::vector<double>& rhs) {

   if(rhs.size() != cols){
      throw std::length_error("Matrix and vector have incompatible sizes for addition.");
   }else{
      for(int i=0; i<rows; i++) 
         for(int j=0; j<cols; j++) 
            mat[j][i] += rhs[j];
   }
   return *this;
}
//Subtraction of this matrix and a vector (subtracts from rows)
Matrix Matrix::operator-(const std::vector<double>& rhs) {

   Matrix result(rows, cols);

   if(rhs.size() != cols){
      throw std::length_error("Matrix and vector have incompatible sizes for subtraction.");
   }else{
      for(int i=0; i<rows; i++) 
         for(int j=0; j<cols; j++) 
            result.mat[j][i] = mat[j][i] - rhs[j];
   }
   return result;
}

//Cumulative subtraction of this matrix and a vector (subtracts from rows)
Matrix& Matrix::operator-=(const std::vector<double>& rhs) {

   if(rhs.size() != cols){
      throw std::length_error("Matrix and vector have incompatible sizes for subtraction.");
   }else{
      for(int i=0; i<rows; i++) 
         for(int j=0; j<cols; j++) 
            mat[j][i] -= rhs[j];
   }
   return *this;
}
//Subtraction of this matrix and another
Matrix Matrix::operator-(const Matrix& rhs) {

   Matrix result(rows, cols);

   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for subtraction.");
   }else{
      for(int i=0; i<rows*cols; i++) 
         result.data[i] = data[i] - rhs.data[i];
   }
   return result;
}

//Cumulative subtraction of this matrix and another
Matrix& Matrix::operator-=(const Matrix& rhs) {

   if(rhs.get_rows() != rows || rhs.get_cols() != cols){
      throw std::length_error("Matrices have incompatible sizes for subtraction.");
   }else{
      for(int i=0; i<rows*cols; i++) 
         data[i] -= rhs.data[i];
   }
   return *this;
}

//Left multiplication of this matrix and another
Matrix Matrix::operator*(const Matrix& rhs) {

   Matrix result(rows, rhs.get_cols(), 0);

   if(rhs.get_rows() != cols){
      throw std::length_error("Matrices have incompatible sizes for multiplication.");
   }else{
      Matrix tran = this->transpose();
      for(int i=0; i<rows; i++) 
         for(int j=0; j<rhs.get_cols(); j++) 
            for(int k=0; k<cols; k++) 
               result.mat[j][i] += tran.mat[i][k] * rhs.mat[j][k];
   }
   return result;
}

//Cumulative left multiplication of this matrix and another
Matrix& Matrix::operator*=(const Matrix& rhs) {
   Matrix result = (*this) * rhs;
   (*this) = result;
   return *this;
}

//Calculuate a transpose of this matrix
Matrix Matrix::transpose() {
   Matrix result(cols, rows, 0.0);
   
   for(int i=0; i<cols; i++) 
      for(int j=0; j<rows; j++) 
         result.mat[j][i] = mat[i][j];

   return result;
}

//Matrix/scalar addition
Matrix Matrix::operator+(const double& rhs) {
   Matrix result(rows, cols);

   for(int i=0; i<rows*cols; i++)
      result.data[i] = data[i] + rhs;

   return result;
}

//Matrix/scalar subtraction
Matrix Matrix::operator-(const double& rhs) {
   Matrix result(rows, cols);

   for(int i=0; i<rows*cols; i++)
      result.data[i] = data[i] - rhs;

   return result;
}

//Matrix/scalar multiplication
Matrix Matrix::operator*(const double& rhs) {
   Matrix result(rows, cols);

   for(int i=0; i<rows*cols; i++)
      result.data[i] = data[i] * rhs;

   return result;
}

//Matrix/scalar division
Matrix Matrix::operator/(const double& rhs) {
   Matrix result(rows, cols);

   for(int i=0; i<rows*cols; i++)
      result.data[i] = data[i] / rhs;

   return result;
}

//Multiply a matrix with a vector<double>
std::vector<double> Matrix::operator*(const std::vector<double>& rhs) {
   std::vector<double> result(rows,0.0);

   if(rhs.size() != cols){
      throw std::length_error("Matrix and vector have incompatible sizes for multiplication.");
   }else{
      for(int i=0; i<rows; i++) 
         for(int j=0; j<cols; j++) 
            result[i] += mat[j][i] * rhs[j];
   }

   return result;
}
//Convert to flattened vector
std::vector<double> Matrix::tovector() {
   std::vector<double> result = data;
   return result;
}
//Obtain a vector<double> of the diagonal elements
std::vector<double> Matrix::diag() {

   std::vector<double> result(rows);
   if(rows != cols){
      throw std::length_error("Requested diagonal of non-square matrix.");
   }else{
      for(int i=0; i<rows; i++) 
         result[i] = mat[i][i];
   }
   return result;
}

//Access the individual elements
double& Matrix::operator()(const int& i, const int& j) {
   return mat[j][i];
}

//Access the individual elements
const double& Matrix::operator()(const int& i, const int& j) const {
   return mat[j][i];
}

//Get the number of rows of the matrix
int Matrix::get_rows() const {
   return rows;
}

std::vector<int> Matrix::shape() const {
   return std::vector<int>{rows,cols};
}

//Get the number of columns of the matrix
int Matrix::get_cols() const {
   return cols;
}

