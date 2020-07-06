/****************************************************************************
 * MeshLab                                                           o o     *
 * A versatile mesh processing toolbox                             o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2005                                                \/)\/    *
 * Visual Computing Lab                                            /\/|      *
 * ISTI - Italian National Research Council                           |      *
 *                                                                    \      *
 * All rights reserved.                                                      *
 *                                                                           *
 * This program is free software; you can redistribute it and/or modify      *
 * it under the terms of the GNU General Public License as published by      *
 * the Free Software Foundation; either version 2 of the License, or         *
 * (at your option) any later version.                                       *
 *                                                                           *
 * This program is distributed in the hope that it will be useful,           *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 * GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
 * for more details.                                                         *
 *                                                                           *
 ****************************************************************************/


#include "stdio.h"
#include "stdbool.h"
#include <math.h>
#include <stdlib.h>
#include <vcg/space/colorspace.h>
#include "spherical_volume_invariant.h"
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/point_outlier.h>

#include <vcg/complex/algorithms/point_outlier.h>

// HERE IS THE EIGENSOLVER!!!!!! :) :) :)
#include <eigenlib/Eigen/Eigenvalues>
#include <complex>

#include <QOpenGLContext>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/normal.h>


#include <iostream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>

using namespace std;
using namespace vcg;

// ERROR CHECKING UTILITY
#define CheckError(x,y); if ((x)) {this->errorMessage = (y); return false;}
///////////////////////////////////////////////////////




// huge swath of SVIPCA code (ends line 1057):

/*Compute minimum distance from edge xy to origin*/
static inline double edge_dist(double *x, double *y){

   double p[3], t;

   sub(y,x,p);
   t = dot(p,y)/norm_squared(p);
   t = MAX(MIN(t,1),0);
   mult(p,t,p);
   return dist_squared(y,p);

}

/*Compute unit outward normal of triangle xyz
 * Also returns triangle area*/
static inline double tri_normal(double *x, double *y, double *z, double *nu){
   
   double p[3], q[3], nnu, t;

   sub(y,x,p);
   sub(z,x,q);
   cross(p,q,nu);
   nnu = norm(nu);
   t = 1/nnu;mult(nu,t,nu);
   return nnu/2;
}

/*Gram-Schmidt othornoralization:
 * Returns normal vector in direction x - proj_y x
 * x and y assumed to be unit normal vectors
 */
static inline void gram_schmidt(double *x, double *y, double *p){

   double a, b;
   b = dot(x,y);
   mult(y,b,p);
   sub(x,p,p);
   a = 1/norm(p);
   mult(p,a,p);
}

/*Compute minimum distance from triangle xyz to origin
 * Works by projecting origin onto triangle
 * x,y,z = vertices of triangle
 * Returns distance to origin
 * Also returns angle theta used in exact integration
 */
double tri_dist(double *x, double *y, double *z){

   double p[3], nxy[3], nyz[3], nzx[3], pxy[3], pyz[3], pzx[3], xs[3], nu[3], txy, tyz, tzx, t;

   sub(y,x,pxy); t=1/norm(pxy); mult(pxy,t,pxy);
   sub(z,y,pyz); t=1/norm(pyz); mult(pyz,t,pyz);
   sub(x,z,pzx); t=1/norm(pzx); mult(pzx,t,pzx);

   gram_schmidt(pyz,pxy,nxy);  //Normal to xy
   gram_schmidt(pzx,pyz,nyz);  //Normal to yz
   gram_schmidt(pxy,pzx,nzx);  //Normal to zx
   
   txy = dot(x,nxy);
   tyz = dot(y,nyz);
   tzx = dot(z,nzx);

   //Compute projection xs of origin (0,0,0) onto triangle xyz
   if(txy > 0){ //projection to plane lies outside of edge xy
      t = MIN(MAX(-dot(x,pxy),0),dist(x,y));
      mult(pxy,t,p);
      add(x,p,xs);
   }else if(tyz > 0){ //projection to plane lies outside of edge yz
      t = MIN(MAX(-dot(y,pyz),0),dist(y,z));
      mult(pyz,t,p);
      add(y,p,xs);
   }else if(tzx > 0){ //projection to plane lies outside of edge zx
      t = MIN(MAX(-dot(z,pzx),0),dist(z,x));
      mult(pzx,t,p);
      add(z,p,xs);
   }else{ //Projection lies inside triangle
      tri_normal(x,y,z,nu);
      t = dot(x,nu);
      mult(nu,t,xs); //projection of origin to triangle plane
   }

   return norm_squared(xs);
}

/*Returns for each vertex the indices of adjacent triangles
 * n = number vertices
 * m = number of triangles
 * T = triangle connectivity list
 * Returns nxk array L[i][j] = index of jth triangle adjacent to vertex i
 *                   L[i][0] = number of adjacent triangles
 */
int** vertex_to_triangle_list(int **T, int n, int m){
  
   int **L;
   int i,j,k,l,p,q;
   int *C = vector_int(n,0);
   for(i=0;i<m;i++){
      C[T[i][0]]++; C[T[i][1]]++; C[T[i][2]]++;
   }
   k = C[0];
   for(i=0;i<n;i++)
      k = MAX(k,C[i]);
   L = array_int(n,k+1,-1);
   for(i=0;i<n;i++)
      L[i][0] = C[i];

   for(i=0;i<m;i++){
      for(j=0;j<3;j++){
         p = T[i][j];
         bool addi = true;
         q = 1;
         for(l=1;l<=L[p][0];l++){
            addi = addi && L[p][l] != i;
            if(L[p][l] != -1)
               q++;
         }
         if(addi)
            L[p][q] = i;
      }
   }

   return L;
}

/*Returns the face graph (triangle adjacency matrix)
 * n = number vertices
 * m = number of triangles,
 * T = triangle connectivity list
 * L = vertex to triangle list
 * Returns mx3 array F[i][j] = index of jth triangle adjacent to triangle i
 * Also updates B to encode mesh boundary points (B[p]=false iff p boundary vertex)
 */
int** face_graph(int **T, int **L, bool *B, int n, int m){
  
   int i,j,k,l,l2,p,q;
   
   //Allocate memory
   int **F = array_int(m,3,-1);

   for(i=0;i<m;i++){
      for(l=0;l<3;l++){
         l2 = (l+1)%3;
         p = T[i][l]; q = T[i][l2];
         for(j=1;j<=L[p][0];j++){
            for(k=1;k<=L[q][0];k++){
               if(L[p][j] == L[q][k] && L[p][j] != i)
                  F[i][l] = L[p][j];
            }
         }
         if(F[i][l] == -1){
            //printf("boundary!!\n");
            B[p] = false; B[q] = false;
         }
      }
   }
  
   return F;

}
/*Non-recursive depth first search to find local patch for integration
 * Replaces expensive range-search
 * T = triangles
 * P = vertices
 * F = face graph adjacency
 * NN = output neighbor indices
 * v = vector of visited triangle indices
 * b = output boundary triangles
 * num = current triangle number
 * ind = index of root triangle
 * r2 = r^2
 * i = current index
 * Returns true if B(x,r) does not interesect boundary of mesh, and false otherwise
 */
bool depth_first_search_nr(int **T, double **P, int **F, int *NN, bool *v, bool *b, int *num, int ind, double r2, int i, int *stack){

   int j,p,q,k,t,stack_size;
   bool bi = true;
   double x[3], y[3], z[3], rmin, rmax;

   stack[0] = ind;
   stack_size = 1;  //Next position in stack

   while(stack_size > 0){

      stack_size--;
      ind = stack[stack_size];
      if(!(v[ind])){ //If not already visited

         //Extract vertices
         tri_vertices(T,P,x,y,z,ind,i);

         //Compute distance from triangle to origin
         rmin = tri_dist(x,y,z);
         rmax = MAX3(norm_squared(x),norm_squared(y),norm_squared(z));
         
         //Check if within B(x,r)
         if(rmin < r2){
            v[ind] = true;
            NN[*num] = ind;
            if(rmax > r2)
               b[*num] = true;
            else
               b[*num] = false;
            ++*num;

            //Add neighbors to stack
            for(j=0;j<3;j++){
               t = F[ind][j];
               if(t == -1){
                  k = (j+1)%3;
                  p = T[ind][j]; q = T[ind][k];
                  x[0] = P[p][0] - P[i][0]; x[1] = P[p][1] - P[i][1]; x[2] = P[p][2] - P[i][2];
                  y[0] = P[q][0] - P[i][0]; y[1] = P[q][1] - P[i][1]; y[2] = P[q][2] - P[i][2];
                  bi = bi && (edge_dist(x,y) >= r2);
               }else if(!(v[t])){
                  stack[stack_size] = t;
                  stack_size++;
               }
            }
         }
      }
   }
   return bi;
}
/*Non-recursive breadth first search to find local patch for integration
 * Replaces expensive range-search
 * T = triangles
 * P = vertices
 * F = face graph adjacency
 * NN = output neighbor indices
 * v = vector of visited triangle indices
 * b = output boundary triangles
 * num = current triangle number
 * ind = index of root triangle
 * r2 = r^2
 * i = current index
 * Returns true if B(x,r) does not interesect boundary of mesh, and false otherwise
 */
bool breadth_first_search(int **T, double **P, int **F, int *NN, bool *v, bool *b, int *num, int ind, double r2, int i, int *stack){

   int j,p,q,k,t,stack_next,stack_curr;
   bool bi = true;
   double x[3], y[3], z[3], rmin, rmax;

   stack[0] = ind;
   stack_next = 1;  //Next position in stack
   stack_curr = 0;  //Current position in stack

   while(stack_curr < stack_next){

      ind = stack[stack_curr];
      stack_curr++;
      if(!(v[ind])){ //If not already visited

         //Extract vertices
         tri_vertices(T,P,x,y,z,ind,i);

         //Compute distance from triangle to origin
         rmin = tri_dist(x,y,z);
         rmax = MAX3(norm_squared(x),norm_squared(y),norm_squared(z));
         
         //Check if within B(x,r)
         if(rmin < r2){
            v[ind] = true;
            NN[*num] = ind;
            if(rmax > r2)
               b[*num] = true;
            else
               b[*num] = false;
            ++*num;

            //Add neighbors to stack
            for(j=0;j<3;j++){
               t = F[ind][j];
               if(t == -1){
                  k = (j+1)%3;
                  p = T[ind][j]; q = T[ind][k];
                  x[0] = P[p][0] - P[i][0]; x[1] = P[p][1] - P[i][1]; x[2] = P[p][2] - P[i][2];
                  y[0] = P[q][0] - P[i][0]; y[1] = P[q][1] - P[i][1]; y[2] = P[q][2] - P[i][2];
                  bi = bi && (edge_dist(x,y) >= r2);
               }else if(!(v[t])){
                  stack[stack_next] = t;
                  stack_next++;
               }
            }
         }
      }
   }
   return bi;
}
/*Depth first search to find local patch for integration
 * Replaces expensive range-search
 * T = triangles
 * P = vertices
 * F = face graph adjacency
 * NN = output neighbor indices
 * v = vector of visited triangle indices
 * b = output boundary triangles
 * num = current triangle number
 * ind = index of current triangle
 * r2 = r^2
 * i = current index
 * Returns true if B(x,r) does not interesect boundary of mesh, and false otherwise
 */
bool depth_first_search(int **T, double **P, int **F, int *NN, bool *v, bool *b, int *num, int ind, double r2, int i){

   int j,p,q,k,t;
   bool bi = true;
   double x[3], y[3], z[3], rmin, rmax;

   if(!(v[ind])){

      //Extract vertices
      tri_vertices(T,P,x,y,z,ind,i);

      //Compute distance from triangle to origin
      rmin = tri_dist(x,y,z);
      rmax = MAX3(norm_squared(x),norm_squared(y),norm_squared(z));

      //Check if within B(x,r)
      if(rmin < r2){
         v[ind] = true;
         NN[*num] = ind;
         if(rmax > r2)
            b[*num] = true;
         else
            b[*num] = false;
         ++*num;

         //Call depth first search at neighbors
         for(j=0;j<3;j++){
            t = F[ind][j];
            if(t == -1){
               k = (j+1)%3;
               p = T[ind][j]; q = T[ind][k];
               x[0] = P[p][0] - P[i][0]; x[1] = P[p][1] - P[i][1]; x[2] = P[p][2] - P[i][2];
               y[0] = P[q][0] - P[i][0]; y[1] = P[q][1] - P[i][1]; y[2] = P[q][2] - P[i][2];
               bi = bi && (edge_dist(x,y) >= r2);
            }else{
               bi = bi && depth_first_search(T,P,F,NN,v,b,num,t,r2,i);
            }
         }
      }
   }

   return bi;
}

/*Numerical integration via triangle refinement. Used at boundary triangles
 * x,y,z = vertices of triangle
 * dx,dy,dz = squared norms of x,y,z
 * A = area of triangle
 * eta = x dot triangle unit normal (same as y or z)
 * eps = tolerance
 * r,r2,r3 = r,r^2,r^3
 * level = current triangle refinement level
 * maxlevel = maximum refinement level
 * num_subtri = running count of number of subtriangles
 * Returns value of integral
 */

   
double integrate(double *x, double *y, double *z, double dx, double dy, double dz, double A, double eta, double eps, double r, double r2, double r3, short level, short *maxlevel, int *num_subtri, int i, int j, double (*integrate_approx) (double *, double *, double *, double, double, double, double, int, int), double (*integrate_exact) (double *, double *, double *, double, double, int, int), bool (*error_computation) (double, double, double, double, double)){

   double p[3], dp, o=0, rmin;
   double A2 = A/2;
   double Lxy = dist_squared(x,y);
   double Lyz = dist_squared(y,z);
   double Lxz = dist_squared(x,z);

   //Update level counters
   if(level == -1)
      *num_subtri = 0;
   level++; *maxlevel = MAX(*maxlevel,level);
   
   //If error not small enough, continue refining
   if(error_computation(MAX3(Lxy,Lyz,Lxz), eta, r, r2, eps)){

      //Split triangle along longest side and integrate two halves separately
      if(Lxy >= MAX(Lyz,Lxz)){ //Split along edge xy, p=midpoint
         average(x,y,p); dp = norm_squared(p);

         //Triangle xpz
         rmin = tri_dist(x,p,z);
         if(MAX3(dx,dp,dz) <= r2){
            o = o + integrate_exact(x,p,z,r2,r3,i,j);
            ++*num_subtri;
         }else if(rmin < r2){
            o = o + integrate(x,p,z,dx,dp,dz,A2,eta,eps,r,r2,r3,level,maxlevel,num_subtri,i,j,integrate_approx,integrate_exact,error_computation);
         }else{
            ++*num_subtri;
         }

         //Triangle yzp
         rmin = tri_dist(y,z,p);
         if(MAX3(dy,dz,dp) <= r2){
            o = o + integrate_exact(y,z,p,r2,r3,i,j);
            ++*num_subtri;
         }else if(rmin < r2){
            o = o + integrate(y,z,p,dy,dz,dp,A2,eta,eps,r,r2,r3,level,maxlevel,num_subtri,i,j,integrate_approx,integrate_exact,error_computation);
         }else{
            ++*num_subtri;
         }
      }else if(Lyz >= MAX(Lxy,Lxz)){  //Split along edge yz, p=midpoint
         average(y,z,p); dp = norm_squared(p);

         //Triangle xyp
         rmin = tri_dist(x,y,p);
         if(MAX3(dx,dy,dp) <= r2){
            o = o + integrate_exact(x,y,p,r2,r3,i,j);
            ++*num_subtri;
         }else if(rmin < r2){
            o = o + integrate(x,y,p,dx,dy,dp,A2,eta,eps,r,r2,r3,level,maxlevel,num_subtri,i,j,integrate_approx,integrate_exact,error_computation);
         }else{
            ++*num_subtri;
         }
         //Triangle xpz
         rmin = tri_dist(x,p,z);
         if(MAX3(dx,dp,dz) <= r2){
            o = o + integrate_exact(x,p,z,r2,r3,i,j);
            ++*num_subtri;
         }else if(rmin < r2){
            o = o + integrate(x,p,z,dx,dp,dz,A2,eta,eps,r,r2,r3,level,maxlevel,num_subtri,i,j,integrate_approx,integrate_exact,error_computation);
         }else{
            ++*num_subtri;
         }
      }else{ //Split along edge xz, p=midpoint
         average(x,z,p); dp = norm_squared(p);

         //Triangle xyp
         rmin = tri_dist(x,y,p);
         if(MAX3(dx,dy,dp) <= r2){
            o = o + integrate_exact(x,y,p,r2,r3,i,j);
            ++*num_subtri;
         }else if(rmin < r2){
            o = o + integrate(x,y,p,dx,dy,dp,A2,eta,eps,r,r2,r3,level,maxlevel,num_subtri,i,j,integrate_approx,integrate_exact,error_computation);
         }else{
            ++*num_subtri;
         }
         //Triangle yzp
         rmin = tri_dist(y,z,p);
         if(MAX3(dy,dz,dp) <= r2){
            o = o + integrate_exact(y,z,p,r2,r3,i,j);
            ++*num_subtri;
         }else if(rmin < r2){
            o = o + integrate(y,z,p,dy,dz,dp,A2,eta,eps,r,r2,r3,level,maxlevel,num_subtri,i,j,integrate_approx,integrate_exact,error_computation);
         }else{
            ++*num_subtri;
         }
      }
   }else{ //Error condition is met; use approximate integration
      o = integrate_approx(x,y,z,A,eta,r2,r3,i,j);
      ++*num_subtri;
   }

   return o;
}


/*Allocate memory for a mxn array of int and initialize to val*/
int** array_int(int m, int n, int val){

   int **ptr = (int**)malloc(m*sizeof(int*));
   ptr[0] = (int*)malloc(m*n*sizeof(int));
   int i,j;
   for(i=0;i<m;i++){
      ptr[i] = ptr[0] + n*i;
      for(j=0;j<n;j++){
         ptr[i][j] = val;
      }
   }
   return ptr;
}

/*Allocate memory for a length m array of ints and initialize to val*/
bool* vector_bool(int m, bool val){

   bool *ptr = (bool*)malloc(m*sizeof(bool));
   int i;
   for(i=0;i<m;i++){
      ptr[i] = val;
   }
   return ptr;
}

/*Allocate memory for a length m array of ints and initialize to val*/
int* vector_int(int m, int val){

   int *ptr = (int*)malloc(m*sizeof(int));
   int i;
   for(i=0;i<m;i++){
      ptr[i] = val;
   }
   return ptr;
}

/*Allocate memory for a length m array of doubles and initialize to val*/
double* vector_double(int m, double val){

   double *ptr = (double*)malloc(m*sizeof(double));
   int i;
   for(i=0;i<m;i++){
      ptr[i] = val;
   }
   return ptr;
}

/*Allocate memory for a mxn array of doubles and initialize to val*/
double** array_double(int m, int n, double val){

   double **ptr = (double**)malloc(m*sizeof(double*));
   ptr[0] = (double*)malloc(m*n*sizeof(double));
   int i,j;
   for(i=0;i<m;i++){
      ptr[i] = ptr[0] + n*i;
      for(j=0;j<n;j++){
         ptr[i][j] = val;
      }
   }
   return ptr;
}

//Main subroutine to compute spherical volume invaraint and PCA on local neighborhoods
void svipca(double *P_ptr, int n, int *T_ptr, int m, bool *ID, double r, double eps_svi, double eps_pca, bool prog, double *S, double *M_ptr){

   double dx,dy,dz,eta,A,sgam, x[3], y[3], z[3], nu[3], mi[3], cij[9], **P, **M;
   double r3 = r*r*r;
   double r2 = r*r;
   double vol = 4*PI*r3/3;
   int i, j, t, q, **L, **F, *NN, *stack, tri_num, num_subtri, max_num_subtri, ki, kj, **T;
   short maxlevel;
   bool *B, *visited, *boundary, bi;

   //Setup 2D array pointers to input
   P = (double**)malloc(n * sizeof(double*));
   T = (int**)malloc(m * sizeof(int*));
   for(i=0;i<n;i++)
      P[i] = P_ptr + 3*i;
   for(i=0;i<m;i++)
      T[i] = T_ptr + 3*i;

   //Setup 2D array pointers to output
   M = (double**)malloc(n * sizeof(double*));
   for(i=0;i<n;i++)
      M[i] = M_ptr + 9*i;

   B = vector_bool(n,true);  //Records boundary points of mesh
   NN = vector_int(m,-1);    //Nearest neighbors
   stack = vector_int(m,-1);    //Stack
   visited = vector_bool(m,false); //Records which triangles have alread been visited in depth first search
   boundary = vector_bool(m,false);//Records which triangles intersect boundary of ball B(x,r)

   //Compute face graph
   L = vertex_to_triangle_list(T, n, m);
   F = face_graph(T, L, B, n, m);

   maxlevel = 0;
   max_num_subtri = 0;
   int interval = n/100;
   if(prog){
      printf("Progress...0%%"); fflush(stdout);
   }
   int perc = 0;
   /*Main loop over vertices*/
   for(i=0;i<n;i++){

      S[i] = -1; //Default value if not computed later
      for(j=0;j<9;j++)
         M[i][j] = j;

      if((i % interval == 1 || i == n-1) && prog){
         if(perc < 10){
            printf("\b\b%d%%",(100*(i+1))/n); fflush(stdout);
         }else{
            printf("\b\b\b%d%%",(100*(i+1))/n); fflush(stdout);
         }
         perc = (100*(i+1))/n;
      }

      if(B[i] && ID[i] && L[i][0] !=0){ //If vertex i is not a boundary vertex, then compute Gamma
         sgam = svigamma(T,P,L,i);
         S[i] = vol*sgam;

         //Depth first search to find neighboring triangles
         tri_num = 0;
         //bi = depth_first_search(T,P,F,NN,visited,boundary,&tri_num,L[i][1],r2,i);
         //bi = breadth_first_search(T,P,F,NN,visited,boundary,&tri_num,L[i][1],r2,i,stack);
         bi = depth_first_search_nr(T,P,F,NN,visited,boundary,&tri_num,L[i][1],r2,i,stack);

         //Zero out PCA variables
         mult(mi,0,mi);
         for(j=0;j<9;j++)
            cij[j]=0;
          
         //Loop over neighboring triangles in B(x,r)
         for(j=0;j<tri_num;j++){ //Compute spherical volume invariant
            t = NN[j];
            visited[t] = false;

            if(bi){ //If B(x,r) does not intersect mesh boundary
               tri_vertices(T,P,x,y,z,t,i); //retrieve vertices
               if(boundary[j]){  //If triangle intersects partial B(x,r), use refinement integration
                  dx = norm_squared(x); dy = norm_squared(y); dz = norm_squared(z);
                  A = tri_normal(x,y,z,nu);
                  eta = dot(x,nu);
                  if(MIN3(dx,dy,dz) != 0 && eta != 0){
                     S[i] = S[i] + integrate(x,y,z,dx,dy,dz,A,eta,eps_svi,r,r2,r3,-1,&maxlevel,
                                             &num_subtri,0,0,svi_integrate_approx,svi_integrate_exact,
                                             svi_error_computation);
                     max_num_subtri = MAX(num_subtri,max_num_subtri);
                  }
                  q=0;
                  for(ki=0;ki<3;ki++){
                     mi[ki] = mi[ki] + integrate(x,y,z,dx,dy,dz,A,eta,eps_pca,r,r2,r3,-1,&maxlevel,
                                                 &num_subtri,ki,0,pcami_integrate_approx,
                                                 pcami_integrate_exact,pcami_error_computation);
                     max_num_subtri = MAX(num_subtri,max_num_subtri);
                     for(kj=0;kj<3;kj++){
                        if(kj >= ki){
                           cij[q] = cij[q] + integrate(x,y,z,dx,dy,dz,A,eta,eps_pca,r,r2,r3,-1,&maxlevel,
                                                    &num_subtri,ki,kj,pcacij_integrate_approx,
                                                    pcacij_integrate_exact,
                                                    pcacij_error_computation);
                           max_num_subtri = MAX(num_subtri,max_num_subtri);
                        }
                        q++;
                     }
                  }
               }else{ //If triangle is interior to ball, use analytic formula
                  S[i] = S[i] + svi_integrate_exact(x,y,z,r2,r3,0,0);
                  q=0;
                  for(ki=0;ki<3;ki++){
                     mi[ki] = mi[ki] + pcami_integrate_exact(x,y,z,r2,r3,ki,kj);
                     for(kj=0;kj<3;kj++){
                        if(kj >= ki)
                           cij[q] = cij[q] + pcacij_integrate_exact(x,y,z,r2,r3,ki,kj);
                        q++;
                     }
                  }
               }
            }
         }
         if(bi){
            //Assemble parts to get matrix M
            q=0;
            for(ki=0;ki<3;ki++){
               for(kj=0;kj<3;kj++){
                  if(kj >= ki)
                     M[i][q] = cij[q] - mi[ki]*mi[kj]/S[i];
                  if(ki == kj)
                     M[i][q] += r2*S[i]/5;
                  q++;
               }
            }
            q=0;
            for(ki=0;ki<3;ki++){
               for(kj=0;kj<3;kj++){
                  if(kj < ki)
                     M[i][q] = M[i][ki+3*kj];
                  q++;
               }
            }
         }
      }
   }
   printf("\n");
   //printf("\nMax number of refinements = %d\n",maxlevel);
   //printf("Max number of subtriangles = %d\n",max_num_subtri);
}


//Main subroutine to compute spherical volume invaraint
void svi(double *P_ptr, int n, int *T_ptr, int m, bool *ID, double r, double eps, bool prog, double *S, double *Q){

   double dx,dy,dz,eta,A,sgam,x[3], y[3], z[3],  nu[3], **P;
   double r3 = r*r*r;
   double r2 = r*r;
   double vol = 4*PI*r3/3;
   int i, j, t, **L, **F, *NN, *stack, tri_num, num_subtri, max_num_subtri, **T;
   short maxlevel;
   bool *B, *visited, *boundary, bi;

   //Setup 2D array pointers to input
   P = (double**)malloc(n * sizeof(double*));
   T = (int**)malloc(m * sizeof(int*));
   for(i=0;i<n;i++)
      P[i] = P_ptr + 3*i;
   for(i=0;i<m;i++)
      T[i] = T_ptr + 3*i;

   B = vector_bool(n,true);  //Records boundary points of mesh
   NN = vector_int(m,-1);    //Nearest neighbors
   stack = vector_int(m,-1);    //Stack
   visited = vector_bool(m,false); //Records which triangles have alread been visited in depth first search
   boundary = vector_bool(m,false);//Records which triangles intersect boundary of ball B(x,r)

   //Compute face graph
   L = vertex_to_triangle_list(T, n, m);
   F = face_graph(T, L, B, n, m);

   maxlevel = 0;
   max_num_subtri = 0;
   int interval = n/100;

   if(prog){
      printf("Progress...0%%"); fflush(stdout);
   }

   int perc = 0;
   /*Main loop over vertices*/
   for(i=0;i<n;i++){

      S[i] = -1.0; //Default value if not computed later
      Q[i] = -1.0;

      if((i % interval == 1 || i == n-1) && prog){
         if(perc < 10){
            printf("\b\b%d%%",(100*(i+1))/n); fflush(stdout);
         }else{
            printf("\b\b\b%d%%",(100*(i+1))/n); fflush(stdout);
         }
         perc = (100*(i+1))/n;
      }

      if(B[i] && ID[i] && L[i][0] !=0){ //If vertex i is not a boundary vertex, then compute Gamma
         sgam = svigamma(T,P,L,i);
         S[i] = vol*sgam;
         Q[i] = sgam;

         //Depth first search to find neighboring triangles
         tri_num = 0;
         //bi = depth_first_search(T,P,F,NN,visited,boundary,&tri_num,L[i][1],r2,i);
         //bi = breadth_first_search(T,P,F,NN,visited,boundary,&tri_num,L[i][1],r2,i,stack);
         bi = depth_first_search_nr(T,P,F,NN,visited,boundary,&tri_num,L[i][1],r2,i,stack);
           
         //Loop over neighboring triangles in B(x,r)

         for(j=0;j<tri_num;j++){
            t = NN[j];
            visited[t] = false;

            if(bi){ //If B(x,r) does not intersect mesh boundary
               tri_vertices(T,P,x,y,z,t,i); //retrieve vertices
               if(boundary[j]){  //If triangle intersects boundary of B(x,r), use refinement integration
                  dx = norm_squared(x); dy = norm_squared(y); dz = norm_squared(z);
                  A = tri_normal(x,y,z,nu);
                  eta = dot(x,nu);
                  if(MIN3(dx,dy,dz) != 0 && eta != 0){
                     num_subtri = 0;
                     S[i] = S[i] + integrate(x,y,z,dx,dy,dz,A,eta,eps,r,r2,r3,-1,&maxlevel,&num_subtri,0,0,svi_integrate_approx,svi_integrate_exact,svi_error_computation);
                     max_num_subtri = MAX(num_subtri,max_num_subtri);
                  }
               }else //If triangle is interior to ball, use analytic formula
                  S[i] = S[i] + svi_integrate_exact(x,y,z,r2,r3,0,0);
            }
         }
      }
   }
   printf("\n");
   //printf("\nMax number of refinements = %d\n",maxlevel);
   //printf("Max number of subtriangles = %d\n",max_num_subtri);

}

/*Change of basis used in exact integration*/
static inline double change_of_basis(double *x, double *y, double *nu, double *xtilde, double *qi, double *pi, double *pi1){

   double e1[3], e2[3], t[3], a;

   /*Compute basis*/
   sub(y,x,e1);
   a=1/sqrt(dot(e1,e1)); mult(e1,a,e1);
   cross(e1,nu,e2);
   mult(e2,-1,e2);

   /*Change of basis formulas*/
   sub(x,xtilde,t);
   *pi = dot(t,e1);
   *qi = dot(t,e2);
   sub(y,xtilde,t);
   *pi1 = dot(t,e1);

   return dot(x,e2);
}

/*Computes stopping condition for triangle refinement
 * lmax = squared max side length of triangle
 * eta = x \cdot \nu on triangle
 * r2 = r^2
 * eps = error tolerance
 *
 * Should return True if error is too large and further refinement necessary
 * Otherwise returns false
 */
bool svi_error_computation(double lmax, double eta, double r, double r2, double eps){

   double err = MAX(r-sqrt(lmax),0);
   err = eps*err*err*err*err/r2;
   return (lmax > err || lmax > r2) && eps < 100;
}


/*Approximate integration over triangle for SVI
 * x,y,z = vertices of triangle
 * r3 = r^3
 * theta = angle used in integration
 */
double svi_integrate_approx(double *x, double *y, double *z, double A, double eta, double r2, double r3, int i, int j){
      double p[3];
      centroid(x,y,z,p);
      return MIN(1 - r3*pow(norm_squared(p),-1.5),0)*A*eta/3;
}

//Returns angle between yx and yz
double angle(double *x, double *y, double *z){
   
   double p[3], q[3], a, b, c;
   sub(x,y,p); sub(z,y,q);
   a = norm(p); b = norm(q); c = dot(p,q);

   return acos(c/(a*b));
}

/*Analytic integration of hypersingular kernel 1/||x||^3$ over triangle xyz
 * x,y,z = vertices of triangle
 * r3 = r^3
 */
double svi_integrate_exact(double *x, double *y, double *z, double r2, double r3, int i, int j){

   double xtilde[3], nu[3];
   double qi, pi, pi1, a, num, den, eta, A, dx, dy, dz, o=0, txy, tyz, tzx, theta;
   bool bxy, byz, bzx;

   dx = norm(x); dy = norm(y); dz = norm(z);
   A = tri_normal(x,y,z,nu);
   eta = dot(x,nu);

   if(MIN3(dx,dy,dz) != 0 && eta != 0){

      mult(nu,eta,xtilde);

      txy = change_of_basis(x,y,nu,xtilde,&qi,&pi,&pi1);
      num = -2*pi*qi*eta*dx;
      den = qi*qi*dx*dx - pi*pi*eta*eta;
      a = atan2(num,den);
      num = -2*pi1*qi*eta*dy;
      den = qi*qi*dy*dy - pi1*pi1*eta*eta;
      a = a - atan2(num,den);
      
      tyz = change_of_basis(y,z,nu,xtilde,&qi,&pi,&pi1);
      num = -2*pi*qi*eta*dy;
      den = qi*qi*dy*dy - pi*pi*eta*eta;
      a = a + atan2(num,den);
      num = -2*pi1*qi*eta*dz;
      den = qi*qi*dz*dz - pi1*pi1*eta*eta;
      a = a - atan2(num,den);
    
      tzx = change_of_basis(z,x,nu,xtilde,&qi,&pi,&pi1);
      num = -2*pi*qi*eta*dz;
      den = qi*qi*dz*dz - pi*pi*eta*eta;
      a = a + atan2(num,den);
      num = -2*pi1*qi*eta*dx;
      den = qi*qi*dx*dx - pi1*pi1*eta*eta;
      a = a - atan2(num,den);

      theta = 0;
      if(txy <=0 && tyz <= 0 && tzx <= 0){
         bxy = txy < 0; byz = tyz < 0; bzx = tzx < 0;
         if(bxy && byz && bzx){
            theta = 2*PI;
         }else if((bxy && byz) || (byz && bzx) || (bzx && bxy)){
            theta = PI;
         }else if(!bxy && !byz){
            theta = angle(x,y,z);
         }else if(!byz && !bzx){
            theta = angle(y,z,x);
         }else if(!bzx && !bxy){
            theta = angle(z,x,y);
         }
      }
      o = A*eta/3 - (r3/3)*(a/2 + theta*SIGN(eta));
   }

   return o;
}

/*Computes gamma
 * T = triangle list
 * P = vertices
 * L = vertex to triangle list
 * i = vertex
 */
double svigamma(int **T, double **P, int **L, int i){

   double p[3], q[3], nu[3], e1[3], e2[3], e3[3];
   double alpha, d, phi1, phi2, Gam, v1, v2, v3, xx[3], yy[3], zz[3], *x, *y, *z, *temp, ne1, ne3;
   int k = L[i][0]; //Number of adjacent triangles
   int j,t;
   x=xx;y=yy;z=zz;

   nu[0] = 0; nu[1] = 0; nu[2] = 0;
  
   //Compute outward unit normal vector
   for(j=1;j<=k;j++){
      t = L[i][j];
      tri_vertices(T,P,x,y,z,t,i);
      tri_normal(x,y,z,p);
      add(nu,p,nu);
   }
   e3[0] = -nu[0]; e3[1] = -nu[1]; e3[2] = -nu[2];
   ne3 = norm(e3);
   d=1/ne3; mult(e3,d,e3);

   //Choose e1 orthogonal to e3
   e1[0] = 0; e1[1] = 0; e1[2] = 0;
   if(e3[0] != 0 || e3[1] !=0){
      e1[0] = -e3[1];
      e1[1] = e3[0];
   }else{
      e1[1] = -e3[2];
      e1[2] = e3[1];
   }
   ne1 = norm(e1);
   d=1/ne1; mult(e1,d,e1);

   //Compute e2 by cross product
   cross(e1,e3,e2);
   mult(e2,-1,e2);
  
   Gam = 0;
   for(j=1;j<=k;j++){
      t = L[i][j];
      tri_vertices(T,P,x,y,z,t,i);

      //vertex i should be first
      if(T[t][1] == i){
         temp = x; x=y; y=z; z=temp;
      }else if(T[t][2] == i){
         temp = z; z=y; y=x; x=temp;
      }

      //Change basis
      new_coordinates(x,e1,e2,e3);
      new_coordinates(y,e1,e2,e3);
      new_coordinates(z,e1,e2,e3);

      tri_normal(x,y,z,nu);
      sub(y,x,q);
      sub(z,x,p);
      alpha = atan2(nu[1],nu[0]);
      d = sqrt(nu[0]*nu[0] + nu[1]*nu[1]);
      phi1 = atan2(p[1],p[0]);
      phi2 = atan2(q[1],q[0]);

      Gam = Gam + asin(d*sin(phi2-alpha)) - asin(d*sin(phi1-alpha));
   }
   Gam = 0.5 - Gam/(4*PI);
   return Gam;
}

/*Computes stopping condition for triangle refinement
 * lmax = squared max side length of triangle
 * eta = x \cdot \nu on triangle
 * r2 = r^2
 * eps = error tolerance
 *
 * Should return True if error is too large and further refinement necessary
 * Otherwise returns false
 */
bool pcami_error_computation(double lmax, double eta, double r, double r2, double eps){
   return (lmax > r2*eps*eps || lmax > r2) && eps < 100;
}
bool pcacij_error_computation(double lmax, double eta, double r, double r2, double eps){
   return (lmax > r2*eps*eps || lmax > r2) && eps < 100;
}
/*Analytic integration for PCA on local neighborhoods
 * x,y,z = vertices of triangle
 * r3 = r^3
 */
double pcami_integrate_exact(double *x, double *y, double *z, double r2, double r3, int i, int j){
   
   double nu[3];
   double A = tri_normal(x, y, z, nu);
   double eta = dot(x,nu);

   return (A/4)*(eta*(x[i] + y[i] + z[i])/3 - r2*nu[i]);
}

double pcacij_integrate_exact(double *x, double *y, double *z, double r2, double r3, int i, int j){
   
   double nu[3];
   double A = tri_normal(x, y, z, nu);
   double eta = dot(x,nu);
   double ai = A*(x[i] + y[i] + z[i])/3;
   double aj = A*(x[j] + y[j] + z[j])/3;

   double bij = (A/12)*(2*x[i]*x[j] + 2*y[i]*y[j] + 2*z[i]*z[j] + x[i]*y[j] + x[j]*y[i] + x[i]*z[j] + x[j]*z[i] + y[i]*z[j] + y[j]*z[i]);

   return (eta/5)*bij - (r2/10)*(nu[i]*aj + nu[j]*ai);
}

/*Approximate integration for PCA on local neighborhoods
 * x,y,z = vertices of triangle
 * r3 = r^3
 * theta = angle used in integration
 */
double pcami_integrate_approx(double *x, double *y, double *z, double A, double eta, double r2, double r3, int i, int j){

   double p[3], nu[3];
   tri_normal(x, y, z, nu);
   centroid(x,y,z,p);

   double integrand = 0;
   if(norm_squared(p) <= r2)
      integrand = eta*x[i] - r2*nu[i];

   return integrand*A/4;
}

double pcacij_integrate_approx(double *x, double *y, double *z, double A, double eta, double r2, double r3, int i, int j){

   double p[3], nu[3];
   tri_normal(x, y, z, nu);
   centroid(x,y,z,p);

   double integrand = 0;
   if(norm_squared(p) <= r2)
      integrand = 2*p[i]*p[j]*eta - r2*(p[j]*nu[i] + p[i]*nu[j]);

   return integrand*A/10;
}

// ACTUAL
// CODE
// STARTS
// IN 5
// 4
// 3
// 2
// 1 HERE:::::::::::::::::::::::::::::::::::::::::::::::::
SphericalVolumeInvariantFilterPlugin::SphericalVolumeInvariantFilterPlugin()
{
  typeList <<
     FP_NORMAL_SPHERICAL_VOLUME_INVARIANT;

  FilterIDType tt;

  foreach(tt , types())
    {
      actionList << new QAction(filterName(tt), this);
      if (tt == FP_NORMAL_SPHERICAL_VOLUME_INVARIANT){
        //If you want a shortcut key, here it is:
          actionList.last()->setShortcut(Qt::CTRL + Qt::SHIFT + Qt::Key_S);
          actionList.last()->setPriority(QAction::HighPriority);
      }
    }
}

QString SphericalVolumeInvariantFilterPlugin::filterName(FilterIDType filter) const
{
 switch(filter)
 {
      //This is the name of the plugin, as it appears in the meshlab menu
       case FP_NORMAL_SPHERICAL_VOLUME_INVARIANT:    return tr("Spherical Volume Invariant");
 }
 assert(0);
 return QString("Unknown filter");
}

QString SphericalVolumeInvariantFilterPlugin::filterInfo(FilterIDType filterId) const
{
 switch(filterId)
 {
   //This is the description of the plugin
    case FP_NORMAL_SPHERICAL_VOLUME_INVARIANT:  return tr("Computes the spherical volume invariant.");
 }
 assert(0);
 return QString("Unknown filter");
}

void SphericalVolumeInvariantFilterPlugin::initParameterSet(QAction *action, MeshModel &md, RichParameterSet &parlst)
{
 switch(ID(action))
 {
    case FP_NORMAL_SPHERICAL_VOLUME_INVARIANT:
     {
         
         parlst.addParam(new RichEnum("Mode", 0, QStringList() << "Spherical Volume Invariant" << "1st Principal Curvature" << "2nd Principal Curvature"<< "Gaussian Curvature"<< "Mean Curvature"<< "Abs Principal Curvature Difference", tr("Plot:"), tr("Select invariant to plot")));
         
         parlst.addParam(new RichFloat("Radius", 1.0f, "Radius", "Radius of patch to use."));
        
         
      //This is to add a popup window with parameters. Other things need to be added as well to enable this, so leave it off for now, and have the radius hard-coded at first.
        //parlst.addParam(new RichBool("allLayers", false, "Apply to all visible Layers", "If selected, the filter will be applied to all visible mesh Layers."));
        
     }
    break;
 }
}




bool SphericalVolumeInvariantFilterPlugin::applyFilter(QAction *action, MeshDocument &md, RichParameterSet & par, vcg::CallBackPos * /*cb*/)
{
   
    MeshModel* mesh = md.mm();
    MeshModel &m=*(md.mm());
    CMeshO::VertexIterator vi;

    //Requirements
    tri::RequirePerVertexNormal(m.cm);
    tri::UpdateNormal<CMeshO>::PerVertexNormalized(m.cm);

    
    
    
    double r;   //radius of computation
    r = static_cast<double>(par.getFloat("Radius"));
    
    int Mode = par.getEnum("Mode");
    
    double eps_svi = .1; //svi error tolerance
    double eps_pca = .1; //pca error tolerance
    
    //Number of vertices in whole mesh
    int num_verts = m.cm.vert.size();
    int num_faces = m.cm.face.size();
  
    //Build connectivity list:
    int *Tri = new int[3*num_faces];

    int i = 0;
    int j = 0;
    int k = 0;
    for(i=0;i<num_faces;i++){
        for(int k=0;k<3;k++){
            Tri[j] = m.cm.face[i].V(k)->Index();
            j++;
        }
    } //end for loop
    
    //Now get points & logical array:
    bool *ID = new bool[num_verts];
    double *Pts = new double[3*num_verts]; //vertices
    double *Ns = new double[3*num_verts]; //surface normals
    i = 0;
    j = 0;
    k = 0;
    
    
    for(CMeshO::VertexIterator vi=m.cm.vert.begin(); vi!=m.cm.vert.end(); ++vi)
    {
        ID[i] = true;
        for(int k=0;k<3;k++){
            Pts[j] = ((*vi).P()[k]);
            Ns[j] =(*vi).N()[k];
            j++;
        }
        i++;
     }
    // Run SVI:
    double *S = new double[num_verts]; //svi
    
    double *V1 = new double[3*num_verts];
    double *V2 = new double[3*num_verts];
    double *V3 = new double[3*num_verts];
    double *K1 = new double[num_verts]; //k1
    double *K2 = new double[num_verts]; //k2
    double *MK = new double[num_verts]; //gaussian
    double *GK = new double[num_verts]; //mean
    double *absK = new double[num_verts]; //absolute difference
    double *Q = new double[num_verts]; //gamma
    
    //double M[3*num_verts][3];
    double *M = new double[9*num_verts];
    
    bool prog = 0;
    
    
    //svi(Pts,num_verts,Tri,num_faces,ID,r,eps_svi,true,S,Q);
    
    svipca(Pts, num_verts, Tri, num_faces, ID, r, eps_svi, eps_pca, prog, S, M);
    
    
    double minc=1e9, maxc=-1e9, minabsc=1e9;
    
    double Kdiff;
    double Ksum;
    double k1t;
    double k2t;
    
    double ev[3];
    double dot;
    double dot2;
    int loc;
    double pi = 3.14159265;
    
    typedef Eigen::Matrix<double, 3, 3> Matrix3d;
    Eigen::Matrix3d A;
    
    i = 0;
    // get max & min for normalization
    for(i=0;i<num_verts;i++){
        
//      eigenvalue stuff begins here:
        A(0,0) = M[9*i];   A(0,1) = M[9*i+1]; A(0,2) = M[9*i+2];
        A(1,0) = M[9*i+3]; A(1,1) = M[9*i+4]; A(1,2) = M[9*i+5];
        A(2,0) = M[9*i+6]; A(2,1) = M[9*i+7]; A(2,2) = M[9*i+8];
        
        Eigen::EigenSolver<Matrix3d> es;
        es.compute(A, /* computeEigenvectors = */ true);
        
        
        // ascertain which is normal to surface:
        j = 0; dot = -1;
        for(j=0;j<3;j++){
            dot2 = abs(Ns[3*i  ]*es.eigenvectors()(0,j).real()+             Ns[3*i+1]*es.eigenvectors()(1,j).real()+ Ns[3*i+2]*es.eigenvectors()(2,j).real() );
            
           if(dot2>dot){
               dot = dot2;
               loc = j;
            }
        }
        
      //throw in loops for building V1,V2,V3
        j = 0;
        switch(loc) {
      case 0 :
                ev[0] = es.eigenvalues()(1,0).real();
                ev[1] = es.eigenvalues()(2,0).real();
                ev[2] = es.eigenvalues()(0,0).real();
                
                for(j=0;j<3;j++){
                    V1[3*i+j] = es.eigenvectors()(j,1).real();
                    V2[3*i+j] = es.eigenvectors()(j,2).real();
                    V3[3*i+j] = es.eigenvectors()(j,0).real();
                }
                break;
      case 1 :
                ev[0] = es.eigenvalues()(0,0).real();
                ev[1] = es.eigenvalues()(2,0).real();
                ev[2] = es.eigenvalues()(1,0).real();
                
                for(j=0;j<3;j++){
                    V1[3*i+j] = es.eigenvectors()(j,0).real();
                    V2[3*i+j] = es.eigenvectors()(j,2).real();
                    V3[3*i+j] = es.eigenvectors()(j,1).real();
                }
                break;
      case 2 :
                ev[0] = es.eigenvalues()(0,0).real();
                ev[1] = es.eigenvalues()(1,0).real();
                ev[2] = es.eigenvalues()(2,0).real();
                
                for(j=0;j<3;j++){
                    V1[3*i+j] = es.eigenvectors()(j,0).real();
                    V2[3*i+j] = es.eigenvectors()(j,1).real();
                    V3[3*i+j] = es.eigenvectors()(j,2).real();
                }
                break;
        }
        
        // SO: don't use ev[2]
        
        Kdiff = (ev[0]-ev[1])*24/(pi*pow(r,6));
        Ksum = (16*pi*pow(r,3)/3 - 8*S[i])/(pi*pow(r,4));
        k1t = (Kdiff + Ksum)/2;
        k2t = (Ksum - Kdiff)/2;
        
        if (k1t>k2t){
            K1[i] = k1t;
            K2[i] = k2t;
        }
        else{
            K1[i] = k2t;
            K2[i] = k1t;
        }
        
        GK[i] = k1t*k2t;
        MK[i] = .5*(k1t+k2t);
        absK[i] =abs(K1[i] - K2[i]);

        
        switch(Mode){
            case 0 :
                m.cm.vert[i].Q() = S[i];
                break;
            case 1 :
                m.cm.vert[i].Q() = K1[i];
                break;
            case 2 :
                m.cm.vert[i].Q() = K2[i];
                break;
            case 3 :
                m.cm.vert[i].Q() = GK[i];
                break;
            case 4 :
                m.cm.vert[i].Q() = MK[i];
                break;
            case 5 :
                m.cm.vert[i].Q() = absK[i];
                break;
        }
        
      
        
        // need to add vector processing here...
        
        //also visualization of vectors
        
        
        
        
        
        
        minc = std::min(S[i],minc);
        maxc = std::max(S[i],maxc);
    } //end for loop


    
           switch(Mode){
            case 0 :
                Log("SVI");
                   break;
            case 1 :
                Log("K1");
                   break;
            case 2 :
                Log("K2");
                   break;
            case 3 :
                Log("Gaussian");
                   break;
            case 4 :
                Log("Mean");
                   break;
            case 5 :
                Log("Absolute");
                   break;
        }
    
    //i = 0;
    //for(i=0;i<num_verts;i++){
    //    S[i] = (S[i]-minc)/(maxc-minc);
        
    //} //end for loop
    
    
    // color the mesh:
    vcg::Histogramf H;
    vcg::tri::Stat<CMeshO>::ComputePerVertexQualityHistogram(m.cm,H);
    vcg::tri::UpdateColor<CMeshO>::PerVertexQualityRamp(m.cm,H.Percentile(0.01f),H.Percentile(0.99f));

          /*
          //Print angle to log
          Log("Break #%d, Angle = %.0f\n", break_number, theta);
          this->RealTimeLog(QString("Virtual Goniometer"),m.shortName(),"Break #%d, Angle = %.0f\n", break_number, theta);

          //Color one region VGcolors[color_i]
          tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
          i=0; j=0;
            for(CMeshO::VertexIterator vi=m.cm.vert.begin(); vi!=m.cm.vert.end(); ++vi){
             if(current_selection[i]){
                if(C[j] == 1)
                   (*vi).SetS();
                j++;
             }
             i++;
          }
            tri::UpdateColor<CMeshO>::PerVertexConstant(m.cm, VGcolors[color_i], TRUE);

          //Color other region VGcolors[color_j]
          tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
          i=0; j=0;
            for(CMeshO::VertexIterator vi=m.cm.vert.begin(); vi!=m.cm.vert.end(); ++vi){
             if(current_selection[i]){
                if(C[j] == 2)
                   (*vi).SetS();
                j++;
             }
             i++;
          }
            tri::UpdateColor<CMeshO>::PerVertexConstant(m.cm, VGcolors[color_j], TRUE);

          //Clear the selection
          tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
         */
    
   Log("%d,%d,%d,%d,%d,$d",S[0],K1[0],K2[0],GK[0],MK[0],absK[0]);
   Log("TEST");
  

   return true;
}

MeshFilterInterface::FilterClass SphericalVolumeInvariantFilterPlugin::getClass(QAction *action)
{
  switch(ID(action))
  {
       //This line controls which menu under Filters the plugin appears
       case FP_NORMAL_SPHERICAL_VOLUME_INVARIANT: return FilterClass(MeshFilterInterface::Normal);
  }
  return MeshFilterInterface::Selection;
}


 int SphericalVolumeInvariantFilterPlugin::getRequirements(QAction *action)
{
 switch(ID(action))
  {
     case FP_NORMAL_SPHERICAL_VOLUME_INVARIANT: return MeshModel::MM_VERTCOLOR;

      default: return MeshModel::MM_NONE;
  }
}

int SphericalVolumeInvariantFilterPlugin::postCondition(QAction *action) const
{
    switch(ID(action))
    {
   }
  return MeshModel::MM_ALL;
}

int SphericalVolumeInvariantFilterPlugin::getPreConditions( QAction * action) const
{
  switch(ID(action))
  {
  }
  return 0;
}
MESHLAB_PLUGIN_NAME_EXPORTER(SphericalVolumeInvariantFilterPlugin)
