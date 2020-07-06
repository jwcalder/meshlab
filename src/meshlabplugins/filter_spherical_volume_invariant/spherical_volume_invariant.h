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

#ifndef FILTER_SELECT_H
#define FILTER_SELECT_H

#include "stdbool.h"
#include <QObject>
#include <common/interfaces.h>

#include "stdbool.h"
#define tri_vertices(T,P,x,y,z,t,i) sub(P[T[t][0]],P[i],x);sub(P[T[t][1]],P[i],y);sub(P[T[t][2]],P[i],z)

#define ABS(a) (((a)<0)?-(a):(a))
#define SIGN(a) (((a)<0)?(-1):(1))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MIN3(a,b,c) (MIN(MIN(a,b),c))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MAX2(a,b) (((a)>(b))?(a):(b))
#define MAX3(a,b,c) (MAX(MAX(a,b),c))
#define PI 3.14159265359

#define dot(x,y) (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])
#define norm(x) (sqrt(x[0]*x[0] + x[1]*x[1] +  x[2]*x[2]))
#define norm_squared(x) (x[0]*x[0] + x[1]*x[1] +  x[2]*x[2])
#define dist(x,y) (sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2])))
#define dist_squared(x,y) ((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]))
#define cross(x,y,z) z[0] = x[1]*y[2] - x[2]*y[1]; z[1] = x[2]*y[0] - x[0]*y[2]; z[2] = x[0]*y[1] - x[1]*y[0]
#define centroid(x,y,z,p) p[0] = (x[0] + y[0] + z[0])/3; p[1] = (x[1] + y[1] + z[1])/3; p[2] = (x[2] + y[2] + z[2])/3
#define average(x,y,z) z[0] = (x[0] + y[0])/2; z[1] = (x[1] + y[1])/2; z[2] = (x[2] + y[2])/2
#define add(x,y,z) z[0] = x[0] + y[0]; z[1] = x[1] + y[1]; z[2] = x[2] + y[2]
#define sub(x,y,z) z[0] = x[0] - y[0]; z[1] = x[1] - y[1]; z[2] = x[2] - y[2]
#define mult(x,a,z) z[0] = a*x[0]; z[1] = a*x[1]; z[2] = a*x[2]
#define new_coordinates(x,a,b,c) v1 = dot(x,e1); v2 = dot(x,e2); v3 = dot(x,e3);x[0] = v1; x[1] = v2; x[2] = v3

// svi_computations.h:
void svipca(double *P_ptr, int n, int *T_ptr, int m, bool *ID, double r, double eps_svi, double eps_pca, bool prog, double *S, double *M_ptr);
void svi(double *P_ptr, int n, int *T_ptr, int m, bool *ID, double r, double eps, bool prog, double *S, double *Q);
bool svi_error_computation(double lmax, double eta, double r, double r2, double eps);
double svi_integrate_approx(double *x, double *y, double *z, double A, double eta, double r2, double r3, int i, int j);
double svi_integrate_exact(double *x, double *y, double *z, double r2, double r3, int i, int j);
double svigamma(int **T, double **P, int **L, int i);
bool pcami_error_computation(double lmax, double eta, double r, double r2, double eps);
bool pcacij_error_computation(double lmax, double eta, double r, double r2, double eps);
double pcami_integrate_exact(double *x, double *y, double *z, double r2, double r3, int i, int j);
double pcacij_integrate_exact(double *x, double *y, double *z, double r2, double r3, int i, int j);
double pcami_integrate_approx(double *x, double *y, double *z, double A, double eta, double r2, double r3, int i, int j);
double pcacij_integrate_approx(double *x, double *y, double *z, double A, double eta, double r2, double r3, int i, int j);

//memory allocation.h:
int** array_int(int m, int n, int val);
bool* vector_bool(int m, bool val);
int* vector_int(int m, int val);
double* vector_double(int m, double val);
double** array_double(int m, int n, double val);

//mesh operations.h
double tri_dist(double *x, double *y, double *z);
int** vertex_to_triangle_list(int **T, int n, int m);
int** face_graph(int **T, int **L, bool *B, int n, int m);
bool depth_first_search(int **T, double **P, int **F, int *NN, bool *v, bool *b, int *num, int ind, double r2, int i);
bool breadth_first_search(int **T, double **P, int **F, int *NN, bool *v, bool *b, int *num, int ind, double r2, int i, int *stack);
bool depth_first_search_nr(int **T, double **P, int **F, int *NN, bool *v, bool *b, int *num, int ind, double r2, int i, int *stack);
double integrate(double *x, double *y, double *z, double dx, double dy, double dz, double A, double eta, double eps, double r, double r2, double r3, short level, short *maxlevel, int *num_subtri, int i, int j, double (*integrate_approx) (double *, double *, double *, double, double, double, double, int, int), double (*integrate_exact) (double *, double *, double *, double, double, int, int), bool (*error_computation) (double, double, double, double, double));



class SphericalVolumeInvariantFilterPlugin : public QObject, public MeshFilterInterface
{
    Q_OBJECT
    MESHLAB_PLUGIN_IID_EXPORTER(MESH_FILTER_INTERFACE_IID)
    Q_INTERFACES(MeshFilterInterface)
        
        public:
    /* naming convention :
         - FP -> Filter Plugin
         - name of the plugin separated by _
    */
    enum {
      FP_NORMAL_SPHERICAL_VOLUME_INVARIANT
    } ;

    SphericalVolumeInvariantFilterPlugin();
   //~SphericalVolumeInvariantFilterPlugin();
   virtual QString filterInfo(FilterIDType filter) const;
   virtual QString filterName(FilterIDType filter) const;
   
   virtual FilterClass getClass(QAction *);
   void initParameterSet(QAction *action, MeshModel &m, RichParameterSet &parlst);
   int getPreConditions(QAction *) const;
   int postCondition( QAction* ) const;
   int getRequirements(QAction *);
   bool applyFilter(QAction *filter, MeshDocument &md, RichParameterSet & /*parent*/, vcg::CallBackPos * cb) ;
   FILTER_ARITY filterArity(QAction *) const {return SINGLE_MESH;}
    
    
};

#endif
