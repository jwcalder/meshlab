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

#include <math.h>
#include <stdlib.h>
#include <vcg/space/colorspace.h>
#include "mesh_segmentation.h"
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/point_outlier.h>

#include <QOpenGLContext>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <random>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include "Matrix.h"

#define TRUE 1
#define FALSE 0

using namespace std;
using namespace vcg;

vector<Color4b> C{Color4b::Red,Color4b::Blue,Color4b::Green,Color4b::Magenta,Color4b::Yellow,Color4b::Cyan,Color4b::LightGreen,Color4b::LightRed,Color4b::DarkGreen,Color4b::DarkRed,Color4b::DarkBlue,Color4b::Black,Color4b::White,Color4b::Gray};

static vector<int> label_ind_verts; //Indices of labeled vertices
static vector<int> label_val;  //Label values

// ERROR CHECKING UTILITY
#define CheckError(x,y); if ((x)) {this->errorMessage = (y); return false;}
///////////////////////////////////////////////////////

MeshSegmentationFilterPlugin::MeshSegmentationFilterPlugin()
{
  typeList <<
     FP_NORMAL_MESH_SEGMENTATION;

  FilterIDType tt;

  foreach(tt , types())
    {
        actionList << new QAction(filterName(tt), this);
        if (tt == FP_NORMAL_MESH_SEGMENTATION){
           //If you want a shortcut key, here it is:
           actionList.last()->setShortcut(Qt::CTRL + Qt::SHIFT + Qt::Key_S);
           actionList.last()->setPriority(QAction::HighPriority);
        }
    }
}

QString MeshSegmentationFilterPlugin::filterName(FilterIDType filter) const
{
 switch(filter)
 {
      //This is the name of the plugin, as it appears in the meshlab menu
	   case FP_NORMAL_MESH_SEGMENTATION:    return tr("Mesh Segmentation");
 }
 assert(0);
 return QString("Unknown filter");
}

QString MeshSegmentationFilterPlugin::filterInfo(FilterIDType filterId) const
{
 switch(filterId)
 {
   //This is the description of the plugin
	case FP_NORMAL_MESH_SEGMENTATION:  return tr("Mesh segmentation using semi-supervised learning.");
 }
 assert(0);
 return QString("Unknown filter");
}

double weight_map(double a){
   double sgn = 1.0;
   if(a < 0)
      sgn = -1.0;
   return sgn*pow(a*sgn,4) + 1;
}

//Removes repeated labels from label_val and label_ind_verts obtained by constantly rereading the .txt file
void remove_repeated_labels(){
   
   vector<bool> is_labeled(max(label_ind_verts)+1,0);
   vector<int> index_label_val(max(label_ind_verts)+1,-1);

   for(int i=0; i<label_val.size(); i++){
      if(is_labeled[label_ind_verts[i]] == 0){
         index_label_val[label_ind_verts[i]] = label_val[i];
         is_labeled[label_ind_verts[i]] = 1;
      }
   }

   label_val.clear();
   label_ind_verts.clear();
   for(int i=0; i<is_labeled.size(); i++){
      if(is_labeled[i]){
         label_val.push_back(index_label_val[i]);
         label_ind_verts.push_back(i);
      }
   }
}
//Returns C-string with mesh name (filename_Mesh.ply returns filename)
void mesh_name(MeshModel &m, char *plyfile){
   
   //Get mesh ply filename
   QString qs_plyfile = m.shortName();
   int len = qs_plyfile.length();
   strcpy(plyfile,(char *)qUtf8Printable(qs_plyfile));

   //Remove .ply extension
   plyfile[len-4] = '\0';

   //Remove _Mesh if found
   if(!strcmp(plyfile+len-9,"_Mesh"))
      plyfile[len-9] = '\0';
}
tuple<float, int, float, vector<double>> load_params(char *out_file){

   float radius, p;
   int num_nodes;
   vector<double> user_weights(C.size(),0.0);
   FILE *pFile;   

   /*label_ind_verts.clear();
   label_val.clear();*/

   pFile = fopen(out_file,"r");

   if(pFile == NULL){  //User defaults if file does not exist
      radius = 3.0;
      p = 1.0;
      num_nodes = 5000;
   }else{

      fscanf(pFile,"Radius,%f\n",&radius);
      fscanf(pFile,"Number of nodes,%d\n",&num_nodes);
      fscanf(pFile,"Weight matrix parameter,%f\n",&p);
      for(int i=0; i<user_weights.size(); i++){
         fscanf(pFile,"Face %*d weight,%lf\n",&user_weights[i]);
      }
      
      //Advance by one line
      char buff[100];
      fscanf(pFile, "%[^\n]\n", buff);

      int a,b;
      int i=0;
      while(fscanf(pFile,"%d,%d\n",&a,&b) == 2){
         label_ind_verts.push_back(a);
         label_val.push_back(b);
      }

      remove_repeated_labels();

      fclose(pFile);
   }

   return make_tuple(radius, num_nodes, p, user_weights);
}
void save_params(char *out_file, float radius, int num_nodes, float p, vector<double> user_weights){

   FILE *pFile;   
   pFile = fopen(out_file,"w");
   
   fprintf(pFile,"Radius,%f\n",radius);
   fprintf(pFile,"Number of nodes,%d\n",num_nodes);
   fprintf(pFile,"Weight matrix parameter,%f\n",p);
   for(int i=0; i<user_weights.size(); i++)
      fprintf(pFile,"Face %d weight,%f\n",i+1,user_weights[i]);
   
   fprintf(pFile,"Vertex index, Label value\n");
   for(int i=0; i<label_val.size(); i++)
      fprintf(pFile,"%d,%d\n",label_ind_verts[i],label_val[i]);

   fclose(pFile);

}

double withness(vector<double> x)
{
   int sum = 0;
   int n = x.size();
   sort(x.begin(), x.end());
   vector<double> v(n-1,0);

   //Compute Variance
   double mean = 0.0;
   for(int i=0; i<n; i++)
      mean += x[i];
   mean = mean/n;
   double nvar = 0.0;
   for(int i = 0; i<n; i++)
      nvar += (x[i] - mean)*(x[i] - mean);

   int minind = 0;
   for(int i = 0; i<n-1; i++)
   {
      double m1 = 0.0;
      double m2 = 0.0;

      for(int j = 0; j<i+1; j++)
         m1 += x[j];
      for(int j = i+1; j<n; j++)
         m2 += x[j];

      m1 = m1/(i+1);
      m2 = m2/(n-i-1);

      v[i] = 0;
      for(int j = 0; j<i+1; j++)
         v[i] += (x[j]-m1)*(x[j]-m1);
      for(int j = i+1; j<n; j++)
         v[i] += (x[j]-m2)*(x[j]-m2);
      v[i] = v[i]/nvar;

      if(v[i] < v[minind])
         minind = i;
   }
   
   print(x);
   printf("minind=%d\n",minind);
   return x[minind];
}
//Returns indices of vertices in mesh that were selected
vector<unsigned int> get_selected_indices(MeshModel &m){
   
   vector<unsigned int> indices;
   int i = 0;
   CMeshO::VertexIterator vi;
   for(vi=m.cm.vert.begin(); vi!=m.cm.vert.end(); ++vi)
   {
      if(!(*vi).IsD() && (*vi).IsS())
         indices.push_back(i);
      i++;
   }
   return indices;
}
void patch_statistics(vector<float> &vecx, vector<float> &vecy, vector<float> &vecz, float *meanx, float *meany, float *meanz, float *surf_meanx, float *surf_meany, float *surf_meanz, float *radius){

   int num_selected_pts = vecx.size();

   //Mean location of patch
   *meanx = accumulate(vecx.begin(), vecx.end(), 0.0)/num_selected_pts;
   *meany = accumulate(vecy.begin(), vecy.end(), 0.0)/num_selected_pts;
   *meanz = accumulate(vecz.begin(), vecz.end(), 0.0)/num_selected_pts;

   int ind = 0;
   float min_dist = (vecx[0] - *meanx)*(vecx[0] - *meanx) + (vecy[0] - *meany)*(vecy[0] - *meany) + (vecz[0] - *meanz)*(vecz[0] - *meanz);
   for(int j=0; j < num_selected_pts; j++){
      float dist = (vecx[j] - *meanx)*(vecx[j] - *meanx) + (vecy[j] - *meany)*(vecy[j] - *meany) + (vecz[j] - *meanz)*(vecz[j] - *meanz);
      if(dist < min_dist){
         min_dist = dist;
         ind = j;
      }
   }

   //Closest point on surface to mean
   *surf_meanx = vecx[ind];
   *surf_meany = vecy[ind];
   *surf_meanz = vecz[ind];

   //Compute radius of patch
   *radius = 0.0;
   for(int j=0; j < num_selected_pts; j++){
      float dist = (vecx[j] - *surf_meanx)*(vecx[j] - *surf_meanx) + (vecy[j] - *surf_meany)*(vecy[j] - *surf_meany) + (vecz[j] - *surf_meanz)*(vecz[j] - *surf_meanz);
      if(dist > *radius)
         *radius = dist;
   }
   *radius = sqrt(*radius);
}
void get_vertices(MeshModel &m, vector<unsigned int> &indices, vector<float> &vecx, vector<float> &vecy, vector<float> &vecz){

   int num_selected_pts = indices.size();
   vecx.resize(num_selected_pts);
   vecy.resize(num_selected_pts);
   vecz.resize(num_selected_pts);
   for(int k=0;k<num_selected_pts;k++){
      vecx[k] = m.cm.vert[indices[k]].P()[0];
      vecy[k] = m.cm.vert[indices[k]].P()[1];
      vecz[k] = m.cm.vert[indices[k]].P()[2];
   }
}
int index_first_selected(MeshModel &m){

   int i=0;
   CMeshO::VertexIterator vi;
   for(vi=m.cm.vert.begin(); vi!=m.cm.vert.end(); ++vi)
   {
      if(!(*vi).IsD() && (*vi).IsS()){
         (*vi).ClearS();
         break;
      }
      i++;
   }   
   return i;
}
//Return elements from indices where C=val
vector <unsigned int>  subset_indices(vector<int> indices, vector<int> C, int val){
   vector<unsigned int> sub_indices;
   for(int k=0;k<indices.size();k++){
      if(C[k] == val)
         sub_indices.push_back(indices[k]); 
   }
   return sub_indices;
}
//Return elements from indices where C=val
vector <unsigned int>  subset_indices(vector<unsigned int> &indices, vector<int> C, int val){
   vector<unsigned int> sub_indices;
   for(int k=0;k<indices.size();k++){
      if(C[k] == val)
         sub_indices.push_back(indices[k]); 
   }
   return sub_indices;
}
void color_patch(MeshModel &m, vector<unsigned int> indices, Color4b color){
   tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
   for(int k=0;k<indices.size();k++)
      m.cm.vert[indices[k]].SetS();
   tri::UpdateColor<CMeshO>::PerVertexConstant(m.cm, color, TRUE);
   tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
}
void color_ball(MeshModel &m, int ind, int k, Color4b color){

   VertexConstDataWrapper<CMeshO> wrapper(m.cm);
   KdTree<typename CMeshO::ScalarType> tree(wrapper);
   KdTree<typename CMeshO::ScalarType>::PriorityQueue queue;

   vector<unsigned int> points;
   tree.doQueryK(m.cm.vert[ind].cP(),k,queue);
   //tree.doQueryDist(m.cm.vert[ind].cP(),rad,points,dists);
   for(int j=0; j < queue.getNofElements(); j++)
      points.push_back(queue.getIndex(j));

   color_patch(m,points,color);
}

void MeshSegmentationFilterPlugin::initParameterSet(QAction *action, MeshModel &m, RichParameterSet &parlst)
{
   //Get *.ply mesh name
   char plyfile[1000];
   char out_file[1000];
   mesh_name(m,plyfile);
   sprintf(out_file,"%s/%s_MeshSegmentation_Parameters.txt",getenv("HOME"),plyfile);

   //Load parameters from file
   float default_radius, p;
   int num_nodes;
   vector<double> user_weights;
   tie(default_radius,num_nodes,p,user_weights) = load_params(out_file);

   static bool first_call = 1;
   if(first_call){
      tri::UpdateColor<CMeshO>::PerVertexConstant(m.cm, Color4b::LightGray, FALSE);
      first_call = 0;
   }
   int num_verts = m.cm.vert.size();
   int num_selected_pts = tri::UpdateSelection<CMeshO>::VertexCount(m.cm);

   if(num_selected_pts > 50){
      vector<unsigned int> indices = get_selected_indices(m);
      tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
      vector<float> vecx, vecy, vecz;
      float meanx, meany, meanz, surf_meanx, surf_meany, surf_meanz, radius;
      get_vertices(m, indices, vecx, vecy, vecz);
      patch_statistics(vecx, vecy, vecz, &meanx, &meany, &meanz, &surf_meanx, &surf_meany, &surf_meanz, &default_radius);
   }

   parlst.addParam(new RichFloat("Radius", default_radius, "Connection radius", "Connection radius for graph construction."));
   parlst.addParam(new RichInt("Nodes", min(num_verts,num_nodes), "Number of Nodes", "Number of nodes to use in graph construction."));
   parlst.addParam(new RichDynamicFloat("p",p,0.1,2,"Weight Matrix parameter", "Parameter controlling how the weight matrix is constructed"));

   QStringList FaceList;
   FaceList.push_back("Face 1:Red");
   FaceList.push_back("Face 2:Blue");
   FaceList.push_back("Face 3:Green");
   FaceList.push_back("Face 4:Magenta");
   FaceList.push_back("Face 5:Yellow");
   FaceList.push_back("Face 6:Cyan");
   FaceList.push_back("Face 7:Light Green");
   FaceList.push_back("Face 8:Light Red");
   FaceList.push_back("Face 9:Dark Green");
   FaceList.push_back("Face 10:Dark Red");
   FaceList.push_back("Face 11:Dark Blue");
   FaceList.push_back("Face 12:Black");
   FaceList.push_back("Face 13:White");
   FaceList.push_back("Face 14:Gray");
   parlst.addParam(new RichEnum("FaceIndex", 0, FaceList, tr("Face index (if adding point):"), QString("Face index, if you are adding new point.")));

   parlst.addParam(new RichBool("RunSeg",0,"Run segmentation", "Check when ready to run segmentation algorithm."));
   parlst.addParam(new RichBool("AdjustWeights",0,"Adjust weights (below)", "Toggles whether to enter weight adjustment mode."));

   parlst.addParam(new RichDynamicFloat("Face1", user_weights[0],-1,1,"Face 1 (Red)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face2", user_weights[1],-1,1,"Face 2 (Blue)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face3", user_weights[2],-1,1,"Face 3 (Green)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face4", user_weights[3],-1,1,"Face 4 (Magenta)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face5", user_weights[4],-1,1,"Face 5 (Yellow)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face6", user_weights[5],-1,1,"Face 6 (Cyan)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face7", user_weights[6],-1,1,"Face 7 (Light Green)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face8", user_weights[7],-1,1,"Face 8 (Light Red)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face9", user_weights[8],-1,1,"Face 9 (Dark Green)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face10",user_weights[9],-1,1,"Face 10 (Dark Red)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face11",user_weights[10],-1,1,"Face 11 (Dark Blue)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face12",user_weights[11],-1,1,"Face 12 (Black)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face13",user_weights[12],-1,1,"Face 13 (White)", "Parameter controlling how heavily to weight the face."));
   parlst.addParam(new RichDynamicFloat("Face14",user_weights[13],-1,1,"Face 14 (Gray)", "Parameter controlling how heavily to weight the face."));
}

void ColorMesh(MeshModel &m, vector<int> v){
   for(int i=0; i<=max(v); i++)
      color_patch(m, subset_indices(arange(m.cm.vert.size()),v,i), C[i]);
}
void ColorMesh(MeshModel &m, vector<double> v){
   vector<int> vv(v.size());
   for(int i=0; i<=v.size(); i++)
      vv[i] = (int)v[i];
   ColorMesh(m,vv);
}
void ColorMesh(MeshModel &m, vector<int> v, vector<int> label_map){
   for(int i=0; i<=max(v); i++)
      color_patch(m, subset_indices(arange(m.cm.vert.size()),v,i), C[label_map[i]]);
}
//Sets up graph on random set of vertices and constructs interpolation matrices, etc.
tuple<SparseMatrix, SparseMatrix, vector<unsigned int>, vector<unsigned int>, int> GraphSetup(MeshModel &m, int num_nodes, float radius){

   int i, j, num_verts = m.cm.vert.size();

   //Random sample of nodes
   vector<unsigned int> subset(num_verts);
   vector<bool> selected(num_verts,FALSE); 
   vector<unsigned int> selected_index(num_verts,-1); 
   for(i=0; i<num_verts; i++)
      subset[i] = i;
   
   random_device rd;
   mt19937 g(rd());
   shuffle(subset.begin(), subset.end(), g);
   subset.resize(num_nodes);
   
   for(i=0; i<num_nodes; i++){
      selected[subset[i]] = TRUE;
      selected_index[subset[i]] = i;
   }

   VertexConstDataWrapper<CMeshO> wrapper(m.cm);
   KdTree<typename CMeshO::ScalarType> tree(wrapper);

   SparseMatrix ND_VertsNodes(num_verts,num_nodes); //Normal distance between vertices and subset nodes
   SparseMatrix ND(num_nodes,num_nodes); //Normal distance between nodes

   vector<vector<unsigned int>> NN(num_nodes);
   vector<unsigned int> min_ind(num_verts,-1);
   vector<double> min_dist(num_verts,0);
   Point3m N[num_nodes];
   int min_nn = num_verts;
   for(i=0; i<num_nodes; i++){
      vector<float> dists;
      vector<unsigned int> points;
      tree.doQueryDist(m.cm.vert[subset[i]].cP(),radius,points,dists);
      ND_VertsNodes.insert(subset[i],i,1E-10);
      min_dist[subset[i]]=0.0;
      min_ind[subset[i]]=i;

      //Compute average normal
      float sum = 1.0;
      N[i] = m.cm.vert[subset[i]].N();
      Point3m P = m.cm.vert[subset[i]].N();
      for(j=0;j<points.size();j++){
         Point3m Nj = m.cm.vert[points[j]].N();
         Point3m Q = P - Nj;
         double nd = Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2];
         double weight = exp(-8*nd);
         sum += weight;
         N[i][0] += weight*Nj[0];
         N[i][1] += weight*Nj[1];
         N[i][2] += weight*Nj[2];
         if(min_ind[points[j]]==-1 || dists[j] < min_dist[points[j]]){
            min_ind[points[j]]=i;
            min_dist[points[j]]=dists[j];
         }
         if(selected[points[j]]){
            int k = selected_index[points[j]];
            NN[i].push_back(k);
         }
      }
      min_nn = min(min_nn,(int)NN[i].size());

      N[i][0]/=sum;
      N[i][1]/=sum;
      N[i][2]/=sum;
      for(j=0;j<points.size();j++){
         Point3m P = N[i] - m.cm.vert[points[j]].N();
         double nd = P[0]*P[0] + P[1]*P[1] + P[2]*P[2];
         ND_VertsNodes.insert(points[j],i,nd/4.0);
      }
   }

   for(i=0; i<num_nodes; i++){
      for(j=0; j<NN[i].size(); j++){
         int k = NN[i][j];
         Point3m Q = N[i] - N[k];
         double nd = Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2];
         ND.insert(i,k,nd/4.0);
      }
   }

   SparseMatrix I = sparse_exp(-32.0*ND_VertsNodes);
   I = I.row_divide(I.sum(1));
   
   return make_tuple(ND,I,subset,min_ind,min_nn);
}
vector<int> get_label_map(vector<int>& label_val){

   vector<int> label_map;
   vector<bool> label_present(C.size(),0);
   for(int i=0; i<label_val.size(); i++)
      label_present[label_val[i]] = 1;
   for(int i=0; i<label_present.size(); i++){
      if(label_present[i])
         label_map.push_back(i);
   }
   return label_map;
}

//Without balancing
Matrix PoissonForcing(int num_nodes, vector<int>& label_ind, vector<int>& label_val){

   int num_clusters = max(label_val)+1;

   Matrix F = Matrix::zeros(num_nodes,num_clusters);
   for(int i=0; i<label_ind.size(); i++)
      F(label_ind[i],label_val[i]) = 1.0;

   vector<double> c = F.sum(0)/label_ind.size();
   for(int i=0; i<label_ind.size(); i++)
      for(int j=0; j<F.get_cols(); j++)
         F(label_ind[i],j) = F(label_ind[i],j) - c[j];

   return F;
}

vector<double> SpectralCluster(SparseMatrix W){
   
   int n = W.get_cols();

   SparseMatrix L = SparseMatrix::spdiags(W.sum(1)) - W;

   double l = L.LargestEigenvalue(sqrt(n)*1E-8);
   SparseMatrix M = L - l*SparseMatrix::speye(n);

   Matrix u = Matrix::ones(n,1)/sqrt(n);

   return M.Largest_Eigenvector(u,sqrt(n)*1E-5).tovector();
}
   

bool MeshSegmentationFilterPlugin::applyFilter(QAction *action, MeshDocument &md, RichParameterSet & par, vcg::CallBackPos * /*cb*/)
{
   if (md.mm() == NULL)
      return false;

   float radius = par.getFloat("Radius");
   int num_nodes = par.getInt("Nodes");
   float p = par.getDynamicFloat("p");
   bool AdjustWeights = par.getBool("AdjustWeights");
   bool RunSeg = par.getBool("RunSeg");

   //Mesh
   MeshModel &m=*(md.mm());
   int num_verts = m.cm.vert.size();
   num_nodes = min(num_nodes,num_verts);

   //Check how many neighbors in radius ball
   VertexConstDataWrapper<CMeshO> wrapper(m.cm);
   KdTree<typename CMeshO::ScalarType> tree(wrapper);
   vector<float> dists;
   vector<unsigned int> points;
   tree.doQueryDist(m.cm.vert[0].cP(),radius,points,dists);
   int num_k = points.size();

  
   //Requirements
   tri::RequirePerVertexNormal(m.cm);
   tri::UpdateNormal<CMeshO>::PerVertexNormalized(m.cm);

   //Get *.ply mesh name
   char plyfile[1000];
   char out_file[1000];
   mesh_name(m,plyfile);
   sprintf(out_file,"%s/%s_MeshSegmentation_Parameters.txt",getenv("HOME"),plyfile);

   static bool first_time = 1;
   static bool RunOnce = 0;
   if(first_time && 0){
      //Test data
      label_ind_verts.push_back(95882);
      label_ind_verts.push_back(103879);
      label_val.push_back(0);
      label_val.push_back(1);
      first_time = 0;
   }

   static Matrix u;
   static SparseMatrix I;
   vector<double> face_weights;
   vector<double> user_weights(C.size(),1.0);

   //Get face weights and recolor if AdjustWeights=True
   user_weights[0] = par.getDynamicFloat("Face1");
   user_weights[1] = par.getDynamicFloat("Face2");
   user_weights[2] = par.getDynamicFloat("Face3");
   user_weights[3] = par.getDynamicFloat("Face4");
   user_weights[4] = par.getDynamicFloat("Face5");
   user_weights[5] = par.getDynamicFloat("Face6");
   user_weights[6] = par.getDynamicFloat("Face7");
   user_weights[7] = par.getDynamicFloat("Face8");
   user_weights[8] = par.getDynamicFloat("Face9");
   user_weights[9] = par.getDynamicFloat("Face10");
   user_weights[10] = par.getDynamicFloat("Face11");
   user_weights[11] = par.getDynamicFloat("Face12");
   user_weights[12] = par.getDynamicFloat("Face13");
   user_weights[13] = par.getDynamicFloat("Face14");

   vector<int> label_map = get_label_map(label_val);
   for(int i=0; i < label_map.size(); i++)
      face_weights.push_back(weight_map(user_weights[label_map[i]]));

   if(AdjustWeights && RunOnce){
      ColorMesh(m,argmax(I*u.col_multiply(face_weights),1),label_map);
      save_params(out_file, radius, num_nodes, p, user_weights);

      for(int i=0; i<label_val.size(); i++){
         color_ball(m,label_ind_verts[i],num_k/8,Color4b::Black);
         color_ball(m,label_ind_verts[i],num_k/16,C[label_val[i]]);
      }
      return true;
   }

   //Number of selected vertices
   int num_selected_pts = tri::UpdateSelection<CMeshO>::VertexCount(m.cm);

   if(num_selected_pts <= 50 && num_selected_pts >= 1)
   {
      //If points were selected, add to label_ind
      int FaceIndex = par.getEnum("FaceIndex");
      int ind = index_first_selected(m);
      color_ball(m,ind,num_k/8,Color4b::Black);
      color_ball(m,ind,num_k/16,C[FaceIndex]);
      label_ind_verts.push_back(ind);
      label_val.push_back(FaceIndex);

      if(!RunSeg){
         //Color all labeled points
         for(int i=0; i<label_val.size(); i++){
            color_ball(m,label_ind_verts[i],num_k/8,Color4b::Black);
            color_ball(m,label_ind_verts[i],num_k/16,C[label_val[i]]);
         }
      }
   }else if(label_map.size() >= 2)
   {
      if(!RunSeg){
         //Color all labeled points
         for(int i=0; i<label_val.size(); i++){
            color_ball(m,label_ind_verts[i],num_k/8,Color4b::Black);
            color_ball(m,label_ind_verts[i],num_k/16,C[label_val[i]]);
         }
         return true;
      }

      //Matrices associated with the graph
      vector<unsigned int> subset; //Indices in {1,num_verts} of selected subset
      vector<unsigned int> min_ind; //Index of closest neighbor in graph
      SparseMatrix ND; //Interpolation matrix I and pairwise normal distance matrix ND
      int min_nn;


      tie(ND,I,subset,min_ind,min_nn) = GraphSetup(m, num_nodes, radius);
      
      //Check if graph was setup correctly
      if(min_nn < 3){
         Log("Radius or number of nodes is too small for connectivity.");
      }else{
         Log("Graph initialized successfully.");

         SparseMatrix W, L; //Weight matrix and Laplacian matrix
         vector<double> deg, sdeg; //Degree vectors

         W = sparse_exp(-32*ND.pow(p/2.0));
         W = W - SparseMatrix::spdiags(W.diag());
         L = SparseMatrix::spdiags(W.sum(1)) - W;
         deg = W.sum(1);
         sdeg = sqrt(deg);
         
         //Transfer label indices to graph
         vector<int> label_ind(label_ind_verts.size());
         for(int i=0; i<label_ind.size(); i++)
            label_ind[i] = min_ind[label_ind_verts[i]];
         
         //Poisson Forcing 
         Matrix F = PoissonForcing(num_nodes, label_ind, label_val);
         vector<int> label_map = get_label_map(label_val);
         F = F(arange(num_nodes),label_map);
        
         //Preconditioning
         SparseMatrix Lp = L.col_divide(sdeg);
         Lp = Lp.row_divide(sdeg);
         F = F.row_divide(sdeg);

         //Conjugate gradient solver
         u = conjgrad(Lp, F, sqrt(num_nodes)*1E-10);
         u = u.row_divide(sdeg);
         
         //Color mesh
         ColorMesh(m,argmax(I*u.col_multiply(face_weights),1),label_map);
         for(int i=0; i<label_val.size(); i++){
            color_ball(m,label_ind_verts[i],num_k/8,Color4b::Black);
            color_ball(m,label_ind_verts[i],num_k/16,C[label_val[i]]);
         }
         RunOnce = 1;
         save_params(out_file, radius, num_nodes, p, user_weights);


//Clutering code
 /*              
               int num_clusters = par.getEnum("Clusters")+2;
               Log("Number of clusters = %d\n",num_clusters);
*/              
               //Laplacian matrix
               /*L = SparseMatrix::spdiags(W.sum(1)) - W;
               double l = L.LargestEigenvalue(sqrt(num_nodes)*1E-8);
               SparseMatrix M = L - l*SparseMatrix::speye(num_nodes);
               Matrix u = Matrix::ones(num_nodes,1)/sqrt(num_nodes);
               u = M.Largest_Eigenvector(u,sqrt(num_nodes)*1E-5);*/
/*               
               vector<double> v = SpectralCluster(W); 
               vector<bool> J = v>0;
               vector<double> w1 = SpectralCluster(W(J,J));
               vector<double> w2 = SpectralCluster(W(!J,!J));
               vector<double> w = splice(w1,w2,J);

               //vector<double> labels = 2*((I*v)>0) + 1*((I*w)>0);
               vector<double> labels = 2*(v>0) + 1*(w>0);
               ColorMesh(m,argmax(I*onehot(labels),1));
               //ColorMesh(m,labels);
               save(v,"vector.dat");

*/
/*               int source = rand() % num_nodes; 
               vector<double> F(num_nodes,-1/(double)num_nodes);
               F[source]++;
               vector<double> v = conjgrad(L, F, sqrt(num_nodes)*1E-8);
               double threshold = withness(v);
               ColorMesh(m,Heaviside(I*v-threshold));
               save(v,"vector.dat");
               print(v);
               printf("Thresh=%f\n",threshold); */

               //vector<double> v = L.Smallest_NonConst_Eigenvector(sqrt(num_nodes)*1E-5);
               //Matrix u = L.Smallest_NonConst_Eigenvectors(num_clusters-1,sqrt(num_nodes)*1E-5);
               //Spectral shift
               //SparseMatrix Lnorm = L.row_divide(sdeg);
               //Lnorm = Lnorm.row_divide(sdeg);
/*               double l = L.LargestEigenvalue(sqrt(num_nodes)*1E-5);
               SparseMatrix M = L - l*SparseMatrix::speye(num_nodes);

               Matrix u = Matrix::ones(num_nodes,1)/sqrt(num_nodes);

               for(int i=1; i<num_clusters+1; i++){
                  Matrix w = M.Largest_Eigenvector(u,sqrt(num_nodes)*1E-5);
                  //w = sign(w);
                  //w = w - w.mean();
                  //w = w/w.norm();
                  u.append_col(w);
               }

               //Normalize rows
               vector<double> V = u.norm(1);
               u = u.row_divide(V);
               u.save("EigenData.dat");
               /u = Heaviside(u);
               //u.print();
               //u.save("EigenData.dat");
               //vector<double> y(num_clusters,1.0);
               //for(int i=1; i<num_clusters; i++)
               //   y[i] = 2*y[i-1];

               //Remove zero columns of u
               //u = onehot(u*y);
               //vector<int> nonz = nonzero(u.norm(0));
               //u = u(arange(num_nodes),nonz);

               //Color mesh
               //ColorMesh(m,argmax(I*u,1));
               vector<int> labels = kmeans(u,num_clusters);
               ColorMesh(m,argmax(I*onehot(labels),1));*/

               //if(num_clusters == 2){

               /*
               int source = rand() % num_nodes; 
               vector<double> F(num_nodes,-1/(double)num_nodes);
               F[source]++;
               //vector<double> v = conjgrad(L, F, sqrt(num_nodes)*1E-8);
               SparseMatrix P = W.row_divide(deg);
               F = F/deg;
               vector<double> v(num_nodes,0.0);
               vector<double> u(num_nodes,0.0);
               for(int i=0; i<100000; i++){
                  u = P*v + F;
                  double err = norm(u-v);
                  printf("i:%d,err:%f\n",i,err);
                  v = u;
               }


               vector<int> labels = kmeans(v,num_clusters);
               ColorMesh(m,argmax(I*onehot(labels),1));*/
               /*double threshold = withness(v);
               ColorMesh(m,Heaviside(I*v-threshold));*/
               //save(v,"vector.dat");
               //}
               /*if(num_clusters >= 3){
                  Matrix F = Matrix::ones(num_nodes,num_clusters);
                  F = F/double(-num_nodes);
                  for(int i=0; i<num_clusters; i++){
                     int source = rand() % num_nodes; 
                     F(source,i) = F(source,i)+1;
                  }
                  Matrix u = conjgrad(L, F, sqrt(num_nodes)*1E-10);
                  vector<int> labels = kmeans(u,num_clusters);
                  ColorMesh(m,argmax(I*onehot(labels),1));
               }*/

        }
   }else{
      Log("Must supply points on at least 2 faces.");
   }

   return true;
}

MeshFilterInterface::FilterClass MeshSegmentationFilterPlugin::getClass(QAction *action)
{
  switch(ID(action))
  {
	   //This line controls which menu under Filters the plugin appears
	   case FP_NORMAL_MESH_SEGMENTATION: return FilterClass(MeshFilterInterface::Normal);
  }
  return MeshFilterInterface::Selection;
}

 int MeshSegmentationFilterPlugin::getRequirements(QAction *action)
{
 switch(ID(action))
  {
     case FP_NORMAL_MESH_SEGMENTATION: return MeshModel::MM_VERTCOLOR;

	  default: return MeshModel::MM_NONE;
  }
}

int MeshSegmentationFilterPlugin::postCondition(QAction *action) const
{
	switch(ID(action))
	{
   }
  return MeshModel::MM_ALL;
}

int MeshSegmentationFilterPlugin::getPreConditions( QAction * action) const
{
  switch(ID(action))
  {
  }
  return 0;
}
MESHLAB_PLUGIN_NAME_EXPORTER(MeshSegmentationFilterPlugin)
