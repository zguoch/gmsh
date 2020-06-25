#include "hxt_vertices.h"


HXTStatus hxtCreateNodalSizeFromFunction(HXTMesh* mesh, double** nodalSizes_ptr,
                                         double (*meshSizeFun)(double x, double y, double z,
                                                               void* meshSizeData),
                                         void* meshSizeData)
{
  HXT_CHECK(hxtAlignedMalloc(nodalSizes_ptr,mesh->vertices.num*sizeof(double)));
  double* nodalSizes = *nodalSizes_ptr;

  #pragma omp parallel for
  for (uint32_t i=0; i<mesh->vertices.num; i++) {
    double* coord = &mesh->vertices.coord[4*i];
    nodalSizes[i] = meshSizeFun(coord[0], coord[1], coord[2], meshSizeData);
  }

  return HXT_STATUS_OK;
}


HXTStatus hxtCreateNodalsizeFromTrianglesAndLines(HXTMesh* mesh, double** nodalSizes_ptr)
{
  HXTVertex* vertices = (HXTVertex*) mesh->vertices.coord;

  HXT_CHECK(hxtAlignedMalloc(nodalSizes_ptr,mesh->vertices.num*sizeof(double)));
  double* nodalSizes = *nodalSizes_ptr;
  
  #pragma omp parallel for
  for (uint32_t i = 0; i<mesh->vertices.num; i++){
    nodalSizes[i] = 0;
    vertices[i].padding.hilbertDist = 0; // we use that as a counter to do the average...
  }

  // only do for triangles
  // we do not take into account hereafter nodalSizes = to DBL_MAX
  // could be changed in another fashion
  for (uint32_t i = 0; i<mesh->triangles.num; i++){
    for (uint32_t j = 0; j<3; j++){  
      for (uint32_t k = j+1; k<3; k++){  
        uint32_t n1 = mesh->triangles.node[3*i+j];
        uint32_t n2 = mesh->triangles.node[3*i+k];
        if (n1 != HXT_GHOST_VERTEX && n2 != HXT_GHOST_VERTEX){
          double *X1 = vertices[n1].coord;
          double *X2 = vertices[n2].coord;
          vertices[n1].padding.hilbertDist++;
          vertices[n2].padding.hilbertDist++;
          double l = sqrt ((X1[0]-X2[0])*(X1[0]-X2[0])+
                           (X1[1]-X2[1])*(X1[1]-X2[1])+
                           (X1[2]-X2[2])*(X1[2]-X2[2]));
          nodalSizes[n1] += l;
          nodalSizes[n2] += l;
        }
      }
    }
  }

  for (uint32_t i = 0; i<mesh->lines.num; i++){
      uint32_t n1 = mesh->lines.node[2*i+0];
      uint32_t n2 = mesh->lines.node[2*i+1];
      if (n1 != HXT_GHOST_VERTEX && n2 != HXT_GHOST_VERTEX && n1!=n2){
        double *X1 = vertices[n1].coord;
        double *X2 = vertices[n2].coord;
        vertices[n1].padding.hilbertDist++;
        vertices[n2].padding.hilbertDist++;
        double l = sqrt ((X1[0]-X2[0])*(X1[0]-X2[0])+
                         (X1[1]-X2[1])*(X1[1]-X2[1])+
                         (X1[2]-X2[2])*(X1[2]-X2[2]));
        nodalSizes[n1] += l;
        nodalSizes[n2] += l;
    }
  }

  #pragma omp parallel for
  for (uint32_t i=0; i<mesh->vertices.num; i++)
  {
    if(vertices[i].padding.hilbertDist == 0) {
      nodalSizes[i] = DBL_MAX;
    }
    else {
      nodalSizes[i] /= (double) vertices[i].padding.hilbertDist;
    }
  }
  return HXT_STATUS_OK;    
}

HXTStatus hxtCreateNodalsizeFromMesh(HXTMesh* mesh, double** nodalSizes_ptr)
{

  HXT_CHECK(hxtAlignedMalloc(nodalSizes_ptr,mesh->vertices.num*sizeof(double)));
  double* nodalSizes = *nodalSizes_ptr;
  
  #pragma omp parallel for
  for (uint32_t i = 0; i<mesh->vertices.num; i++){
    nodalSizes[i] = DBL_MAX;
  }

  // only do for triangles
  // we do not take into account hereafter nodalSizes = to DBL_MAX
  // could be changed in another fashion
  for (uint32_t i = 0; i<mesh->tetrahedra.num; i++){
    for (uint32_t j = 0; j<4; j++){  
      for (uint32_t k = j+1; k<4; k++){  
        uint32_t n1 = mesh->tetrahedra.node[4*i+j];
        uint32_t n2 = mesh->tetrahedra.node[4*i+k];
        if (n1 != HXT_GHOST_VERTEX && n2 != HXT_GHOST_VERTEX){
          double *X1 = mesh->vertices.coord + (size_t) 4*n1;
          double *X2 = mesh->vertices.coord + (size_t) 4*n2;
          double l = sqrt ((X1[0]-X2[0])*(X1[0]-X2[0])+
               (X1[1]-X2[1])*(X1[1]-X2[1])+
               (X1[2]-X2[2])*(X1[2]-X2[2]));
          if(l<nodalSizes[n1]) nodalSizes[n1] = l;
          if(l<nodalSizes[n2]) nodalSizes[n2] = l;
        }
      }
    }
  }
  return HXT_STATUS_OK;
}

HXTStatus hxtDestroyNodalsize(double** nodalSizes_ptr)
{
  HXT_CHECK( hxtAlignedFree(nodalSizes_ptr) );
  return HXT_STATUS_OK;
}