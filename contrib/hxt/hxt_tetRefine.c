#include "hxt_tetRefine.h"
#include "predicates.h"
#include "hxt_tetFlag.h"


HXTStatus hxtEmptyMesh(HXTMesh* mesh, HXTDelaunayOptions* delOptions)
{
// we assume that the input is a surface mesh
  if (mesh->tetrahedra.num)  
    return HXT_ERROR_MSG(HXT_STATUS_FAILED, "The input mesh should only contain triangles");
  if (mesh->triangles.num == 0)  
    return HXT_ERROR_MSG(HXT_STATUS_FAILED, "The input mesh should contain triangles");

  double minDist2 = DBL_MAX;
  #pragma omp parallel for reduction(min:minDist2)
  for (uint64_t i=0; i<mesh->triangles.num; i++){
    uint32_t* node = mesh->triangles.node + 3*i;
    for (int j=0; j<3; j++) {
      double* n1 = mesh->vertices.coord + (size_t) 4*node[j];
      double* n2 = mesh->vertices.coord + (size_t) 4*node[(j+1)%3];

      double dist2 = (n1[0]-n2[0])*(n1[0]-n2[0])
                   + (n1[1]-n2[1])*(n1[1]-n2[1])
                   + (n1[2]-n2[2])*(n1[2]-n2[2]);

      if(dist2<minDist2)
        minDist2 = dist2;
    }
  }

  hxtNodeInfo* nodeInfo;
  HXT_CHECK( hxtAlignedMalloc(&nodeInfo, sizeof(hxtNodeInfo)*mesh->vertices.num) );

  #pragma omp parallel for simd aligned(nodeInfo:SIMD_ALIGN)
  for (uint32_t i=0; i<mesh->vertices.num; i++) {
    nodeInfo[i].node = i;
    nodeInfo[i].status = HXT_STATUS_TRYAGAIN;
  }

  HXT_CHECK( hxtDelaunaySteadyVertices(mesh, delOptions, nodeInfo, mesh->vertices.num) );
  delOptions->numVerticesInMesh = mesh->vertices.num;

#ifdef DEBUG
  #pragma omp parallel for simd aligned(nodeInfo:SIMD_ALIGN)
  for (uint32_t i=0; i<mesh->vertices.num; i++) {
    if(nodeInfo[i].status!=HXT_STATUS_TRUE){
      HXT_WARNING("vertex %u of the empty mesh was not inserted\n", nodeInfo[i].node);
    }
  }
#endif

  HXT_CHECK( hxtAlignedFree(&nodeInfo) );
  
  return HXT_STATUS_OK;
}



// refine 


double hxtTetCircumcenter(double a[3], double b[3], double c[3], double d[3],
                            double circumcenter[3], double *xi, double *eta, double *zeta)
{
  double xba, yba, zba, xca, yca, zca, xda, yda, zda;
  double balength, calength, dalength;
  double xcrosscd, ycrosscd, zcrosscd;
  double xcrossdb, ycrossdb, zcrossdb;
  double xcrossbc, ycrossbc, zcrossbc;
  double denominator;
  double xcirca, ycirca, zcirca;

  /* Use coordinates relative to point `a' of the tetrahedron. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  zba = b[2] - a[2];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  zca = c[2] - a[2];
  xda = d[0] - a[0];
  yda = d[1] - a[1];
  zda = d[2] - a[2];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba + zba * zba;
  calength = xca * xca + yca * yca + zca * zca;
  dalength = xda * xda + yda * yda + zda * zda;
  /* Cross products of these edges. */
  xcrosscd = yca * zda - yda * zca;
  ycrosscd = zca * xda - zda * xca;
  zcrosscd = xca * yda - xda * yca;
  xcrossdb = yda * zba - yba * zda;
  ycrossdb = zda * xba - zba * xda;
  zcrossdb = xda * yba - xba * yda;
  xcrossbc = yba * zca - yca * zba;
  ycrossbc = zba * xca - zca * xba;
  zcrossbc = xba * yca - xca * yba;

  /* Calculate the denominator of the formulae. */
  /* Use orient3d() from http://www.cs.cmu.edu/~quake/robust.html     */
  /*   to ensure a correctly signed (and reasonably accurate) result, */
  /*   avoiding any possibility of division by zero.                  */
  const double xxx =  orient3d(b, c, d, a);
  denominator = 0.5 / xxx;

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) *
    denominator;
  ycirca = (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) *
    denominator;
  zcirca = (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) *
    denominator;
  circumcenter[0] =  xcirca + a[0];
  circumcenter[1] =  ycirca + a[1];
  circumcenter[2] =  zcirca + a[2];

  /*
 printf(" %g %g %g %g\n",
         sqrt((a[0]-xcirca)*(a[0]-xcirca)+(a[1]-ycirca)*(a[1]-ycirca)+(a[2]-zcirca)*(a[2]-zcirca)),
         sqrt((b[0]-xcirca)*(b[0]-xcirca)+(b[1]-ycirca)*(b[1]-ycirca)+(b[2]-zcirca)*(b[2]-zcirca)),
         sqrt((c[0]-xcirca)*(c[0]-xcirca)+(c[1]-ycirca)*(c[1]-ycirca)+(c[2]-zcirca)*(c[2]-zcirca)),
         sqrt((d[0]-xcirca)*(d[0]-xcirca)+(d[1]-ycirca)*(d[1]-ycirca)+(d[2]-zcirca)*(d[2]-zcirca)) );
  */

  if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a    */
    /*   coordinate system with a xi-axis directed from `a' to `b',      */
    /*   an eta-axis directed from `a' to `c', and a zeta-axis directed  */
    /*   from `a' to `d'.  The values for xi, eta, and zeta are computed */
     /*   by Cramer's Rule for solving systems of linear equations.       */
    *xi = (xcirca * xcrosscd + ycirca * ycrosscd + zcirca * zcrosscd) *
      (2.0 * denominator);
    *eta = (xcirca * xcrossdb + ycirca * ycrossdb + zcirca * zcrossdb) *
      (2.0 * denominator);
    *zeta = (xcirca * xcrossbc + ycirca * ycrossbc + zcirca * zcrossbc) *
      (2.0 * denominator);
  }
  return xxx;
}


static inline double square_dist(double v0[3], double v1[3])
{
  return (v1[0] - v0[0])*(v1[0] - v0[0]) +
         (v1[1] - v0[1])*(v1[1] - v0[1]) +
         (v1[2] - v0[2])*(v1[2] - v0[2]);
}


static int is_too_close(HXTMesh* mesh, double tetNodeSize, double vtaSize,
                        double squareDist)
{
  if(tetNodeSize!=DBL_MAX) {
    double meanSize = 0.5*(tetNodeSize+vtaSize);
    if(squareDist < /*(0.94*0.94) * */meanSize * meanSize) {
      return 1;
    }
  }

  return 0;
}


HXTStatus hxtRefineTetrahedra(HXTMesh* mesh, HXTDelaunayOptions* delOptions,
                              double (*meshSizeFun)(double x, double y, double z,
                                                    void* meshSizeData),
                              void* meshSizeData)
{
  int maxThreads = omp_get_max_threads();
  uint32_t* numCreated;
  uint64_t* startIndex;
  HXT_CHECK( hxtMalloc(&numCreated, maxThreads*sizeof(uint32_t)) );
  HXT_CHECK( hxtMalloc(&startIndex, (maxThreads+1)*sizeof(uint64_t)) );

  for(int iter=0; iter<42; iter++) {
    double *newVertices;

    /******************************************************************
     * share the creation of new vertices between threads
     *****************************************************************/
    uint64_t numToProcessPerThread;
    {
      // get the number of unprocessed tet to avid wasting memory ressources...
      uint64_t hxtDeclareAligned64 numToProcessSIMD[4] = {0};
      for(uint64_t i=0; i<mesh->tetrahedra.num/4; i++) {
        numToProcessSIMD[0] += getProcessedFlag(mesh, i*4+0)==0;
        numToProcessSIMD[1] += getProcessedFlag(mesh, i*4+1)==0;
        numToProcessSIMD[2] += getProcessedFlag(mesh, i*4+2)==0;
        numToProcessSIMD[3] += getProcessedFlag(mesh, i*4+3)==0;
      }

      for(uint64_t i=mesh->tetrahedra.num/4*4; i<mesh->tetrahedra.num; i++) {
        numToProcessSIMD[0] += getProcessedFlag(mesh, i)==0;
      }

      // each thread will have to process a certain number of tet
      numToProcessPerThread =((numToProcessSIMD[0]+numToProcessSIMD[1]) +
                              (numToProcessSIMD[2]+numToProcessSIMD[3]))/maxThreads + 1;

      uint64_t numToProcess = 0;
      uint64_t step = 0;
      int threadID = 0;
      for(uint64_t i=0; i<mesh->tetrahedra.num; i++) {
        if(numToProcess==step) {
          startIndex[threadID++] = i;
          step += numToProcessPerThread;
        }
        numToProcess += getProcessedFlag(mesh, i)==0;
      }

      while(threadID <= maxThreads) {
        startIndex[threadID++] = mesh->tetrahedra.num;
      }
      
      HXT_CHECK( hxtAlignedMalloc(&newVertices, sizeof(double)*4*numToProcess) );
    }



    
    uint32_t add = 0;
    HXTStatus status = HXT_STATUS_OK;
    #pragma omp parallel reduction(+:add)
    {
      int threadID = omp_get_thread_num();

      uint64_t vertexStart = numToProcessPerThread * threadID;
      uint32_t localAdd = 0;

      for(uint64_t i=startIndex[threadID]; i<startIndex[threadID+1]; i++) {

        if (getProcessedFlag(mesh, i)==0){
          int pointAdded = 0;
          double* p[4];
          p[0] = mesh->vertices.coord + (size_t) 4*mesh->tetrahedra.node[4*i+0];
          p[1] = mesh->vertices.coord + (size_t) 4*mesh->tetrahedra.node[4*i+1];
          p[2] = mesh->vertices.coord + (size_t) 4*mesh->tetrahedra.node[4*i+2];
          p[3] = mesh->vertices.coord + (size_t) 4*mesh->tetrahedra.node[4*i+3];
          double newCoord [3];
          double bary[4];
          double vtaSize;
          setProcessedFlag(mesh, i); // we do not need to refine that tetrahedra anymore

          double tetNodeSize[4];
          for(int j=0; j<4; j++) {
            tetNodeSize[j] = delOptions->nodalSizes[mesh->tetrahedra.node[4*i+j]];
          }

          // try to add the point at the circumcenter
          hxtTetCircumcenter(p[0], p[1], p[2], p[3], newCoord, &bary[1], &bary[2], &bary[3]);
          bary[0] = 1.0 - bary[1] - bary[2] - bary[3];

          // we don't insert at circumcenter if it is outside the tet.
          if(bary[0] > 0 && bary[1] > 0 && bary[2] > 0 && bary[3] > 0) {

            // computing the mesh size
            if(meshSizeFun!=NULL) {
              vtaSize = meshSizeFun(newCoord[0], newCoord[1], newCoord[2], meshSizeData);
            }
            else { // we suppose delOptions->nodalSize!=NULL
              double denom = 0.0;
              double num = 0.0;
              for (int j=0;j<4;j++){
                double size = delOptions->nodalSizes[mesh->tetrahedra.node[4*i+j]];
                if (size != DBL_MAX){
                  if(bary[j]>0.999) {
                    num = size;
                    denom = 1;
                    break;
                  }
                  double weight = bary[j];
                  denom += weight;
                  num += weight * size;
                }
              }

              vtaSize = num/denom;
            }

            double circumradius2 = square_dist(p[0], newCoord);
            
            if(!is_too_close(mesh, tetNodeSize[0], vtaSize, circumradius2) &&
               !is_too_close(mesh, tetNodeSize[1], vtaSize, circumradius2) &&
               !is_too_close(mesh, tetNodeSize[2], vtaSize, circumradius2) &&
               !is_too_close(mesh, tetNodeSize[3], vtaSize, circumradius2)) {
              size_t vertexIndex = vertexStart + localAdd;
              newVertices[vertexIndex*4  ] = newCoord[0];
              newVertices[vertexIndex*4+1] = newCoord[1];
              newVertices[vertexIndex*4+2] = newCoord[2];
              newVertices[vertexIndex*4+3] = vtaSize;
              pointAdded = 1;
              localAdd++;
            }
          }
          
          // try to add the point at barycentric coordinates
          // ponderated by the inverse of the mesh size
          if(!pointAdded){
            double denom = 0.0;
            int ndef = 0;
            for (int j=0; j<4; j++)
            {
              double size = tetNodeSize[j];
              if(size!=DBL_MAX && size!=0.0) {
                bary[j] = 1./size;
                denom += bary[j];
              }
              else {
                bary[j] = -1.0;
                ndef++;
              }
            }

            if(ndef && ndef!=4) {
              denom *= 4.0/(4.0-ndef);
            }

            for (int j=0; j<4; j++)
            {
              if(bary[j]==-1.0) {
                bary[j] = 0.25;
              }
              else {
                bary[j] /= denom;
              }
            }

            for (int j=0; j<3; j++)
            {
              newCoord[j] = bary[0]*p[0][j] +
                            bary[1]*p[1][j] +
                            bary[2]*p[2][j] +
                            bary[3]*p[3][j];
            }

            // this does not work with ffast-math, GCC completely ignore this...
  #if !defined(__FAST_MATH__) && !defined(NDEBUG)
            if(!isfinite(newCoord[0]) || !isfinite(newCoord[1]) || !isfinite(newCoord[2])){
              HXT_ERROR_MSG(HXT_STATUS_ERROR, "new coordinates are not finite");
              exit(EXIT_FAILURE);
            }
  #endif

            if(meshSizeFun!=NULL) {
              vtaSize = meshSizeFun(newCoord[0], newCoord[1], newCoord[2], meshSizeData);
            }
            else {
              if(ndef==4)
                vtaSize = DBL_MAX;
              else {
                vtaSize = 4.0/denom;
              }
            }

            if(!is_too_close(mesh, tetNodeSize[0], vtaSize, square_dist(p[0], newCoord)) &&
               !is_too_close(mesh, tetNodeSize[1], vtaSize, square_dist(p[1], newCoord)) &&
               !is_too_close(mesh, tetNodeSize[2], vtaSize, square_dist(p[2], newCoord)) &&
               !is_too_close(mesh, tetNodeSize[3], vtaSize, square_dist(p[3], newCoord))){
              size_t vertexIndex = vertexStart + localAdd;
              newVertices[vertexIndex*4  ] = newCoord[0];
              newVertices[vertexIndex*4+1] = newCoord[1];
              newVertices[vertexIndex*4+2] = newCoord[2];
              newVertices[vertexIndex*4+3] = vtaSize;
              pointAdded = 1;
              localAdd++;
            }
          }
        }
      }


      numCreated[threadID] = localAdd;

      #pragma omp barrier
      #pragma omp single
      {
        int nthreads = omp_get_num_threads();
        add = 0;
        for (int i=0; i<nthreads; i++) {
          uint32_t tsum = add + numCreated[i];
          numCreated[i] = add;
          add = tsum;
        }

        if(mesh->vertices.num + add>mesh->vertices.size){
          status=hxtAlignedRealloc(&mesh->vertices.coord, sizeof(double)*4*(mesh->vertices.num + add));
          if(status==HXT_STATUS_OK){
            status=hxtAlignedRealloc(&delOptions->nodalSizes, (mesh->vertices.num + add)*sizeof(double));
            mesh->vertices.size = mesh->vertices.num + add;
          }
        }
      }

      size_t v = numCreated[threadID] + mesh->vertices.num;

      if(status==HXT_STATUS_OK){
        for (uint64_t i=0; i<localAdd; i++){
          size_t vertexIndex = vertexStart + i;
          mesh->vertices.coord[v*4  ] = newVertices [4*vertexIndex  ];
          mesh->vertices.coord[v*4+1] = newVertices [4*vertexIndex+1];
          mesh->vertices.coord[v*4+2] = newVertices [4*vertexIndex+2];
          delOptions->nodalSizes[v] = newVertices [4*vertexIndex+3];
          v++;
        }
      }
    }

    if(status!=HXT_STATUS_OK){
      HXT_TRACE(status);
      return status;
    }

    mesh->vertices.num += add;

    HXT_CHECK(hxtAlignedFree(&newVertices));
    
    HXT_CHECK(hxtDelaunay(mesh, delOptions));

    uint32_t numAdd = mesh->vertices.num - delOptions->numVerticesInMesh;
    delOptions->numVerticesInMesh = mesh->vertices.num;

    if (numAdd == 0) break;
  }

  HXT_CHECK( hxtFree(&startIndex) );
  HXT_CHECK( hxtFree(&numCreated) );

  return HXT_STATUS_OK;
}

