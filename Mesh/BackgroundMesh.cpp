// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include "GmshMessage.h"
#include "BackgroundMesh.h"
#include "Numeric.h"
#include "Context.h"
#include "GVertex.h"
#include "GEdge.h"
#include "GFace.h"
#include "GModel.h"
#include "OS.h"
#include "Field.h"
#include "MElement.h"
#include "MElementOctree.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MQuadrangle.h"
#include "MVertex.h"

#if defined(HAVE_SOLVER)
#include "dofManager.h"
#include "laplaceTerm.h"
#include "linearSystemGMM.h"
#include "linearSystemCSR.h"
#include "linearSystemFull.h"
#include "linearSystemPETSc.h"
#endif

#if defined(HAVE_ANN)
static int _NBANN = 2;
#endif

void backgroundMesh::set(GFace *gf)
{
  if(_current) delete _current;
  _current = new backgroundMesh(gf);
}

void backgroundMesh::setCrossFieldsByDistance(GFace *gf)
{
  if(_current) delete _current;
  _current = new backgroundMesh(gf, true);
}

void backgroundMesh::unset()
{
  if(_current) delete _current;
  _current = 0;
}

backgroundMesh::backgroundMesh(GFace *_gf, bool cfd)
#if defined(HAVE_ANN)
  : _octree(0), uv_kdtree(0), nodes(0), angle_nodes(0), angle_kdtree(0)
#endif
{
  if(cfd) {
    Msg::Info("Building A Cross Field Using Closest Distance");
    propagateCrossFieldByDistance(_gf);
    return;
  }

  // create a bunch of triangles on the parametric space
  // those triangles are local to the backgroundMesh so that
  // they do not depend on the actual mesh that can be deleted

  std::set<SPoint2> myBCNodes;
  for(unsigned int i = 0; i < _gf->triangles.size(); i++) {
    MTriangle *e = _gf->triangles[i];
    MVertex *news[3];
    for(int j = 0; j < 3; j++) {
      MVertex *v = e->getVertex(j);
      std::map<MVertex *, MVertex *>::iterator it = _3Dto2D.find(v);
      MVertex *newv = 0;
      if(it == _3Dto2D.end()) {
        SPoint2 p;
        reparamMeshVertexOnFace(v, _gf, p);
        newv = new MVertex(p.x(), p.y(), 0.0);
        _vertices.push_back(newv);
        _3Dto2D[v] = newv;
        _2Dto3D[newv] = v;
        if(v->onWhat()->dim() < 2) myBCNodes.insert(p);
      }
      else
        newv = it->second;
      news[j] = newv;
    }
    MTriangle *T2D = new MTriangle(news[0], news[1], news[2]);
    _triangles.push_back(T2D);
  }

#if defined(HAVE_ANN)
  // printf("creating uv kdtree %d \n", myBCNodes.size());
  index = new ANNidx[2];
  dist = new ANNdist[2];
  nodes = annAllocPts(myBCNodes.size(), 3);
  std::set<SPoint2>::iterator itp = myBCNodes.begin();
  int ind = 0;
  while(itp != myBCNodes.end()) {
    SPoint2 pt = *itp;
    // fprintf(of, "SP(%g,%g,%g){%g};\n", pt.x(), pt.y(), 0.0, 10000);
    nodes[ind][0] = pt.x();
    nodes[ind][1] = pt.y();
    nodes[ind][2] = 0.0;
    itp++;
    ind++;
  }
  uv_kdtree = new ANNkd_tree(nodes, myBCNodes.size(), 3);
#endif

  // build a search structure
  _octree = new MElementOctree(_triangles);

  // compute the mesh sizes at nodes
  if(CTX::instance()->mesh.lcFromPoints) {
    propagate1dMesh(_gf);
  }
  else {
    std::map<MVertex *, MVertex *>::iterator itv2 = _2Dto3D.begin();
    for(; itv2 != _2Dto3D.end(); ++itv2) {
      _sizes[itv2->first] = CTX::instance()->mesh.lcMax;
    }
  }
  // ensure that other criteria are fullfilled
  updateSizes(_gf);

  // compute optimal mesh orientations
  propagateCrossField(_gf);

  _3Dto2D.clear();
  _2Dto3D.clear();
}

backgroundMesh::~backgroundMesh()
{
  for(unsigned int i = 0; i < _vertices.size(); i++) delete _vertices[i];
  for(unsigned int i = 0; i < _triangles.size(); i++) delete _triangles[i];
  if(_octree) delete _octree;
#if defined(HAVE_ANN)
  if(uv_kdtree) delete uv_kdtree;
  if(angle_kdtree) delete angle_kdtree;
  if(nodes) annDeallocPts(nodes);
  if(angle_nodes) annDeallocPts(angle_nodes);
  delete[] index;
  delete[] dist;
#endif
}

static void propagateValuesOnFace(GFace *_gf,
                                  std::map<MVertex *, double> &dirichlet,
                                  simpleFunction<double> *ONE,
                                  bool in_parametric_plane = false)
{
#if defined(HAVE_SOLVER)
  linearSystem<double> *_lsys = 0;
#if defined(HAVE_PETSC)
  _lsys = new linearSystemPETSc<double>;
#elif defined(HAVE_GMM)
  linearSystemGmm<double> *_lsysb = new linearSystemGmm<double>;
  _lsysb->setGmres(1);
  _lsys = _lsysb;
#else
  _lsys = new linearSystemFull<double>;
#endif

  dofManager<double> myAssembler(_lsys);

  // fix boundary conditions
  std::map<MVertex *, double>::iterator itv = dirichlet.begin();
  for(; itv != dirichlet.end(); ++itv) {
    myAssembler.fixVertex(itv->first, 0, 1, itv->second);
  }

  // Number vertices
  std::set<MVertex *> vs;
  for(unsigned int k = 0; k < _gf->triangles.size(); k++)
    for(int j = 0; j < 3; j++) vs.insert(_gf->triangles[k]->getVertex(j));
  for(unsigned int k = 0; k < _gf->quadrangles.size(); k++)
    for(int j = 0; j < 4; j++) vs.insert(_gf->quadrangles[k]->getVertex(j));

  std::map<MVertex *, SPoint3> theMap;
  if(in_parametric_plane) {
    for(std::set<MVertex *>::iterator it = vs.begin(); it != vs.end(); ++it) {
      SPoint2 p;
      reparamMeshVertexOnFace(*it, _gf, p);
      theMap[*it] = SPoint3((*it)->x(), (*it)->y(), (*it)->z());
      (*it)->setXYZ(p.x(), p.y(), 0.0);
    }
  }

  for(std::set<MVertex *>::iterator it = vs.begin(); it != vs.end(); ++it)
    myAssembler.numberVertex(*it, 0, 1);

  // Assemble
  laplaceTerm l(0, 1, ONE);
  for(unsigned int k = 0; k < _gf->triangles.size(); k++) {
    MTriangle *t = _gf->triangles[k];
    SElement se(t);
    l.addToMatrix(myAssembler, &se);
  }

  // Solve
  if(myAssembler.sizeOfR()) {
    _lsys->systemSolve();
  }

  // save solution
  for(std::set<MVertex *>::iterator it = vs.begin(); it != vs.end(); ++it) {
    myAssembler.getDofValue(*it, 0, 1, dirichlet[*it]);
  }

  if(in_parametric_plane) {
    for(std::set<MVertex *>::iterator it = vs.begin(); it != vs.end(); ++it) {
      SPoint3 p = theMap[(*it)];
      (*it)->setXYZ(p.x(), p.y(), p.z());
    }
  }
  delete _lsys;
#endif
}

void backgroundMesh::propagate1dMesh(GFace *_gf)
{
  std::vector<GEdge *> const &e = _gf->edges();
  std::vector<GEdge *>::const_iterator it = e.begin();
  std::map<MVertex *, double> sizes;

  for(; it != e.end(); ++it) {
    if(!(*it)->isSeam(_gf)) {
      for(unsigned int i = 0; i < (*it)->lines.size(); i++) {
        MVertex *v1 = (*it)->lines[i]->getVertex(0);
        MVertex *v2 = (*it)->lines[i]->getVertex(1);
        if(v1 != v2) {
          double d = sqrt((v1->x() - v2->x()) * (v1->x() - v2->x()) +
                          (v1->y() - v2->y()) * (v1->y() - v2->y()) +
                          (v1->z() - v2->z()) * (v1->z() - v2->z()));
          for(int k = 0; k < 2; k++) {
            MVertex *v = (*it)->lines[i]->getVertex(k);
            std::map<MVertex *, double>::iterator itv = sizes.find(v);
            if(itv == sizes.end())
              sizes[v] = log(d);
            else
              itv->second = 0.5 * (itv->second + log(d));
          }
        }
      }
    }
  }

  simpleFunction<double> ONE(1.0);
  propagateValuesOnFace(_gf, sizes, &ONE);

  std::map<MVertex *, MVertex *>::iterator itv2 = _2Dto3D.begin();
  for(; itv2 != _2Dto3D.end(); ++itv2) {
    MVertex *v_2D = itv2->first;
    MVertex *v_3D = itv2->second;
    _sizes[v_2D] = exp(sizes[v_3D]);
  }
}

crossField2d::crossField2d(MVertex *v, GEdge *ge)
{
  double p;
  bool success = reparamMeshVertexOnEdge(v, ge, p);
  if(!success) {
    Msg::Warning("cannot reparametrize a point in crossField");
    _angle = 0;
    return;
  }
  SVector3 t = ge->firstDer(p);
  t.normalize();
  _angle = atan2(t.y(), t.x());
  crossField2d::normalizeAngle(_angle);
}

void backgroundMesh::propagateCrossFieldByDistance(GFace *_gf)
{
  std::vector<GEdge *> const &e = _gf->edges();
  std::vector<GEdge *>::const_iterator it = e.begin();
  std::map<MVertex *, double> _cosines4, _sines4;
  std::map<MVertex *, SPoint2> _param;

  for(; it != e.end(); ++it) {
    if(!(*it)->isSeam(_gf)) {
      for(unsigned int i = 0; i < (*it)->lines.size(); i++) {
        MVertex *v[2];
        v[0] = (*it)->lines[i]->getVertex(0);
        v[1] = (*it)->lines[i]->getVertex(1);
        SPoint2 p1, p2;
        reparamMeshEdgeOnFace(v[0], v[1], _gf, p1, p2);
        /* a correct way of computing angles  */
        Pair<SVector3, SVector3> der = _gf->firstDer((p1 + p2) * .5);
        SVector3 t1 = der.first();
        SVector3 t2(v[1]->x() - v[0]->x(), v[1]->y() - v[0]->y(),
                    v[1]->z() - v[0]->z());
        t1.normalize();
        t2.normalize();
        double _angle = angle(t1, t2);
        //        double angle = atan2 ( p1.y()-p2.y() , p1.x()-p2.x() );
        crossField2d::normalizeAngle(_angle);
        for(int i = 0; i < 2; i++) {
          std::map<MVertex *, double>::iterator itc = _cosines4.find(v[i]);
          std::map<MVertex *, double>::iterator its = _sines4.find(v[i]);
          if(itc != _cosines4.end()) {
            itc->second = 0.5 * (itc->second + cos(4 * _angle));
            its->second = 0.5 * (its->second + sin(4 * _angle));
          }
          else {
            _param[v[i]] = (i == 0) ? p1 : p2;
            _cosines4[v[i]] = cos(4 * _angle);
            _sines4[v[i]] = sin(4 * _angle);
          }
        }
      }
    }
  }

#if defined(HAVE_ANN)
  index = new ANNidx[_NBANN];
  dist = new ANNdist[_NBANN];
  angle_nodes = annAllocPts(_cosines4.size(), 3);
  std::map<MVertex *, double>::iterator itp = _cosines4.begin();
  int ind = 0;
  _sin.clear();
  _cos.clear();
  while(itp != _cosines4.end()) {
    MVertex *v = itp->first;
    double c = itp->second;
    SPoint2 pt = _param[v];
    double s = _sines4[v];
    angle_nodes[ind][0] = pt.x();
    angle_nodes[ind][1] = pt.y();
    angle_nodes[ind][2] = 0.0;
    _cos.push_back(c);
    _sin.push_back(s);
    itp++;
    ind++;
  }
  angle_kdtree = new ANNkd_tree(angle_nodes, _cosines4.size(), 3);
#endif
}

inline double myAngle(const SVector3 &a, const SVector3 &b, const SVector3 &d)
{
  double cosTheta = dot(a, b);
  double sinTheta = dot(crossprod(a, b), d);
  return atan2(sinTheta, cosTheta);
}

// smoothness = h * (|grad (cos 4 a)| + |grad (sin 4 a)|)
// smoothness is of order 1 if not smooth
// smoothness is of order h/L if smooth
// h --> mesh size
// L --> domain size
double backgroundMesh::getSmoothness(MElement *e)
{
  MVertex *v0 = _3Dto2D[e->getVertex(0)];
  MVertex *v1 = _3Dto2D[e->getVertex(1)];
  MVertex *v2 = _3Dto2D[e->getVertex(2)];
  std::map<MVertex *, double>::const_iterator i0 = _angles.find(v0);
  std::map<MVertex *, double>::const_iterator i1 = _angles.find(v1);
  std::map<MVertex *, double>::const_iterator i2 = _angles.find(v2);
  double a[3] = {cos(4 * i0->second), cos(4 * i1->second), cos(4 * i2->second)};
  double b[3] = {sin(4 * i0->second), sin(4 * i1->second), sin(4 * i2->second)};
  //      printf("coucou\n");
  double f[3];
  e->interpolateGrad(a, 0, 0, 0, f);
  const double gradcos = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
  e->interpolateGrad(b, 0, 0, 0, f);
  // const double gradsin = sqrt (f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
  const double h = e->maxEdge();
  return (gradcos /*+ gradsin*/) * h;
}

double backgroundMesh::getSmoothness(double u, double v, double w)
{
  MElement *e = _octree->find(u, v, w, 2, true);
  if(!e) return -1.0;
  MVertex *v0 = e->getVertex(0);
  MVertex *v1 = e->getVertex(1);
  MVertex *v2 = e->getVertex(2);
  std::map<MVertex *, double>::const_iterator i0 = _angles.find(v0);
  std::map<MVertex *, double>::const_iterator i1 = _angles.find(v1);
  std::map<MVertex *, double>::const_iterator i2 = _angles.find(v2);
  double a[3] = {cos(4 * i0->second), cos(4 * i1->second), cos(4 * i2->second)};
  double b[3] = {sin(4 * i0->second), sin(4 * i1->second), sin(4 * i2->second)};
  //      printf("coucou\n");
  double f[3];
  e->interpolateGrad(a, 0, 0, 0, f);
  const double gradcos = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
  e->interpolateGrad(b, 0, 0, 0, f);
  // const double gradsin = sqrt (f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
  const double h = e->maxEdge();
  return (gradcos /*+ gradsin*/) * h;
}

void backgroundMesh::propagateCrossField(GFace *_gf)
{
  propagateCrossFieldHJ(_gf);
  // solve the non liear problem
  constantPerElement<double> C;
  int ITER = 0;
  //  int NSMOOTH = _gf->triangles.size();
  while(0) {
    //    int NSMOOTH_NOW = 0;
    for(unsigned int i = 0; i < _gf->triangles.size(); i++) {
      double smoothness = getSmoothness(_gf->triangles[i]);
      double val = smoothness < .5 ? 1.0 : 1.e-3; // exp(-absf/10);
      C.set(_gf->triangles[i], val);
    }
    //    if (NSMOOTH_NOW == NSMOOTH) break;
    //    NSMOOTH = NSMOOTH_NOW;
    //    break;
    _angles.clear();
    propagateCrossField(_gf, &C);
    if(++ITER > 0) break;
  }
    //  printf("converged in %d iterations\n", ITER);
#if 0 // debug print
  char name[256];
  sprintf(name, "cross-%d-%d.pos", _gf->tag(), ITER);
  print(name, 0, 1);
  sprintf(name, "smooth-%d-%d.pos", _gf->tag(), ITER);
  print(name, _gf, 2);
#endif
}

void backgroundMesh::propagateCrossFieldHJ(GFace *_gf)
{
  simpleFunction<double> ONE(1.0);
  propagateCrossField(_gf, &ONE);
}

void backgroundMesh::propagateCrossField(GFace *_gf,
                                         simpleFunction<double> *ONE)
{
  std::map<MVertex *, double> _cosines4, _sines4;
  std::vector<GEdge *> const &e = _gf->edges();
  std::vector<GEdge *>::const_iterator it = e.begin();
  for(; it != e.end(); ++it) {
    if(!(*it)->isSeam(_gf)) {
      for(unsigned int i = 0; i < (*it)->lines.size(); i++) {
        MVertex *v[2];
        v[0] = (*it)->lines[i]->getVertex(0);
        v[1] = (*it)->lines[i]->getVertex(1);
        SPoint2 p1, p2;
        reparamMeshEdgeOnFace(v[0], v[1], _gf, p1, p2);
        Pair<SVector3, SVector3> der = _gf->firstDer((p1 + p2) * .5);
        SVector3 t1 = der.first();
        SVector3 t2 = der.second();
        SVector3 n = crossprod(t1, t2);
        n.normalize();
        SVector3 d1(v[1]->x() - v[0]->x(), v[1]->y() - v[0]->y(),
                    v[1]->z() - v[0]->z());
        t1.normalize();
        d1.normalize();
        double _angle = myAngle(t1, d1, n);
        crossField2d::normalizeAngle(_angle);
        for(int i = 0; i < 2; i++) {
          std::map<MVertex *, double>::iterator itc = _cosines4.find(v[i]);
          std::map<MVertex *, double>::iterator its = _sines4.find(v[i]);
          if(itc != _cosines4.end()) {
            itc->second = 0.5 * (itc->second + cos(4 * _angle));
            its->second = 0.5 * (its->second + sin(4 * _angle));
          }
          else {
            _cosines4[v[i]] = cos(4 * _angle);
            _sines4[v[i]] = sin(4 * _angle);
          }
        }
      }
    }
  }

  propagateValuesOnFace(_gf, _cosines4, ONE, false);
  propagateValuesOnFace(_gf, _sines4, ONE, false);

  //    print("cos4.pos",0,_cosines4,0);
  //    print("sin4.pos",0,_sines4,0);

  std::map<MVertex *, MVertex *>::iterator itv2 = _2Dto3D.begin();
  for(; itv2 != _2Dto3D.end(); ++itv2) {
    MVertex *v_2D = itv2->first;
    MVertex *v_3D = itv2->second;
    double angle = atan2(_sines4[v_3D], _cosines4[v_3D]) / 4.0;
    crossField2d::normalizeAngle(angle);
    _angles[v_2D] = angle;
  }
}

void backgroundMesh::updateSizes(GFace *_gf)
{
  std::map<MVertex *, double>::iterator itv = _sizes.begin();
  for(; itv != _sizes.end(); ++itv) {
    SPoint2 p;
    MVertex *v = _2Dto3D[itv->first];
    double lc;
    if(v->onWhat()->dim() == 0) {
      lc = BGM_MeshSize(v->onWhat(), 0, 0, v->x(), v->y(), v->z());
    }
    else if(v->onWhat()->dim() == 1) {
      double u;
      v->getParameter(0, u);
      lc = BGM_MeshSize(v->onWhat(), u, 0, v->x(), v->y(), v->z());
    }
    else {
      reparamMeshVertexOnFace(v, _gf, p);
      lc = BGM_MeshSize(_gf, p.x(), p.y(), v->x(), v->y(), v->z());
    }
    // printf("2D -- %g %g 3D -- %g %g\n",p.x(),p.y(),v->x(),v->y());
    itv->second = std::min(lc, itv->second);
    itv->second = std::max(itv->second, CTX::instance()->mesh.lcMin);
    itv->second = std::min(itv->second, CTX::instance()->mesh.lcMax);
  }
  // do not allow large variations in the size field
  // (Int. J. Numer. Meth. Engng. 43, 1143-1165 (1998) MESH GRADATION
  // CONTROL, BOROUCHAKI, HECHT, FREY)
  std::set<MEdge, Less_Edge> edges;
  for(unsigned int i = 0; i < _triangles.size(); i++) {
    for(int j = 0; j < _triangles[i]->getNumEdges(); j++) {
      edges.insert(_triangles[i]->getEdge(j));
    }
  }
  const double _beta = 1.3;
  for(int i = 0; i < 3; i++) {
    std::set<MEdge, Less_Edge>::iterator it = edges.begin();
    for(; it != edges.end(); ++it) {
      MVertex *v0 = it->getVertex(0);
      MVertex *v1 = it->getVertex(1);
      MVertex *V0 = _2Dto3D[v0];
      MVertex *V1 = _2Dto3D[v1];
      std::map<MVertex *, double>::iterator s0 = _sizes.find(V0);
      std::map<MVertex *, double>::iterator s1 = _sizes.find(V1);
      if(s0->second < s1->second)
        s1->second = std::min(s1->second, _beta * s0->second);
      else
        s0->second = std::min(s0->second, _beta * s1->second);
    }
  }
}

bool backgroundMesh::inDomain(double u, double v, double w) const
{
  return _octree->find(u, v, w, 2, true) != 0;
}

double backgroundMesh::operator()(double u, double v, double w) const
{
  double uv[3] = {u, v, w};
  double uv2[3];
  MElement *e = _octree->find(u, v, w, 2, true);
  if(!e) {
#if defined(HAVE_ANN)
    // printf("BGM octree not found --> find in kdtree \n");
    if(uv_kdtree->nPoints() < 2) return -1000.;
    double pt[3] = {u, v, 0.0};
    uv_kdtree->annkSearch(pt, 2, index, dist);
    SPoint3 p1(nodes[index[0]][0], nodes[index[0]][1], nodes[index[0]][2]);
    SPoint3 p2(nodes[index[1]][0], nodes[index[1]][1], nodes[index[1]][2]);
    SPoint3 pnew;
    double d;
    signedDistancePointLine(p1, p2, SPoint3(u, v, 0.), d, pnew);
    e = _octree->find(pnew.x(), pnew.y(), 0.0, 2, true);
#endif
    if(!e) {
      Msg::Error("BGM octree: cannot find UVW=%g %g %g", u, v, w);
      return -1000.0; // 0.4;
    }
  }
  e->xyz2uvw(uv, uv2);
  std::map<MVertex *, double>::const_iterator itv1 =
    _sizes.find(e->getVertex(0));
  std::map<MVertex *, double>::const_iterator itv2 =
    _sizes.find(e->getVertex(1));
  std::map<MVertex *, double>::const_iterator itv3 =
    _sizes.find(e->getVertex(2));
  return itv1->second * (1 - uv2[0] - uv2[1]) + itv2->second * uv2[0] +
         itv3->second * uv2[1];
}

double backgroundMesh::getAngle(double u, double v, double w) const
{
  // JFR :
  // we can use closest point for computing
  // cross field angles : this allow NOT to
  // generate a spurious mesh and solve a PDE
  if(!_octree) {
#if defined(HAVE_ANN)
    double angle = 0.;
    if(angle_kdtree->nPoints() >= _NBANN) {
      double pt[3] = {u, v, 0.0};
      angle_kdtree->annkSearch(pt, _NBANN, index, dist);
      double SINE = 0.0, COSINE = 0.0;
      for(int i = 0; i < _NBANN; i++) {
        SINE += _sin[index[i]];
        COSINE += _cos[index[i]];
        //      printf("%2d %2d %12.5E
        //      %12.5E\n",i,index[i],_sin[index[i]],_cos[index[i]]);
      }
      angle = atan2(SINE, COSINE) / 4.0;
    }
    crossField2d::normalizeAngle(angle);
    return angle;
#endif
  }

  // HACK FOR LEWIS
  // h = 1+30(y-x^2)^2  + (1-x)^2
  //  double x = u;
  //  double y = v;
  //  double dhdx = 30 * 2 * (y-x*x) * (-2*x) - 2 * (1-x);
  //  double dhdy = 30 * 2 * (y-x*x);
  //  double angles = atan2(y,x*x);
  //  crossField2d::normalizeAngle (angles);
  //  return angles;

  double uv[3] = {u, v, w};
  double uv2[3];
  MElement *e = _octree->find(u, v, w, 2, true);
  if(!e) {
#if defined(HAVE_ANN)
    // printf("BGM octree not found --> find in kdtree \n");
    if(uv_kdtree->nPoints() < 2) return -1000.0;
    double pt[3] = {u, v, 0.0};
    uv_kdtree->annkSearch(pt, 2, index, dist);
    SPoint3 p1(nodes[index[0]][0], nodes[index[0]][1], nodes[index[0]][2]);
    SPoint3 p2(nodes[index[1]][0], nodes[index[1]][1], nodes[index[1]][2]);
    SPoint3 pnew;
    double d;
    signedDistancePointLine(p1, p2, SPoint3(u, v, 0.), d, pnew);
    e = _octree->find(pnew.x(), pnew.y(), 0., 2, true);
#endif
    if(!e) {
      Msg::Error("BGM octree angle: cannot find UVW=%g %g %g", u, v, w);
      return -1000.0;
    }
  }
  e->xyz2uvw(uv, uv2);
  std::map<MVertex *, double>::const_iterator itv1 =
    _angles.find(e->getVertex(0));
  std::map<MVertex *, double>::const_iterator itv2 =
    _angles.find(e->getVertex(1));
  std::map<MVertex *, double>::const_iterator itv3 =
    _angles.find(e->getVertex(2));

  double cos4 = cos(4 * itv1->second) * (1 - uv2[0] - uv2[1]) +
                cos(4 * itv2->second) * uv2[0] + cos(4 * itv3->second) * uv2[1];
  double sin4 = sin(4 * itv1->second) * (1 - uv2[0] - uv2[1]) +
                sin(4 * itv2->second) * uv2[0] + sin(4 * itv3->second) * uv2[1];
  double angle = atan2(sin4, cos4) / 4.0;
  crossField2d::normalizeAngle(angle);

  return angle;
}

void backgroundMesh::print(const std::string &filename, GFace *gf,
                           const std::map<MVertex *, double> &_whatToPrint,
                           int smooth)
{
  FILE *f = Fopen(filename.c_str(), "w");
  if(!f) {
    Msg::Error("Could not open file '%s'", filename.c_str());
    return;
  }
  fprintf(f, "View \"Background Mesh\"{\n");
  if(smooth) {
    for(unsigned int i = 0; i < gf->triangles.size(); i++) {
      MVertex *v1 = gf->triangles[i]->getVertex(0);
      MVertex *v2 = gf->triangles[i]->getVertex(1);
      MVertex *v3 = gf->triangles[i]->getVertex(2);
      double x = getSmoothness(gf->triangles[i]);
      fprintf(f, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n", v1->x(),
              v1->y(), v1->z(), v2->x(), v2->y(), v2->z(), v3->x(), v3->y(),
              v3->z(), x, x, x);
    }
  }
  else {
    for(unsigned int i = 0; i < _triangles.size(); i++) {
      MVertex *v1 = _triangles[i]->getVertex(0);
      MVertex *v2 = _triangles[i]->getVertex(1);
      MVertex *v3 = _triangles[i]->getVertex(2);
      std::map<MVertex *, double>::const_iterator itv1 = _whatToPrint.find(v1);
      std::map<MVertex *, double>::const_iterator itv2 = _whatToPrint.find(v2);
      std::map<MVertex *, double>::const_iterator itv3 = _whatToPrint.find(v3);
      if(!gf) {
        fprintf(f, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n", v1->x(),
                v1->y(), v1->z(), v2->x(), v2->y(), v2->z(), v3->x(), v3->y(),
                v3->z(), itv1->second, itv2->second, itv3->second);
      }
      else {
        GPoint p1 = gf->point(SPoint2(v1->x(), v1->y()));
        GPoint p2 = gf->point(SPoint2(v2->x(), v2->y()));
        GPoint p3 = gf->point(SPoint2(v3->x(), v3->y()));
        fprintf(f, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n", p1.x(),
                p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), p3.x(), p3.y(), p3.z(),
                itv1->second, itv2->second, itv3->second);
      }
    }
  }
  fprintf(f, "};\n");
  fclose(f);
}

MElementOctree *backgroundMesh::get_octree() { return _octree; }

MElement *backgroundMesh::getMeshElementByCoord(double u, double v, double w,
                                                bool strict)
{
  if(!_octree) {
    Msg::Debug("Rebuilding BackgroundMesh element octree");
    _octree = new MElementOctree(_triangles);
  }
  return _octree->find(u, v, w, 2, strict);
}

backgroundMesh *backgroundMesh::_current = 0;
