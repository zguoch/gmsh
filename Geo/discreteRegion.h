// Gmsh - Copyright (C) 1997-2020 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef DISCRETE_REGION_H
#define DISCRETE_REGION_H

#include "GModel.h"
#include "GRegion.h"
#include "MFace.h"

class discreteRegion : public GRegion {
public:
  discreteRegion(GModel *model, int num);
  discreteRegion(GModel *model);
  virtual ~discreteRegion() {}
  virtual GeomType geomType() const { return DiscreteVolume; }
  void setBoundFaces(const std::set<int> &tagFaces);
  void setBoundFaces(const std::vector<int> &tagFaces,
                     const std::vector<int> &signFaces);
  void findFaces(std::map<MFace, std::vector<int>, MFaceLessThan> &map_faces);
  virtual void remesh();
};

#endif
