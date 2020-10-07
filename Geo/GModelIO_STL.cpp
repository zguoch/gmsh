// Gmsh - Copyright (C) 1997-2020 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include <stdio.h>
#include <string>
#include <algorithm>
#include <sstream>
#include "GModel.h"
#include "OS.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MQuadrangle.h"
#include "MVertexRTree.h"
#include "discreteFace.h"
#include "StringUtils.h"
#include "Context.h"

static bool invalidChar(char c) { return !(c >= 32 && c <= 126); }

int GModel::readSTL(const std::string &name, double tolerance)
{
  FILE *fp = Fopen(name.c_str(), "rb");
  if(!fp) {
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  // store triplets of points for each solid found in the file
  std::vector<std::vector<SPoint3> > points;
  SBoundingBox3d bbox;
  std::vector<std::string> names;

  // "solid", or binary data header
  char buffer[256];
  if(!fgets(buffer, sizeof(buffer), fp)) {
    fclose(fp);
    return 0;
  }

  //SPoint3 p0(1.9e6, 4e6, 0);

  bool binary = strncmp(buffer, "solid", 5) && strncmp(buffer, "SOLID", 5);

  // ASCII STL
  if(!binary) {
    if(strlen(buffer) > 6)
      names.push_back(&buffer[6]);
    else
      names.push_back("");
    points.resize(1);
    while(!feof(fp)) {
      // "facet normal x y z" or "endsolid"
      if(!fgets(buffer, sizeof(buffer), fp)) break;
      if(!strncmp(buffer, "endsolid", 8) || !strncmp(buffer, "ENDSOLID", 8)) {
        // "solid"
        if(!fgets(buffer, sizeof(buffer), fp)) break;
        if(!strncmp(buffer, "solid", 5) || !strncmp(buffer, "SOLID", 5)) {
          if(strlen(buffer) > 6)
            names.push_back(&buffer[6]);
          else
            names.push_back("");
          points.resize(points.size() + 1);
          // "facet normal x y z"
          if(!fgets(buffer, sizeof(buffer), fp)) break;
        }
      }
      // "outer loop"
      if(!fgets(buffer, sizeof(buffer), fp)) break;
      // "vertex x y z"
      for(int i = 0; i < 3; i++) {
        if(!fgets(buffer, sizeof(buffer), fp)) break;
        char s1[256];
        double x, y, z;
        if(sscanf(buffer, "%s %lf %lf %lf", s1, &x, &y, &z) != 4) break;
        SPoint3 p(x, y, z);
        //p -= p0;
        points.back().push_back(p);
        bbox += p;
      }
      // "endloop"
      if(!fgets(buffer, sizeof(buffer), fp)) break;
      // "endfacet"
      if(!fgets(buffer, sizeof(buffer), fp)) break;
    }
  }

  // check if we could parse something
  bool empty = true;
  for(std::size_t i = 0; i < points.size(); i++) {
    if(points[i].size()) {
      empty = false;
      break;
    }
  }
  if(empty) points.clear();

  // binary STL (we also try to read in binary mode if the header told
  // us the format was ASCII but we could not read any vertices)
  if(binary || empty) {
    if(binary)
      Msg::Info("Mesh is in binary format");
    else
      Msg::Info("Wrong ASCII header or empty file: trying binary read");
    rewind(fp);
    while(!feof(fp)) {
      char header[80];
      if(!fread(header, sizeof(char), 80, fp)) break;
      unsigned int nfacets = 0;
      size_t ret = fread(&nfacets, sizeof(unsigned int), 1, fp);
      bool swap = false;
      if(nfacets > 100000000) {
        Msg::Info("Swapping bytes from binary file");
        swap = true;
        SwapBytes((char *)&nfacets, sizeof(unsigned int), 1);
      }
      if(ret && nfacets) {
        names.push_back(header);
        points.resize(points.size() + 1);
        char *data = new char[nfacets * 50 * sizeof(char)];
        ret = fread(data, sizeof(char), nfacets * 50, fp);
        if(ret == nfacets * 50) {
          for(std::size_t i = 0; i < nfacets; i++) {
            float *xyz = (float *)&data[i * 50 * sizeof(char)];
            if(swap) SwapBytes((char *)xyz, sizeof(float), 12);
            for(int j = 0; j < 3; j++) {
              SPoint3 p(xyz[3 + 3 * j], xyz[3 + 3 * j + 1], xyz[3 + 3 * j + 2]);
              //p -= p0;
              points.back().push_back(p);
              bbox += p;
            }
          }
        }
        delete[] data;
      }
    }
  }

  // cleanup names
  if(names.size() != points.size()){
    Msg::Debug("Invalid number of names in STL file - should never happen");
    names.resize(points.size());
  }
  for(std::size_t i = 0; i < names.size(); i++){
    names[i].erase(remove_if(names[i].begin(), names[i].end(), invalidChar),
                   names[i].end());
  }

  std::vector<GFace *> faces;
  for(std::size_t i = 0; i < points.size(); i++) {
    if(points[i].empty()) {
      Msg::Error("No facets found in STL file for solid %d %s", i,
                 names[i].c_str());
      fclose(fp);
      return 0;
    }
    if(points[i].size() % 3) {
      Msg::Error("Wrong number of points (%d) in STL file for solid %d %s",
                 points[i].size(), i, names[i].c_str());
      fclose(fp);
      return 0;
    }
    Msg::Info("%d facets in solid %d %s", points[i].size() / 3, i,
              names[i].c_str());
    // create face
    GFace *face = new discreteFace(this, getMaxElementaryNumber(2) + 1);
    faces.push_back(face);
    add(face);
    if(!names[i].empty()) setElementaryName(2, face->tag(), names[i]);
  }

  // create triangles using unique vertices
  double eps = norm(SVector3(bbox.max(), bbox.min())) * tolerance;
  std::vector<MVertex *> vertices;
  for(std::size_t i = 0; i < points.size(); i++)
    for(std::size_t j = 0; j < points[i].size(); j++)
      vertices.push_back(
        new MVertex(points[i][j].x(), points[i][j].y(), points[i][j].z()));
  MVertexRTree pos(eps);
  pos.insert(vertices);

  std::set<MFace, MFaceLessThan> unique;
  int nbDuplic = 0, nbDegen = 0;
  for(std::size_t i = 0; i < points.size(); i++) {
    for(std::size_t j = 0; j < points[i].size(); j += 3) {
      MVertex *v[3];
      for(int k = 0; k < 3; k++) {
        double x = points[i][j + k].x();
        double y = points[i][j + k].y();
        double z = points[i][j + k].z();
        v[k] = pos.find(x, y, z);
        if(!v[k])
          Msg::Error("Could not find node at position (%g, %g, %g) with tol=%g",
                     x, y, z, eps);
      }
      if(!v[0] || !v[1] || !v[2]) {
        // error
      }
      else if(v[0] == v[1] || v[0] == v[2] || v[1] == v[2]) {
        Msg::Debug("Skipping degenerated triangle %lu %lu %lu",
                   v[0]->getNum(), v[1]->getNum(), v[2]->getNum());
        nbDegen++;
      }
      else if(CTX::instance()->mesh.stlRemoveDuplicateTriangles) {
        MFace mf(v[0], v[1], v[2]);
        if(unique.find(mf) == unique.end()) {
          faces[i]->triangles.push_back(new MTriangle(v[0], v[1], v[2]));
          unique.insert(mf);
        }
        else {
          nbDuplic++;
        }
      }
      else{
        faces[i]->triangles.push_back(new MTriangle(v[0], v[1], v[2]));
      }
    }
  }
  if(nbDuplic || nbDegen)
    Msg::Warning("%d duplicate/%d degenerate triangles in STL file",
                 nbDuplic, nbDegen);

  _associateEntityWithMeshVertices();

  _storeVerticesInEntities(vertices); // will delete unused vertices

  fclose(fp);
  return 1;
}

static void writeSTLfaces(FILE *fp, std::vector<GFace*> &faces, bool binary,
                          double scalingFactor, const std::string &name)
{
  bool useGeoSTL = false;
  unsigned int nfacets = 0;
  for(std::vector<GFace*>::iterator it = faces.begin(); it != faces.end(); ++it) {
    nfacets += (*it)->triangles.size() + 2 * (*it)->quadrangles.size();
  }
  if(!nfacets) { // use CAD STL if there is no mesh
    useGeoSTL = true;
    for(std::vector<GFace*>::iterator it = faces.begin(); it != faces.end(); ++it) {
      (*it)->buildSTLTriangulation();
      nfacets += (*it)->stl_triangles.size() / 3;
    }
  }

  if(!binary) {
    fprintf(fp, "solid %s\n", name.c_str());
  }
  else {
    char header[80];
    strncpy(header, name.c_str(), 79);
    header[79] = '\0';
    fwrite(header, sizeof(char), 80, fp);
    fwrite(&nfacets, sizeof(unsigned int), 1, fp);
  }

  for(std::vector<GFace*>::iterator it = faces.begin(); it != faces.end(); ++it) {
    if(useGeoSTL && (*it)->stl_vertices_uv.size()) {
      for(std::size_t i = 0; i < (*it)->stl_triangles.size(); i += 3) {
        SPoint2 &p1((*it)->stl_vertices_uv[(*it)->stl_triangles[i]]);
        SPoint2 &p2((*it)->stl_vertices_uv[(*it)->stl_triangles[i + 1]]);
        SPoint2 &p3((*it)->stl_vertices_uv[(*it)->stl_triangles[i + 2]]);
        GPoint gp1 = (*it)->point(p1);
        GPoint gp2 = (*it)->point(p2);
        GPoint gp3 = (*it)->point(p3);
        double x[3] = {gp1.x(), gp2.x(), gp3.x()};
        double y[3] = {gp1.y(), gp2.y(), gp3.y()};
        double z[3] = {gp1.z(), gp2.z(), gp3.z()};
        double n[3];
        normal3points(x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2], n);
        if(!binary) {
          fprintf(fp, "\tfacet normal %g %g %g\n", n[0], n[1], n[2]);
          fprintf(fp, "  \touter loop\n");
          for(int j = 0; j < 3; j++)
            fprintf(fp, "    \tvertex %g %g %g\n", x[j] * scalingFactor,
                    y[j] * scalingFactor, z[j] * scalingFactor);
          fprintf(fp, "  \tendloop\n");
          fprintf(fp, "\tendfacet\n");
        }
        else {
          char data[50];
          float *coords = (float *)data;
          coords[0] = (float)n[0];
          coords[1] = (float)n[1];
          coords[2] = (float)n[2];
          for(int j = 0; j < 3; j++) {
            coords[3 + 3 * j] = (float)(x[j] * scalingFactor);
            coords[3 + 3 * j + 1] = (float)(y[j] * scalingFactor);
            coords[3 + 3 * j + 2] = (float)(z[j] * scalingFactor);
          }
          data[48] = data[49] = 0;
          fwrite(data, sizeof(char), 50, fp);
        }
      }
    }
    else {
      for(std::size_t i = 0; i < (*it)->triangles.size(); i++)
        (*it)->triangles[i]->writeSTL(fp, binary, scalingFactor);
      for(std::size_t i = 0; i < (*it)->quadrangles.size(); i++)
        (*it)->quadrangles[i]->writeSTL(fp, binary, scalingFactor);
    }
  }

  if(!binary) fprintf(fp, "endsolid %s\n", name.c_str());
}

int GModel::writeSTL(const std::string &name, bool binary, bool saveAll,
                     double scalingFactor, int oneSolidPerSurface)
{
  FILE *fp = Fopen(name.c_str(), binary ? "wb" : "w");
  if(!fp) {
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }
  bool writePhysicalVolumes_Surfaces = false; //是否使用自己定义的方式输出stl文件
  if(noPhysicalGroups()) saveAll = true;
  if(oneSolidPerSurface == 1){ // one solid per surface
    
    for(fiter it = firstFace(); it != lastFace(); ++it) {
      
      if(saveAll || (*it)->physicals.size()) {
        std::vector<GFace*> faces(1, *it);
        std::string name = getElementaryName(2, (*it)->tag());
        if(name.empty()){
          std::ostringstream s;
          s << "Gmsh Surface " << (*it)->tag();
          name = s.str();
        }
        writeSTLfaces(fp, faces, binary, scalingFactor, name);
      }
    }
  }
  else if(oneSolidPerSurface == 2){ // one solid per physical surface
    std::map<int, std::vector<GEntity *> > phys;
    getPhysicalGroups(2, phys);
    for(std::map<int, std::vector<GEntity *> >::iterator it = phys.begin();
        it != phys.end(); it++){
      std::vector<GFace*> faces;
      for(std::size_t i = 0; i < it->second.size(); i++){
        faces.push_back(static_cast<GFace*>(it->second[i]));
      }
      std::string name = getPhysicalName(2, it->first);
      if(name.empty()){
        std::ostringstream s;
        s << "Gmsh Physical Surface " << it->first;
        name = s.str();
      }
      writeSTLfaces(fp, faces, binary, scalingFactor, name);
    }
  }
  else{ // one solid
    std::vector<GFace*> faces;
    for(fiter it = firstFace(); it != lastFace(); ++it) {
      if(saveAll || (*it)->physicals.size()) {
        faces.push_back(*it);
      }
    }
    std::string name = "Created by Gmsh";
    writeSTLfaces(fp, faces, binary, scalingFactor, name);
    writePhysicalVolumes_Surfaces=true; //如果没有指定Mesh.oneSolidPerSurface，则按照自己的方式输出stl文件
  }
  fclose(fp);
  if(writePhysicalVolumes_Surfaces)
  {
      const std::string path = SplitFileName(name)[0];
      const std::string baseName = SplitFileName(name)[1];
      const std::string extName = SplitFileName(name)[2];  // with dot, e.g. .stl
      std::vector<std::string> surfaceName, volumeName; 
      // 1. 获取所有的Physical Groups
      std::map<int, std::vector<GEntity *> > physicals[4]; //0，1，2，3分别表示点，线，面，体
      getPhysicalGroups(physicals);
      std::map<int, std::vector<GEntity *> >::const_iterator it;
      // 2. 处理Physical Volume
      const int dimVolume = 3;
      // open file, save all physical volumes into a single stl file
      std::string fnameRegion=path+baseName+"_Volume"+extName;
      FILE *fp_volume = Fopen(fnameRegion.c_str(), binary ? "wb" : "w");
      if(!fp_volume) {
        Msg::Error("Unable to open file '%s'", fnameRegion.c_str());
        return 0;
      }
      for(it=physicals[dimVolume].begin();it!=physicals[dimVolume].end();++it)
      {
        std::string nameRegion = getPhysicalName(dimVolume, it->first);       //获取某个physical Volume 的name
        std::vector<int> tags=getTagsForPhysicalName(dimVolume,nameRegion);   //获取给定physical name的所有entries的tag值
        std::vector<GFace *> faces;
        for (size_t i = 0; i < tags.size(); i++)
        {
          GRegion* region=getRegionByTag(tags[i]);
          std::vector<GFace *> faces_region=region->faces();
          for (size_t iface = 0; iface < faces_region.size(); iface++)
          {
            faces.push_back(faces_region[iface]);
          }
          Msg::Info("Physical Volume: %s , the %dth, faces %d",nameRegion.c_str(),i,faces_region.size());
        }
        // 将每个physical volume都写入一个stl文件
          std::string fnameRegion=path+baseName+"_"+nameRegion+extName;
          FILE *fp_region = Fopen(fnameRegion.c_str(), binary ? "wb" : "w");
          if(!fp_region) {
            Msg::Error("Unable to open file '%s'", fnameRegion.c_str());
            return 0;
          }
          
          writeSTLfaces(fp_volume, faces, binary, scalingFactor, nameRegion);
          writeSTLfaces(fp_region, faces, binary, scalingFactor, nameRegion);
          Msg::Info("Volume(%s) has %d faces, done writing %s",nameRegion.c_str(),faces.size(),fnameRegion.c_str());
          fclose(fp_region);
          // volume name
          volumeName.push_back(nameRegion);
      }
      fclose(fp_volume);
      // 3. 处理Physical Surfaces
      const int dimSurface = 2;
      std::string fnameSurface=path+baseName+"_Surface"+extName;
      FILE *fp_surface = Fopen(fnameSurface.c_str(), binary ? "wb" : "w");
      if(!fp_surface) {
        Msg::Error("Unable to open file '%s'", fnameSurface.c_str());
        return 0;
      }
      for(it=physicals[dimSurface].begin();it!=physicals[dimSurface].end();++it)
      {
        std::string nameSurface = getPhysicalName(dimSurface, it->first);
        std::vector<int> tags=getTagsForPhysicalName(dimSurface,nameSurface); //某个physical name对应的surface可能不止一个
        std::vector<GFace *> faces;
        for (size_t i = 0; i < tags.size(); i++)
        {
          faces.push_back(getFaceByTag(tags[i]));
        }
        // Msg::Info("Physical Surface: %s , faces %d",nameSurface.c_str(),i,faces.size());
          // std::string fnameSurface=path+baseName+"_"+nameSurface+extName;
          // FILE *fp_surface = Fopen(fnameSurface.c_str(), binary ? "wb" : "w");
          // if(!fp) {
          //   Msg::Error("Unable to open file '%s'", fnameSurface.c_str());
          //   return 0;
          // }
          writeSTLfaces(fp_surface, faces, binary, scalingFactor, nameSurface);
          Msg::Info("Surface(%s) has %d face(s), done writing %s",nameSurface.c_str(),faces.size(),fnameSurface.c_str());
          // fclose(fp_surface);
          // surface name
          surfaceName.push_back(nameSurface);
      }
      fclose(fp_surface);

      // write snappyHexMesh sub-dicts 
      std::string fnameSnappyHexMeshDict=path+baseName+".snappyHexMeshDict";
      FILE *fp_SnappyHexMeshDict = Fopen(fnameSnappyHexMeshDict.c_str(), "w");
      if(!fp_SnappyHexMeshDict) {
        Msg::Error("Unable to open file '%s'", fnameSnappyHexMeshDict.c_str());
        return 0;
      }
      // 1. geometry
      fprintf(fp_SnappyHexMeshDict,"geometry\n{\n");
      fprintf(fp_SnappyHexMeshDict,"\t%s_Surface.stl\n\t{\n\t\ttype triSurfaceMesh;\n\t\tregions\n\t\t{\n",baseName.c_str());
      for(int i=0;i<surfaceName.size();i++)
      {
        fprintf(fp_SnappyHexMeshDict,"\t\t\t%s {name %s;}\n",surfaceName[i].c_str(), surfaceName[i].c_str());
      }
      fprintf(fp_SnappyHexMeshDict,"\t\t}\n\t}\n");
      for(int i=0;i<volumeName.size();i++)
      {
        fprintf(fp_SnappyHexMeshDict,"\t%s_%s.stl {type triSurfaceMesh; name %s;}\n", 
                baseName.c_str(), volumeName[i].c_str(), volumeName[i].c_str());
      }
      fprintf(fp_SnappyHexMeshDict,"};\n\n");
      // 2. features
      fprintf(fp_SnappyHexMeshDict,"features\n(\n\t{file \"%s_Surface.eMesh\"; level 1;}\n", baseName.c_str());
      for(int i=0;i<volumeName.size();i++)
      {
        fprintf(fp_SnappyHexMeshDict,"\t{file \"%s_%s.eMesh\"; level 1;}\n",baseName.c_str(),volumeName[i].c_str());
      }
      fprintf(fp_SnappyHexMeshDict,");\n");
      // 3. refinementSurfaces
      fprintf(fp_SnappyHexMeshDict,"refinementSurfaces\n{\n\t%s_Surface.stl\n\t{\n\t\tlevel (0,0);\n\t\tregions\n\t\t{\n",baseName.c_str());
      for(int i=0;i<surfaceName.size();i++)
      {
        fprintf(fp_SnappyHexMeshDict,"\t\t\t%s {level (1 1); patchInfo { type patch; }}\n",surfaceName[i].c_str());
      }
      fprintf(fp_SnappyHexMeshDict,"\t\t}\n\t}\n");
      for(int i=0;i<volumeName.size();i++)
      {
        fprintf(fp_SnappyHexMeshDict,"\t%s\n\t{\n\t\tlevel (1 1);\n",volumeName[i].c_str());
        fprintf(fp_SnappyHexMeshDict,"\t\tfaceZone %s;\n\t\tcellZone %s;\n\t\tcellZoneInside inside;\n\t\tboundary    internal;\n\t}\n",
                volumeName[i].c_str(), volumeName[i].c_str());
      }
      // 4. refinementRegions
      fprintf(fp_SnappyHexMeshDict,"refinementRegions\n{\n");
      for(int i=0;i<volumeName.size();i++)
      {
        fprintf(fp_SnappyHexMeshDict,"\t%s {levels ((0 3)); mode inside;}\n",volumeName[i].c_str());
      }
      fprintf(fp_SnappyHexMeshDict,"}\n");
      // close file
      fclose(fp_SnappyHexMeshDict);
      Msg::Info("Snippet of snappyHexMeshDict has been saved: %s",fnameSnappyHexMeshDict.c_str());
  }
  
  return 1;
}
