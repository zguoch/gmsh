// Gmsh - Copyright (C) 1997-2020 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <map>
#include <vector>
#include <list>
#include "GmshConfig.h"
#include "STensor3.h"
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <algorithm>

#if defined(HAVE_POST)
class PView;
#endif

class Field;
class GEntity;

typedef enum {
  FIELD_OPTION_DOUBLE = 0,
  FIELD_OPTION_INT,
  FIELD_OPTION_STRING,
  FIELD_OPTION_PATH,
  FIELD_OPTION_BOOL,
  FIELD_OPTION_LIST,
  FIELD_OPTION_LIST_DOUBLE
} FieldOptionType;

class FieldCallback {
private:
  std::string _help;

public:
  virtual void run() = 0;
  FieldCallback(const std::string &help) { _help = help; }
  virtual ~FieldCallback(){};
  virtual std::string getDescription() { return _help; }
};

class FieldOption {
private:
  std::string _help;

protected:
  bool *status;
  inline void modified()
  {
    if(status) *status = true;
  }

public:
  FieldOption(const std::string &help, bool *_status)
    : _help(help), status(_status)
  {
  }
  virtual ~FieldOption() {}
  virtual FieldOptionType getType() = 0;
  virtual void getTextRepresentation(std::string &v_str) = 0;
  virtual std::string getDescription() { return _help; }
  std::string getTypeName()
  {
    switch(getType()) {
    case FIELD_OPTION_INT: return "integer";
    case FIELD_OPTION_DOUBLE: return "float";
    case FIELD_OPTION_BOOL: return "boolean";
    case FIELD_OPTION_PATH: return "path";
    case FIELD_OPTION_STRING: return "string";
    case FIELD_OPTION_LIST: return "list";
    case FIELD_OPTION_LIST_DOUBLE: return "list_double";
    default: return "unknown";
    }
  }
  virtual void numericalValue(double val) {}
  virtual double numericalValue() const { return 0.; }
  virtual const std::list<int> &list() const
  {
    static std::list<int> l;
    return l;
  }
  virtual const std::list<double> &listdouble() const
  {
    static std::list<double> l;
    return l;
  }
  virtual void list(std::list<int> value) {}
  virtual void listdouble(std::list<double> value) {}
  virtual std::string string() const { return ""; }
  virtual void string(const std::string value) {}
};

class Field {
public:
  Field() : update_needed(false) {}
  virtual ~Field();
  virtual void update() {}
  int id;
  std::map<std::string, FieldOption *> options;
  std::map<std::string, FieldCallback *> callbacks;
  virtual bool isotropic() const { return true; }
  // isotropic
  virtual double operator()(double x, double y, double z, GEntity *ge = 0) = 0;
  // anisotropic
  virtual void operator()(double x, double y, double z, SMetric3 &,
                          GEntity *ge = 0)
  {
  }
  bool update_needed;
  virtual const char *getName() = 0;
#if defined(HAVE_POST)
  void putOnView(PView *view, int comp = -1);
#endif
  void putOnNewView(int viewTag = -1);
  virtual std::string getDescription() { return ""; }
  FieldOption *getOption(const std::string &optionName);
};

class FieldFactory {
public:
  virtual ~FieldFactory() {}
  virtual Field *operator()() = 0;
};

class FieldManager : public std::map<int, Field *> {
private:
  int _background_field;
  std::vector<int> _boundaryLayer_fields;

public:
  std::map<std::string, FieldFactory *> map_type_name;
  void initialize();
  void reset();
  Field *get(int id);
  Field *newField(int id, const std::string &type_name);
  void deleteField(int id);
  int newId();
  int maxId();
  FieldManager();
  ~FieldManager();
  // compatibility with -bgm
  void setBackgroundMesh(int iView);
  // set and get background field
  void setBackgroundField(Field *BGF);
  inline void setBackgroundFieldId(int id) { _background_field = id; };
  inline void addBoundaryLayerFieldId(int id)
  {
    for(std::size_t i = 0; i < _boundaryLayer_fields.size(); ++i) {
      if(_boundaryLayer_fields[i] == id) return;
    }
    _boundaryLayer_fields.push_back(id);
  }
  inline void addBoundaryLayerFieldId(std::vector<int> &tags)
  {
    for(std::size_t i = 0; i < tags.size(); ++i)
      addBoundaryLayerFieldId(tags[i]);
  }
  inline int getBackgroundField() { return _background_field; }
  inline int getNumBoundaryLayerFields()
  {
    return (int)_boundaryLayer_fields.size();
  }
  inline int getBoundaryLayerField(int i) { return _boundaryLayer_fields[i]; }
};

// Boundary Layer Field (used both for anisotropic meshing and BL
// extrusion)

class DistanceField;

class BoundaryLayerField : public Field {
private:
  std::list<DistanceField *> _att_fields;
  std::list<double> hwall_n_nodes;
  std::list<int> nodes_id, edges_id;
  std::list<int> edges_id_saved, nodes_id_saved, fan_nodes_id;
  std::list<int> excluded_faces_id;
  void operator()(DistanceField *cc, double dist, double x, double y, double z,
                  SMetric3 &metr, GEntity *ge);

public:
  double hwall_n, ratio, hfar, thickness;
  double current_distance, tgt_aniso_ratio;
  SPoint3 _closest_point;
  int iRecombine, iIntersect;
  DistanceField *current_closest;
  virtual bool isotropic() const { return false; }
  virtual const char *getName();
  virtual std::string getDescription();
  BoundaryLayerField();
  ~BoundaryLayerField() { removeAttractors(); }
  virtual double operator()(double x, double y, double z, GEntity *ge = 0);
  virtual void operator()(double x, double y, double z, SMetric3 &metr,
                          GEntity *ge = 0);
  bool isEdgeBL(int iE) const
  {
    return std::find(edges_id.begin(), edges_id.end(), iE) != edges_id.end();
  }
  bool isEdgeBLSaved(int iE) const
  {
    return std::find(edges_id_saved.begin(), edges_id_saved.end(), iE) !=
           edges_id_saved.end();
  }
  bool isFanNode(int iV) const
  {
    return std::find(fan_nodes_id.begin(), fan_nodes_id.end(), iV) !=
           fan_nodes_id.end();
  }
  bool isEndNode(int iV) const
  {
    return std::find(nodes_id.begin(), nodes_id.end(), iV) != nodes_id.end();
  }
  double hwall(int iV)
  {
    for(std::list<double>::iterator it = hwall_n_nodes.begin();
        it != hwall_n_nodes.end(); ++it) {
      int i = (int)*it;
      ++it;
      double h = *it;
      if(i == iV) return h;
    }
    return hwall_n;
  }

  void computeFor1dMesh(double x, double y, double z, SMetric3 &metr);
  void setupFor1d(int iE);
  void setupFor2d(int iF);
  void removeAttractors();
};

class FieldOptionString : public FieldOption {
public:
  std::string &val;
  virtual FieldOptionType getType() { return FIELD_OPTION_STRING; }
  FieldOptionString(std::string &_val, const std::string &_help,
                    bool *_status = 0)
    : FieldOption(_help, _status), val(_val)
  {
  }
  void string(const std::string value)
  {
    modified();
    val = value;
  }
  std::string string() const { return val; }
  void getTextRepresentation(std::string &v_str)
  {
    std::ostringstream sstream;
    sstream << "\"" << val << "\"";
    v_str = sstream.str();
  }
};

class FieldOptionDouble : public FieldOption {
public:
  double &val;
  FieldOptionType getType() { return FIELD_OPTION_DOUBLE; }
  FieldOptionDouble(double &_val, const std::string &_help, bool *_status = 0)
    : FieldOption(_help, _status), val(_val)
  {
  }
  double numericalValue() const { return val; }
  void numericalValue(double v)
  {
    modified();
    val = v;
  }
  void getTextRepresentation(std::string &v_str)
  {
    std::ostringstream sstream;
    sstream.precision(16);
    sstream << val;
    v_str = sstream.str();
  }
};

class FieldOptionInt : public FieldOption {
public:
  int &val;
  FieldOptionType getType() { return FIELD_OPTION_INT; }
  FieldOptionInt(int &_val, const std::string &_help, bool *_status = 0)
    : FieldOption(_help, _status), val(_val)
  {
  }
  double numericalValue() const { return val; }
  void numericalValue(double v)
  {
    modified();
    val = (int)v;
  }
  void getTextRepresentation(std::string &v_str)
  {
    std::ostringstream sstream;
    sstream << val;
    v_str = sstream.str();
  }
};

class FieldOptionList : public FieldOption {
public:
  std::list<int> &val;
  FieldOptionType getType() { return FIELD_OPTION_LIST; }
  FieldOptionList(std::list<int> &_val, const std::string &_help,
                  bool *_status = 0)
    : FieldOption(_help, _status), val(_val)
  {
  }
  void list(std::list<int> value)
  {
    modified();
    val = value;
  }
  const std::list<int> &list() const { return val; }
  void getTextRepresentation(std::string &v_str)
  {
    std::ostringstream sstream;
    sstream << "{";
    for(std::list<int>::iterator it = val.begin(); it != val.end(); it++) {
      if(it != val.begin()) sstream << ", ";
      sstream << *it;
    }
    sstream << "}";
    v_str = sstream.str();
  }
};

class FieldOptionListDouble : public FieldOption {
public:
  std::list<double> &val;
  FieldOptionType getType() { return FIELD_OPTION_LIST_DOUBLE; }
  FieldOptionListDouble(std::list<double> &_val, const std::string &_help,
                        bool *_status = 0)
    : FieldOption(_help, _status), val(_val)
  {
  }
  void listdouble(std::list<double> value)
  {
    modified();
    val = value;
  }
  const std::list<double> &listdouble() const { return val; }
  void getTextRepresentation(std::string &v_str)
  {
    std::ostringstream sstream;
    sstream.precision(16);
    sstream << "{";
    for(std::list<double>::iterator it = val.begin(); it != val.end(); it++) {
      if(it != val.begin()) sstream << ", ";
      sstream << *it;
    }
    sstream << "}";
    v_str = sstream.str();
  }
};

class FieldOptionPath : public FieldOptionString {
public:
  virtual FieldOptionType getType() { return FIELD_OPTION_PATH; }
  FieldOptionPath(std::string &_val, const std::string &_help,
                  bool *_status = 0)
    : FieldOptionString(_val, _help, _status)
  {
  }
};

class FieldOptionBool : public FieldOption {
public:
  bool &val;
  FieldOptionType getType() { return FIELD_OPTION_BOOL; }
  FieldOptionBool(bool &_val, const std::string &_help, bool *_status = 0)
    : FieldOption(_help, _status), val(_val)
  {
  }
  double numericalValue() const { return val; }
  void numericalValue(double v)
  {
    modified();
    val = v;
  }
  void getTextRepresentation(std::string &v_str)
  {
    std::ostringstream sstream;
    sstream << val;
    v_str = sstream.str();
  }
};

template <class t> class FieldCallbackGeneric : public FieldCallback {
  t *_field;
  void (t::*_callback)();

public:
  void run() { (_field->*_callback)(); }
  FieldCallbackGeneric(t *field, void (t::*callback)(),
                       const std::string description)
    : FieldCallback(description)
  {
    _field = field;
    _callback = callback;
  }
};

template <class F> class FieldFactoryT : public FieldFactory {
public:
  Field *operator()() { return new F; }
};

// the class GenericField contains a set of void* functions, which give a mesh
// size All these functions are called when calling operator() ; then, the
// minimum size is returned.
class GenericField : public Field {
public:
  // callback prototypes:
  // this callback is called with a void* previously given to the GenericField !
  typedef bool (*ptrfunction)(double, double, double, void *, double &);
  // this callback also takes the GEntity object into account
  typedef bool (*ptrfunctionextended)(double, double, double, void*, void*, double&);

  GenericField();
  ~GenericField();
  using Field::operator();
  virtual double operator()(double x, double y, double z, GEntity *ge = 0);
  virtual const char *getName() { return "GenericField"; };

  // sets the callbacks
  void setCallbackWithData(ptrfunction fct, void *data);
  void setCallbackWithData(ptrfunctionextended fct, void *data);

private:
  std::vector<std::pair<ptrfunction, void*> > cbs_with_data; // the callbacks with the data to be sent to them
  std::vector<std::pair<ptrfunctionextended, void*> > cbs_extended_with_data; // the extended callbacks with the data to be sent to them
};

#endif
