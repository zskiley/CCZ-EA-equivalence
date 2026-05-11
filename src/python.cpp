#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstdint>
#include <exception>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "algorithms.h"
#include "ccz_equivalence.h"
#include "dfs_auto.h"
#include "ea_equivalence.h"
#include "gf2_linear.h"
#include "graph.h"
#include "groups/graph_point_perm.h"
#include "groups/group_ops.h"
#include "partial_map.h"

namespace {

bool PyToLongLong(PyObject* obj, long long* out) {
  if (out == nullptr) return false;
  PyObject* idx = PyNumber_Index(obj);
  if (idx == nullptr) return false;
  const long long value = PyLong_AsLongLong(idx);
  Py_DECREF(idx);
  if (PyErr_Occurred()) return false;
  *out = value;
  return true;
}

bool PyToUInt32(PyObject* obj, uint32_t* out) {
  if (out == nullptr) return false;
  PyObject* idx = PyNumber_Index(obj);
  if (idx == nullptr) return false;
  const unsigned long long value = PyLong_AsUnsignedLongLong(idx);
  Py_DECREF(idx);
  if (PyErr_Occurred()) return false;
  if (value > std::numeric_limits<uint32_t>::max()) {
    PyErr_SetString(PyExc_OverflowError, "value does not fit in uint32");
    return false;
  }
  *out = static_cast<uint32_t>(value);
  return true;
}

bool ParseBits(PyObject* obj, const char* name, int* out) {
  if (out == nullptr) return false;
  long long v = 0;
  if (!PyToLongLong(obj, &v)) {
    PyErr_Format(PyExc_TypeError, "%s must be an integer", name);
    return false;
  }
  if (v <= 0 || v >= 31) {
    PyErr_Format(PyExc_ValueError, "%s must be in [1, 30]", name);
    return false;
  }
  *out = static_cast<int>(v);
  return true;
}

bool ParseOptionalBits(PyObject* obj, const char* name, int default_value,
                       int* out) {
  if (out == nullptr) return false;
  if (obj == Py_None) {
    *out = default_value;
    return true;
  }
  return ParseBits(obj, name, out);
}

bool ParseOptionalTimeLimit(PyObject* obj, double default_value,
                            double* out) {
  if (out == nullptr) return false;
  if (obj == Py_None) {
    *out = default_value;
    return true;
  }
  const double value = PyFloat_AsDouble(obj);
  if (PyErr_Occurred()) {
    PyErr_SetString(PyExc_TypeError,
                    "time_limit_seconds must be a float or None");
    return false;
  }
  *out = value;
  return true;
}

bool ParseOptionalSizeT(PyObject* obj, const char* name,
                        std::size_t default_value, std::size_t* out) {
  if (out == nullptr) return false;
  if (obj == Py_None) {
    *out = default_value;
    return true;
  }
  PyObject* idx = PyNumber_Index(obj);
  if (idx == nullptr) {
    PyErr_Format(PyExc_TypeError, "%s must be an integer", name);
    return false;
  }
  const unsigned long long value = PyLong_AsUnsignedLongLong(idx);
  Py_DECREF(idx);
  if (PyErr_Occurred()) {
    PyErr_Format(PyExc_TypeError, "%s must be a non-negative integer", name);
    return false;
  }
  if (value > std::numeric_limits<std::size_t>::max()) {
    PyErr_Format(PyExc_OverflowError, "%s is too large", name);
    return false;
  }
  *out = static_cast<std::size_t>(value);
  return true;
}

bool ParseOptionalTimeLimitAsOptional(PyObject* obj, std::optional<double>* out) {
  if (out == nullptr) return false;
  if (obj == Py_None) {
    *out = std::nullopt;
    return true;
  }
  const double value = PyFloat_AsDouble(obj);
  if (PyErr_Occurred()) {
    PyErr_SetString(PyExc_TypeError,
                    "time_limit_seconds must be a float or None");
    return false;
  }
  *out = value;
  return true;
}

bool ParseTruthTable(PyObject* obj, std::vector<uint32_t>* out) {
  if (out == nullptr) return false;
  PyObject* seq = PySequence_Fast(obj, "truth_table must be a sequence");
  if (seq == nullptr) return false;

  const Py_ssize_t n = PySequence_Fast_GET_SIZE(seq);
  out->clear();
  out->reserve(static_cast<std::size_t>(n));
  PyObject** items = PySequence_Fast_ITEMS(seq);
  for (Py_ssize_t i = 0; i < n; ++i) {
    uint32_t value = 0u;
    if (!PyToUInt32(items[i], &value)) {
      Py_DECREF(seq);
      PyErr_Format(PyExc_TypeError,
                   "truth_table[%zd] is not a valid non-negative integer", i);
      return false;
    }
    out->push_back(value);
  }

  Py_DECREF(seq);
  return true;
}

bool BuildGraphFromTruthTableArgs(PyObject* truth_table_obj, int n_bits,
                                  int m_bits, GraphData* out) {
  if (out == nullptr) return false;
  if (n_bits + m_bits >= 31) {
    PyErr_SetString(PyExc_ValueError, "n_bits + m_bits must be <= 30");
    return false;
  }

  std::vector<uint32_t> truth_table;
  if (!ParseTruthTable(truth_table_obj, &truth_table)) return false;

  const std::size_t expected_size =
      static_cast<std::size_t>(1u) << static_cast<std::size_t>(n_bits);
  if (truth_table.size() != expected_size) {
    PyErr_Format(PyExc_ValueError,
                 "truth_table must have exactly 2^n_bits elements (expected %zu, got %zu)",
                 expected_size, truth_table.size());
    return false;
  }

  *out = BuildGraph(truth_table, n_bits, m_bits);
  return true;
}

template <typename PointMap>
PyObject* PointMapToPyDict(const PointMap& map) {
  PyObject* dict = PyDict_New();
  if (dict == nullptr) return nullptr;

  for (const auto& kv : map) {
    PyObject* key = PyLong_FromUnsignedLong(static_cast<unsigned long>(kv.first));
    if (key == nullptr) {
      Py_DECREF(dict);
      return nullptr;
    }
    PyObject* value =
        PyLong_FromUnsignedLong(static_cast<unsigned long>(kv.second));
    if (value == nullptr) {
      Py_DECREF(key);
      Py_DECREF(dict);
      return nullptr;
    }
    if (PyDict_SetItem(dict, key, value) != 0) {
      Py_DECREF(key);
      Py_DECREF(value);
      Py_DECREF(dict);
      return nullptr;
    }
    Py_DECREF(key);
    Py_DECREF(value);
  }

  return dict;
}

PyObject* AutoListToPyList(const std::vector<GraphPointMap>& autos) {
  PyObject* out = PyList_New(static_cast<Py_ssize_t>(autos.size()));
  if (out == nullptr) return nullptr;

  for (std::size_t i = 0; i < autos.size(); ++i) {
    PyObject* map_dict = PointMapToPyDict(autos[i]);
    if (map_dict == nullptr) {
      Py_DECREF(out);
      return nullptr;
    }
    PyList_SET_ITEM(out, static_cast<Py_ssize_t>(i), map_dict);
  }
  return out;
}

PyObject* GeneratorToGraphPointDict(const GraphData& graph,
                                    const groups::Permutation& generator) {
  const std::size_t n = graph.points.size();
  if (generator.Degree() != n) {
    PyErr_SetString(PyExc_RuntimeError,
                    "generator degree does not match graph size");
    return nullptr;
  }

  PyObject* dict = PyDict_New();
  if (dict == nullptr) return nullptr;

  for (std::size_t i = 0; i < n; ++i) {
    const uint32_t image_index = generator.Apply(static_cast<uint32_t>(i));
    if (image_index >= n) {
      Py_DECREF(dict);
      PyErr_SetString(PyExc_RuntimeError,
                      "generator image index is out of bounds");
      return nullptr;
    }
    const uint32_t x = graph.points[i];
    const uint32_t y = graph.points[image_index];

    PyObject* key = PyLong_FromUnsignedLong(static_cast<unsigned long>(x));
    if (key == nullptr) {
      Py_DECREF(dict);
      return nullptr;
    }
    PyObject* value = PyLong_FromUnsignedLong(static_cast<unsigned long>(y));
    if (value == nullptr) {
      Py_DECREF(key);
      Py_DECREF(dict);
      return nullptr;
    }
    if (PyDict_SetItem(dict, key, value) != 0) {
      Py_DECREF(key);
      Py_DECREF(value);
      Py_DECREF(dict);
      return nullptr;
    }
    Py_DECREF(key);
    Py_DECREF(value);
  }

  return dict;
}

PyObject* GeneratorListToPyList(const GraphData& graph,
                                const std::vector<groups::Permutation>& generators) {
  PyObject* out = PyList_New(0);
  if (out == nullptr) return nullptr;

  for (const auto& generator : generators) {
    if (generator.IsIdentity()) continue;
    PyObject* gen_dict = GeneratorToGraphPointDict(graph, generator);
    if (gen_dict == nullptr) {
      Py_DECREF(out);
      return nullptr;
    }
    if (PyList_Append(out, gen_dict) != 0) {
      Py_DECREF(gen_dict);
      Py_DECREF(out);
      return nullptr;
    }
    Py_DECREF(gen_dict);
  }

  return out;
}

PyObject* AmbientGeneratorToPyDict(const AffineMapData& generator) {
  PyObject* dict = PyDict_New();
  if (dict == nullptr) return nullptr;

  PyObject* translation_obj = PyLong_FromUnsignedLong(
      static_cast<unsigned long>(generator.translation));
  if (translation_obj == nullptr) {
    Py_DECREF(dict);
    return nullptr;
  }
  if (PyDict_SetItemString(dict, "translation", translation_obj) != 0) {
    Py_DECREF(translation_obj);
    Py_DECREF(dict);
    return nullptr;
  }
  Py_DECREF(translation_obj);

  PyObject* cols_obj =
      PyList_New(static_cast<Py_ssize_t>(generator.linear_cols.size()));
  if (cols_obj == nullptr) {
    Py_DECREF(dict);
    return nullptr;
  }
  for (std::size_t i = 0; i < generator.linear_cols.size(); ++i) {
    PyObject* col_obj = PyLong_FromUnsignedLong(
        static_cast<unsigned long>(generator.linear_cols[i]));
    if (col_obj == nullptr) {
      Py_DECREF(cols_obj);
      Py_DECREF(dict);
      return nullptr;
    }
    PyList_SET_ITEM(cols_obj, static_cast<Py_ssize_t>(i), col_obj);
  }
  if (PyDict_SetItemString(dict, "linear_cols", cols_obj) != 0) {
    Py_DECREF(cols_obj);
    Py_DECREF(dict);
    return nullptr;
  }
  Py_DECREF(cols_obj);

  return dict;
}

PyObject* AmbientGeneratorListToPyList(
    const std::vector<AffineMapData>& generators) {
  PyObject* out = PyList_New(0);
  if (out == nullptr) return nullptr;

  for (const AffineMapData& generator : generators) {
    if (generator.IsIdentity()) continue;
    PyObject* gen_dict = AmbientGeneratorToPyDict(generator);
    if (gen_dict == nullptr) {
      Py_DECREF(out);
      return nullptr;
    }
    if (PyList_Append(out, gen_dict) != 0) {
      Py_DECREF(gen_dict);
      Py_DECREF(out);
      return nullptr;
    }
    Py_DECREF(gen_dict);
  }

  return out;
}

PyObject* BuildAutoGroupResult(const GraphData& graph) {
  PyObject* out = PyDict_New();
  if (out == nullptr) return nullptr;
  const std::vector<groups::Permutation>& graph_generators =
      GetAutoGroupGenerators();
  const std::vector<AffineMapData>& affine_generators =
      GetAffineAutoGenerators();

  PyObject* order_obj = PyLong_FromUnsignedLongLong(
      static_cast<unsigned long long>(GetTotalAutoGroup()));
  if (order_obj == nullptr) {
    Py_DECREF(out);
    return nullptr;
  }
  if (PyDict_SetItemString(out, "order", order_obj) != 0) {
    Py_DECREF(order_obj);
    Py_DECREF(out);
    return nullptr;
  }
  Py_DECREF(order_obj);

  PyObject* complete_obj = FoundEntireAutoGroup() ? Py_True : Py_False;
  Py_INCREF(complete_obj);
  if (PyDict_SetItemString(out, "found_entire_group", complete_obj) != 0) {
    Py_DECREF(complete_obj);
    Py_DECREF(out);
    return nullptr;
  }
  Py_DECREF(complete_obj);

  PyObject* generators_obj = AmbientGeneratorListToPyList(affine_generators);
  if (generators_obj == nullptr) {
    Py_DECREF(out);
    return nullptr;
  }
  if (PyDict_SetItemString(out, "generators", generators_obj) != 0) {
    Py_DECREF(generators_obj);
    Py_DECREF(out);
    return nullptr;
  }
  Py_DECREF(generators_obj);

  PyObject* graph_generators_obj =
      GeneratorListToPyList(graph, graph_generators);
  if (graph_generators_obj == nullptr) {
    Py_DECREF(out);
    return nullptr;
  }
  if (PyDict_SetItemString(out, "graph_generators", graph_generators_obj) != 0) {
    Py_DECREF(graph_generators_obj);
    Py_DECREF(out);
    return nullptr;
  }
  Py_DECREF(graph_generators_obj);

  return out;
}

bool GraphGeneratorDictToPermutation(PyObject* item, const GraphData& graph,
                                     groups::Permutation* out) {
  if (out == nullptr) return false;
  if (!PyDict_Check(item)) return false;

  groups::GraphPointMap map;
  PyObject *key = nullptr, *value = nullptr;
  Py_ssize_t pos = 0;
  while (PyDict_Next(item, &pos, &key, &value)) {
    uint32_t x = 0u;
    uint32_t y = 0u;
    if (!PyToUInt32(key, &x) || !PyToUInt32(value, &y)) {
      return false;
    }
    map.emplace_back(x, y);
  }

  const groups::GraphPointIndex index(graph);
  return groups::GraphPointMapToPermutation(map, index, out);
}

bool ParseAmbientGeneratorDict(PyObject* item, const GraphData& graph,
                               AffineMapData* out) {
  if (out == nullptr) return false;
  if (!PyDict_Check(item)) return false;

  PyObject* translation_obj = PyDict_GetItemString(item, "translation");
  PyObject* cols_obj = PyDict_GetItemString(item, "linear_cols");
  if (translation_obj == nullptr || cols_obj == nullptr) return false;

  uint32_t translation = 0u;
  if (!PyToUInt32(translation_obj, &translation)) {
    PyErr_SetString(PyExc_TypeError,
                    "ambient generator translation must be an integer");
    return false;
  }

  PyObject* cols_seq = PySequence_Fast(
      cols_obj, "ambient generator linear_cols must be a sequence");
  if (cols_seq == nullptr) return false;

  const Py_ssize_t n_cols = PySequence_Fast_GET_SIZE(cols_seq);
  if (n_cols != graph.d_bits) {
    Py_DECREF(cols_seq);
    PyErr_SetString(
        PyExc_ValueError,
        "ambient generator linear_cols length must equal n_bits + m_bits");
    return false;
  }

  if (graph.d_bits <= 0 || graph.d_bits >= 31) {
    Py_DECREF(cols_seq);
    PyErr_SetString(PyExc_ValueError, "graph dimension must be in [1, 30]");
    return false;
  }
  const uint32_t ambient_size = (1u << graph.d_bits);
  if (translation >= ambient_size) {
    Py_DECREF(cols_seq);
    PyErr_SetString(PyExc_ValueError,
                    "ambient generator translation is out of range");
    return false;
  }

  AffineMapData map;
  map.dimension_bits = graph.d_bits;
  map.translation = translation;
  map.linear_cols.reserve(static_cast<std::size_t>(n_cols));

  PyObject** items = PySequence_Fast_ITEMS(cols_seq);
  for (Py_ssize_t i = 0; i < n_cols; ++i) {
    uint32_t col = 0u;
    if (!PyToUInt32(items[i], &col)) {
      Py_DECREF(cols_seq);
      PyErr_SetString(PyExc_TypeError,
                      "ambient generator linear_cols entries must be integers");
      return false;
    }
    if (col >= ambient_size) {
      Py_DECREF(cols_seq);
      PyErr_SetString(PyExc_ValueError,
                      "ambient generator linear_cols entry is out of range");
      return false;
    }
    map.linear_cols.push_back(col);
  }

  Py_DECREF(cols_seq);
  *out = std::move(map);
  return true;
}

bool AmbientGeneratorToPermutation(const GraphData& graph,
                                   const AffineMapData& map,
                                   groups::Permutation* out) {
  if (out == nullptr) return false;
  if (map.dimension_bits != graph.d_bits) return false;
  if (map.linear_cols.size() != static_cast<std::size_t>(graph.d_bits)) {
    return false;
  }

  std::vector<uint32_t> inverse_cols;
  if (!gf2_linear::InvertLinearColumns(map.linear_cols, graph.d_bits,
                                       &inverse_cols)) {
    return false;
  }

  const std::size_t degree = graph.points.size();
  if (degree > static_cast<std::size_t>(
                   std::numeric_limits<uint16_t>::max()) +
                   1u) {
    return false;
  }

  const groups::GraphPointIndex index(graph);
  std::vector<uint16_t> images(degree, 0u);
  for (std::size_t i = 0; i < degree; ++i) {
    const uint32_t image_point = map.Apply(graph.points[i]);
    const auto image_idx = index.IndexOf(image_point);
    if (!image_idx.has_value() || *image_idx >= degree) return false;
    images[i] = static_cast<uint16_t>(*image_idx);
  }

  try {
    *out = groups::Permutation(std::move(images));
  } catch (...) {
    return false;
  }
  return true;
}

bool IsValidGeneratorAuto(const GraphData& graph, const groups::Permutation& g,
                          bool require_ea) {
  if (g.Degree() != graph.points.size()) return false;
  PartialAffineMap A(graph.d_bits);
  for (std::size_t i = 0; i < graph.points.size(); ++i) {
    const uint32_t x = graph.points[i];
    const uint32_t y = graph.points[g.Apply(static_cast<uint32_t>(i))];
    if (!A.Update(x, y)) return false;
  }
  return require_ea ? A.valid_ea(graph) : A.valid_ccz(graph);
}

bool EquivalencePointMapToAffineMap(const EquivalencePointMap& point_map,
                                    const GraphData& graph,
                                    bool require_ea,
                                    AffineMapData* out) {
  if (out == nullptr) return false;
  PartialAffineMap A(graph.d_bits);
  for (const auto& kv : point_map) {
    if (!A.Update(kv.first, kv.second)) return false;
  }
  if (require_ea) {
    return A.ExtractRepresentativeAffineMap(out, graph.n_bits, graph.m_bits);
  }
  return A.ExtractRepresentativeAffineMap(out);
}

bool ParseProvidedAutoGroupGenerators(PyObject* auto_group_obj,
                                      const GraphData& graph,
                                      bool require_ea,
                                      std::vector<groups::Permutation>* out) {
  if (out == nullptr) return false;
  out->clear();

  std::vector<groups::Permutation> generators;
  auto parse_generator_sequence = [&](PyObject* generators_obj,
                                      const char* error_context) -> bool {
    PyObject* seq = PySequence_Fast(generators_obj, error_context);
    if (seq == nullptr) return false;

    const Py_ssize_t n = PySequence_Fast_GET_SIZE(seq);
    PyObject** items = PySequence_Fast_ITEMS(seq);
    generators.reserve(generators.size() + static_cast<std::size_t>(n));

    for (Py_ssize_t i = 0; i < n; ++i) {
      groups::Permutation perm;
      const bool has_ambient_marker =
          PyDict_Check(items[i]) &&
          (PyDict_GetItemString(items[i], "translation") != nullptr ||
           PyDict_GetItemString(items[i], "linear_cols") != nullptr);
      if (has_ambient_marker) {
        AffineMapData ambient_map;
        if (!ParseAmbientGeneratorDict(items[i], graph, &ambient_map)) {
          Py_DECREF(seq);
          if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError,
                            "auto_group contains invalid ambient generator");
          }
          return false;
        }
        if (!AmbientGeneratorToPermutation(graph, ambient_map, &perm)) {
          Py_DECREF(seq);
          PyErr_SetString(
              PyExc_ValueError,
              "auto_group contains ambient generator that does not act as a "
              "graph automorphism");
          return false;
        }
      } else {
        if (!GraphGeneratorDictToPermutation(items[i], graph, &perm)) {
          Py_DECREF(seq);
          PyErr_SetString(PyExc_ValueError,
                          "auto_group contains invalid graph generator");
          return false;
        }
      }
      if (!IsValidGeneratorAuto(graph, perm, require_ea)) {
        Py_DECREF(seq);
        PyErr_SetString(
            PyExc_ValueError,
            "auto_group contains generator that is not a valid automorphism");
        return false;
      }
      if (!perm.IsIdentity()) generators.push_back(std::move(perm));
    }

    Py_DECREF(seq);
    return true;
  };

  if (PyDict_Check(auto_group_obj)) {
    PyObject* generators_obj = PyDict_GetItemString(auto_group_obj, "generators");
    PyObject* graph_generators_obj =
        PyDict_GetItemString(auto_group_obj, "graph_generators");
    if (generators_obj == nullptr && graph_generators_obj == nullptr) {
      PyErr_SetString(
          PyExc_ValueError,
          "auto_group dict must contain key 'generators' or 'graph_generators'");
      return false;
    }
    if (generators_obj != nullptr &&
        !parse_generator_sequence(
            generators_obj,
            "auto_group['generators'] must be a sequence of generators")) {
      return false;
    }
    if (graph_generators_obj != nullptr &&
        !parse_generator_sequence(
            graph_generators_obj,
            "auto_group['graph_generators'] must be a sequence of generators")) {
      return false;
    }
  } else {
    if (!parse_generator_sequence(
            auto_group_obj,
            "auto_group must be a generator list or a dict with "
            "'generators'/'graph_generators'")) {
      return false;
    }
  }

  *out = groups::DeduplicateGenerators(std::move(generators));
  return true;
}

PyObject* PyCCZAuto(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* truth_table_obj = nullptr;
  PyObject* m_bits_obj = Py_None;
  PyObject* time_limit_obj = Py_None;
  PyObject* min_active_hyperplanes_obj = Py_None;
  int n_bits = 0;
  static const char* kKwlist[] = {"truth_table", "n_bits", "m_bits",
                                   "time_limit_seconds",
                                   "min_active_hyperplanes", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi|OOO",
                                   const_cast<char**>(kKwlist),
                                   &truth_table_obj, &n_bits, &m_bits_obj,
                                   &time_limit_obj,
                                   &min_active_hyperplanes_obj)) {
    return nullptr;
  }

  int m_bits = 0;
  if (!ParseOptionalBits(m_bits_obj, "m_bits", n_bits, &m_bits)) return nullptr;
  double time_limit_seconds = 0.0;
  if (!ParseOptionalTimeLimit(time_limit_obj, 0.0, &time_limit_seconds)) {
    return nullptr;
  }
  std::size_t min_active_hyperplanes = 0;
  if (!ParseOptionalSizeT(min_active_hyperplanes_obj, "min_active_hyperplanes",
                          0, &min_active_hyperplanes)) {
    return nullptr;
  }

  GraphData F;
  if (!BuildGraphFromTruthTableArgs(truth_table_obj, n_bits, m_bits, &F)) {
    return nullptr;
  }

  try {
    (void)algorithms::ccz_auto(F, time_limit_seconds, min_active_hyperplanes);
    return BuildAutoGroupResult(F);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "ccz_auto failed");
    return nullptr;
  }
}

PyObject* PyEAAuto(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* truth_table_obj = nullptr;
  PyObject* m_bits_obj = Py_None;
  PyObject* time_limit_obj = Py_None;
  PyObject* min_active_hyperplanes_obj = Py_None;
  int n_bits = 0;
  static const char* kKwlist[] = {"truth_table", "n_bits", "m_bits",
                                   "time_limit_seconds",
                                   "min_active_hyperplanes", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi|OOO",
                                   const_cast<char**>(kKwlist),
                                   &truth_table_obj, &n_bits, &m_bits_obj,
                                   &time_limit_obj,
                                   &min_active_hyperplanes_obj)) {
    return nullptr;
  }

  int m_bits = 0;
  if (!ParseOptionalBits(m_bits_obj, "m_bits", n_bits, &m_bits)) return nullptr;
  double time_limit_seconds = 0.0;
  if (!ParseOptionalTimeLimit(time_limit_obj, 0.0, &time_limit_seconds)) {
    return nullptr;
  }
  std::size_t min_active_hyperplanes = 0;
  if (!ParseOptionalSizeT(min_active_hyperplanes_obj, "min_active_hyperplanes",
                          0, &min_active_hyperplanes)) {
    return nullptr;
  }

  GraphData F;
  if (!BuildGraphFromTruthTableArgs(truth_table_obj, n_bits, m_bits, &F)) {
    return nullptr;
  }

  try {
    (void)algorithms::ea_auto(F, time_limit_seconds, min_active_hyperplanes);
    return BuildAutoGroupResult(F);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "ea_auto failed");
    return nullptr;
  }
}

PyObject* PyCCZEquivalence(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* truth_table_f_obj = nullptr;
  PyObject* truth_table_g_obj = nullptr;
  PyObject* m_bits_obj = Py_None;
  PyObject* time_limit_obj = Py_None;
  PyObject* min_active_hyperplanes_obj = Py_None;
  PyObject* auto_group_obj = Py_None;
  int n_bits = 0;
  static const char* kKwlist[] = {"truth_table_f", "truth_table_g", "n_bits",
                                   "m_bits", "time_limit_seconds",
                                   "min_active_hyperplanes", "auto_group", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOi|OOOO",
                                   const_cast<char**>(kKwlist),
                                   &truth_table_f_obj, &truth_table_g_obj,
                                   &n_bits, &m_bits_obj, &time_limit_obj,
                                   &min_active_hyperplanes_obj, &auto_group_obj)) {
    return nullptr;
  }

  int m_bits = 0;
  if (!ParseOptionalBits(m_bits_obj, "m_bits", n_bits, &m_bits)) return nullptr;
  std::optional<double> time_limit_seconds;
  if (!ParseOptionalTimeLimitAsOptional(time_limit_obj, &time_limit_seconds)) {
    return nullptr;
  }
  std::size_t min_active_hyperplanes = 0;
  if (!ParseOptionalSizeT(min_active_hyperplanes_obj, "min_active_hyperplanes",
                          0, &min_active_hyperplanes)) {
    return nullptr;
  }

  GraphData F;
  GraphData G;
  if (!BuildGraphFromTruthTableArgs(truth_table_f_obj, n_bits, m_bits, &F)) {
    return nullptr;
  }
  if (!BuildGraphFromTruthTableArgs(truth_table_g_obj, n_bits, m_bits, &G)) {
    return nullptr;
  }

  try {
    std::optional<EquivalencePointMap> eq;
    if (auto_group_obj != Py_None) {
      std::vector<groups::Permutation> seed_generators;
      if (!ParseProvidedAutoGroupGenerators(auto_group_obj, G,
                                            /*require_ea=*/false,
                                            &seed_generators)) {
        return nullptr;
      }
      eq = RunCCZEquivalence(F, G, seed_generators, min_active_hyperplanes);
    } else {
      eq = algorithms::ccz_equivalence(F, G, time_limit_seconds,
                                       min_active_hyperplanes);
    }
    if (!eq.has_value()) Py_RETURN_NONE;
    AffineMapData affine_map;
    if (!EquivalencePointMapToAffineMap(*eq, F, /*require_ea=*/false,
                                        &affine_map)) {
      PyErr_SetString(PyExc_RuntimeError,
                      "failed to extract affine equivalence map");
      return nullptr;
    }
    return AmbientGeneratorToPyDict(affine_map);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "ccz_equivalence failed");
    return nullptr;
  }
}

PyObject* PyEAEquivalence(PyObject*, PyObject* args, PyObject* kwargs) {
  PyObject* truth_table_f_obj = nullptr;
  PyObject* truth_table_g_obj = nullptr;
  PyObject* m_bits_obj = Py_None;
  PyObject* time_limit_obj = Py_None;
  PyObject* min_active_hyperplanes_obj = Py_None;
  PyObject* auto_group_obj = Py_None;
  int n_bits = 0;
  static const char* kKwlist[] = {"truth_table_f", "truth_table_g", "n_bits",
                                   "m_bits", "time_limit_seconds",
                                   "min_active_hyperplanes", "auto_group", nullptr};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOi|OOOO",
                                   const_cast<char**>(kKwlist),
                                   &truth_table_f_obj, &truth_table_g_obj,
                                   &n_bits, &m_bits_obj, &time_limit_obj,
                                   &min_active_hyperplanes_obj, &auto_group_obj)) {
    return nullptr;
  }

  int m_bits = 0;
  if (!ParseOptionalBits(m_bits_obj, "m_bits", n_bits, &m_bits)) return nullptr;
  std::optional<double> time_limit_seconds;
  if (!ParseOptionalTimeLimitAsOptional(time_limit_obj, &time_limit_seconds)) {
    return nullptr;
  }
  std::size_t min_active_hyperplanes = 0;
  if (!ParseOptionalSizeT(min_active_hyperplanes_obj, "min_active_hyperplanes",
                          0, &min_active_hyperplanes)) {
    return nullptr;
  }

  GraphData F;
  GraphData G;
  if (!BuildGraphFromTruthTableArgs(truth_table_f_obj, n_bits, m_bits, &F)) {
    return nullptr;
  }
  if (!BuildGraphFromTruthTableArgs(truth_table_g_obj, n_bits, m_bits, &G)) {
    return nullptr;
  }

  try {
    std::optional<EquivalencePointMap> eq;
    if (auto_group_obj != Py_None) {
      std::vector<groups::Permutation> seed_generators;
      if (!ParseProvidedAutoGroupGenerators(auto_group_obj, G,
                                            /*require_ea=*/true,
                                            &seed_generators)) {
        return nullptr;
      }
      eq = RunEAEquivalence(F, G, seed_generators, min_active_hyperplanes);
    } else {
      eq = algorithms::ea_equivalence(F, G, time_limit_seconds,
                                      min_active_hyperplanes);
    }
    if (!eq.has_value()) Py_RETURN_NONE;
    AffineMapData affine_map;
    if (!EquivalencePointMapToAffineMap(*eq, F, /*require_ea=*/true,
                                        &affine_map)) {
      PyErr_SetString(PyExc_RuntimeError,
                      "failed to extract affine equivalence map");
      return nullptr;
    }
    return AmbientGeneratorToPyDict(affine_map);
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "ea_equivalence failed");
    return nullptr;
  }
}

PyMethodDef kMethods[] = {
    {"ccz_auto", reinterpret_cast<PyCFunction>(PyCCZAuto),
     METH_VARARGS | METH_KEYWORDS,
     "ccz_auto(truth_table, n_bits, m_bits=None, time_limit_seconds=None, min_active_hyperplanes=None) -> dict"},
    {"ea_auto", reinterpret_cast<PyCFunction>(PyEAAuto),
     METH_VARARGS | METH_KEYWORDS,
     "ea_auto(truth_table, n_bits, m_bits=None, time_limit_seconds=None, min_active_hyperplanes=None) -> dict"},
    {"ccz_equivalence", reinterpret_cast<PyCFunction>(PyCCZEquivalence),
     METH_VARARGS | METH_KEYWORDS,
     "ccz_equivalence(truth_table_f, truth_table_g, n_bits, m_bits=None, time_limit_seconds=None, min_active_hyperplanes=None, auto_group=None) -> affine dict | None"},
    {"ea_equivalence", reinterpret_cast<PyCFunction>(PyEAEquivalence),
     METH_VARARGS | METH_KEYWORDS,
     "ea_equivalence(truth_table_f, truth_table_g, n_bits, m_bits=None, time_limit_seconds=None, min_active_hyperplanes=None, auto_group=None) -> affine dict | None"},
    {nullptr, nullptr, 0, nullptr}};

PyModuleDef kModuleDef = {
    PyModuleDef_HEAD_INIT,
    "ccz_bindings",
    "Python bindings for CCZ/EA auto and equivalence algorithms.",
    -1,
    kMethods};

}  // namespace

PyMODINIT_FUNC PyInit_ccz_bindings(void) { return PyModule_Create(&kModuleDef); }
