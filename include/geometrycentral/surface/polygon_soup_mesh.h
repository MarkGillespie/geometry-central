#pragma once

#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
#include <geometrycentral/utilities/utilities.h>

namespace geometrycentral {

class PolygonSoupMesh {
public:
  PolygonSoupMesh();
  PolygonSoupMesh(std::string meshFilename, bool loadTexture = false);
  PolygonSoupMesh(const std::vector<std::vector<size_t>>& polygons_, const std::vector<Vector3>& vertexCoordinates_);

  // Mutate this mesh and by naively triangulating polygons
  void triangulate();

  // Mesh data
  std::vector<std::vector<size_t>> polygons;
  std::vector<Vector3> vertexCoordinates;

  std::vector<std::vector<Vector2>> textureCoordinates;

private:
  void readMeshFromFile(std::string filename, bool loadTexture = false);
};

} // namespace geometrycentral
