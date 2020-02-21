#ifndef _chi_region_h
#define _chi_region_h

#include "../MeshContinuum/chi_meshcontinuum.h"
#include "../VolumeMesher/chi_volumemesher.h"



//######################################################### Class definition
class chi_mesh::Region
{
  friend void VolumeMesher::AddContinuumToRegion(MeshContinuum* grid,
                                                 Region& region);
public:
  std::vector<chi_mesh::Boundary*> boundaries;

private:
  std::vector<chi_mesh::MeshContinuum*> volume_mesh_continua;
public:
  chi_mesh::MeshContinuum* GetGrid();
};

#endif