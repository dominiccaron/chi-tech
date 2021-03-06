#include "../k_eigenvalue_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

using namespace LinearBoltzmann;

//###################################################################
/**Compute the total fission production in the problem.*/
double KEigenvalue::Solver::ComputeProduction()
{
  int first_grp = groups.front().id;
  int last_grp = groups.back().id;

  auto pwl = std::static_pointer_cast<SpatialDiscretization_PWLD>(discretization);

  // ----- Loop over cells
  double local_F = 0.0;
  for (auto& cell : grid->local_cells)
  {
    // ----- Cell information
    int cell_matid = cell.material_id;
    int xs_id = matid_to_xs_map[cell_matid];
    auto xs = material_xs[xs_id];
    auto& cell_fe_view = pwl->GetUnitIntegrals(cell);
    auto& transport_view = cell_transport_views[cell.local_id];

    // ----- Loop over cell dofs
    for (int i = 0; i < cell_fe_view.NumNodes(); ++i)
    {
      size_t ir = transport_view.MapDOF(i, 0, 0);

      double intV_shapeI = cell_fe_view.IntV_shapeI(i);

      // ----- Contribute fission source
      if (xs->is_fissile)
      {
        for (size_t g = first_grp; g <= last_grp; ++g)
        {
          double nu_sigma_f = 0.0;
          if (options.use_precursors)
            nu_sigma_f = xs->nu_prompt_sigma_f[g] +
                         xs->nu_delayed_sigma_f[g];
          else
            nu_sigma_f = xs->nu_prompt_sigma_f[g];

          local_F += nu_sigma_f * phi_new_local[ir + g] * intV_shapeI;
        }
      }
    }//for i
  }//for c

  double global_F = 0.0;
  MPI_Allreduce(&local_F, &global_F, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return global_F;
}