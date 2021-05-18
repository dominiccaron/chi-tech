-- 1D transport test in pointsymmetric spherical geometry with
-- vacuum boundary condition - monoenergetic.
-- SDM: PWLD
-- Test: Max-value=1.00000
num_procs = 2
--Structured mesh




--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--------------------------------------------------------------------------------
--  mesh
--------------------------------------------------------------------------------
chiMeshHandlerCreate()
dim = 1
length = {1, }
ncells = {50, }
nodes = {}
for d = 1, dim do
  delta = length[d]/ncells[d]
  nodes[d] = {}
  for i = 0, ncells[d] do
    nodes[d][i+1] = i*delta
  end
end
surf_mesh, region0 = chiMeshCreateUnpartitioned1DOrthoMesh(nodes[1])
chiVolumeMesherSetProperty(PARTITION_TYPE, PARMETIS)
chiVolumeMesherExecute()

vol0 = chiLogicalVolumeCreate(RPP, 0, 0, 0, 0, 0, length[1])
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--------------------------------------------------------------------------------
--  materials
--------------------------------------------------------------------------------
ngrp = 1
sigmat = 25.0
ratioc = 0.1
source = sigmat * (1 - ratioc)

material0 = chiPhysicsAddMaterial("Material_0");
chiPhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(material0,ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialSetProperty(material0, TRANSPORT_XSECTIONS,
                              SIMPLEXS1, ngrp, sigmat, ratioc)
chiPhysicsMaterialSetProperty(material0, ISOTROPIC_MG_SOURCE,
                              SINGLE_VALUE, source)

--------------------------------------------------------------------------------
--  physics
--------------------------------------------------------------------------------
phys0 = chiLBSCurvilinearCreateSolver(LBSCurvilinear.SPHERICAL)
chiSolverAddRegion(phys0, region0)

--  angular quadrature
pquad = chiCreateSphericalProductQuadrature(GAUSS_LEGENDRE, 8)

--  groups
groups = {}
for g = 1, ngrp do
  groups[g] = chiLBSCreateGroup(phys0)
end

--  groupsets
gs0 = chiLBSCreateGroupset(phys0)
chiLBSGroupsetAddGroups(phys0, gs0, 0, ngrp-1)
chiLBSGroupsetSetQuadrature(phys0, gs0, pquad)
chiLBSGroupsetSetAngleAggregationType(phys0, gs0, LBSGroupset.ANGLE_AGG_AZIMUTHAL)
chiLBSGroupsetSetIterativeMethod(phys0, gs0, NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys0, gs0, 1.0e-12)
chiLBSGroupsetSetMaxIterations(phys0, gs0, 100)
chiLBSGroupsetSetGMRESRestartIntvl(phys0, gs0, 30)

--  spatial discretisation
chiLBSSetProperty(phys0, DISCRETIZATION_METHOD, PWLD)

--  scattering order
chiLBSSetProperty(phys0, SCATTERING_ORDER, 0)

--------------------------------------------------------------------------------
--  boundary conditions
--------------------------------------------------------------------------------
dirichlet_value = {}
for g = 1, ngrp do
  dirichlet_value[g] = 0
end
chiLBSSetProperty(phys0, BOUNDARY_CONDITION,
                  ZMIN, LBSBoundaryTypes.REFLECTING)
chiLBSSetProperty(phys0, BOUNDARY_CONDITION,
                  ZMAX, LBSBoundaryTypes.INCIDENT_ISOTROPIC, dirichlet_value)

--------------------------------------------------------------------------------
--  solvers
--------------------------------------------------------------------------------
chiLBSInitialize(phys0)
chiLBSExecute(phys0)

--------------------------------------------------------------------------------
--  output
--------------------------------------------------------------------------------
--  field functions
fflist, count = chiLBSGetScalarFieldFunctionList(phys0)
if master_export == nil then
  chiExportFieldFunctionToVTKG(fflist,
                               "Transport2DSph1-scalar_flux",
                               "scalar_flux")
end

--  volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi, OPERATION, OP_MAX)
chiFFInterpolationSetProperty(curffi, LOGICAL_VOLUME, vol0)
chiFFInterpolationSetProperty(curffi, ADD_FIELDFUNCTION, fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0, string.format("Max-value=%.5f", maxval))
