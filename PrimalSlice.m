
%% run from TetMeshFunctions
path2HMeshFuncs='./../HexMeshSingularitySheetingTest';
addpath(path2HMeshFuncs);

HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/bunny_ours.vtk']);
TMesh = HexToTet(HMesh);
X = TMesh.V2P;
T = TMesh.T2V;

%[X,T]=paul_loadTetGenMesh('./meshes/sphere_61k');
data = paul_getTetData(T,X);
