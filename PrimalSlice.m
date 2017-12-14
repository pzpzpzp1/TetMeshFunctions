
%% run from TetMeshFunctions
path2HMeshFuncs='./../HexMeshSingularitySheetingTest';
addpath(path2HMeshFuncs);

% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/bunny_ours.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/bunny.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/rocker_arm.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/polycut2013/fertility.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/sculpture-B_ours.vtk']);

HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/rod_ours.vtk']);
HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/fertility.vtk']);
HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/double_torus.vtk']);
HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/ellipsoid-B.vtk']);

TMesh = HexToTet(HMesh);
X = TMesh.V2P;
T = TMesh.T2V;

%[X,T]=paul_loadTetGenMesh('./meshes/sphere_61k');
%data = paul_getTetData(T,X,1);
data = paul_getTetData(T,X,0);

dtree = DualVolumeVertexSpanningTree(data);











