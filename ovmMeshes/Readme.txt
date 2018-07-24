This folder contains a set of input files for the SCOF command line tool in order to generate an octahedral field for hexmeshing(visit https://.../libSCOF for more information).
The input file is in OpenVolumeMesh file format(see https://www.openvolumemesh.org). It stores the tetrahedral mesh as well as the constraints saved as properties. The file starts with the mesh information which are vertex coordinates, edges, faces, and tetrahedra. Then follows the input constraints extracted from a valid singularity graph. In principle, there are four types of constraints including edge type, corner constraint, sector constraint and tangent continuity constraint.
Edge type is stored as an integer in edge property. For interior edges, there are three different integers -1, 0, 1, meaning the edge index -1/4, 0, 1/4; for boundary edges, we use -1, 0, 1, 2 to represent the edge index -1/4, 0, 1/4 and 1/2 respectively.
The corner has the data structure as follows: 
typedef corner{
VertexHandle center_vh, vh0, vh1, vh2; 
vector<HalffaceHandle> dualPath0, dualPath1, dualPath2;
}. The "vhi" indicates the outgoing Halfedge starting from the center vertex. The "dualPathi" means the dual path from the tetrahedron which is incident to the corresponding Halfedge to the common tetrahedron in the corner. The vertex handles of corners are stored as integers while the dual paths are stored as vectors of Halfface handles sequentially in mesh property.
The sector's data structure is: 
typedef sector{
VertexHandle center_vh, vh0, vh1;  
double angle}. It encodes two outgoing Halfedges starting from the center vertex to "vhi" which span an angle.
The data structure of the tangent continuity is the following: 
typedef tangentContinuity{
VertexHandle center_vh, vh0, vh1; 
vector<HalffaceHandle> dualPath;}. These variables basically share the same meaning as in the corner. The two Halfedges implied by the three vertices are connected via the "dualPath".
To read the property, e.g. edge type, simply call the function:
OpenVolumeMesh::EdgePropertyT<int> valance =  mesh_.template request_edge_property<int>("edge_valance”);

