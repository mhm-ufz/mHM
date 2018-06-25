#ifndef misc_INC
#define misc_INC
#include <vector>
#include <string>

using namespace std;

namespace  Mesh_Group{class CFEMesh; }
using  Mesh_Group::CFEMesh;
extern vector<Mesh_Group::CFEMesh*> fem_msh_vector;

extern bool FEMRead(std::string);

void AssignMatID_by_Ply(CFEMesh* a_msh, const vector<int> ly_idx, int a_type);
void mappingPointData2ARCGISShape(std::string data_file, std::string file_path);
#endif
