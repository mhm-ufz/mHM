// expre_new_Operator.cpp
// compile with: /EHsc
#include <cmath>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include <float.h>

//#include <time.h>

#include "msh_mesh.h"
#include "geo.h"
#include "misc.h"


std::string FileName;
std::string FilePath; //WW



 using namespace std;
int main(int argc, char* argv[])
{
 

  const int max_size = 1028; 
  char str1[max_size];
  char str2[max_size];
  vector<int> exclude_ly;


  int option = 1;
  CFEMesh *a_mesh = NULL;
  fstream ofile;
  string ofname;

 
  cout<<"\t|======================================================| "<<endl;
  cout<<"\t|                                                      | "<<endl;
  cout<<"\t|        Toolkit for hydraulic modeling                | "<<endl;
  cout<<"\t|                                                      | "<<endl;
  cout<<"\t|                 By    WW                             | "<<endl;
  cout<<"\t|                                                      | "<<endl;
  cout<<"\t|    Argument:                                         | "<<endl;
  cout<<"\t|      option and file name (no extension)             | "<<endl;
  cout<<"\t|    Option:                                           | "<<endl;
  cout<<"\t|     1.  Generation OGS Neumman BC from               | "<<endl;
  cout<<"\t|         raster files of GIS (3D)                     | "<<endl;
  cout<<"\t|     2.  Generation OGS Neumman BC from               | "<<endl;
  cout<<"\t|         raster files of GIS (2D)                     | "<<endl;
  cout<<"\t|     3.  Top surface integration for 3D mesh          | "<<endl;
  cout<<"\t|     31. Find top/bottom surface nodes of 3D mesh     | "<<endl;
  cout<<"\t|     32. Ouput elavation of nodes on the top surface  | "<<endl;
  cout<<"\t|     33. Ouput elavation of the selected nodes on top | "<<endl;
  cout<<"\t|     4.  Convert GIS raster cells into FE mesh        | "<<endl;
  cout<<"\t|     51. Specify meterial groups in a mesh.           | "<<endl;
  cout<<"\t|         Assign ID for differet levels in the areas   | "<<endl;
  cout<<"\t|         defined by polylines.                        | "<<endl;
  cout<<"\t|     52. Specify meterial groups in a mesh.           | "<<endl;
  cout<<"\t|         Assign the same ID for the same area         | "<<endl;
  cout<<"\t|     53. Specify meterial groups in a mesh            | "<<endl;
  cout<<"\t|         Change IDs of the spefcified material domain | "<<endl;
  cout<<"\t|     6.  Write the top surface                        | "<<endl;
  cout<<"\t|     7.  2D finite element wise recharge data to      | "<<endl;
  cout<<"\t|         Neumann BC of the top surface                | "<<endl;
  cout<<"\t|     8.  Point-wise recharge data to GIS Shape file   | "<<endl;
  cout<<"\t|                                                      | "<<endl;
  cout<<"\t|======================================================| "<<endl;
  cout<<"\tInput file name: ";

  if(argc>1) 
  {
    std::strcpy(str2,argv[1]);
    std::strcpy(str1,argv[2]);
  }
  else 
      scanf(" %s %s%*[^\n]%*c",str2, str1);

  sscanf(str2, "%d", &option);
 
  double  time_cpu = -clock();

  FileName = str1;
  basic_string <char>::size_type indexChWin, indexChLinux; 
  indexChWin = indexChLinux = 0;
  indexChWin = FileName.find_last_of('\\');
  indexChLinux = FileName.find_last_of('/');
  //
  if(indexChWin!=string::npos)
     FilePath = FileName.substr(0,indexChWin)+"\\";
  else if(indexChLinux!=string::npos)
     FilePath = FileName.substr(0,indexChLinux)+"/";

  if(option!=4 && option!=7 && option!=8)
  {
     FEMRead(FileName);
     a_mesh = fem_msh_vector[0];
  }
 

  std::string aline;
  std::stringstream ss;
  std::string key;
  std::ifstream ins;
  std::string infiltration_files;

  switch(option)
  {
    case 1:
      a_mesh->mHM2NeumannBC(); 
      break;
    case 2:

      infiltration_files = FileName+".pcp";
      ins.open(infiltration_files.c_str());
      a_mesh->ConstructGrid();
      while(!ins.eof())
      {
         getline(ins, aline);
         ss.str(aline);
         ss>> key>>key;
         ss.clear();

         if(key.size()==0)                        // An empty line
            continue;

         if(key.find("#STOP")!=std::string::npos)
            break;

         ofname = FilePath+key+"_node_source.asc";

         key = FilePath+key;
         a_mesh->Precipitation2NeumannBC(key, ofname, true, 1.0);
      }
      break;
    case 3:
      a_mesh->TopSurfaceIntegration(); 
      break;
    case 31:
      a_mesh->ConstructGrid(); 
      //a_mesh->GenerateHighOrderNodes();
      a_mesh->MarkInterface_mHM_Hydro_3D(false, true); 
      break;
    case 32:
      a_mesh->Output_Z_TopSurface(); 
      break;
    case 33:
      a_mesh->Output_Z_TopSurface(true); //Selected points
      break;
    case 4:
      a_mesh = new CFEMesh();
      a_mesh->ConvertShapeCells(FileName+".asc");
      ofname = FileName+".msh";
      ofile.open(ofname.c_str(), ios::out|ios::trunc);
      a_mesh->Write(&ofile);
      break;

    case 51:      
      aline = FileName+".ecl";
      ins.open(aline.c_str());
      if(ins.good())
      {
         while(!ins.eof())
         {
            getline(ins, aline);
            if(aline.find("...")!=string::npos)
              break;
            if(aline.find("--- Exclude layer")!=string::npos)
              continue;

			 ss.str(aline);
             int layer_idx;
             ss>> layer_idx;
             ss.clear();
             exclude_ly.push_back(layer_idx);
         }   
      }        
      
      
      GeoRead(FileName);
      a_mesh->ConstructGrid(); 
      AssignMatID_by_Ply(a_mesh, exclude_ly, 51);
      ofname = FileName+"_new.msh";
      ofile.open(ofname.c_str(), ios::out|ios::trunc);
      a_mesh->Write(&ofile);
      break;

    case 52:      
      GeoRead(FileName);
      a_mesh->ConstructGrid(); 
      AssignMatID_by_Ply(a_mesh, exclude_ly, 52);
      ofname = FileName+"_new.msh";
      ofile.open(ofname.c_str(), ios::out|ios::trunc);
      a_mesh->Write(&ofile);
      break;

    case 53:      
      aline = FileName+".cgr";
      ins.open(aline.c_str());
      if(ins.good())
      {
         while(!ins.eof())
         {
            getline(ins, aline);
            if(aline.find("...")!=string::npos)
              break;
            if(aline.find("--- ID to be changed")!=string::npos)
              continue;

             ss.str(aline);
             int layer_idx;
             ss>> layer_idx;
             ss.clear();
             exclude_ly.push_back(layer_idx);
         }   
      }        
      
      
      GeoRead(FileName);
      a_mesh->ConstructGrid(); 
      AssignMatID_by_Ply(a_mesh, exclude_ly, 53);
      ofname = FileName+"_new.msh";
      ofile.open(ofname.c_str(), ios::out|ios::trunc);
      a_mesh->Write(&ofile);
      break;

	case 6: 
      a_mesh->ConstructGrid();
      a_mesh-> MarkInterface_mHM_Hydro_3D();
      a_mesh->Write_Surface_GMSH(FileName);
      break;
	case 7:
       ins.open(FileName.c_str());
       
	   double ratio;
	   ratio = 1.;
	   getline(ins, aline);  
	   ss.str(aline);
	   ss >> aline >> ratio;
	   ss.clear();

	   int unit_type;
	   unit_type = 1;
	   getline(ins, aline);  
       if(aline.find("year") != string::npos)
         unit_type = 2; 

	   getline(ins, aline);  
       aline = FilePath +  aline;
	   FEMRead(aline);
       a_mesh = fem_msh_vector[0];

	   getline(ins, aline);  
       aline = FilePath +  aline;
	   a_mesh->TwoDRechargeDatato3DSurface(aline, FilePath, ratio, unit_type);
       break;  
	case 8:
       mappingPointData2ARCGISShape(FileName, FilePath);
       break; 

  }

  if(a_mesh)
    delete a_mesh;

  time_cpu += clock();
  cout<<"\n\tCPU time elapsed: "  <<(double)time_cpu / CLOCKS_PER_SEC<<"s"<<endl;

  return 0;
}
