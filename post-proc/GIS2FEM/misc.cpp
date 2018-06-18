# include <cfloat>

#include "misc.h"
#include "geo.h"

#include "msh_mesh.h"
vector<Mesh_Group::CFEMesh*> fem_msh_vector;

/**************************************************************************
FEMLib-Method: 
Task:
Programing:
03/2010 WW  Based on FEMRead in MSHLib
**************************************************************************/
bool FEMRead(string file_base_name)
{
 
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  CFEMesh *m_fem_msh = NULL;
  char line[1024];
  string sub_line;
  string line_string;
  ios::pos_type position;
  //========================================================================
  // File handling
  bool msh_file_binary = false;
  string msh_file_name_bin = file_base_name + "_binary" + ".msh";
  string msh_file_name_ascii = file_base_name + ".msh";
  ifstream msh_file_bin;
  ifstream msh_file_ascii;

  msh_file_bin.open(msh_file_name_bin.c_str(),ios::binary|ios::in);
  if(msh_file_bin.good()){ 
    msh_file_binary = true;
  }

  //----------------------------------------------------------------------
  cout << "MSHRead: ";

  cout << "ASCII file" << endl;
  msh_file_ascii.open(msh_file_name_ascii.data(),ios::in);
  if (!msh_file_ascii.good()){
    cout<<"Opening MSH file failed"<<endl;
    return false;
  }
  //----------------------------------------------------------------------
  //========================================================================
  // Keyword loop
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  while(!msh_file_ascii.eof()){
    msh_file_ascii.getline(line,1024);
    line_string = line;

    if(line_string.find("$MeshFormat")!=string::npos) 
	{
       m_fem_msh = new CFEMesh();
       fem_msh_vector.push_back(m_fem_msh) ;
	   m_fem_msh->ReadGmsh(&msh_file_ascii);
       break;
	} 

    if(line_string.find("#STOP")!=string::npos)
      return true;
    //..................................................................
    if(line_string.find("#FEM_MSH")!=string::npos) { // keyword found
      m_fem_msh = new CFEMesh();
      position = m_fem_msh->Read(&msh_file_ascii);
      fem_msh_vector.push_back(m_fem_msh);
      msh_file_ascii.seekg(position,ios::beg);
    } // keyword found
  } // eof

  msh_file_ascii.close();
 

  //========================================================================
  return true;
}

using namespace Mesh_Group;
using namespace _GEO_MISC;
void AssignMatID_by_Ply(CFEMesh* a_msh, const vector<int> ly_idx, int a_type)
{
   size_t j, k, l; 
   long i;
   double *e_center;
   CElem *ele;
   bool exclude = false;
   
   int max_matD0 = a_msh->max_mmp_groups+1;
   int max_matD1;
   int mat_idx;  

   vector<int> mat_id;

   for(j=0; j<polylines.size(); j++)
   {
     max_matD1 = 0;
     for(i=0; i<(long)a_msh->ele_vector.size(); i++)
     {
         ele = a_msh->ele_vector[i];

		 if(a_type == 51)
		 {
            exclude = false;  

            mat_idx = ele->GetPatchIndex();
           
            for(k=0; k<ly_idx.size(); k++)
            {
               if(mat_idx==ly_idx[k])
               {
                  exclude = true;
                  break; 
               } 
            }
            if(exclude)
              continue;

            for(k=0; k<ly_idx.size(); k++)
            {
               if(mat_idx>ly_idx[k])
                 mat_idx--;
            }
		 }
         if(a_type == 53)
		 {
			 exclude = true;  
             mat_idx = ele->GetPatchIndex();
             for(k=0; k<ly_idx.size(); k++)
			 {
                if(mat_idx==ly_idx[k])
				{
                   exclude = false;
                   break; 

				}
			 }
             if(exclude)
               continue;

		 }

         e_center = ele->GetGravityCenter();

         if(polylines[j]->PointInDomain(e_center[0], e_center[1]))
         {
             if(a_type == 51)
			 {
                mat_idx += max_matD0; 
                ele->SetPatchIndex(mat_idx); 
                if(max_matD1<mat_idx) 
                   max_matD1 = mat_idx;
			 }
			 else if(a_type == 52)
               ele->SetPatchIndex(max_matD0); 
			 else if(a_type == 53)
			 {
               for(k=0; k<ly_idx.size(); k++)
			   {
                  if(mat_idx==ly_idx[k])
			      {
                     ele->SetPatchIndex(max_matD0+k); 
                     break; 
	              }
		       }

			 }

         }
        
         if(a_type == 51)
		 {
            /// Check whether this material index is saved at  mat_id
            exclude = false;
            for(l=0; l<mat_id.size(); l++)
            {
               if(mat_id[l] == ele->GetPatchIndex())
               {
                  exclude = true;
                  break;
               }              
            }
            if(!exclude)
              mat_id.push_back(ele->GetPatchIndex());
		 }
          
      }
	   
	  if(a_type == 51)
         max_matD0 = max_matD1;
      else if(a_type == 52)
         max_matD0++;
      else if(a_type == 53)
          max_matD0 += static_cast<int>(ly_idx.size());

   }

   if(a_type != 51)
      return;

   /// Re-order material ID
   max_matD0 = 0;
   for(i=0; i<(long)a_msh->ele_vector.size(); i++)
   {
      ele = a_msh->ele_vector[i];
      for(l=0; l<mat_id.size(); l++)
      {
         if(mat_id[l] == ele->GetPatchIndex())
         {
            ele->SetPatchIndex(l);
            break;
         }              
      }
   }
}

// 05.2013. WW
void mappingPointData2ARCGISShape(std::string data_file,  std::string file_path)
{
  std::string aline;
  std::stringstream ss;
  std::string key, sub_key;
  std::ifstream ins;
  CNode *node;

  size_t i, ncols, nrows;
  double *x;
  vector<string> names;


  double x_min, y_min, x_max, y_max;
  x_min = DBL_MAX;
  y_min = DBL_MAX;
  x_max = DBL_MIN;
  y_max = DBL_MIN;


  ins.open(data_file.c_str());
  getline(ins, aline);  

  ss.str(aline);

  ss >> key >> key ;
  while(!ss.eof())
  {
     ss >> key >> sub_key;
	 if (!key.empty() && !sub_key.empty()  )
	 {
         key = key + sub_key + ".asc";
		 names.push_back(key);
	 }		  
  }
  ss.clear();
  ncols = names.size();



  // Output file names
  string ofname = file_path + "3Dmesh.pcp";
  ofstream of(ofname.c_str(), ios::trunc);
  of << "Unit: Day" << endl;
  of << "Ratio: 0.8" << endl;
  for(size_t k=0; k<ncols; k++)
  {
	  of << names[k] << endl;
	  names[k] = file_path + names[k];
  }
  of.clear();
  of.close();




  getline(ins, aline);  
  ss.str(aline);
  ss >> nrows;
  ss.clear();
 
  x = new double[2*nrows]; 

  //double *rchg_data =  new double[nrows * ncols];   
  vector<double*> rchg_data(ncols);
  for(i=0; i<ncols; i++)
  {
     double *new_row = new double[nrows];
     rchg_data[i] = new_row; 
  }


  for(i=0; i<nrows; i++)
  {
     double xx, yy;
     getline(ins, aline);  
     ss.str(aline);
     ss >> xx;
	 ss >> yy;

	  if(xx < x_min)
        x_min = xx; 
	  if(xx > x_max)
        x_max = xx; 
	  if(yy < y_min)
        y_min = yy; 
	  if(yy > x_max)
        y_max = yy; 

     x[2*i] = xx;
     x[2*i+1] = yy;

     for(size_t k=0; k<ncols; k++)
	 {
        double *col_dat = rchg_data[k];
        ss >> col_dat[i]; //rchg_data[i*ncols+k];  
	 }
	 ss.clear();
  }


  //
  const double d_val = -9999.0;
  const double cell_size = 100.0;
  size_t cols = static_cast<size_t>((x_max-x_min)/cell_size) ;
  size_t rows = static_cast<size_t>((y_max-y_min)/cell_size) ;
  const size_t ext_cells = 5;

  double *gdata = new double [cols*rows];
  const size_t size =  cols * rows;

  for(i=0; i<size; i++)
  {
	  gdata[i] = d_val;
  }


  ofstream os;
  for(i=0; i<ncols; i++)
  {

	 os.open(names[i].c_str(), ios::trunc);
	 os << "ncols "<<cols << endl;
	 os << "nrows "<<rows << endl;
     os << "xllcorner " <<  x_min << endl;
     os << "yllcorner " <<  y_min << endl;
     os << "cellsize " << cell_size <<endl;
     os << "NODATA_value "<< d_val <<endl;

     double *col_dat = rchg_data[i];

	 for (size_t k=0; k<nrows; k++)
	 {
        size_t i_col = static_cast<size_t>((x[2*k]-x_min)/cell_size);
        size_t i_row = rows - static_cast<size_t>((x[2*k+1]-y_min)/cell_size);
        if(i_row <0 )
           i_row = 0;
		if(i_row > rows-1)
           i_row = rows-1;
        if(i_col <0 )
           i_col = 0;
		if(i_col > cols-1)
           i_col = cols-1;
        
        gdata[cols*i_row + i_col] =  col_dat[k]; //rchg_data[k*ncols+i]; 
	 }


     //smear
 	 for (size_t k=0; k<rows; k++)
	 {
         for (size_t j=0; j<cols; j++) 
         {
            if (fabs(gdata[cols*k + j] - d_val) < DBL_EPSILON) 
			{
                double v1, v2, v3, v4;
				v1 =v1 = v2 =v3 = v4 = 0.;
                for(size_t m=1; m<ext_cells; m++)
				{
                   size_t kr = k-m;
                   size_t jc = j-m;
                   if(kr < 0 || kr >= rows)
                      continue; 
                   if(jc < 0 || jc >= cols)
                      continue; 
                  
				   if((fabs(gdata[cols*kr + jc] - d_val) > DBL_EPSILON))
				   {
                      v1 =  gdata[cols*kr + jc];

				   }                      
				}
                


			}
		 }
	 }



	 for (size_t k=0; k<rows; k++)
	 {
         for (size_t j=0; j<cols; j++) 
         {
            os << gdata[cols*k + j] << " ";
		 }
		 os << endl;
	 }






	 os.clear();
	 os.close();

  } 



  for(i=0; i<ncols; i++)
  {
    delete [] rchg_data[i];   
  }

  delete [] x;
}
