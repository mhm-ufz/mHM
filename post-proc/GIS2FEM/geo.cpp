/*!
  \file geo.cpp

   Definition of member fucntions of class Point and Polyline, and the function 
   used to read geometrical data

   14.04.2011. WW   
*/

#include <vector>
#include <sstream>
#include <cfloat>
#include <cstdlib>

#include "geo.h"


namespace _GEO_MISC
{

/*!
\fn  DeleteArray(num *an_array)
    
  Release the memory of arrary allocated by using 'new'. 

  10.2010
*/
/*!
\fn  DeleteArray(num *an_array)
    
  Release the memory of arrary allocated by using 'new'. 

  10.2010
*/
template<class num> void  DeleteArray(num *an_array)
{
   if(an_array) delete [] an_array;
   an_array = NULL; 
}
template<class num> void  DeleteVector(vector<num*> &a_vec)
{
   
   while (a_vec.size()>0)
   {
      delete a_vec[(int)a_vec.size()-1];
      a_vec.pop_back();
   }
   
}

/*!
     Computer area of a triangle
*/
double ComputeDetTri(const double *x1, const double *x2,
                                const double *x3)
{
    double u[3], v[3], z[3];
    
    u[0] = x3[0] - x1[0];	
    u[1] = x3[1] - x1[1];
    u[2] = x3[2] - x1[2];

    v[0] = x2[0] - x1[0];	
    v[1] = x2[1] - x1[1];
    v[2] = x2[2] - x1[2];
 
    z[0] = u[1]*v[2] - u[2]*v[1];
    z[1] = u[2]*v[0] - u[0]*v[2];
    z[2] = u[0]*v[1] - u[1]*v[0];

    return 0.5*sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2] );   
} 
 string string_To_lower(string strToConvert)
{
   for(unsigned int i=0;i<strToConvert.length();i++)
   {
      strToConvert[i] = tolower(strToConvert[i]);
   }
   return strToConvert;
}; 


   /*!
      \fn Point::Write(ostream &os)
        
       output the coordinate of a point
   */
   void Point::Write(ostream &os)
   {
      os<<index<<" "<< coordinates[0]<<" "<<coordinates[1]<<"  0."<<endl;      
   }
   void Point::Write_VTK(ostream &os)
   {
      os<< coordinates[0]<<" "<<coordinates[1]<<"  0."<<endl;      
   }


   /*!
      \fn constructor of class Polyline
      
   */
   Polyline::Polyline(ifstream &ins, string ply_name): name(ply_name)
   {
      long id; 
      string aline;
      std::stringstream ss;
  
      for(;;)
      {
         getline(ins, aline); 
         if(aline.find("...")!=string::npos)
           break; 
         
         ss.str(aline);
         ss>> id ;
         ss.clear();

         points.push_back(GetPointByID(id)); 
      }
         
   }
    
   /// Destructor
   Polyline::~Polyline()
   {
      points.clear();
   }

   /*!
      \fn Point::Write(ostream &os)
        
       output the coordinate of a polyline
   */
   void Polyline::Write(ostream &os)
   {
      os<<"--- polyline "<<name<<endl;
      for(int i=0; i<(int)points.size(); i++)
        os<<points[i]->Index()<<endl;  
      os<<"...\n"<<endl;    
   }

   /*!
     \fn  PointInDomain(double x, double y)

     Determine whether a point is in the domain 

     04.2011. WW
   */ 
   bool Polyline::PointInDomain(double x, double y)
   {
      Point *pnt_i, *pnt_j; 
      double xi, yi, xj, yj;
      int size = (int)points.size();
      int   i, j;
      j = size-1 ;
      bool  inside = false;

      for (i=0; i<size; i++) 
      {
         pnt_i = points[i];
         pnt_j = points[j];
         xi =  pnt_i->X();
         yi =  pnt_i->Y();
         xj =  pnt_j->X();
         yj =  pnt_j->Y();

         if ( yi<y && yj>=y || yj<y && yi>=y) 
         {
            if (xi+(y-yi)/(yj-yi)*(xj-xi)<x)
            {
              inside = !inside; 
            }
         }
         j=i;
       }

       return inside;  
   }
   /*!
     \fn real MinDisttanceTo_a_Point(Point *pnt);

      Calculate the distance to a point 
      
      05.2011 WW   
   */ 
   double Polyline::MinDisttanceTo_a_Point(const Point *pnt)
   {
      int i, j;
      double a[3], b[3], c[3], dist; 
      double min_dist = DBL_MAX;
      for(i=0; i<(int)points.size()-1; i++) 
      {
         for(j=0; j<3; j++) 
         {
            a[j] = pnt->coordinates[j];
            b[j] = points[0]->coordinates[j];
            c[j] = points[1]->coordinates[j];
         }
         dist = 2.0*ComputeDetTri(a,b,c)/points[0]->GetDistanceTo(points[1]);
         if(dist<min_dist)
           min_dist = dist;
      }
      return min_dist; 
   } 
}

using namespace _GEO_MISC;
vector<_GEO_MISC::Point*> points;
vector<_GEO_MISC::Polyline*> polylines;
/*!
  \fn  ReadPolyline()
  Read geometrical data

  14.04.2011. WW
*/
void GeoRead(string file_name)
{   
   long id; 
   float xy[2];
   string aline;
   std::stringstream ss;
     
   string geo_fname = file_name+".geo"; 
   ifstream ins(geo_fname.c_str());

   if(!ins.good()) 
   {
      cout<<"Could not find file "<<geo_fname<<". Stop now!"<<endl;
      exit(1);
   } 

   while(!ins.eof())
   {
      getline(ins, aline); 

      aline = string_To_lower(aline);
      if(aline.find("point")!=string::npos)  
      {   
         for(;;)
         {
            getline(ins, aline); 
            if(aline.find("...")!=string::npos)
              break; 
         
            ss.str(aline);
            ss>> id >> xy[0] >> xy[1];
            ss.clear();

            points.push_back(new _GEO_MISC::Point(id, xy[0], xy[1])); 
         }
      }

      if(aline.find("polyline")!=string::npos) 
      {
         ss.str(aline);
         // Skip "---" and "polyline"
         ss>>aline>>aline>>aline;
         ss.clear();
         polylines.push_back(new _GEO_MISC::Polyline(ins, aline)); 

      }
       

   }

} 

/*!
    \fn Polyline *GetPolylineByName(string name); 
   
     Find a polyline by name 
      
     04.2011. WW 
*/
_GEO_MISC::Polyline *GetPolylineByName(string name)
{
   int i;
  
   for(i=0; i<(int)polylines.size(); i++)
   {
      if(polylines[i]->Name().find(name)!=string::npos)
      {
         return polylines[i]; 
      }
   }          
   return NULL;
}
/*!
    \fn FDM::Point *GetPointByID(long ID); 
   
     Find a polyline by name 
      
     04.2011. WW 
*/

_GEO_MISC::Point *GetPointByID(long ID)
{
   int i;
  
   for(i=0; i<(int)points.size(); i++)
   {
      if(points[i]->Index() == ID)
      {
         return points[i]; 
      }
   }          
   return NULL;
}
//-------------------------------------------
/*!
   \fn  WriteGeoData(ostream &os = cout);
   
    Output Geometrical Data
*/
void WriteGeoData(ostream &os)
{
   size_t i;

   os<<"--- point"<<endl;
   for(i=0; i<points.size(); i++)
     points[i]->Write(os);
   os<<"...\n"<<endl;

  
   for(i=0; i<polylines.size(); i++)
      polylines[i]->Write(os);

}

void GeoReleaseMemory()
{
   DeleteVector(points);
   DeleteVector(polylines);
}

