/*!
  \file geo.h
   
   Declaration of class Point and Polyline, which are used
   to assign  boundary conditions
   
   13.04.2011. WW
*/
#ifndef geo_INC
#define geo_INC

#include "misc.h"
#include<iostream>
#include<fstream>
#include<cmath>

namespace _GEO_MISC
{
   
   /*!
     \class  Point
         
   */
  
   /// class Geo_Root;
   class Point
   {
      public:
        Point(long id, float x, float y) 
        {
           index = id;
           coordinates = new float[3];
           coordinates[0] = x;
           coordinates[1] = y;
           coordinates[2] = 0.;
        }
        
        ~Point()
        {
           delete [] coordinates;            
        }  

        float X() const {return coordinates[0];}
        float Y() const {return coordinates[1];}
        float GetDistanceTo(const Point *a_p)
           {
              float *a_coord = a_p->coordinates; 
              return sqrt((coordinates[0]-a_coord[0])
                        *(coordinates[0]-a_coord[0])
                      +  (coordinates[1]-a_coord[1])
                        *(coordinates[1]-a_coord[1]));
            }
 
        void Write(ostream &os = cout);
        void Write_VTK(ostream &os = cout);
        long Index() const {return index;}        

      private:
        long index;
        long grid_i;
        long grid_j;

        float *coordinates; 

        friend class Polyline;
   };

   /*!
      \class Polyline
          
      Define a polyine that consists of points
   */
    class Polyline
    {
        public:
          Polyline(ifstream &ins, string ply_name);
          ~Polyline();

        string Name() const {return name;}

        void Write(ostream &os = cout);
        
        bool PointInDomain(double x, double y);
        double MinDisttanceTo_a_Point(const Point *pnt);

        private:
          vector<Point*> points;
          string name;
          friend class Polyline;
          friend class FiniteDifference;
          friend class ConditionData;           
    };
   
}

/// Contains all points of the domain
extern vector<_GEO_MISC::Point*> points;
extern vector<_GEO_MISC::Polyline*> polylines;
extern void GeoRead(string file_name); 
extern void GeoReleaseMemory(); 
extern _GEO_MISC::Polyline *GetPolylineByName(string name); 
extern _GEO_MISC::Point *GetPointByID(long ID); 
extern void WriteGeoData(ostream &os = cout);



#endif