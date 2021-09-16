/**
  @file water.surf.transf.SU_modified.cpp
*/


/*
<sim2doc>

</sim2doc>
*/


#include <openfluid/ware/PluggableSimulator.hpp>
#include <openfluid/tools/DataHelpers.hpp>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("water.surf.transf.RS.weir.id");

 // Informations
  DECLARE_NAME("Water transfer on surface units using the weir equation");
  DECLARE_DESCRIPTION("Combination of submerged and not submerged weir equation with 4th order Runge Kutta solution of ODE");
  DECLARE_VERSION("");
  DECLARE_STATUS(openfluid::ware::EXPERIMENTAL);

  DECLARE_AUTHOR("Christina Orieschnig", "christina.orieschnig@ird.fr");

  // Required and used variables

  //DECLARE_REQUIRED_VARIABLE("water.surf.H", "SU", "water height on surface of SU at the previous time step", "m");
  DECLARE_USED_VARIABLE("water_level_input", "RS", "water level in the Bassac", "m") ;
  DECLARE_USED_VARIABLE("z_new_flow", "RS", "water level in the RS", "m") ;

  // Produced variable

  DECLARE_PRODUCED_VARIABLE("z_new_flow", "RS", "water level in the RS", "m");


  // Required attributes for the test simulator:

  DECLARE_REQUIRED_ATTRIBUTE("length", "RS", "length of RS", "m");
  DECLARE_REQUIRED_ATTRIBUTE("conn_len", "RS", "intersection between the two RS and other RS and SU - basic weir length", "m");
  DECLARE_REQUIRED_ATTRIBUTE("OFLD_TO", "RS", "list of connections between the RS and other RS and SU", "");


  DECLARE_REQUIRED_ATTRIBUTE("elev", "RS", "mean elev of RS", "m"); /// this should really be pulled from the elevation-area curve

  // Simulator parameters (for later on)

  //DECLARE_USED_PARAMETER("coeff", "basic weir coefficient", "-");
  //DECLARE_USED_PARAMETER("length_modify", "coefficient to modify the basic intersection / weir length", "-");
  //DECLARE_USED_PARAMETER("weir_height", "basic weir height - to be modified", "-");


END_SIMULATOR_SIGNATURE;


// =====================================================================
// =====================================================================


/**
/////////////////////////// function to calculate the new water level using the runge kutta scheme - this works perfectly according to the sample simulation!
*/
class Weir_modified : public openfluid::ware::PluggableSimulator
{
 private:

   double Q_weir(double Z1, double Z2, double Z_weir, double length, double coeff){

            /// declare necessary variables
            double Z_up{}, Z_down{};
            double H_up{}, H_down{}; // water level upstream and downstream, calculated from Z1, Z2, and Z_weir
            double Q{};
            double sign{};
            int i, n; // interval and number of time steps

            /// choose direction
            Z_up = fmax(Z1, Z2);
            Z_down = fmin (Z1, Z2);

            /// determine flow direction by sign of flow
            if (Z1 >= Z2){
                sign = 1;
                }
            else {
                sign = -1;
                }

            /// convert water elevations into heights above weir
            H_up = Z_up - Z_weir;
            H_down = Z_down - Z_weir;

            double coeff2 = 1.5*sqrt(3)*coeff; /// define the second coefficient

            if (H_up <= 0 && H_up <= H_down){ //// no flow
              Q = 0;
              }

            else if (H_up > (1.5*H_down)){ //// not submerged
                Q = coeff * sqrt(2*9.81) * pow(H_up, 1.5);
              }
            else{ //// submerged
                Q = coeff2 * sqrt(2*9.81) * H_up * sqrt(H_up - H_down);
              }

             Q = Q * sign;

            return Q;
        }


    //// Define Runge Kutta explicit scheme to calculate new water level in the SU at the end of the time step

    double water_level_new(double Z1, double Z2, double Z_weir, double length, double coeff, double elev, double A, double dts){

        double k1 = Q_weir(Z1, Z2, Z_weir, length, coeff)/A; // divided by A gives height 
        double z2 = elev + dts / 2 * k1;
        double k2 = Q_weir (Z1, z2, Z_weir, length, coeff)/A;
        double z3 = elev + dts/2 * k2;
        double k3 = Q_weir (Z1, z3, Z_weir, length, coeff)/A;
        double z4 = elev + dts/2 * k3;
        double k4 = Q_weir (Z1, z4, Z_weir, length, coeff)/A;
        double km = (k1 + 2 * k2 + 2* k3 + k4)/6;

        double  Q_new = dts*km*A; // <------------ get flow in the time interval  Question : positive / negative! 

        return  Q_new;
      }

  public:




    Weir_modified(): PluggableSimulator()
    {


    }


    // =====================================================================
    // =====================================================================


    ~Weir_modified()
    {


    }


    // =====================================================================
    // =====================================================================


    void initParams(const openfluid::ware::WareParams_t& /*Params*/)
    {


    }


    // =====================================================================
    // =====================================================================


    void prepareData()
    {


    }


    // =====================================================================
    // =====================================================================


    void checkConsistency()
    {


    }


    // =====================================================================
    // =====================================================================

   openfluid::base::SchedulingRequest initializeRun()
    {
    openfluid::core::SpatialUnit* RS;
    OPENFLUID_UNITS_ORDERED_LOOP("RS",RS)
  {
    OPENFLUID_InitializeVariable(RS,"z_new_flow",0.0); //// initialize new water level height to 0 at the beginning

  }

      return DefaultDeltaT(); /// leave at default deltaT
    }


    // =====================================================================
    // =====================================================================

  ///

    openfluid::base::SchedulingRequest runStep()
    {

      openfluid::core::SpatialUnit* RS; /// define spatial units
      openfluid::core::UnitID_t ID;


      OPENFLUID_UNITS_ORDERED_LOOP("RS",RS) // run this loop over each unit for this time step!
      {

       int  ID =  RS->getID(); // get the ID of the UNIT

       if (ID==1 || ID==28) { // if  the RS is #1 (Bassac river) or #28 (Stung)

          openfluid::core::IndexedValue water_level_input;  // define water level input variable

          OPENFLUID_GetLatestVariable(RS,"water_level_input",water_level_input);  /// get the water level input

          openfluid::core::DoubleValue water_level_Bassac = (water_level_input.value()->asDoubleValue().get()); /// get water level in the river at time step!

          openfluid::core::DoubleValue z_new_flow = water_level_Bassac; /// save the water level measured at this time step as the water level

          OPENFLUID_AppendVariable(RS,"z_new_flow",z_new_flow);

      }

       else {    // if it's any other timestep, calculate the flow into/from neighbouring units

          openfluid::core::DoubleValue length; // define length as double value type OF object
          OPENFLUID_GetAttribute(RS,"length",length); // get the length of the RS

          openfluid::core::DoubleValue width = 10.0; // define the  width of the Prek - standard as of now, later according to type

          double dtd = OPENFLUID_GetDefaultDeltaT(); /// timestep in days
          openfluid::core::TimeIndex_t CurrentTimeIndex = OPENFLUID_GetCurrentTimeIndex();
          openfluid::core::TimeIndex_t PreviousTimeIndex = OPENFLUID_GetPreviousRunTimeIndex();

           openfluid::core::DoubleValue elev;  /// elevation of RS
           OPENFLUID_GetAttribute(RS,"elev",elev);

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          /// this is where the loop starts!!
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          /////// 1 - access OFLD_TO and conn_lengt fields

          openfluid::core::StringValue conn_lengt; // define connection length as double value type OF object
          OPENFLUID_GetAttribute(RS,"conn_lengt",conn_lengt); // get the connection lengths of ther RS to the SU/RS in question

          openfluid::core::StringValue OFLD_TO; // define connection length as double value type OF object
          OPENFLUID_GetAttribute(RS,"OFLD_TO",OFLD_TO); // get the connection lengths of ther RS to the SU/RS in question

          /////// 2 - save them to vector

          // above


          std::vector<std::string> conn_lengt_list = openfluid::tools::splitString(conn_lengt, ";");

          std::vector<std::string> OFLD_TO_list = openfluid::tools::splitString(OFLD_TO, ";");



          /////// 3 - get length of vector

          int vector_length;
          vector_length = OFLD_TO_list.size();

          // add later: error if the lenght of the two vectors is different


          /////// 4 - loop through each element in the vector

          std::vector<double> values; // create an empty vector for the results of the loop

          for (int i = 0; i < vector_length ; i++) { /// for however many elements the vector has (= or < ?)


              openfluid::core::DoubleValue z_weir = 1.5; // define weir height length as double value type OF object - will later be adjusted to connection type

              std::string  ID_TO = OFLD_TO_list[i]; /// get the ID in question

              std::string  ID_TO_number = ID_TO.erase(0,3); /// erase the first three  characters in the string  (RS#)  to get just the number
             
              int neighbour_ID = stoi(ID_TO_number); /// convert string to integer
             
              std::string ID_TO_type = ID_TO.substr(0, 2); /// check whether the neighbour is a RS or a SU 

              openfluid::core::DoubleValue water_level_neighbour;

              openfluid::core::SpatialUnit* RS_neighbour = OPENFLUID_GetUnit (ID_TO_type, neighbour_ID);  // use the neighbour type + number identified here

              OPENFLUID_GetAttribute(RS_neighbour,"z_new_flow",water_level_neighbour); /// not sure if this is in the right order  

              double conn_lengt_ati = stold(conn_lengt_list[i]); /// <--------- convert to double  value

              double Z1 = water_level_neighbour ; /// absolute water level OUTSIDE the SU at the beginning of the time step, detract gauge elev if it's the Bassac!

              openfluid::core::DoubleValue Z2;
             
              OPENFLUID_GetVariable(RS,"z_new_flow",PreviousTimeIndex, Z2); /// absolute water level INSIDE the SU at the beginning of the time step

              double Z_weir = z_weir; /// absolute elevation of the weir

              double weir_length = conn_lengt_ati; /// length of the weir (= length of the intersection SU/RS, or SU/SU)

              double coeff = 0.4; /// coefficient to be calibrated => start with 0.4

              double dts = 60*dtd; /// timestep in seconds (for calculation)

              double A = length*width ; /// area of the RS

              //// Use Runge Kutta explicit scheme and weir equation (defined in separate flow formula to calculate new flow)

              double Q_new = water_level_new(Z1, Z2, Z_weir, weir_length, coeff, elev, A, dts);  // new water level as result of interaction with this unit
             
              values.push_back(Q_new); // add the new  value to the vector

          }

           /// append new variable for export!

          openfluid::core::DoubleValue z_new_avg ; // get the average of the water levels calculated in the loop above
          
         //double flow_sum = std::accumulate(values.begin(), values.end(), 0.0);
        
        double flow_sum;
        
        for (auto& n : values)
               flow_sum += n;
         
        double A = length*width ; /// area of the RS
         
        openfluid::core::DoubleValue z_new_flow = flow_sum / A;

        OPENFLUID_AppendVariable(RS,"z_new_flow",z_new_flow);
        }

        }
    return DefaultDeltaT(); /// leave at default deltaT

    }



    // =====================================================================
    // =====================================================================


    void finalizeRun()
    {


    }

};


// =====================================================================
// =====================================================================


DEFINE_SIMULATOR_CLASS(Weir_modified);


DEFINE_WARE_LINKUID(WARE_LINKUID)
