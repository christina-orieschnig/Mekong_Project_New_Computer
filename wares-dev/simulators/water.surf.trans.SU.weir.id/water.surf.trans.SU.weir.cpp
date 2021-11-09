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
#include <iostream>

// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("water.surf.trans.SU.weir.id");

 // Informations
  DECLARE_NAME("Water transfer on surface units using the weir equation");
  DECLARE_DESCRIPTION("Combination of submerged and not submerged weir equation with 4th order Runge Kutta solution of ODE");
  DECLARE_VERSION("");
  DECLARE_STATUS(openfluid::ware::EXPERIMENTAL);

  DECLARE_AUTHOR("Christina Orieschnig", "christina.orieschnig@ird.fr");

  // Required and used variables

  //DECLARE_REQUIRED_VARIABLE("water.surf.H", "SU", "water height on surface of SU at the previous time step", "m");
  DECLARE_USED_VARIABLE("z_new_flow_SU", "SU", "water level in the SU", "m") ;
  DECLARE_USED_VARIABLE("z_new_flow_RS", "RS", "water level in the RS", "m") ;

  // Produced variable

  DECLARE_PRODUCED_VARIABLE("z_new_flow_SU", "SU", "water level in the SU", "m");


  // Required ys for the test simulator:

  DECLARE_REQUIRED_ATTRIBUTE("area", "SU", "area of SU", "m");
  DECLARE_REQUIRED_ATTRIBUTE("conn_lengt", "SU", "intersection between the two SU and other RS and SU - basic weir length", "m");
  DECLARE_REQUIRED_ATTRIBUTE("OFLD_TO", "SU", "list of connections between the SU and other RS and SU - basic weir length", "m");
  DECLARE_REQUIRED_ATTRIBUTE("elev", "SU", "mean elev of SU", "m"); /// this should really be pulled from the elevation-area curve

  // Simulator parameters (for later on)

  DECLARE_USED_PARAMETER("coefficient", "basic weir coefficient", "-");
  //DECLARE_USED_PARAMETER("length_modify", "coefficient to modify the basic intersection / weir length", "-");
  DECLARE_USED_PARAMETER("weir_height_rel", "basic weir height - to be modified", "-");


END_SIMULATOR_SIGNATURE;


// =====================================================================
// =====================================================================


/**
/////////////////////////// function to calculate the new water level using the runge kutta scheme - this works perfectly according to the sample simulation!
*/
class Weir_modified : public openfluid::ware::PluggableSimulator
{
 private:
   double coefficient;

   double weir_height_rel;

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

    double water_level_new(double Z1, double Z2, double Z_weir, double length, double coeff, double A, double dts){

        /*

        double k1 = Q_weir(Z1, Z2, Z_weir, length, coeff)///A; // divided by A gives height 
        double z2 = Z2 + dts / 2 * k1;
        double k2 = Q_weir (Z1, z2, Z_weir, length, coeff)/A;
        double z3 = Z2 + dts/2 * k2;
        double k3 = Q_weir (Z1, z3, Z_weir, length, coeff)/A;
        double z4 = Z2 + dts/2 * k3;
        double k4 = Q_weir (Z1, z4, Z_weir, length, coeff)/A;
        double km = (k1 + 2 * k2 + 2* k3 + k4)/6;

        double  Q_new = dts*km*A; // <------------ get flow in the time interval  Question : positive / negative! 
        */

        double Q_new = Q_weir(Z1, Z2, Z_weir, length, coeff);
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


    void initParams(const openfluid::ware::WareParams_t& Params)
    {
      openfluid::core::DoubleValue coefficient;
      openfluid::core::DoubleValue weir_height_rel;

      coefficient = 0.0;
      if (!OPENFLUID_GetSimulatorParameter(Params,"coefficient",coefficient))
        OPENFLUID_RaiseError("coefficinet ???");

      weir_height_rel = 0.0;
      if (!OPENFLUID_GetSimulatorParameter(Params,"weir_height_rel",weir_height_rel))
        OPENFLUID_RaiseError("weir_height_rel ???");


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
    openfluid::core::SpatialUnit* SU;
    OPENFLUID_UNITS_ORDERED_LOOP("SU",SU)
  {
    OPENFLUID_InitializeVariable(SU,"z_new_flow_SU",0.0); //// initialize new water level height to 0 at the beginning

  }

      return MultipliedDefaultDeltaT(60); /// leave at default deltaT
    }


    // =====================================================================
    // =====================================================================

  ///

    openfluid::base::SchedulingRequest runStep()
    {

      openfluid::core::SpatialUnit* SU; /// define spatial units
      openfluid::core::UnitID_t ID;


      OPENFLUID_UNITS_ORDERED_LOOP("SU",SU) // run this loop over each unit for this time step!
      {

        int  ID =  SU->getID(); // get the ID of the UNIT

        openfluid::core::DoubleValue elevation;  /// elevation of RS
        OPENFLUID_GetAttribute(SU,"elev",elevation);
        double elev_at_unit = elevation.get();

        openfluid::core::DoubleValue area; // define length as double value type OF object
        OPENFLUID_GetAttribute(SU,"area",area); // get the area of the SU
        double area_at_unit = area.get();

        double dtd = OPENFLUID_GetDefaultDeltaT(); /// timestep in days
        openfluid::core::TimeIndex_t CurrentTimeIndex = OPENFLUID_GetCurrentTimeIndex();
        openfluid::core::TimeIndex_t PreviousTimeIndex = CurrentTimeIndex-dtd*60;

        openfluid::core::DoubleValue previous_value; /// previous water level
        OPENFLUID_GetVariable(SU,"z_new_flow_SU",PreviousTimeIndex, previous_value); /// absolute water level INSIDE the SU at the beginning of the time step
        double previous_water_level_at_unit = previous_value.get()+elev_at_unit;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// this is where the loop starts!!
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /////// 1 - access OFLD_TO and conn_lengt fields

        openfluid::core::StringValue conn_lengt; // define connection length as double value type OF object
        
        try  {
         OPENFLUID_GetAttribute(SU,"conn_lengt",conn_lengt);
        }

        catch (...) {
           openfluid::core::IntegerValue temp_length;
           OPENFLUID_GetAttribute(SU,"conn_lengt",temp_length);
          int length = temp_length.get();
          std::string length_str = std::to_string(length) ;
          conn_lengt.set(length_str);
        } // get the connection lengths of ther SU to the SU/RS in question

        openfluid::core::StringValue OFLD_TO; // define connection length as double value type OF object
        OPENFLUID_GetAttribute(SU,"OFLD_TO",OFLD_TO); // get the connections lengths of the SU to the SU/RS in question

        /////// 2 - save them to vector

        std::vector<std::string> conn_lengt_list = openfluid::tools::splitString(conn_lengt, ";");

        std::vector<std::string> OFLD_TO_list = openfluid::tools::splitString(OFLD_TO, ";");

        /////// 3 - get length of vector

        int vector_length;
        vector_length = OFLD_TO_list.size();

        // add later: error if the lenght of the two vectors is different


        /////// 4 - loop through each element in the vector

        std::vector<double> values; // create an empty vector for the results of the loop

        for (int i = 0; i < vector_length ; i++) { /// for however many elements the vector has (= or < ?)

            std::string  ID_TO = OFLD_TO_list[i]; /// get the ID in question

            std::string  ID_TO_number = ID_TO.substr(3,3); /// erase the first three  characters in the string  (SU#)  to get just the number

            int neighbour_ID = stoi(ID_TO_number); /// convert string to integer

            std::string ID_TO_type = ID_TO.substr(0, 2); /// check whether the neighbour is a RS or a SU

            std::string neighbour_variable;

            if (ID_TO_type == "SU"){
             neighbour_variable = "z_new_flow_SU";
           }
            else {
             neighbour_variable = "z_new_flow_RS";
           }

           /// define the neighbouring unit that we're dealing with
            openfluid::core::SpatialUnit* SU_neighbour = OPENFLUID_GetUnit (ID_TO_type, neighbour_ID);  // use the neighbour type + number identified here

            /// get water level of neighbour 
            openfluid::core::DoubleValue  water_level_neighbour; /// define variable for previous water level in neighbouring element
            OPENFLUID_GetVariable(SU_neighbour,neighbour_variable, PreviousTimeIndex,water_level_neighbour);

            /// get elevation of neighbour 
            openfluid::core::DoubleValue elevation_neighbour;  /// get elevation of neighbour unit as OFLD double
            OPENFLUID_GetAttribute(SU_neighbour,"elev",elevation_neighbour); /// get attribute 
            double elev_at_neighbour = elevation_neighbour.get(); /// convert to regular double 

            /// Get ABSOLUTE water level of neighbour
            double  previous_water_level_at_neighbour = water_level_neighbour.get()+elev_at_neighbour; /// absolute water level OUTSIDE the SU at the beginning of the time step, detract gauge elev if it's the Bassac!

            /// Get Connection Length 
            double conn_lengt_ati = stold(conn_lengt_list[i]); /// <--------- convert to double  value

            /// Get ABSOLUTE elevation of the weir 
            double larger_elev;

            if (elev_at_unit >= elev_at_neighbour){
              larger_elev = elev_at_unit;
            }

            else {
              larger_elev = elev_at_neighbour;
            }
            
            double z_weir = weir_height_rel; // define weir height length as double value type OF object - will later be adjusted to connection type
            double Z_weir = z_weir+larger_elev; /// absolute elevation of the weir

            double weir_length = conn_lengt_ati; /// length of the weir (= length of the intersection SU/RS, or SU/SU)

            double coeff = coefficient; /// coefficient to be calibrated => start with 0.4

            double dts = 60*dtd; /// timestep in seconds (for calculation)
           
            double Q_new;

            if (Z_weir > previous_water_level_at_unit && Z_weir > previous_water_level_at_neighbour){
              Q_new = 0;
            } 
            else {
            //// Use Runge Kutta explicit scheme and weir equation (defined in separate flow formula to calculate new flow)
             Q_new = water_level_new(previous_water_level_at_unit, previous_water_level_at_neighbour, Z_weir, weir_length, coeff, area_at_unit, dts);  // new water level as result of interaction with this unit
            }
/*
            OPENFLUID_LogAndDisplayInfo("ID");
            OPENFLUID_LogAndDisplayInfo(ID);
            OPENFLUID_LogAndDisplayInfo("ID_TO_type");
            OPENFLUID_LogAndDisplayInfo(neighbour_variable);
            OPENFLUID_LogAndDisplayInfo("neighbour_ID");
            OPENFLUID_LogAndDisplayInfo(neighbour_ID);
            OPENFLUID_LogAndDisplayInfo("weir height");
            OPENFLUID_LogAndDisplayInfo(Z_weir);
            OPENFLUID_LogAndDisplayInfo("elevation of neighbour");
            OPENFLUID_LogAndDisplayInfo(elev_at_neighbour);
            OPENFLUID_LogAndDisplayInfo("water level neighbour");
            OPENFLUID_LogAndDisplayInfo(previous_water_level_at_neighbour);
            OPENFLUID_LogAndDisplayInfo("water level at unit");
            OPENFLUID_LogAndDisplayInfo(previous_water_level_at_unit);
            OPENFLUID_LogAndDisplayInfo("flow");
            OPENFLUID_LogAndDisplayInfo(Q_new);
              */
            values.push_back(Q_new); // add the new  value to the vector

        }

         /// append new variable for export!

        //double flow_sum = std::accumulate(values.begin(), values.end(), 0.0);

        double flow_sum;

        for (auto& n : values)
               flow_sum += n;

        double z_new_flow = flow_sum/ area_at_unit+ previous_water_level_at_unit;

        double z_new_flow_final;

        if (z_new_flow < elev_at_unit) {
          z_new_flow_final = 0;
        }
        else {
          z_new_flow_final = z_new_flow-elev_at_unit;
        }

        openfluid::core::DoubleValue z_new_flow_SU;

        z_new_flow_SU.set(z_new_flow_final);

        OPENFLUID_AppendVariable(SU,"z_new_flow_SU",z_new_flow_SU);


        }
    return MultipliedDefaultDeltaT(60); /// leave at default deltaT

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
