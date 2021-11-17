/**
  @file water.surf.transf.SU_modified.cpp
*/


/*
<sim2doc>

</sim2doc>
*/


#include <openfluid/ware/PluggableSimulator.hpp>


// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("water.surf.transf.SU_modified.id")

 // Informations
  DECLARE_NAME("Water transfer on surface units using the weir equation"); 
  DECLARE_DESCRIPTION("Combination of submerged and not submerged weir equation with 4th order Runge Kutta solution of ODE"); 
  DECLARE_VERSION(""); 
  DECLARE_STATUS(openfluid::ware::EXPERIMENTAL); 
 
  DECLARE_AUTHOR("Christina Orieschnig", "christina.orieschnig@ird.fr"); 
  
  // Required and used variables 
  
  //DECLARE_REQUIRED_VARIABLE("water.surf.H", "SU", "water height on surface of SU at the previous time step", "m"); 
  DECLARE_USED_VARIABLE("water_level_input", "SU", "water level in the Bassac", "m") 
  
  // Produced variable 
  
  DECLARE_PRODUCED_VARIABLE("z_new_flow", "SU", "new water level", "m3/s"); 
  
  
  // Required attributes for the test simulator: 
  
  DECLARE_REQUIRED_ATTRIBUTE("area", "SU", "area of SU", "m2"); 
  DECLARE_REQUIRED_ATTRIBUTE("conn_len", "SU", "intersection between the two SU - basic weir lenght", "m"); 
  
    // Required attributes for the actual simulator later on:  
  //DECLARE_REQUIRED_ATTRIBUTE("weir_height", "SU", "height of the assumed weir", "m");
  //DECLARE_REQUIRED_ATTRIBUTE("slope", "SU", "slope of SU", "m/m"); 
  //DECLARE_REQUIRED_ATTRIBUTE("SU_elev", "SU", "mean elev of SU", "m"); /// this should really be pulled from the elevation-area curve 
  //DECLARE_REQUIRED_ATTRIBUTE("nmanning", "SU", "manning roughness coefficient", "s/m(-1/3)");
  
  
  // Simulator parameters (for later on) 
  
  //DECLARE_USED_PARAMETER("coeff", "basic weir coefficient", "-");
  //DECLARE_USED_PARAMETER("length_modify", "coefficient to modify the basic intersection / weir length", "-"); 
  //DECLARE_USED_PARAMETER("weir_height", "basic weir height - to be modified", "-");
  

END_SIMULATOR_SIGNATURE


// =====================================================================
// =====================================================================


/**

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

        double k1 = Q_weir(Z1, Z2, Z_weir, length, coeff)/A/364; // divided by A gives height but there's an issue with the units of A, I think! 
        double z2 = elev + dts / 2 * k1;  
        double k2 = Q_weir (Z1, z2, Z_weir, length, coeff)/A/864;
        double z3 = elev + dts/2 * k2; 
        double k3 = Q_weir (Z1, z3, Z_weir, length, coeff)/A/864;
        double z4 = elev + dts/2 * k3; 
        double k4 = Q_weir (Z1, z4, Z_weir, length, coeff)/A/864;
        double km = (k1 + 2 * k2 + 2* k3 + k4)/6;  

        double  z_new = elev + dts*km; // 

        return z_new;
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
    openfluid::core::SpatialUnit* SU; 
    OPENFLUID_UNITS_ORDERED_LOOP("SU",SU) 
  { 
    OPENFLUID_InitializeVariable(SU,"z_new_flow",0.0); //// initialize new water level height to 0 at the beginning 
    
  }
      
      return DefaultDeltaT(); /// leave at default deltaT
    }


    // =====================================================================
    // =====================================================================
  
  
    openfluid::base::SchedulingRequest runStep()
    {

      openfluid::core::SpatialUnit* SU; /// define spatial units 

      OPENFLUID_UNITS_ORDERED_LOOP("SU",SU) // run this loop over each unit for this time step! 
      {

        openfluid::core::DoubleValue area; // define area as double value type OF object
        OPENFLUID_GetAttribute(SU,"area",area); // get the area of the SU
        
        openfluid::core::DoubleValue conn_len; // define connection length as double value type OF object
        OPENFLUID_GetAttribute(SU,"conn_len",conn_len); // get the area of the SU
        
        openfluid::core::DoubleValue z_weir = 1.5; // define weir height length as double value type OF object
        //OPENFLUID_GetAttribute(SU,"z_weir",z_weirlen); // get the area of the SU (for later! => according to type )
        
       
        openfluid::core::IndexedValue water_level_input; 
        
        OPENFLUID_GetLatestVariable(SU,"water_level_input",water_level_input); 

        
        openfluid::core::DoubleValue water_level_Bassac = (water_level_input.value()->asDoubleValue().get()); /// get water level in the river at time step! 

        
        // 
        
        double dtd = OPENFLUID_GetDefaultDeltaT(); /// timestep in days 
        openfluid::core::TimeIndex_t CurrentTimeIndex = OPENFLUID_GetCurrentTimeIndex();
        openfluid::core::TimeIndex_t PreviousTimeIndex = OPENFLUID_GetPreviousRunTimeIndex();
        
        
        double Z1 = water_level_Bassac-1; /// absolute water level OUTSIDE the SU at the beginning of the time step, detract gauge elev
        openfluid::core::DoubleValue Z2;
        OPENFLUID_GetVariable(SU,"z_new_flow",PreviousTimeIndex, Z2); /// absolute water level INSIDE the SU at the beginning of the time step 
        double Z_weir = z_weir; /// absolute elevation of the weir 
        double length = conn_len; /// length of the weir (= length of the intersection SU/RS, or SU/SU) 
        double coeff = 0.6; /// coefficient to be calibrated => start with 0.4 
        
        double dts = 86400*dtd; /// timestep in seconds (for calculation)

        double A = area ; /// area of the SU 

        double elev = 0; /// elevation of the SU (get from elevatin curve later on) 


        //// Use Runge Kutta explicit scheme and weir equation (defined in separate flow formula to calculate new flow) 
        
        double z_new = water_level_new(Z1, Z2, Z_weir, length, coeff, elev, A, dts); 
        
        
        /// append new variable for export! 

        openfluid::core::DoubleValue z_new_flow = z_new;
        
        OPENFLUID_AppendVariable(SU,"z_new_flow",z_new_flow); 

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


