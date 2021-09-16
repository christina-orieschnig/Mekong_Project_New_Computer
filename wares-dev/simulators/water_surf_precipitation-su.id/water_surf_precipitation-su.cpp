/**
  @file water_surf_precipitation-su.cpp
*/


/*
<sim2doc>

</sim2doc>
*/


#include <openfluid/ware/PluggableSimulator.hpp>


// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("water_surf_precipitation-su.id")

  // Informations
  DECLARE_NAME("Conversion of precipitation data to component of water height on surface unig")
  DECLARE_DESCRIPTION("Takes input of rainfall data, and calculates rainfall height during time step for each SU")
  //DECLARE_VERSION("")
  DECLARE_STATUS(openfluid::ware::EXPERIMENTAL)
  
  
  DECLARE_REQUIRED_VARIABLE("water_rain_rate", "SU", "precipitation rate at time step x", "mm/d"); 
  
  DECLARE_REQUIRED_ATTRIBUTE("area", "SU", "area of SU", "m2"); 
  
  //DECLARE_USED_PARAMETER("dtd", "time step in days"); 
  
  DECLARE_PRODUCED_VARIABLE("rainfall_m", "SU", "rainfall height on Su at time step x", "m"); 
  DECLARE_PRODUCED_VARIABLE("rainfall_vol", "SU", "rainfall volume on Su at time step x", "m3"); 

END_SIMULATOR_SIGNATURE


// =====================================================================
// =====================================================================


/**

*/
class Precipitation : public openfluid::ware::PluggableSimulator
{
  private:

  
  public:
    
    


  
    Precipitation(): PluggableSimulator()
    {
  
  
    }
  
  
    // =====================================================================
    // =====================================================================
  
  
    ~Precipitation()
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
    OPENFLUID_InitializeVariable(SU,"rainfall_m",0.0);
    OPENFLUID_InitializeVariable(SU,"rainfall_vol",0.0);
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

        openfluid::core::IndexedValue water_rain_rate; 
        
        OPENFLUID_GetLatestVariable(SU,"water_rain_rate",water_rain_rate); 

        
        openfluid::core::DoubleValue rainfall_mm = (water_rain_rate.value()->asDoubleValue().get());


        // compute rainfall height in m and rainfall volume per SU 
        openfluid::core::DoubleValue  rainfall_m = rainfall_mm/1000;// daily rainfall in m per m2
        openfluid::core::DoubleValue  rainfall_vol = rainfall_m*area; // volume of rainfall over SU 

        OPENFLUID_AppendVariable(SU,"rainfall_m",rainfall_m); // produce for export 
        OPENFLUID_AppendVariable(SU,"rainfall_vol",rainfall_vol); // produce for export


      return DefaultDeltaT(); /// leave at default deltaT
    }
    
    }


    // =====================================================================
    // =====================================================================
  
  
  /// what this bit has to do: get the area of each SU, get the water_rain_rate for this time step, and calculate the rainfall height in m and the volume for the SU, what it does currently: nothing but break 
  
  
    void finalizeRun()
    {
  
  
    }

};





// =====================================================================
// =====================================================================


DEFINE_SIMULATOR_CLASS(Precipitation); 


DEFINE_WARE_LINKUID(WARE_LINKUID);

