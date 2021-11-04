/**
  @file observe_flow_simulator.cpp
*/


/*
<sim2doc>

</sim2doc>
*/


#include <openfluid/ware/PluggableSimulator.hpp>


// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("observe_flow_simulator.id")

  // Informations
  DECLARE_NAME("observe_flow_simulator")
  DECLARE_DESCRIPTION("for daily instead of minute output of flow simulation")
  DECLARE_VERSION("")
  DECLARE_STATUS(openfluid::ware::EXPERIMENTAL)
  DECLARE_AUTHOR("Christina Orieschnig", "christina.orieschnig@ird.fr");

  DECLARE_USED_VARIABLE("z_new_flow_SU", "SU", "water level in the SU", "m") ;
  DECLARE_PRODUCED_VARIABLE("z_new_flow_SU_daily", "SU", "daily water level in the SU", "m") ;

END_SIMULATOR_SIGNATURE


// =====================================================================
// =====================================================================


/**

*/
class observe_flow_simulator : public openfluid::ware::PluggableSimulator
{
  private:


  public:


    observe_flow_simulator(): PluggableSimulator()
    {


    }


    // =====================================================================
    // =====================================================================


    ~observe_flow_simulator()
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

    ////////// initialize the variable to 0 for each SU
    openfluid::base::SchedulingRequest initializeRun()
    {

       openfluid::core::SpatialUnit* SU;
    OPENFLUID_UNITS_ORDERED_LOOP("SU",SU)
      {
        
 
        OPENFLUID_InitializeVariable(SU,"z_new_flow_SU_daily",0.0); //// initialize new water level height to 0 at the beginning

      }

      return MultipliedDefaultDeltaT(86400);
    }


    // =====================================================================
    // =====================================================================

    /// run loop over the SUs each day 
    
    openfluid::base::SchedulingRequest runStep()
      {
       //// define SUs
        openfluid::core::SpatialUnit* SU;

       //// spatial loop over SUs
        OPENFLUID_UNITS_ORDERED_LOOP("SU",SU) 
          {
            // get time index 
            openfluid::core::TimeIndex_t CurrentTimeIndex = OPENFLUID_GetCurrentTimeIndex();

            // OFLD variable for flow value 
            openfluid::core::DoubleValue flow_value; /// previous water level

            // get OFLD variable from flow simulator 
            OPENFLUID_GetVariable(SU,"z_new_flow_SU",CurrentTimeIndex, flow_value); /// absolute water level INSIDE the SU at the beginning of the time step

            // save it to regular double variable
            double z_new_flow = flow_value.get();

            // declare daily output variable 
            openfluid::core::DoubleValue z_new_flow_SU_daily;

            // save the value of the flow simulator at the daily time step to the daily output variable 
            z_new_flow_SU_daily.set(z_new_flow);

            // append the daily output variable 
            OPENFLUID_AppendVariable(SU,"z_new_flow_SU_daily",z_new_flow_SU_daily);
          }

       // repeat every day, 24 h, 86400 seconds 
       return MultipliedDefaultDeltaT(86400);
      }


    // =====================================================================
    // =====================================================================


    void finalizeRun()
    {


    }

};


// =====================================================================
// =====================================================================


DEFINE_SIMULATOR_CLASS(observe_flow_simulator);


DEFINE_WARE_LINKUID(WARE_LINKUID)
