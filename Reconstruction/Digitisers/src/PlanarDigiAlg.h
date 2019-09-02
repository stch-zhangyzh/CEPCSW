#ifndef PLANAR_DIGI_ALG_H
#define PLANAR_DIGI_ALG_H

#include "FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include <gsl/gsl_rng.h>
#include "plcio/EventHeaderCollection.h"
#include "plcio/SimTrackerHitCollection.h"
#include "plcio/TrackerHitCollection.h"
#include "plcio/LCRelationCollection.h"

/** ======= PlanarDigiProcessor / PlanarDigiAlg ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 * The SimTrackerHits should come from a planar detector like VXD, SIT, SET or FTD.
 * 
 * WARNING: this processor depends on correctly set CellID0s and is NOT backwards compatible to 
 * SimTrackerHit output with wrong CellID0s!!!
 * 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits in u and v direction. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param SimTrkHitRelCollection The name of the TrackerHit SimTrackerHit relation collection <br>
 * (default name VTXTrackerHitRelations) <br>
 * @param ResolutionU resolution in direction of u (in mm) <br>
 * (default value 0.004) <br>
 * @param ResolutionV Resolution in direction of v (in mm) <br>
 * (default value 0.004) <br>
 * @param IsStrip whether the hits are 1 dimensional strip measurements <br>
 * (default value false)<br>
 * <br>
 * 
 */

namespace gear { class GearMgr; }

class IEventSeeder;

class PlanarDigiAlg : public GaudiAlgorithm
{
  friend class AlgFactory<PlanarDigiAlg>;
 
public:
 
  PlanarDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;
 
  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ; 
 
  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;
 
 
protected:
 
  typedef std::vector<float> FloatVec;

  int _nEvt ;

  gear::GearMgr* _GEAR;
  IEventSeeder * _SEEDER;
  gsl_rng* _rng ;
 
  // resolution in direction of u - either one per layer or one for all layers
  Gaudi::Property<FloatVec> _resU{ this, "ResolutionU", {0.0040} };
  // resolution in direction of v - either one per layer or one for all layers
  Gaudi::Property<FloatVec> _resV{ this, "ResolutionV", {0.0040} };
  // whether hits are 1D strip hits
  Gaudi::Property<bool> _isStrip{ this, "IsStrip", false };

  // Input collections
  DataHandle<plcio::EventHeaderCollection> _headerCol{"EventHeaderCol", Gaudi::DataHandle::Reader, this};
  DataHandle<plcio::SimTrackerHitCollection> _inColHdl{"VXDCollection", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<plcio::TrackerHitCollection> _outColHdl{"VXDTrackerHits", Gaudi::DataHandle::Writer, this};
  DataHandle<plcio::LCRelationCollection> _outRelColHdl{"VTXTrackerHitRelations", Gaudi::DataHandle::Reader, this};
};

#endif
