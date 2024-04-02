#ifndef TRACCCALG_H
#define TRACCCALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "DD4hep/Detector.h"
#include <DDRec/DetectorData.h>

#include "GaudiKernel/NTuple.h"
#include "DetInterface/IGeomSvc.h"
#include "syclevent.h"

namespace edm4hep {
    class SimTrackerHitCollection;
}

namespace dd4hep {
    namespace DDSegmentation {
        class BitFieldCoder;
    }
}

class TracccAlg : public GaudiAlgorithm
{

    public :

        TracccAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private :
        DataHandle<edm4hep::SimTrackerHitCollection> m_hitCol{"VXDCollection", Gaudi::DataHandle::Reader, this};
        SmartIF<IGeomSvc> m_geosvc;

        dd4hep::DDSegmentation::BitFieldCoder *m_decoder;

        uint64_t m_layer;
        uint64_t m_module;
        int m_barrelside;
        int m_x;
        int m_y;
        uint64_t gid;

        uint64_t event_num = 0;

        std::vector<my_cell> cells{};
        
        CEPCAlg* my_acts = factory();
};

uint64_t generate_gid(uint64_t layer, uint64_t module, int barrelside) {
    uint64_t acts_module = (module + 1) * 2;
    if (barrelside == 1) acts_module -= 1;
    uint64_t acts_layer = (layer + 1) * 2;
    return 0x0300000000000000 | (acts_layer << 36) | (acts_module << 8);
}

#endif  // PIXEL_H
