// read csv
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

// dependence
#include "TracccAlg.h"
#include "edm4hep/SimTrackerHitCollection.h"

// root
#include <TFile.h> 
#include <TTree.h>

DECLARE_COMPONENT(TracccAlg)

TracccAlg::TracccAlg(const std::string& name, ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc)
{
    declareProperty("SimTrackerHitCollection", m_hitCol);
}

StatusCode TracccAlg::initialize()
{
    info() << "begin initialize TracccAlg" << endmsg;

    m_geosvc = service<IGeomSvc>("GeomSvc");
    info() << "getDecoder VXDCollection" << endmsg;
    
    m_decoder = m_geosvc->getDecoder("VXDCollection");

    if(!m_decoder){
        error() << "Failed to get the decoder. " << endmsg;
        return StatusCode::FAILURE;
    }

    std::string detect = "CEPCSW_sim/tml_detector/detectors.csv";
    std::string digit  = "CEPCSW_sim/tml_detector/double_VXD_config.json";

    my_acts->initial(&detect, &digit);

    return GaudiAlgorithm::initialize();
}

StatusCode TracccAlg::execute()
{
    auto hits = m_hitCol.get();

    // system:5,side:-2,layer:9,module:8,sensor:8,barrelside:-2,x:-11,y:-14
    for (auto hit: *hits) {
        auto cellid = hit.getCellID();
        m_layer   = m_decoder->get(cellid, "layer");
        m_module  = m_decoder->get(cellid, "module");
        m_barrelside  = m_decoder->get(cellid, "barrelside");
        m_x       = m_decoder->get(cellid, "x");
        m_y       = m_decoder->get(cellid, "y");
        // info()  << " x: " << m_x << " y: " << m_y << endmsg;

        gid = generate_gid(m_layer, m_module, m_barrelside);
        if (m_layer < 2){
            m_x += 220;
            m_y += 1250;
        }
        else{
            m_x += 440;
            m_y += 2500;
        }

        // geometry_id,hit_id,channel0,channel1,timestamp,value,particle_id
        my_cell cell_col = {gid, m_layer, (uint32_t)m_x, (uint32_t)m_y,
                         0, 0.02, event_num};
        
        cells.push_back(cell_col);
    }

    if ((event_num % 50) == 49){
        int en = event_num / 50;
        my_acts->run(&cells, en);
        cells.clear();
    }

    event_num ++;
    return StatusCode::SUCCESS;
}

StatusCode TracccAlg::finalize()
{
    info() << "finalize TracccAlg" << endmsg;
    return GaudiAlgorithm::finalize();
}
