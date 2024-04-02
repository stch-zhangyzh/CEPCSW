#!/usr/bin/env python

from Gaudi.Configuration import *
import os
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc",
                 input="CRD-oi-v0j-Sim00.root"
)

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections=[
    "VXDCollection",
    ])

##############################################################################
# Geometry Svc
##############################################################################

geometry_option = "CRD_o1_v01/CRD_o1_v01.xml"

if not os.getenv("DETCRDROOT"):
    print("Can't find the geometry. Please setup envvar DETCRDROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCRDROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geosvc = GeomSvc("GeomSvc")
geosvc.compact = geometry_path


##############################################################################
# TracccAlg
##############################################################################

from Configurables import TracccAlg
alg = TracccAlg("TracccAlg")
alg.SimTrackerHitCollection.Path = "VXDCollection"

##############################################################################
# TracccAlg
##############################################################################
# from Configurables import PodioOutput
# out = PodioOutput("out")
# out.filename = "sim_RawTimeSeries.root"
# out.outputCommands = ["keep *"]

# VXDCollection.cellID
# system:5,side:-2,layer:9,module:8,sensor:8,barrelside:-2,x:-10,y:-10

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [podioinput, alg],
                EvtSel = 'NONE',
                EvtMax = 5000,
                ExtSvc = [dsvc, geosvc],
                # OutputLevel=DEBUG
)
