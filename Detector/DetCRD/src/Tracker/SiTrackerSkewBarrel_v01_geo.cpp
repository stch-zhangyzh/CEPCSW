//====================================================================
//  cepcvxdgeo - CEPC vertex detector models in DD4hep 
//--------------------------------------------------------------------
//  Hao Zeng, IHEP
//  email: zenghao@ihep.ac.cn
//  $Id$
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"
#include "XMLHandlerDB.h"
#include <cmath>

using namespace std;

using dd4hep::Box;
using dd4hep::DetElement;
using dd4hep::Material;
using dd4hep::Position;
using dd4hep::RotationY;
using dd4hep::RotationZYX;
using dd4hep::Transform3D;
using dd4hep::Rotation3D;
using dd4hep::Volume;
using dd4hep::_toString;
using dd4hep::rec::volSurfaceList;
using dd4hep::rec::ZPlanarData;
using dd4hep::mm;

/** helper struct */
struct VXD_Layer {
  int     n_ladders;
  int     n_sensors_per_ladder;
  double  sensor_length;
  double  half_z;
  double  sensitive_inner_radius ;
  double  support_inner_radius ;
  double  ladder_width ;
  double  ladder_dphi ;
};    

//std::vector<VXD_Layer> _VXD_Layers;
  
// /** helper struct */
// struct extended_reconstruction_parameters {
//   double sensor_length_mm;
//   double strip_width_mm;
//   double strip_length_mm;
//   double strip_pitch_mm;
//   double strip_angle_deg;
// };

//extended_reconstruction_parameters _e_r_p;


/** Construction of the VXD detector, ported from Mokka driver SIT_Simple_Pixel.cc
 *
 *  Mokka History:
 *  Feb 7th 2011, Steve Aplin - original version
 *  F.Gaede, DESY, Jan 2014   - dd4hep SIT pixel

 *  @author Hao Zeng, IHEP, July 2021
 */
static dd4hep::Ref_t create_element(dd4hep::Detector& theDetector, xml_h e, dd4hep::SensitiveDetector sens)  {

  xml_det_t  x_det    = e;
  Material   air      = theDetector.air();
  int        det_id   = x_det.id();
  string     name     = x_det.nameStr();
  DetElement vxd(name, det_id);

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, e, vxd);
  dd4hep::xml::setDetectorTypeFlag(e, vxd) ;
  if(theDetector.buildType()==dd4hep::BUILD_ENVELOPE) return vxd;
  envelope.setVisAttributes(theDetector.visAttributes("SeeThrough"));

  sens.setType("tracker");
  std::cout << " ** building SiTrackerSkewBarrel_v01 ..." << std::endl ;

  dd4hep::rec::ZPlanarData* zPlanarData = new dd4hep::rec::ZPlanarData;

   // fetch the global parameters
   XMLHandlerDB db0 = XMLHandlerDB(  x_det.child( _Unicode( global ) ) ) ;
   //Material support_mat         = theDetector.material(db0->fetchString("support_mat"));
   
   //fetch the display parameters
   db0 = XMLHandlerDB(  x_det.child( _Unicode( display ) ) ) ;
   std::string ladderVis      = db0->fetchString("ladder");
   std::string supportVis     = db0->fetchString("support");
   std::string flexVis        = db0->fetchString("flex");
   std::string sensEnvVis     = db0->fetchString("sens_env");
   std::string sensVis        = db0->fetchString("sens");
   std::string deadsensVis    = db0->fetchString("deadsensor");
   std::string deadwireVis    = db0->fetchString("deadwire");


 for(xml_coll_t layer_i(x_det,_U(layer)); layer_i; ++layer_i){
   xml_comp_t x_layer(layer_i);
   XMLHandlerDB db = XMLHandlerDB( x_layer );
   
   dd4hep::PlacedVolume pv;
   int layer_id                 = db->fetchInt("layer_id");

   std::cout << "layer_id: " << layer_id << endl;

   double sensitive_radius      = db->fetchDouble("ladder_radius");
   int n_sensors_per_ladder     = db->fetchInt("n_sensors_per_side");
   int n_ladders                = db->fetchInt("n_ladders") ;
   double ladder_offset         = db->fetchDouble("ladder_offset");
   double ladder_radius         = sqrt(ladder_offset*ladder_offset + sensitive_radius*sensitive_radius); 
   double ladder_phi0           = -atan(ladder_offset/sensitive_radius);
   std::cout << "ladder_radius: " << ladder_radius/mm <<" mm" << endl;

   std::cout << "sensitive_radius: " << sensitive_radius/mm << " mm" << endl;
   std::cout << "n_sensors_per_ladder: " << n_sensors_per_ladder << endl;

   std::string layerName = dd4hep::_toString( layer_id , "layer_%d"  );
   dd4hep::Assembly layer_assembly( layerName ) ;
   pv = envelope.placeVolume( layer_assembly ) ;
   dd4hep::DetElement layerDE( vxd , layerName  , x_det.id() );
   layerDE.setPlacement( pv ) ;
   
   const double ladder_dphi = ( dd4hep::twopi / n_ladders ) ;
   std::cout << "ladder_dphi: " << ladder_dphi << endl;

   //fetch the ladder parameters
   xml_comp_t x_ladder(x_layer.child(_Unicode(ladder)));
   db = XMLHandlerDB(x_ladder);
   
   //fetch the ladder support parameters
   db = XMLHandlerDB(x_ladder.child(_Unicode(ladderSupport)));
   double support_length        = db->fetchDouble("length");
   double support_thickness     = db->fetchDouble("thickness");
   double support_height        = db->fetchDouble("height");
   double support_width         = db->fetchDouble("width");
   Material support_mat         = theDetector.material(db->fetchString("mat"));
   std::cout << "support_length: " << support_length/mm << " mm" << endl;
   std::cout << "support_thickness: " << support_thickness/mm << " mm" << endl;
   std::cout << "support_width: " << support_width/mm << " mm" << endl;
   
   //fetch the flex parameters
   double flex_thickness(0);
   double flex_width(0);
   double flex_length(0);
   xml_comp_t x_flex(x_ladder.child(_Unicode(flex)));
   for(xml_coll_t flex_i(x_flex,_U(layer)); flex_i; ++flex_i){
     xml_comp_t x_flex_layer(flex_i);
     db = XMLHandlerDB(x_flex_layer);
     double x_flex_layer_thickness = db->fetchDouble("thickness");
     double x_flex_layer_width = db->fetchDouble("width");
     double x_flex_layer_length = db->fetchDouble("length");
     flex_thickness += x_flex_layer_thickness;
     if (x_flex_layer_width > flex_width) flex_width = x_flex_layer_width;
     if (x_flex_layer_length > flex_length) flex_length = x_flex_layer_length;
     std::cout << "x_flex_layer_thickness: " << x_flex_layer_thickness/mm << " mm" << endl;
   }
   std::cout << "flex_thickness: " << flex_thickness/mm << " mm" << endl;
   std::cout << "flex_width: " << flex_width/mm << " mm" << endl;
   std::cout << "flex_length: " << flex_length/mm << " mm" << endl;
   
   //fetch the sensor parameters
   db = XMLHandlerDB(x_ladder.child(_Unicode(sensor)));
   int n_sensors_per_side                  = db->fetchInt("n_sensors");
   double dead_gap                         = db->fetchDouble("gap");
   double sensor_thickness                 = db->fetchDouble("thickness");
   double sensor_active_len                = db->fetchDouble("active_length");
   double sensor_active_width              = db->fetchDouble("active_width");
   double sensor_dead_width                = db->fetchDouble("dead_width");
   double sensor_deadwire_length           = db->fetchDouble("deadwire_length");
   double sensor_deadwire_width            = db->fetchDouble("deadwire_width");
   double sensor_deadwire_thickness        = db->fetchDouble("deadwire_thickness");
   Material sensor_mat                     = theDetector.material(db->fetchString("sensor_mat"));
   Material sensor_deadwire_mat            = theDetector.material(db->fetchString("deadwire_mat"));

   std::cout << "n_sensors_per_side: " << n_sensors_per_side << endl;
   std::cout << "dead_gap: " << dead_gap/mm << " mm" << endl;
   std::cout << "sensor_thickness: " << sensor_thickness/mm << " mm" << endl;
   std::cout << "sensor_active_len: " << sensor_active_len/mm << " mm" << endl;
   std::cout << "sensor_active_width: " << sensor_active_width/mm << " mm" << endl;
   std::cout << "sensor_dead_width: " << sensor_dead_width/mm << " mm" << endl;
  
  //create ladder logical volume
  Box LadderSolid((support_height+2*sensor_thickness+2*flex_thickness)/2.0, 
                    support_width / 2.0, support_length / 2.0);
  Volume LadderLogical(name + dd4hep::_toString( layer_id, "_LadderLogical_%02d"),
                      LadderSolid, air);
  // create flex envelope logical volume
  Box FlexEnvelopeSolid(flex_thickness / 2.0, flex_width / 2.0, flex_length / 2.0);
  Volume FlexEnvelopeLogical(name + dd4hep::_toString( layer_id, "_FlexEnvelopeLogical_%02d"), FlexEnvelopeSolid, air);
  FlexEnvelopeLogical.setVisAttributes(theDetector.visAttributes("SeeThrough"));
  //vxd.setVisAttributes(theDetector, flexVis, FlexEnvelopeLogical);

  //create the flex layers inside the flex envelope
  double flex_height(0); 
  int index = 0;
  for(xml_coll_t flex_i(x_flex,_U(layer)); flex_i; ++flex_i){
    xml_comp_t x_flex_layer(flex_i);
    db = XMLHandlerDB(x_flex_layer);
    double x_flex_layer_thickness = db->fetchDouble("thickness");
    double x_flex_layer_width = db->fetchDouble("width");
    double x_flex_layer_length = db->fetchDouble("length");
    Material x_flex_layer_mat = theDetector.material(db->fetchString("mat"));
    flex_height += x_flex_layer_thickness;
    Box FlexLayerSolid(x_flex_layer_thickness/2.0, x_flex_layer_width/2.0, x_flex_layer_length/2.0);
    Volume FlexLayerLogical(name + dd4hep::_toString( layer_id, "_FlexLayerLogical_%02d") + dd4hep::_toString( index, "index_%02d"), FlexLayerSolid, x_flex_layer_mat);
    FlexLayerLogical.setVisAttributes(theDetector.visAttributes(flexVis));
    pv = FlexEnvelopeLogical.placeVolume(FlexLayerLogical, Position(flex_height/2.0, 0., 0.));
    std::cout << "flex thickness = " << x_flex_layer_thickness << std::endl;
    std::cout << "flex width = " << x_flex_layer_width << std::endl;
    std::cout << "flex length = " << x_flex_layer_length << std::endl;
    // std::cout << "flex material: " << x_flex_layer_mat << std::endl;
    index++;
  }

  //place the flex envelope inside the ladder envelope
  pv = LadderLogical.placeVolume(FlexEnvelopeLogical, Position((support_height + flex_thickness) / 2.0, 0., 0.)); //top side
  //define the transformation3D(only need a combination of translation and rotation)
  Transform3D tran_mirro(RotationZYX(0., dd4hep::twopi/2.0, 0.), Position(-(support_height + flex_thickness) / 2.0, 0., 0.));
  pv = LadderLogical.placeVolume(FlexEnvelopeLogical, tran_mirro); //bottom side
  
  //create sensor envelope logical volume
  Box SensorTopEnvelopeSolid(sensor_thickness / 2.0, support_width / 2.0, support_length / 2.0);
  Volume SensorTopEnvelopeLogical(name + dd4hep::_toString( layer_id, "_SensorEnvelopeLogical_%02d"), SensorTopEnvelopeSolid, air);
  Box SensorBottomEnvelopeSolid(sensor_thickness / 2.0, support_width / 2.0, support_length / 2.0);
  Volume SensorBottomEnvelopeLogical(name + dd4hep::_toString( layer_id, "_SensorEnvelopeLogical_%02d"), SensorBottomEnvelopeSolid, air);
  SensorTopEnvelopeLogical.setVisAttributes(theDetector.visAttributes(sensEnvVis));

  //create sensor logical volume
  Box SensorSolid(sensor_thickness / 2.0, sensor_active_width / 2.0, sensor_active_len / 2.0);
  Volume SensorLogical(name + dd4hep::_toString( layer_id, "_SensorLogical_%02d"), SensorSolid, sensor_mat);
  SensorLogical.setSensitiveDetector(sens);
  //vxd.setVisAttributes(theDetector, deadsensVis, SensorDeadLogical);
  SensorLogical.setVisAttributes(theDetector.visAttributes(sensVis));

  //create dead sensor logical volume
  Box SensorDeadSolid(sensor_thickness / 2.0, sensor_dead_width / 2.0, sensor_active_len / 2.0);
  Volume SensorDeadLogical(name + dd4hep::_toString( layer_id, "_SensorDeadLogical_%02d"), SensorDeadSolid, sensor_mat);
  SensorDeadLogical.setVisAttributes(theDetector.visAttributes(deadsensVis));

  //create dead wire logical volume
  Box SensorDeadWireSolid(sensor_deadwire_thickness / 2.0, sensor_deadwire_width / 2.0, sensor_deadwire_length / 2.0);
  Volume SensorDeadWireLogical(name + dd4hep::_toString( layer_id, "_SensorDeadWireLogical_%02d"), SensorDeadWireSolid, sensor_deadwire_mat);
  SensorDeadWireLogical.setVisAttributes(theDetector.visAttributes(deadwireVis));
  
  //place the dead wire in the sensor envelope
  // pv = SensorTopEnvelopeLogical.placeVolume(SensorDeadWireLogical, Position(0.0, (sensor_active_width-support_width/2.0) + sensor_dead_width/2.0 + sensor_deadwire_width/2.0, 0.0));
  // pv = SensorBottomEnvelopeLogical.placeVolume(SensorDeadWireLogical, Position(0.0, (sensor_active_width-support_width/2.0) + sensor_dead_width/2.0 + sensor_deadwire_width/2.0, 0.0));
  pv = SensorTopEnvelopeLogical.placeVolume(SensorDeadWireLogical, Position(0.0, (-support_width/2.0) + (sensor_deadwire_width/2.0), 0.0));
  pv = SensorBottomEnvelopeLogical.placeVolume(SensorDeadWireLogical, Position(0.0, (-support_width/2.0) + (sensor_deadwire_width/2.0), 0.0));

  // place the active sensor and dead sensor inside the sensor envelope
  std::vector<dd4hep::PlacedVolume> TopSensor_pv;
  std::vector<dd4hep::PlacedVolume> BottomSensor_pv;
  for(int isensor=0; isensor < n_sensors_per_side; ++isensor){
     double sensor_total_z = n_sensors_per_side*sensor_active_len + dead_gap*(n_sensors_per_side-1);
     double xpos = 0.0;
     double ypos_active = (support_width/2.0) - (sensor_active_width/2.0);
     double ypos_dead = (-support_width/2.0) + sensor_deadwire_width + (sensor_dead_width/2.0);
     double zpos = -sensor_total_z/2.0 + sensor_active_len/2.0 + isensor*(sensor_active_len + dead_gap);
     pv = SensorTopEnvelopeLogical.placeVolume(SensorLogical, Position(xpos,ypos_active,zpos));
     //pv.addPhysVolID("topsensor",  isensor ) ;
     pv.addPhysVolID("sensor",  isensor ).addPhysVolID("barrelside", 1) ;
     TopSensor_pv.push_back(pv);
     pv = SensorBottomEnvelopeLogical.placeVolume(SensorLogical, Position(xpos,ypos_active,zpos));
     //pv.addPhysVolID("bottomsensor",  isensor ) ;
     pv.addPhysVolID("sensor",  isensor ).addPhysVolID("barrelside", -1) ;
     BottomSensor_pv.push_back(pv);
     pv = SensorTopEnvelopeLogical.placeVolume(SensorDeadLogical, Position(xpos,ypos_dead,zpos));
     pv = SensorBottomEnvelopeLogical.placeVolume(SensorDeadLogical, Position(xpos,ypos_dead,zpos));

  }
  //place the sensor envelope inside the ladder envelope
  pv = LadderLogical.placeVolume(SensorTopEnvelopeLogical,
                                Position(support_height/2.0 + flex_thickness + sensor_thickness/2.0, 0., 0.));//top-side sensors
  Transform3D tran_sen(RotationZYX(0., dd4hep::twopi/2.0, 0.), Position(-(support_height/2.0 + flex_thickness + sensor_thickness/2.0), 0., 0.));
  pv = LadderLogical.placeVolume(SensorBottomEnvelopeLogical,tran_sen);//bottom-side sensors

  //create the ladder support envelope
  Box LadderSupportEnvelopeSolid(support_height/2.0, support_width/2.0, support_length/2.0);
  Volume LadderSupportEnvelopeLogical(name + _toString( layer_id,"_SupEnvLogical_%02d"), LadderSupportEnvelopeSolid, air);
  vxd.setVisAttributes(theDetector, "seeThrough", LadderSupportEnvelopeLogical);

  //create ladder support volume
  Box LadderSupportSolid(support_thickness / 2.0 , support_width / 2.0 , support_length / 2.0);
  Volume LadderSupportLogical(name + _toString( layer_id,"_SupLogical_%02d"), LadderSupportSolid, support_mat);
  LadderSupportLogical.setVisAttributes(theDetector.visAttributes(supportVis));
  
  //vxd.setVisAttributes(theDetector, sensVis, SensorLogical);
  // vxd.setVisAttributes(theDetector, sensEnvVis, SensorBottomEnvelopeLogical);
  // vxd.setVisAttributes(theDetector, ladderVis, LadderLogical);
 
  pv = LadderSupportEnvelopeLogical.placeVolume(LadderSupportLogical);
  pv = LadderLogical.placeVolume(LadderSupportEnvelopeLogical);

  for(int i = 0; i < n_ladders; i++){
    std::stringstream ladder_enum; 
    ladder_enum << "vxt_ladder_" << layer_id << "_" << i;
    DetElement ladderDE(layerDE, ladder_enum.str(), x_det.id());
    std::cout << "start building " << ladder_enum.str() << ":" << endl;
    for(int isensor=0; isensor < n_sensors_per_side; ++isensor){
      std::stringstream topsensor_str;
      std::stringstream bottomsensor_str;
      topsensor_str << ladder_enum.str() << "_top_" << isensor;
      // std::cout << "\tstart building " << topsensor_str.str() << ":" << endl;
      bottomsensor_str << ladder_enum.str() << "_bottom_" << isensor;
      // std::cout << "\tstart building " << bottomsensor_str.str() << ":" << endl;
      DetElement topsensorDE(ladderDE, topsensor_str.str(), x_det.id());
      DetElement bottomsensorDE(ladderDE, bottomsensor_str.str(), x_det.id());
      topsensorDE.setPlacement(TopSensor_pv[isensor]);
      // std::cout << "\t" << topsensor_str.str() << " done." << endl;
      bottomsensorDE.setPlacement(BottomSensor_pv[isensor]);
      // std::cout << "\t" << bottomsensor_str.str() << " done." << endl;
    }
    Transform3D tr (RotationZYX(ladder_dphi*i,0.,0.),Position(ladder_radius*cos(ladder_phi0+ladder_dphi*i), ladder_radius*sin(ladder_phi0+ladder_dphi*i), 0.));
    pv = layer_assembly.placeVolume(LadderLogical,tr);
    pv.addPhysVolID("layer", layer_id ).addPhysVolID("module", i ) ; 
    ladderDE.setPlacement(pv);
    std::cout << ladder_enum.str() << " done." << endl;

  }
  
  // package the reconstruction data
  dd4hep::rec::ZPlanarData::LayerLayout topLayer;
  dd4hep::rec::ZPlanarData::LayerLayout bottomLayer;

  topLayer.ladderNumber         = n_ladders;
  topLayer.phi0                 = 0.;
  topLayer.sensorsPerLadder     = n_sensors_per_side;
  topLayer.lengthSensor         = sensor_active_len;
  topLayer.distanceSupport      = sensitive_radius;
  topLayer.thicknessSupport     = support_thickness / 2.0;
  topLayer.offsetSupport        = ladder_offset;
  topLayer.widthSupport         = support_width;
  topLayer.zHalfSupport         = support_length / 2.0;
  topLayer.distanceSensitive    = sensitive_radius + support_thickness / 2.0;
  topLayer.thicknessSensitive   = sensor_thickness;
  topLayer.offsetSensitive      = ladder_offset;
  topLayer.widthSensitive       = sensor_active_width;
  topLayer.zHalfSensitive       = (n_sensors_per_side*(sensor_active_len + dead_gap) - dead_gap) / 2.0;

  bottomLayer.ladderNumber         = n_ladders;
  bottomLayer.phi0                 = 0.;
  bottomLayer.sensorsPerLadder     = n_sensors_per_side;
  bottomLayer.lengthSensor         = sensor_active_len;
  bottomLayer.distanceSupport      = sensitive_radius - support_thickness / 2.0;
  bottomLayer.thicknessSupport     = support_thickness / 2.0;
  bottomLayer.offsetSupport        = ladder_offset;
  bottomLayer.widthSupport         = support_width;
  bottomLayer.zHalfSupport         = support_length / 2.0;
  bottomLayer.distanceSensitive    = sensitive_radius - support_thickness / 2.0 - sensor_thickness;
  bottomLayer.thicknessSensitive   = sensor_thickness;
  bottomLayer.offsetSensitive      = ladder_offset;
  bottomLayer.widthSensitive       = sensor_active_width;
  bottomLayer.zHalfSensitive       = (n_sensors_per_side*(sensor_active_len + dead_gap) - dead_gap) / 2.0;

  zPlanarData->layers.push_back(topLayer);
  zPlanarData->layers.push_back(bottomLayer);
 }
 std::cout << (*zPlanarData) << endl;
 vxd.addExtension< ZPlanarData >(zPlanarData);
 if ( x_det.hasAttr(_U(combineHits)) ) {
    vxd.setCombineHits(x_det.attr<bool>(_U(combineHits)),sens);
 }
 std::cout << "vxd done." << endl; 
  return vxd;
}
DECLARE_DETELEMENT(SiTrackerSkewBarrel_v01,create_element)
