#include <PHGarfield.h>

#include "PHGarfield.h"

#include "Garfield/ComponentUser.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"

#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/caen_correction.h>
#include <Event/packet.h>

#include <TSystem.h>

#include <algorithm>
#include <cstdint>                             // for uint16_t
#include <cstdlib>                             // for exit, size_t
#include <iostream>                             // for basic_ostream, operat...
#include <filesystem>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/TpcRawHitContainer.h>
#include <ffarawobjects/TpcRawHit.h>
#include <ffaobjects/SyncObject.h>
#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>
#include <ffaobjects/FlagSave.h>
#include <ffaobjects/CdbUrlSave.h>

#include <phfield/PHField3DCartesian.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include <TH1I.h>
#include <TH2I.h>
#include <TNtuple.h>
#include <TPolyLine3D.h>

using namespace std;
using namespace findNode;
using namespace Garfield;


bool         PHGarfield::search_complete=false;
unsigned int PHGarfield::calls=0;

PHGarfield::PHGarfield(const std::string &name)
  : SubsysReco(name)
{
 
  PHI_MIN = -M_PI;  // Local handling of Phi valued that wrap around.

  /*
  if (Verbosity() >= VERBOSITY_SOME) cout << "You have constructed :" << name << endl;
  SingleWaveForms = new TH2I("SingleWaveForms","SingleWaveForms",10,-0.5,9.5,1024,-0.5,1023.5);
  FeeHisto     = new TH1I("FeeHisto","FeeHisto",626,-1.5,625.5);
  ChannelHisto = new TH1I("ChannelHisto","ChannelHisto",258,-1.5,257.5);
  MOOLtiplicity= new TH1I("MOOLtiplicity","MOOLtiplicity",2000,-0.5,200000.5);
  dR       = new TH1I("dR","dR",100,-0.1,0.1);
  dPhi     = new TH1I("dPhi","dPhi",1000,-2.0*M_PI*2.0/800.0,2.0*M_PI*2.0/800.0);
  RdPhi     = new TH1I("RdPhi","RdPhi",1000,-2.0,2.0);
  PhiVsPHI =  new TH2I("PhiVsPHI","PhiVsPHI",1000,PHI_MIN,PHI_MIN+2.0*M_PI,1000,PHI_MIN,PHI_MIN+2.0*M_PI);
  PhiVsPHIII =  new TH2I("PhiVsPHIII","PhiVsPHIII",1000,PHI_MIN,PHI_MIN+2.0*M_PI,1000,PHI_MIN,PHI_MIN+2.0*M_PI);
  PhiVsPHIII_0 =  new TH2I("PhiVsPHIII_0","PhiVsPHIII_0",1000,PHI_MIN,PHI_MIN+2.0*M_PI,1000,PHI_MIN,PHI_MIN+2.0*M_PI);
  PhiVsPHIII_1 =  new TH2I("PhiVsPHIII_1","PhiVsPHIII_1",1000,PHI_MIN,PHI_MIN+2.0*M_PI,1000,PHI_MIN,PHI_MIN+2.0*M_PI);

  //MAPPY = new TNtuple("MAPPY","MAPPY","sector:side:phi:R:LAY:RAD:PHIII");
  
  
  Fun4AllServer *se = Fun4AllServer::instance();
  se->registerHisto(SingleWaveForms);
  se->registerHisto(FeeHisto);
  se->registerHisto(ChannelHisto);
  se->registerHisto(MOOLtiplicity);
  se->registerHisto(dR);
  se->registerHisto(dPhi);
  se->registerHisto(RdPhi);
  se->registerHisto(PhiVsPHI);
  se->registerHisto(PhiVsPHIII);
  se->registerHisto(PhiVsPHIII_0);
  se->registerHisto(PhiVsPHIII_1);
  //se->registerHisto(MAPPY);
  */
}

int PHGarfield::InitRun(PHCompositeNode *topNode)
{
  std::cout << "PHGarfield::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;
  m_cdb = CDBInterface::instance();
  
  //  Here we use the CDBInterface to set up the magnetic field map:
  std::string url = m_cdb->getUrl("FIELDMAP_TRACKING");
  m_field = new PHField3DCartesian(url, 1.0);
  
  //  Here we use the CDBInterface to set up the channel making of the TPC:
  std::string geofile = m_cdb->getUrl("Tracking_Geometry");
  std::cout << "===========================================Dude!  " << geofile << std::endl;
  
  std::string text = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");
  m_cdbTPCMAPttree = new CDBTTree(text.c_str());
  m_cdbTPCMAPttree->LoadCalibrations();

  //  Make the Garfield Component and register the methods that will interface to our fields...
  m_component = new Garfield::ComponentUser();
  m_component->SetMagneticField([this](double x, double y, double z,double& bx, double& by, double& bz) { GetMagneticFieldTesla(x, y, z, bx, by, bz); });
  m_component->SetElectricField([this](double x, double y, double z,double& ex, double& ey, double& ez) { GetElectricFieldVcm  (x, y, z, ex, ey, ez); });
  InitializeGas("/direct/phenix+u/workarea/hemmick/code.sphenix/tkh/gas/gasfiles/");
  
  //  Diagnostic during code development...
  PrintMaps();
  
  // Find all the interesting places...DST
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    { 
      if (Verbosity() >= VERBOSITY_MORE) cout << PHWHERE << "Dude, the DST node was not there, so I made one." << endl;
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
    }
  
  // Find all the interesting places...DST/TPC/
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "TPC"));
  if (!detNode)
    {
      if (Verbosity() >= VERBOSITY_MORE) cout << PHWHERE << "Dude, the TPC node was not there, so I made one." << endl;
      detNode = new PHCompositeNode("TPC");
      dstNode->addNode(detNode);
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGarfield::PrintMaps()
{
  //  Print out a few test points of the Garfield information
  double ex, ey, ez, bx, by, bz, vx, vy, vz;
  double x=0.0;
  double y=0.0;
  double z=0.0;
  GetMagneticFieldTesla( x, y, z, bx, by, bz);
  GetElectricFieldVcm  ( x, y, z, ex, ey, ez);
  m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  cout << " x:" << x
       << " y:" << y
       << " z:" << z
       << " ex:" << ex
       << " ey:" << ey
       << " ez:" << ez
       << " bx:" << bx
       << " by:" << by
       << " bz:" << bz
       << " vx:" << vx
       << " vy:" << vy
       << " vz:" << vz
       << endl;
  
  x=0.0;
  y=0.0;
  z=100.0;
  GetMagneticFieldTesla( x, y, z, bx, by, bz);
  GetElectricFieldVcm  ( x, y, z, ex, ey, ez);
  m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  cout << " x:" << x
       << " y:" << y
       << " z:" << z
       << " ex:" << ex
       << " ey:" << ey
       << " ez:" << ez
       << " bx:" << bx
       << " by:" << by
       << " bz:" << bz
       << " vx:" << vx
       << " vy:" << vy
       << " vz:" << vz
       << endl;
  
  x=0.0;
  y=40.0;
  z=100.0;
  GetMagneticFieldTesla( x, y, z, bx, by, bz);
  GetElectricFieldVcm  ( x, y, z, ex, ey, ez);
  m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  cout << " x:" << x
       << " y:" << y
       << " z:" << z
       << " ex:" << ex
       << " ey:" << ey
       << " ez:" << ez
       << " bx:" << bx
       << " by:" << by
       << " bz:" << bz
       << " vx:" << vx
       << " vy:" << vy
       << " vz:" << vz
       << endl;
  
  x=0.0;
  y=78.0;
  z=100.0;
  GetMagneticFieldTesla( x, y, z, bx, by, bz);
  GetElectricFieldVcm  ( x, y, z, ex, ey, ez);
  m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  cout << " x:" << x
       << " y:" << y
       << " z:" << z
       << " ex:" << ex
       << " ey:" << ey
       << " ez:" << ez
       << " bx:" << bx
       << " by:" << by
       << " bz:" << bz
       << " vx:" << vx
       << " vy:" << vy
       << " vz:" << vz
       << endl;
  
  //  Print out the pad coordinate map:
  int MAX = 10;
  int prints = 0;
  for (unsigned int side=0; side<2; side++)
    {
      for (unsigned int sector=0; sector<12; sector++)
	{
	  for (unsigned int fee=0; fee<26; fee++)
	    {
	      for (unsigned int channel=0; channel<256; channel++)
		{
		  unsigned int key = (256 * (fee)) + channel;
		  int    layer  = m_cdbTPCMAPttree->GetIntValue   (key, "layer");
		  double phi    = ((side == 1 ? 1 : -1) * (m_cdbTPCMAPttree->GetDoubleValue(key, "phi") - M_PI / 2.)) + ((sector % 12) * M_PI / 6);
		  double r      = m_cdbTPCMAPttree->GetDoubleValue(key, "R")/CLHEP::cm;

		  phi = bounder(phi, PHI_MIN);
	      
		  if (layer > 6)
		    {
		      radii[layer-7] = r;
		      if (prints < MAX)
			{
			  prints++;
			  cout << " side: "    << side;
			  cout << " sector: "  << sector;
			  cout << " fee: "     << fee;
			  cout << " channel: " << channel;
			  cout << " layer: "   << layer;
			  cout << " phi: "     << phi;
			  cout << " r: "       << r;
			  cout << endl;
			  /*
			    dR->Fill(R-RAD);
			    dPhi->Fill((phi-PHIII));
			    RdPhi->Fill(R*(phi-PHIII));
			    PhiVsPHI->Fill(phi,PHI);
			    PhiVsPHIII->Fill(phi,PHIII);
			    if (side == 0) PhiVsPHIII_0->Fill(phi,PHIII);
			    if (side == 1) PhiVsPHIII_1->Fill(phi,PHIII);
			  */
			}
		    }
		}
	    }
	}
    }

  for (int i=0; i<48; i++)
    {
      cout << " " << i<<":" <<radii[i];
    }
  cout << endl;
  
}

void PHGarfield::GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double& bx_t, double& by_t, double& bz_t)
{
  // NOTE:  Garfield uses  cm, V/cm, and Tesla.
  //        CLHEP    uses  mm, V/mm, and kiloTesla
  //        PHField3DCartesian follows the CLHEP conventions for magnetic fields.
  
  double point[4] =
  {
    x_cm * CLHEP::cm,
    y_cm * CLHEP::cm,
    z_cm * CLHEP::cm,
    0.0
  };

  double bfield[3] = {0.0, 0.0, 0.0};

  //  Get the magnetic field via the PHField3DCartesian object constructed usinf the CDB url reference.
  m_field->GetFieldValue(point, bfield);

  bx_t = bfield[0] / CLHEP::tesla;
  by_t = bfield[1] / CLHEP::tesla;
  bz_t = bfield[2] / CLHEP::tesla;
}

void PHGarfield::GetElectricFieldVcm(double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm)
{
  // NOTE:  Garfield uses  cm, V/cm, and Tesla.
  (void) x_cm;
  (void) y_cm;

  ex_vcm = 0.0;
  ey_vcm = 0.0;
  ez_vcm = z_cm > 0 ? -400.0 : 400.0;
  
}

void PHGarfield::InitializeGas(std::string dir)
{
  //  Create and fill the gas object so that we can trace particles through the gas...
  m_gas = new Garfield::MediumMagboltz();
  
  auto filename = [&](const int i) { return dir + "/PART_" + std::to_string(i) + ".gas"; };
  
  const std::string first = filename(0);
  if (!std::filesystem::exists(first))
    {
      std::cerr << "Missing first gas file: " << first << std::endl;
      return;
    }
  
  if (!m_gas->LoadGasFile(first))
    {
      std::cerr << "Failed to load " << first << std::endl;
      return;
    }
  
  int nFiles = 1;
  for (int i = 1; ; ++i)
    {
      const std::string file = filename(i);
      
      if (!std::filesystem::exists(file))
	{
	  std::cout << "Stopping at first missing file: " << file << std::endl;
	  break;
	}
      
      std::cout << "Merging " << file << std::endl;
      
      if (!m_gas->MergeGasFile(file, true))
	{
	  std::cerr << "Failed to merge " << file << std::endl;
	  return;
	}
      
      ++nFiles;
    }
}

int PHGarfield::process_event(PHCompositeNode *topNode)
{  
  calls++;
  if (calls%2 == 0)
    {
      cout << PHWHERE << "Calls: " << calls << endl;
      topNode->print();
    }
  //if (!search_complete) {SearchNodeTree(topNode);}

  if (Verbosity() >= VERBOSITY_A_LOT)
    {
      cout << "Dude, print the CDB" << endl;
      m_cdb->Print();
      cout << "Did it work?" << endl;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHGarfield::bounder(double phi, double phi_min)
{

  double phi_max = phi_min + 2.0*M_PI;
  while (phi <  phi_min) phi = phi + 2.0*M_PI;
  while (phi >= phi_max) phi = phi - 2.0*M_PI;

  return phi;
}


TPolyLine3D *PHGarfield::ReverseDrift (double x, double y, double z, double step_ns)
{
  vector<double> xlist;
  vector<double> ylist;
  vector<double> zlist;

  xlist.push_back(x);
  ylist.push_back(y);
  zlist.push_back(z);
  
  double ex, ey, ez, bx, by, bz, vx, vy, vz;
  
  double zPrevious = z;
  while (!StopHere(x,y,z,zPrevious))
    {
      zPrevious = z;
      GetMagneticFieldTesla( x, y, z, bx, by, bz);
      GetElectricFieldVcm  ( x, y, z, ex, ey, ez);
      m_gas->ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);

      x = x - vx*step_ns;
      y = y - vy*step_ns;
      z = z - vz*step_ns;
      
      xlist.push_back(x);
      ylist.push_back(y);
      zlist.push_back(z);
    }

  TPolyLine3D *poly = new TPolyLine3D(xlist.size() - 1);
  for (unsigned int i=0; i<xlist.size() -1; i++)
    {
      poly->SetPoint(i,xlist[i], ylist[i], zlist[i]);
    }

  return poly;
}

bool PHGarfield::StopHere(const double x, const double y, const double z,
              const double zPrevious)
{
  const double r = std::hypot(x, y);
  
  if (r < 18.0) return true;
  if (r > 82.0) return true;
  if (z > 120.0) return true;
  if (z < -120.0) return true;
  
  // z crossed the central membrane.
  if (z * zPrevious < 0.0) return true;
  
  return false;
}

/*
int PHGarfield() {
  // ------------------------------------------------------------
  // Gas: Ar/CF4/isobutane = 75/20/5.
  // ------------------------------------------------------------
  MediumMagboltz gas;
  gas.SetComposition("ar", 75., "cf4", 20., "isobutane", 5.);
  gas.SetTemperature(293.15);  // K
  gas.SetPressure(760.);       // Torr

  // Try to load an existing gas table.
  if (!gas.LoadGasFile("ar_cf4_iso_75_20_5.gas")) {
    std::cout << "No gas file found; generating gas table...\n";
    
    gas.SetFieldGrid(300., 500., 5, false,
		     0., 0., 1,
		     0., 0., 1);
    
    gas.GenerateGasTable(10);
    
    gas.WriteGasFile("ar_cf4_iso_75_20_5.gas");
  }
  
  // Precompute electron transport data.
  // For serious work, you may want finer E/B grids and save/load gas files.
  gas.Initialise(true);

  // ------------------------------------------------------------
  // User electric field.
  // Coordinates in cm, E in V/cm.
  // ------------------------------------------------------------
  ComponentUser field;

  field.SetElectricField(
    [](const double x, const double y, const double z,
       double& ex, double& ey, double& ez) {

      (void)x;
      (void)y;

      ex = 0.;
      ey = 0.;
      ez = (z > 0.) ? -400. : +400.;
    }
  );

  //PHField* phfield = PHFieldUtility::BuildFieldMap(PHFieldUtility::DefaultFieldConfig());
  //PHField* phfield = BuildFieldMap(DefaultFieldConfig());
  field.SetMagneticField([](const double x, const double y, const double z, double& bx, double& by, double& bz) {
    
    double point[4] = {x, y, z, 0.0};
    double b[3] = {0.0, 0.0, 0.0};

    cout << point << endl;
    
    //phfield->GetFieldValue(point, b);
    
    bx = b[0];
    by = b[1];
    bz = b[2];
  }
    );

  // Active region, in cm.
  field.SetArea(-78., -78., -100., 78., 78., 100.);
  field.SetMedium(&gas);

  // ------------------------------------------------------------
  // Sensor holds the field component.
  // ------------------------------------------------------------
  Sensor sensor;
  sensor.AddComponent(&field);
  sensor.SetArea(-78., -78., -100., 78., 78., 100.);

  // ------------------------------------------------------------
  // Drift engine.
  // ------------------------------------------------------------
  DriftLineRKF drift;
  drift.SetSensor(&sensor);

  // This is the important line for your case.
  // drift.EnableDiffusion(false);  // Bad ChatGPT

  // Launch point.
  const double x0 = 30.0;   // cm
  const double y0 = 0.0;    // cm
  const double z0 = 50.0;   // cm
  const double t0 = 0.0;    // ns

  drift.DriftElectron(x0, y0, z0, t0);

  // ------------------------------------------------------------
  // Inspect the drift line and apply your stopping criteria.
  // ------------------------------------------------------------
  // const auto& path = drift.GetDriftLine(); //wrong again albert
  const size_t nPoints = drift.GetNumberOfDriftLinePoints();

  double zPrevious = z0;
  
  for (size_t i = 0; i < nPoints; ++i) {
    double x = 0., y = 0., z = 0., t = 0.;
    drift.GetDriftLinePoint(i, x, y, z, t);
    
    std::cout << x << "  " << y << "  " << z << "  " << t << "\n";
    
    if (StopHere(x, y, z, zPrevious)) {
      std::cout << "Stopping condition reached at: "
		<< "x=" << x << " y=" << y << " z=" << z
		<< " r=" << std::hypot(x, y)
		<< " t=" << t << "\n";
      break;
    }
    
    zPrevious = z;
  }

  return 0;
}
*/
