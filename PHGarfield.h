#ifndef PHGARFIELD__H
#define PHGARFIELD__H

#include <cmath>
#include <iostream>

#include <fun4all/SubsysReco.h>

//class Gl1Packet;
//class TpcRawHitContainer;
//class SyncObject;
//class EventHeader;
//class RunHeader;
//class FlagSave;
class CdbUrlSave;

#include <string>
#include <vector>

class CDBInterface;
class CDBTTree;
class PHCompositeNode;
class TpcMap;
class PHField3DCartesian;
class ComponentUser;
class TPolyLine3D;

namespace Garfield
{
  class ComponentUser;
  class MediumMagboltz;
}

//class TH1;
//class TH2;
//class TNtuple;


class PHGarfield : public SubsysReco
{
 public:
  PHGarfield(const std::string &name = "PHGarfield");
  ~PHGarfield() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  bool StopHere(const double x, const double y, const double z, const double zPrevious);
  void PrintMaps();

  //  These are left in public namespace for easy plotting macros...
  TPolyLine3D *ReverseDrift (double x_cm, double y_cm, double z_cm, double step_ns=50.0); // Drifts electrons from some initial point until they hit a detector boundary...
  double radii[48];  // Radius on each layer just for test purposes...need to be cm!
  
 private:
  CDBInterface             *m_cdb            {nullptr}; // Access to all thiungs CDB...
  CDBTTree                 *m_cdbTPCMAPttree {nullptr}; // Locations of the pads from CDB...
  PHField3DCartesian       *m_field          {nullptr}; // The stanards sPHENIX field holding container.
  Garfield::ComponentUser  *m_component      {nullptr}; // This handles the interface of the electric and magnetic fields as handed to Garfield
  Garfield::MediumMagboltz *m_gas            {nullptr}; // This is the pre-tabulated gas properties required by Garfield...

  void GetMagneticFieldTesla(double x_cm, double y_cm, double z_cm, double& bx_t,   double& by_t,   double& bz_t) ;   // Feeds magnetic field to Garfield
  void GetElectricFieldVcm  (double x_cm, double y_cm, double z_cm, double& ex_vcm, double& ey_vcm, double& ez_vcm) ; // Feeds electric field to Garfield
  void InitializeGas        (std::string dir);
  
  //  These are utilities for a spot check of the overall routine:
  
  std::string calibdir;
  std::string m_DiodeContainerName;

  double bounder(double phi, double phi_min);
  double PHI_MIN;

  static bool search_complete;
  static unsigned int calls;

  /*
  void SearchNodeTree(PHCompositeNode *);
  Gl1Packet          *m_gl1;
  TpcRawHitContainer *m_tpc[24][2];
  SyncObject         *m_sync;
  EventHeader        *m_event;
  RunHeader          *m_runheader;
  FlagSave           *m_flagsave;
  //CdbUrlSave         *m_cdb;


  // Repository copy of a mapping...
  int mc_sectors[12] {5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
  int FEE_map[26]    {4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};
  int FEE_R[26]      {2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
  
  TH2 *SingleWaveForms;
  TH1 *FeeHisto;
  TH1 *ChannelHisto;

  TH1 *dR;
  TH1 *dPhi;
  TH1 *RdPhi;
  TH2 *PhiVsPHI;
  TH2 *PhiVsPHIII;
  TH2 *PhiVsPHIII_0;
  TH2 *PhiVsPHIII_1;
  TH1 *MOOLtiplicity;
  
  TNtuple *MAPPY;
  */
};

#endif
