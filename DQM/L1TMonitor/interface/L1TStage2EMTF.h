#ifndef DQM_L1TMonitor_L1TStage2EMTF_h
#define DQM_L1TMonitor_L1TStage2EMTF_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/L1TMuon/interface/EMTFDaqOut.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"


class L1TStage2EMTF : public DQMEDAnalyzer {

 public:

  L1TStage2EMTF(const edm::ParameterSet& ps);
  virtual ~L1TStage2EMTF();

 protected:

  virtual void dqmBeginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);
  virtual void bookHistograms(DQMStore::IBooker&, const edm::Run&, const edm::EventSetup&) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

 private:

  edm::EDGetTokenT<l1t::EMTFDaqOutCollection> daqToken;
  edm::EDGetTokenT<l1t::EMTFHitCollection> hitToken;
  edm::EDGetTokenT<l1t::EMTFTrackCollection> trackToken;
  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> muonToken;
  std::string monitorDir;
  bool isData;
  bool filterBX;
  bool verbose;

  MonitorElement* emtfErrors;

  MonitorElement* emtfHitBX;
  MonitorElement* emtfHitStrip[20];
  MonitorElement* emtfHitWire[20];
  MonitorElement* emtfChamberStrip[20];
  MonitorElement* emtfChamberWire[20];
  MonitorElement* emtfHitOccupancy;

  MonitorElement* emtfnTracks;
  MonitorElement* emtfTracknHits;
  MonitorElement* emtfTrackBX;
  MonitorElement* emtfTrackPt;
  MonitorElement* emtfTrackPtCoarse;
  MonitorElement* emtfTrackBX1D;
  MonitorElement* emtfTrackEta;
  MonitorElement* emtfTrackEtaCoarse;
  MonitorElement* emtfTrackPhi;
  MonitorElement* emtfTrackPhiCoarse;
  MonitorElement* emtfTrackPhiHighQuality;
  MonitorElement* emtfTrackOccupancy;
  MonitorElement* emtfTrackMode;
  MonitorElement* emtfTrackQuality;
  MonitorElement* emtfTrackQualityVsMode;
  MonitorElement* emtfTrackSectorIndex;

  MonitorElement* emtfMuonBX;
  MonitorElement* emtfMuonhwPt;
  MonitorElement* emtfMuonhwPtCoarse;
  MonitorElement* emtfMuonhwEta;
  MonitorElement* emtfMuonhwEtaCoarse;
  MonitorElement* emtfMuonhwPhi;
  MonitorElement* emtfMuonhwPhiCoarse;
  MonitorElement* emtfMuonhwQual;

//=========QUALITY TESTER=================
  MonitorElement* emtfChamberStrip_QT[20]; //have to reallocate for two additional histograms
  MonitorElement* emtfChamberStrip_QT_sqrt[20];
//  MonitorElement* emtfChamberStrip_QT_mean[20];
  MonitorElement* emtfChamberStrip_QT_hot[20];

  MonitorElement* emtfChamberStrip_QT1D[20];
  MonitorElement* emtfChamberStrip_QT_sqrt1D[20];
  MonitorElement* emtfChamberStrip_QT_dead[20];

  MonitorElement* emtfTrackBX_QT;
  MonitorElement* emtfTrackBX_QT_sqrt;
  MonitorElement* emtfTrackBX_QT_hot;
  MonitorElement* emtfTrackBX_QT_dead;

  MonitorElement* emtfTrackPhiCoarse_QT;
  MonitorElement* emtfTrackPhiCoarse_QT_sqrt;
  MonitorElement* emtfTrackPhiCoarse_QT_hot;
  MonitorElement* emtfTrackPhiCoarse_QT_dead;

  MonitorElement* emtfTrackPhiHighQuality_QT;
  MonitorElement* emtfTrackPhiHighQuality_QT_sqrt;
  MonitorElement* emtfTrackPhiHighQuality_QT_hot;
  MonitorElement* emtfTrackPhiHighQuality_QT_dead;
 //========================================
};

#endif
