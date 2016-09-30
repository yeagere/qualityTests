#include <string>
#include <vector>

#include <iostream>


#include "TMath.h"
Double_t pi=TMath::Pi(); //for trackPhi quality tests

#include "DQM/L1TMonitor/interface/L1TStage2EMTF.h"

L1TStage2EMTF::L1TStage2EMTF(const edm::ParameterSet& ps)
    : daqToken(consumes<l1t::EMTFDaqOutCollection>(ps.getParameter<edm::InputTag>("emtfSource"))),
      hitToken(consumes<l1t::EMTFHitCollection>(ps.getParameter<edm::InputTag>("emtfSource"))),
      trackToken(consumes<l1t::EMTFTrackCollection>(ps.getParameter<edm::InputTag>("emtfSource"))),
      muonToken(consumes<l1t::RegionalMuonCandBxCollection>(ps.getParameter<edm::InputTag>("emtfSource"))),
      monitorDir(ps.getUntrackedParameter<std::string>("monitorDir", "")),
      isData(ps.getUntrackedParameter<bool>("isData", false)),
      filterBX(ps.getUntrackedParameter<bool>("filterBX", false)),
      verbose(ps.getUntrackedParameter<bool>("verbose", false)) {}

L1TStage2EMTF::~L1TStage2EMTF() {}

void L1TStage2EMTF::dqmBeginRun(const edm::Run& r, const edm::EventSetup& c) {}

void L1TStage2EMTF::beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) {}

void L1TStage2EMTF::bookHistograms(DQMStore::IBooker& ibooker, const edm::Run&, const edm::EventSetup&) {

  ibooker.setCurrentFolder(monitorDir);
  
  // DAQ Output Monitor Elements
  emtfErrors = ibooker.book1D("emtfErrors", "EMTF Errors", 6, 0, 6);
  emtfErrors->setAxisTitle("Error Type (Corruptions Not Implemented)", 1);
  emtfErrors->setAxisTitle("Number of Errors", 2);
  emtfErrors->setBinLabel(1, "Corruptions", 1);
  emtfErrors->setBinLabel(2, "Synch. Err.", 1);
  emtfErrors->setBinLabel(3, "Synch. Mod.", 1);
  emtfErrors->setBinLabel(4, "BX Mismatch", 1);
  emtfErrors->setBinLabel(5, "Time Misalign.", 1);
  emtfErrors->setBinLabel(6, "FMM != Ready", 1);

  // Hit (LCT) Monitor Elements
  int n_xbins, n_ybins;
  std::string name, label;
  std::vector<std::string> suffix_name = {"42", "41", "32", "31", "22", "21", "13", "12", "11b", "11a"};
  std::vector<std::string> suffix_label = {"4/2", "4/1", "3/2", "3/1", " 2/2", "2/1", "1/3", "1/2", "1/1b", "1/1a"}; 

  emtfHitBX = ibooker.book2D("emtfHitBX", "EMTF Hit BX", 8, -3, 5, 18, 0, 18);
  emtfHitBX->setAxisTitle("BX", 1);
  for (int xbin = 1, xbin_label = -3; xbin <= 8; ++xbin, ++xbin_label) {
    emtfHitBX->setBinLabel(xbin, std::to_string(xbin_label), 1);
  }
  for (int ybin = 1; ybin <= 9; ++ybin) {
    emtfHitBX->setBinLabel(ybin, "ME-" + suffix_label[ybin - 1], 2);
    emtfHitBX->setBinLabel(19 - ybin, "ME+" + suffix_label[ybin - 1], 2);
  }

  for (int hist = 0, i = 0; hist < 20; ++hist, i = hist % 10) {

    int yMin;
    int yMax;

    if (hist < 10) {
      name = "MENeg" + suffix_name[i];
      label = "ME-" + suffix_label[i];
    } else {
      name = "MEPos" + suffix_name[9 - i];
      label = "ME+" + suffix_label[9 - i];
    }

    if (hist < 6) { 
      n_xbins = (i % 2) ? 18 : 36;
    } else if (hist > 13) {
      n_xbins = ((i+1) % 2) ? 18 : 36;
    } else {
      n_xbins = 36;
    }

    if (hist == 6 || hist == 9 || hist == 10 || hist == 13) {
      yMin = 0;
      yMax = 128; //8 bins
    } else if (hist == 8 || hist == 11) {
      yMin = 128;
      yMax = 224; //6 bins 
    } else {
      yMin = 0;
      yMax = 160; //10 bins
    }

    n_ybins = (yMax-yMin)/16;
     
    emtfHitStrip[hist] = ibooker.book1D("emtfHitStrip" + name, "EMTF Halfstrip " + label, 256, 0, 256);
    emtfHitStrip[hist]->setAxisTitle("Cathode Halfstrip, " + label, 1);

    emtfHitWire[hist] = ibooker.book1D("emtfHitWire" + name, "EMTF Wiregroup " + label, 128, 0, 128);
    emtfHitWire[hist]->setAxisTitle("Anode Wiregroup, " + label, 1);

    emtfChamberStrip[hist] = ibooker.book2D("emtfChamberStrip" + name, "EMTF Halfstrip " + label, n_xbins, 1, 1+n_xbins, 256, 0, 256);
    emtfChamberStrip[hist]->setAxisTitle("Chamber, " + label, 1);
    emtfChamberStrip[hist]->setAxisTitle("Cathode Halfstrip", 2);

    emtfChamberWire[hist] = ibooker.book2D("emtfChamberWire" + name, "EMTF Wiregroup " + label, n_xbins, 1, 1+n_xbins, 128, 0, 128);
    emtfChamberWire[hist]->setAxisTitle("Chamber, " + label, 1);
    emtfChamberWire[hist]->setAxisTitle("Anode Wiregroup", 2);

    for (int bin = 1; bin <= n_xbins; ++bin) {
      emtfChamberStrip[hist]->setBinLabel(bin, std::to_string(bin), 1);
      emtfChamberWire[hist]->setBinLabel(bin, std::to_string(bin), 1);
    }

//================QUALITY TESTER===============================================================================
//================CHAMBER STRIP================================================================================
    emtfChamberStrip_QT[hist] = ibooker.book2D("emtfChamberStrip_QT" + name, "EMTF Halfstrip QT " + label, n_xbins, 1, 1+n_xbins, n_ybins, yMin, yMax);
    emtfChamberStrip_QT[hist]->setAxisTitle("Chamber, " + label, 1);
    emtfChamberStrip_QT[hist]->setAxisTitle("Cathode Halfstrip", 2);

    emtfChamberStrip_QT_sqrt[hist] = ibooker.book2D("emtfChamberStrip_QT_sqrt" + name, "EMTF Halfstrip QT " + label, n_xbins, 1, 1+n_xbins, n_ybins, yMin, yMax);
    emtfChamberStrip_QT_sqrt[hist]->setAxisTitle("Chamber, " + label, 1);
    emtfChamberStrip_QT_sqrt[hist]->setAxisTitle("Cathode Halfstrip", 2);

    // if we want to calculate the average by line
/*    emtfChamberStrip_QT_mean[hist] = ibooker.book1D("emtfChamberStrip_QT_mean" + name, "EMTF Halfstrip QT Mean " + label, 16, yMin, yMax);
    emtfChamberStrip_QT_mean[hist]->setAxisTitle("Cathode Halfstrip", 1); */

    emtfChamberStrip_QT_hot[hist] = ibooker.book2D("emtfChamberStrip_QT_hot" + name, "EMTF Halfstrip QT Hot " + label, n_xbins, 1, 1+n_xbins, n_ybins, yMin, yMax);
    emtfChamberStrip_QT_hot[hist]->setAxisTitle("Chamber, " + label, 1);
    emtfChamberStrip_QT_hot[hist]->setAxisTitle("Cathode Halfstrip", 2);

    //for all noisychannel inverted histograms with 1 bin for every column
    emtfChamberStrip_QT1D[hist] = ibooker.book1D("emtfChamberStrip_QT1D" + name, "EMTF Halfstrip QT 1D " + label,  n_xbins, 1, 1+n_xbins);
    emtfChamberStrip_QT1D[hist]->setAxisTitle("Chamber, " + label, 1);

    emtfChamberStrip_QT_sqrt1D[hist] = ibooker.book1D("emtfChamberStrip_QT_sqrt1D" + name, "EMTF Halfstrip QT 1D " + label, n_xbins, 1, 1+n_xbins);
    emtfChamberStrip_QT_sqrt1D[hist]->setAxisTitle("Chamber, " + label, 1);

    emtfChamberStrip_QT_dead[hist] = ibooker.book1D("emtfChamberStrip_QT_dead" + name, "EMTF Halfstrip QT Dead " + label, n_xbins, 1, 1+n_xbins);
    emtfChamberStrip_QT_dead[hist]->setAxisTitle("Chamber, " + label, 1);

    for (int bin = 1; bin <= n_xbins; ++bin) {
      emtfChamberStrip_QT[hist]->setBinLabel(bin, std::to_string(bin), 1);
      emtfChamberStrip_QT_sqrt[hist]->setBinLabel(bin, std::to_string(bin), 1);
      emtfChamberStrip_QT_hot[hist]->setBinLabel(bin, std::to_string(bin), 1);

      emtfChamberStrip_QT1D[hist]->setBinLabel(bin, std::to_string(bin), 1);
      emtfChamberStrip_QT_sqrt1D[hist]->setBinLabel(bin, std::to_string(bin), 1);
      emtfChamberStrip_QT_dead[hist]->setBinLabel(bin, std::to_string(bin), 1);
    }
//=============================================================================================================
//=============================================================================================================
}

  emtfHitOccupancy = ibooker.book2D("emtfHitOccupancy", "EMTF Chamber Occupancy", 54, 1, 55, 10, -5, 5);
  emtfHitOccupancy->setAxisTitle("Sector (CSCID 1-9 Unlabelled)", 1);
  for (int bin = 1; bin <= 46; bin += 9) {
    emtfHitOccupancy->setBinLabel(bin, std::to_string(bin % 8), 1);
  }
  emtfHitOccupancy->setBinLabel(1, "ME-4", 2);
  emtfHitOccupancy->setBinLabel(2, "ME-3", 2);
  emtfHitOccupancy->setBinLabel(3, "ME-2", 2);
  emtfHitOccupancy->setBinLabel(4, "ME-1b", 2);
  emtfHitOccupancy->setBinLabel(5, "ME-1a", 2);
  emtfHitOccupancy->setBinLabel(6, "ME+1a", 2);
  emtfHitOccupancy->setBinLabel(7, "ME+1b", 2);
  emtfHitOccupancy->setBinLabel(8, "ME+2", 2);
  emtfHitOccupancy->setBinLabel(9, "ME+3", 2);
  emtfHitOccupancy->setBinLabel(10, "ME+4", 2);

  // Track Monitor Elements
  emtfnTracks = ibooker.book1D("emtfnTracks", "Number of EMTF Tracks per Event", 11, 0, 11);
  for (int bin = 1; bin <= 10; ++bin) {
    emtfnTracks->setBinLabel(bin, std::to_string(bin - 1), 1);
  }
  emtfnTracks->setBinLabel(11, "Overflow", 1);

  emtfTracknHits = ibooker.book1D("emtfTracknHits", "Number of Hits per EMTF Track", 5, 0, 5);
  for (int bin = 1; bin <= 5; ++bin) {
    emtfTracknHits->setBinLabel(bin, std::to_string(bin - 1), 1);
  }

  emtfTrackBX = ibooker.book2D("emtfTrackBX", "EMTF Track Bunch Crossing", 12, -6, 6, 8, -3, 5);
  emtfTrackBX->setAxisTitle("Sector (Endcap)", 1);
  for (int i = 0; i < 6; ++i) {
    emtfTrackBX->setBinLabel(i + 1, std::to_string(6 - i) + " (-)", 1);
    emtfTrackBX->setBinLabel(12 - i, std::to_string(6 - i) + " (+)", 1);
  }

  emtfTrackBX->setAxisTitle("Track BX", 2);
  emtfTrackBX1D = ibooker.book1D("emtfTrackBX1D", "EMTF Track BX", 8, -3, 5);

  emtfTrackBX->setAxisTitle("Track BX", 1);
  for (int bin = 1, i = -3; bin <= 8; ++bin, ++i) {
    emtfTrackBX->setBinLabel(bin, std::to_string(i), 2); 
    emtfTrackBX1D->setBinLabel(bin, std::to_string(i), 1);
  }

//==================QUALITY TESTER=============================================================================
//===================TRACK BX==================================================================================

  emtfTrackBX_QT = ibooker.book1D("emtfTrackBX_QT", "EMTF Track BX QT", 12, -6, 6);
  emtfTrackBX_QT->setAxisTitle("Sector (Endcap)", 1);
  for (int i = 0; i < 6; ++i) {
    emtfTrackBX_QT->setBinLabel(i + 1, std::to_string(6 - i) + " (-)", 1);
    emtfTrackBX_QT->setBinLabel(12 - i, std::to_string(6 - i) + " (+)", 1);
  }

  emtfTrackBX_QT_sqrt = ibooker.book1D("emtfTrackBX_QT_sqrt", "EMTF Track BX QT", 12, -6, 6);
  emtfTrackBX_QT_sqrt->setAxisTitle("Sector (Endcap)", 1);
  for (int i = 0; i < 6; ++i) {
    emtfTrackBX_QT_sqrt->setBinLabel(i + 1, std::to_string(6 - i) + " (-)", 1);
    emtfTrackBX_QT_sqrt->setBinLabel(12 - i, std::to_string(6 - i) + " (+)", 1);
  }

  emtfTrackBX_QT_hot = ibooker.book1D("emtfTrackBX_QT_hot", "EMTF Track BX QT Hot", 12, -6, 6);
  emtfTrackBX_QT_hot->setAxisTitle("Sector (Endcap)", 1);
  for (int i = 0; i < 6; ++i) {
    emtfTrackBX_QT_hot->setBinLabel(i + 1, std::to_string(6 - i) + " (-)", 1);
    emtfTrackBX_QT_hot->setBinLabel(12 - i, std::to_string(6 - i) + " (+)", 1);
  }

  emtfTrackBX_QT_dead = ibooker.book1D("emtfTrackBX_QT_dead", "EMTF Track BX QT Dead", 12, -6, 6);
  emtfTrackBX_QT_dead->setAxisTitle("Sector (Endcap)", 1);
  for (int i = 0; i < 6; ++i) {
    emtfTrackBX_QT_dead->setBinLabel(i + 1, std::to_string(6 - i) + " (-)", 1);
    emtfTrackBX_QT_dead->setBinLabel(12 - i, std::to_string(6 - i) + " (+)", 1);
  }

//=============================================================================================================
//=============================================================================================================

  emtfTrackPt = ibooker.book1D("emtfTrackPt", "EMTF Track p_{T}", 256, 1, 257);
  emtfTrackPt->setAxisTitle("Track p_{T} [GeV]", 1);

  emtfTrackPtCoarse = ibooker.book1D("emtfTrackPtCoarse", "EMTF Track p_{T} (coarse)", 32, 1, 257);
  emtfTrackPtCoarse->setAxisTitle("Track p_{T} [GeV]", 1);
  //emtfTrackPtCoarse->getTH1F()->Sumw2();

  emtfTrackEta = ibooker.book1D("emtfTrackEta", "EMTF Track #eta", 100, -2.5, 2.5);
  emtfTrackEta->setAxisTitle("Track #eta", 1);

  emtfTrackEtaCoarse = ibooker.book1D("emtfTrackEtaCoarse", "EMTF Track #eta (coarse)", 20, -2.5, 2.5);
  emtfTrackEtaCoarse->setAxisTitle("Track #eta", 1);
  //emtfTrackEtaCoarse->getTH1F()->Sumw2();

  emtfTrackPhi = ibooker.book1D("emtfTrackPhi", "EMTF Track #phi", 128, -3.2, 3.2);
  emtfTrackPhi->setAxisTitle("Track #phi", 1);

  emtfTrackPhiCoarse = ibooker.book1D("emtfTrackPhiCoarse", "EMTF Track #phi (coarse)", 32, -3.2, 3.2);
  emtfTrackPhiCoarse->setAxisTitle("Track #phi", 1);
  //emtfTrackPhiCoarse->getTH1F()->Sumw2();

  emtfTrackPhiHighQuality = ibooker.book1D("emtfTrackPhiHighQuality", "EMTF High Quality #phi", 128, -3.2, 3.2);
  emtfTrackPhiHighQuality->setAxisTitle("Track #phi (High Quality)", 1);

//==================QUALITY TESTER=============================================================================
//=================TRACK PHI===================================================================================
  Double_t xMin = -(63.0/64.0)*pi;
  Double_t xMax = (65.0/64.0)*pi;

  emtfTrackPhiCoarse_QT = ibooker.book1D("emtfTrackPhiCoarse_QT", "EMTF Track #phi (coarse) QT", 32, xMin, xMax);
  emtfTrackPhiCoarse_QT->setAxisTitle("Track #phi", 1);

  emtfTrackPhiCoarse_QT_sqrt = ibooker.book1D("emtfTrackPhiCoarse_QT_sqrt", "EMTF Track #phi (coarse) QT", 32, xMin, xMax);
  emtfTrackPhiCoarse_QT_sqrt->setAxisTitle("Track #phi", 1);

  emtfTrackPhiCoarse_QT_hot = ibooker.book1D("emtfTrackPhiCoarse_QT_hot", "EMTF Track #phi (coarse) QT Hot", 32, xMin, xMax);
  emtfTrackPhiCoarse_QT_hot->setAxisTitle("Track #phi", 1);

  emtfTrackPhiCoarse_QT_dead = ibooker.book1D("emtfTrackPhiCoarse_QT_dead", "EMTF Track #phi (coarse) QT Dead", 32, xMin, xMax);
  emtfTrackPhiCoarse_QT_dead->setAxisTitle("Track #phi", 1);

  emtfTrackPhiHighQuality_QT = ibooker.book1D("emtfTrackPhiHighQuality_QT", "EMTF High Quality #phi QT", 128, xMin, xMax);
  emtfTrackPhiHighQuality_QT->setAxisTitle("Track #phi (High Quality)", 1);

  emtfTrackPhiHighQuality_QT_sqrt = ibooker.book1D("emtfTrackPhiHighQuality_QT_sqrt", "EMTF High Quality #phi QT", 128, xMin, xMax);
  emtfTrackPhiHighQuality_QT_sqrt->setAxisTitle("Track #phi (High Quality)", 1);

  emtfTrackPhiHighQuality_QT_hot = ibooker.book1D("emtfTrackPhiHighQuality_QT_hot", "EMTF High Quality #phi QT Hot", 128, xMin, xMax);
  emtfTrackPhiHighQuality_QT_hot->setAxisTitle("Track #phi (High Quality)", 1);

  emtfTrackPhiHighQuality_QT_dead = ibooker.book1D("emtfTrackPhiHighQuality_QT_dead", "EMTF High Quality #phi QT Dead", 128, xMin, xMax);
  emtfTrackPhiHighQuality_QT_dead->setAxisTitle("Track #phi (High Quality)", 1);
//=============================================================================================================
//=============================================================================================================


  emtfTrackOccupancy = ibooker.book2D("emtfTrackOccupancy", "EMTF Track Occupancy", 100, -2.5, 2.5, 126, -3.15, 3.15);
  emtfTrackOccupancy->setAxisTitle("#eta", 1);
  emtfTrackOccupancy->setAxisTitle("#phi", 2);

  emtfTrackSectorIndex = ibooker.book1D("emtfTrackSectorIndex", "EMTF Track Sector Index", 13, -6.5, 6.5);

  emtfTrackMode = ibooker.book1D("emtfTrackMode", "EMTF Track Mode", 16, 0, 16);
  emtfTrackMode->setAxisTitle("Mode", 1);


  emtfTrackQuality = ibooker.book1D("emtfTrackQuality", "EMTF Track Quality", 16, 0, 16);
  emtfTrackQuality->setAxisTitle("Quality", 1);


  emtfTrackQualityVsMode = ibooker.book2D("emtfTrackQualityVsMode", "EMTF Track Quality vs Mode", 16, 0, 16, 16, 0, 16);
  emtfTrackQualityVsMode->setAxisTitle("Mode", 1);
  emtfTrackQualityVsMode->setAxisTitle("Quality", 2);

  for (int bin = 1; bin <= 16; ++bin) {
    emtfTrackMode->setBinLabel(bin, std::to_string(bin - 1), 1);
    emtfTrackQuality->setBinLabel(bin, std::to_string(bin - 1), 1);
    emtfTrackQualityVsMode->setBinLabel(bin, std::to_string(bin - 1), 1);
    emtfTrackQualityVsMode->setBinLabel(bin, std::to_string(bin - 1), 2);
  }

  // Regional Muon Candidate Monitor Elements
  ibooker.setCurrentFolder(monitorDir + "/MuonCand");

  emtfMuonBX = ibooker.book1D("emtfMuonBX", "EMTF Muon Cand BX", 7, -3, 4);
  emtfMuonBX->setAxisTitle("BX", 1);
  for (int bin = 1, bin_label = -3; bin <= 7; ++bin, ++bin_label) {
    emtfMuonBX->setBinLabel(bin, std::to_string(bin_label), 1);
  }

  emtfMuonhwPt = ibooker.book1D("emtfMuonhwPt", "EMTF Muon Cand p_{T}", 512, 0, 512);
  emtfMuonhwPt->setAxisTitle("Hardware p_{T}", 1);

  emtfMuonhwPtCoarse = ibooker.book1D("emtfMuonhwPtCoarse", "EMTF Muon Cand p_{T} (coarse)", 32, 0, 512);
  emtfMuonhwPtCoarse->setAxisTitle("Hardware p_{T}", 1);
  //emtfMuonhwPtCoarse->getTH1F()->Sumw2();


  emtfMuonhwEta = ibooker.book1D("emtfMuonhwEta", "EMTF Muon Cand #eta", 460, -230, 230);
  emtfMuonhwEta->setAxisTitle("Hardware #eta", 1);
  //emtfMuonhwEta->getTH1F()->Sumw2();

  emtfMuonhwEtaCoarse = ibooker.book1D("emtfMuonhwEtaCoarse", "EMTF Muon Cand #eta (coarse)", 46, -230, 230);
  emtfMuonhwEtaCoarse->setAxisTitle("Hardware #eta", 1);
  //emtfMuonhwEtaCoarse->getTH1F()->Sumw2();


  emtfMuonhwPhi = ibooker.book1D("emtfMuonhwPhi", "EMTF Muon Cand #phi", 125, -20, 105);
  emtfMuonhwPhi->setAxisTitle("Hardware #phi", 1);

  emtfMuonhwPhiCoarse = ibooker.book1D("emtfMuonhwPhiCoarse", "EMTF Muon Cand #phi (coarse)", 25, -20, 105);
  emtfMuonhwPhiCoarse->setAxisTitle("Hardware #phi", 1);
  //emtfMuonhwPhiCoarse->getTH1F()->Sumw2();

  emtfMuonhwQual = ibooker.book1D("emtfMuonhwQual", "EMTF Muon Cand Quality", 16, 0, 16);
  emtfMuonhwQual->setAxisTitle("Quality", 1);
  for (int bin = 1; bin <= 16; ++bin) {
    emtfMuonhwQual->setBinLabel(bin, std::to_string(bin - 1), 1);
  }

}

void L1TStage2EMTF::analyze(const edm::Event& e, const edm::EventSetup& c) {

  if (verbose) edm::LogInfo("L1TStage2EMTF") << "L1TStage2EMTF: analyze..." << std::endl;
  if (isData){
  // DAQ Output
    edm::Handle<l1t::EMTFDaqOutCollection> DaqOutCollection;
    e.getByToken(daqToken, DaqOutCollection);
  
    for (std::vector<l1t::EMTFDaqOut>::const_iterator DaqOut = DaqOutCollection->begin(); DaqOut != DaqOutCollection->end(); ++DaqOut) {
      const l1t::emtf::MECollection* MECollection = DaqOut->PtrMECollection();
      for (std::vector<l1t::emtf::ME>::const_iterator ME = MECollection->begin(); ME != MECollection->end(); ++ME) {
        if (ME->SE()) emtfErrors->Fill(1);
        if (ME->SM()) emtfErrors->Fill(2);
        if (ME->BXE()) emtfErrors->Fill(3);
        if (ME->AF()) emtfErrors->Fill(4);
      }
  
      const l1t::emtf::EventHeader* EventHeader = DaqOut->PtrEventHeader();
      if (!EventHeader->Rdy()) emtfErrors->Fill(5);
    }
  }
  // Hits (LCTs)
  edm::Handle<l1t::EMTFHitCollection> HitCollection;
  e.getByToken(hitToken, HitCollection);

  for (std::vector<l1t::EMTFHit>::const_iterator Hit = HitCollection->begin(); Hit != HitCollection->end(); ++Hit) {
    if (filterBX && (Hit->BX() > 1 || Hit->BX() < -1)) continue;//restricts BX to -1,0,1 for emulator comparisons, filterBX only true when emulator is run
    int endcap = Hit->Endcap();
    int sector = Hit->Sector();
    int station = Hit->Station();
    int ring = Hit->Ring();
    int cscid = Hit->CSC_ID();
    int chamber = Hit->Chamber();
    int strip = Hit->Strip();
    int wire = Hit->Wire();

    int neighbor = Hit->Neighbor(); //prevent overlap for QT

    int hist_index = 0;
    int cscid_offset = (sector - 1) * 9;

    // The following logic determines the index of the monitor element
    // to which a hit belongs, exploiting the symmetry of the endcaps.
    if (station == 1) {
      if (ring == 4) {
        hist_index = 8;
      } else if (ring == 1) {
        hist_index = 9;
      } else if (ring == 2) {
        hist_index = 7;
      } else if (ring == 3) {
        hist_index = 6;
      }
    } else if (ring == 1) {
      if (station == 2) {
        hist_index = 5;
      } else if (station == 3) {
        hist_index = 3;
      } else if (station == 4) {
        hist_index = 1;
      }
    } else if (ring == 2) {
      if (station == 2) {
        hist_index = 4;
      } else if (station == 3) {
        hist_index = 2;
      } else if (station == 4) {
        hist_index = 0;
      }
    }

    if (endcap > 0) hist_index = 19 - hist_index;

    emtfHitBX->Fill(Hit->BX(), hist_index);

    emtfHitStrip[hist_index]->Fill(strip);
    emtfHitWire[hist_index]->Fill(wire);

    emtfChamberStrip[hist_index]->Fill(chamber, strip);
    emtfChamberWire[hist_index]->Fill(chamber, wire);

//=====================QUALITY TESTER==========================================================================
//====================CHAMBER STRIP============================================================================

    Double_t sumChamberStrip[hist_index];
    Double_t meanChamberStrip[hist_index];
    Double_t sumChamberStrip1D[hist_index];
    Double_t meanChamberStrip1D[hist_index];

    if (neighbor == 0) {
        emtfChamberStrip_QT[hist_index]->Fill(chamber, strip);
        emtfChamberStrip_QT1D[hist_index]->Fill(chamber);
    }
	
/*    for (int binY=1, binX=1; binY <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsY(); binY++, binX=binY) {
        sumChamberStrip[hist_index] = 0;
        for (int binX=1; binX <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsX(); binX++) {
	    emtfChamberStrip_QT_sqrt[hist_index]->getTH1()->SetBinContent(binX, binY, sqrt(emtfChamberStrip_QT[hist_index]->getTH1()->GetBinContent(binX, binY)));

	    sumChamberStrip[hist_index] += emtfChamberStrip_QT_sqrt[hist_index]->getTH1()->GetBinContent(binX, binY);
	}
        emtfChamberStrip_QT_mean[hist_index]->getTH1()->SetBinContent(binX, sumChamberStrip[hist_index]/(emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsX()));
    } */

    for (int binX=1; binX <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsX(); binX++) {
          for (int binY=1; binY <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsY(); binY++) {	
	    emtfChamberStrip_QT_sqrt[hist_index]->getTH1()->SetBinContent(binX, binY, sqrt(emtfChamberStrip_QT[hist_index]->getTH1()->GetBinContent(binX, binY)));

	    sumChamberStrip[hist_index] += emtfChamberStrip_QT_sqrt[hist_index]->getTH1()->GetBinContent(binX, binY);
      }
    }

    for (int binX=1; binX <= emtfChamberStrip_QT1D[hist_index]->getTH1()->GetNbinsX(); binX++) {	
	    emtfChamberStrip_QT_sqrt1D[hist_index]->getTH1()->SetBinContent(binX, sqrt(emtfChamberStrip_QT1D[hist_index]->getTH1()->GetBinContent(binX)));

	    sumChamberStrip1D[hist_index] += emtfChamberStrip_QT_sqrt1D[hist_index]->getTH1()->GetBinContent(binX);
    }

    //find average of each histogram to reset average to 1
    meanChamberStrip[hist_index] = sumChamberStrip[hist_index] /((emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsX())*(emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsY()));

    meanChamberStrip1D[hist_index] = sumChamberStrip1D[hist_index] /(emtfChamberStrip_QT1D[hist_index]->getTH1()->GetNbinsX());

/*    for (int binX=1; binX <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsX(); binX++) {
        for (int binY=1; binY <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsY(); binY++) {
	  emtfChamberStrip_QT_hot[hist_index]->getTH1()->SetBinContent(binX, binY, 1 - emtfChamberStrip_QT_mean[hist_index]->getTH1()->GetBinContent(binX) + emtfChamberStrip_QT_sqrt[hist_index]->getTH1()->GetBinContent(binX, binY));
	}
    } */

    for (int binX=1; binX <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsX(); binX++) {
      for (int binY=1; binY <= emtfChamberStrip_QT[hist_index]->getTH1()->GetNbinsY(); binY++) {
	emtfChamberStrip_QT_hot[hist_index]->getTH1()->SetBinContent(binX, binY, 1 - meanChamberStrip[hist_index] + emtfChamberStrip_QT_sqrt[hist_index]->getTH1()->GetBinContent(binX, binY));
      }
    }

    for (int binX=1; binX <= emtfChamberStrip_QT1D[hist_index]->getTH1()->GetNbinsX(); binX++) {
	emtfChamberStrip_QT_dead[hist_index]->getTH1()->SetBinContent(binX, 1 + meanChamberStrip1D[hist_index] - emtfChamberStrip_QT_sqrt1D[hist_index]->getTH1()->GetBinContent(binX));
    }

//=============================================================================================================
//=============================================================================================================

    if (Hit->Subsector() == 1) {
      emtfHitOccupancy->Fill(cscid + cscid_offset, endcap * (station - 0.5));
    } else {
      emtfHitOccupancy->Fill(cscid + cscid_offset, endcap * (station + 0.5));
    }
  }

  // Tracks
  edm::Handle<l1t::EMTFTrackCollection> TrackCollection;
  e.getByToken(trackToken, TrackCollection);

  int nTracks = TrackCollection->size();

  if (nTracks <= 10) {
    emtfnTracks->Fill(nTracks);
  } else {
    emtfnTracks->Fill(10);
  }

  for (std::vector<l1t::EMTFTrack>::const_iterator Track = TrackCollection->begin(); Track != TrackCollection->end(); ++Track) { 
    if ( filterBX && (Track->BX() > 1 || Track->BX() < -1)) continue;
    int endcap = Track->Endcap();
    int sector = Track->Sector();
    float eta = Track->Eta();
    float phi_glob_rad = Track->Phi_glob_rad();;
    int mode = Track->Mode();
    int quality = Track->Quality();
    int BX = Track->BX();

    int has_neighbor = Track->Has_neighbor(); //for quality tests

    emtfTracknHits->Fill(Track->NumHits());
    emtfTrackBX->Fill(endcap * (sector - 0.5), BX);
    emtfTrackBX1D->Fill(BX);
//==================QUALITY TESTER=============================================================================
//===================TRACK BX==================================================================================
    Double_t sumTrackBX = 0;
    Double_t meanTrackBX;

    if (has_neighbor == 0) {
      if (Track->BX() == 0) {
        emtfTrackBX_QT->Fill(endcap * (sector - 0.5));
      }
    }

    for (int binX=1; binX <= emtfTrackBX_QT->getTH1()->GetNbinsX(); binX++) 
    {
	emtfTrackBX_QT_sqrt->getTH1()->SetBinContent(binX, sqrt(emtfTrackBX_QT->getTH1()->GetBinContent(binX)));	
	sumTrackBX += emtfTrackBX_QT_sqrt->getTH1()->GetBinContent(binX);
    }

    meanTrackBX = sumTrackBX / emtfTrackBX_QT->getTH1()->GetNbinsX();
    for (int binX=1; binX <= emtfTrackBX_QT->getTH1()->GetNbinsX(); binX++) {
	emtfTrackBX_QT_hot->getTH1()->SetBinContent(binX, 1 - meanTrackBX + emtfTrackBX_QT_sqrt->getTH1()->GetBinContent(binX));
	emtfTrackBX_QT_dead->getTH1()->SetBinContent(binX, 1 + meanTrackBX - emtfTrackBX_QT_sqrt->getTH1()->GetBinContent(binX));
    }

//=============================================================================================================
//=============================================================================================================
    emtfTrackPt->Fill(Track->Pt());
    emtfTrackPtCoarse->Fill(Track->Pt());
    emtfTrackEta->Fill(eta);
    emtfTrackEtaCoarse->Fill(eta);
    emtfTrackPhi->Fill(phi_glob_rad);
//================QUALITY TESTER===============================================================================
//================TRACK PHI====================================================================================

    Double_t shift = pi/64.0;
    Double_t sumPhiHighQuality = 0;
    Double_t meanPhiHighQuality;
    Double_t sumPhiCoarse = 0;
    Double_t meanPhiCoarse;

    if (has_neighbor == 0) 
    {
	emtfTrackPhiCoarse_QT->Fill(phi_glob_rad + shift);
        if (mode == 15) emtfTrackPhiHighQuality_QT->Fill(phi_glob_rad + shift);
    }

    for (int binX=1; binX <= emtfTrackPhiCoarse_QT->getTH1()->GetNbinsX(); binX++) 
    {
	emtfTrackPhiCoarse_QT_sqrt->getTH1()->SetBinContent(binX, sqrt(emtfTrackPhiCoarse_QT->getTH1()->GetBinContent(binX)));	
	sumPhiCoarse += emtfTrackPhiCoarse_QT_sqrt->getTH1()->GetBinContent(binX);
    }

    meanPhiCoarse = sumPhiCoarse / emtfTrackPhiCoarse_QT->getTH1()->GetNbinsX();
    for (int binX=1; binX <= emtfTrackPhiCoarse_QT->getTH1()->GetNbinsX(); binX++) {
	emtfTrackPhiCoarse_QT_hot->getTH1()->SetBinContent(binX, 1 - meanPhiCoarse + emtfTrackPhiCoarse_QT_sqrt->getTH1()->GetBinContent(binX));
	emtfTrackPhiCoarse_QT_dead->getTH1()->SetBinContent(binX, 1 + meanPhiCoarse - emtfTrackPhiCoarse_QT_sqrt->getTH1()->GetBinContent(binX));
    }
    
    for (int binX=1; binX <= emtfTrackPhiHighQuality_QT->getTH1()->GetNbinsX(); binX++) 
    {
	emtfTrackPhiHighQuality_QT_sqrt->getTH1()->SetBinContent(binX, sqrt(emtfTrackPhiHighQuality_QT->getTH1()->GetBinContent(binX)));
	sumPhiHighQuality += emtfTrackPhiHighQuality_QT_sqrt->getTH1()->GetBinContent(binX); 	
    }

    meanPhiHighQuality = sumPhiHighQuality / emtfTrackPhiHighQuality_QT->getTH1()->GetNbinsX();
    for (int binX=1; binX <= emtfTrackPhiHighQuality_QT_sqrt->getTH1()->GetNbinsX(); binX++) {
	emtfTrackPhiHighQuality_QT_hot->getTH1()->SetBinContent(binX, 1 - meanPhiHighQuality + emtfTrackPhiHighQuality_QT_sqrt->getTH1()->GetBinContent(binX));
	emtfTrackPhiHighQuality_QT_dead->getTH1()->SetBinContent(binX, 1 + meanPhiHighQuality - emtfTrackPhiHighQuality_QT_sqrt->getTH1()->GetBinContent(binX));
    }
//=============================================================================================================
//=============================================================================================================
    emtfTrackPhiCoarse->Fill(phi_glob_rad);
    emtfTrackOccupancy->Fill(eta, phi_glob_rad);
    emtfTrackMode->Fill(mode);
    emtfTrackQuality->Fill(quality);
    emtfTrackQualityVsMode->Fill(mode, quality);
    emtfTrackSectorIndex->Fill(endcap*sector);
    if (mode == 15) emtfTrackPhiHighQuality->Fill(phi_glob_rad);

  } //end of track collection

  // Regional Muon Candidates
  edm::Handle<l1t::RegionalMuonCandBxCollection> MuonBxCollection;
  e.getByToken(muonToken, MuonBxCollection);

  for (int itBX = MuonBxCollection->getFirstBX(); itBX <= MuonBxCollection->getLastBX(); ++itBX){
    if (filterBX && (itBX > 1 || itBX < -1)) continue;
    for (l1t::RegionalMuonCandBxCollection::const_iterator Muon = MuonBxCollection->begin(itBX); Muon != MuonBxCollection->end(itBX); ++Muon){
      emtfMuonBX->Fill(itBX);
      emtfMuonhwPt->Fill(Muon->hwPt());
      emtfMuonhwPtCoarse->Fill(Muon->hwPt());
      emtfMuonhwEta->Fill(Muon->hwEta());
      emtfMuonhwEtaCoarse->Fill(Muon->hwEta());
      emtfMuonhwPhi->Fill(Muon->hwPhi());
      emtfMuonhwPhiCoarse->Fill(Muon->hwPhi());
      emtfMuonhwQual->Fill(Muon->hwQual());
    }
  }
}




