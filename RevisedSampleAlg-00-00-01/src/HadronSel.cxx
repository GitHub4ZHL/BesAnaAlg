//###########################################//
//                                           //
//    Code for hadronic events selection     //
//             10th version                  //
//                                           //
//         WANG Weiping (at USTC)            //
//             Jul. 05th, 2018               //
//                                           //
//                                           //
//                                           //
//                                           //
//###########################################//

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "McDecayModeSvc/McDecayModeSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "TrigEvent/TrigEvent.h"
#include "TrigEvent/TrigData.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "EvtRecEvent/EvtRecEtaToGG.h"
#include "DstEvent/TofHitStatus.h"

#include "McTruth/McParticle.h"

#include "VertexFit/Helix.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "ParticleID/ParticleID.h"

#include "TMath.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "DTagTool/DTagTool.h"
#include "HadronSelAlg/HadronSel.h"
#include "RscanDQ/RscanDQ.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"
#include "TString.h"
#include <cstdlib>
#include "TList.h"

typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<double> Vdb;
typedef Vdb::iterator VdbIt ;
typedef std::vector<HepLorentzVector> Vp4;

int NcutIni=0, NcutTot=0, NcutVetobb=0, NcutNGoodgeZero=0, NcutNGoodgeTwo=0, NonDDEvts=0, NcutFinal=0;
int VetoDoubleCount=0;
int NcutNGoodeqZero=0, NcutNGoodeqOne=0, NcutNGoodeqTwo=0, NcutNGoodeqThree=0, NcutNGoodgtThree=0;
int NcutNGoodeqZero_EvtECut=0, NcutNGoodeqZero_NPi0Cut=0, NcutNGoodeqZero_NIsogamCut=0;
int NcutNGoodeqOne_EvtECut=0, NcutNGoodeqOne_NPi0Cut=0, NcutNGoodeqOne_BalanceCut=0;
int NcutNGoodeqTwo_AngCut=0, NcutNGoodeqTwo_Vetobb=0;
int NcutNGoodeqThree_AngCut=0, NcutNGoodeqThree_EOPCut=0, NcutNGoodeqThree_Vetobb=0;

const double velc = 299.792458;
const double mpion0 = 0.1349766;
const double mD0 = 1.86483;  // +/- 0.05
const double mDp = 1.86958;  // +/- 0.09
const double mDs = 1.96827;  // +/- 0.10
const double mDStar0 = 2.00685;  // +/- 0.05
const double mDStarp = 2.01026;  // +/- 0.05
const double mDStars = 2.11210;  // +/- 0.40

const double xmass[5] = {0.000511,  0.105658, 0.139570, 0.493677, 0.938272};

const double me  = 0.000511;
const double mmu = 0.105658;
const double mpi = 0.13957;
const double mk  = 0.49368;
const double mp  = 0.93827;

///////////////////////////////////////////////////////////////

HadronSel::HadronSel(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator){

		declareProperty("m_RealData",      m_RealData = false);  // The switch of the RealData
		declareProperty("m_QQbarMC",       m_QQbarMC = false);  // The switch of the QQbarMC
		declareProperty("m_QEDBKGMC",      m_QEDBKGMC = false);  // The switch of the QEDBKGMC
		declareProperty("m_InputEcm",      m_InputEcm = false);  // The switch of the Ecm assignment
		declareProperty("m_CheckBhaBha",   m_CheckBhaBha = false);  // The switch of the BhaBha MCtruth check
		declareProperty("m_CheckLUARLW",   m_CheckLUARLW = false);  // The switch of the MCtruth check
		declareProperty("m_CheckHybrid",   m_CheckHybrid = false);  // The switch of the MCtruth check
		declareProperty("m_CheckPhokha",   m_CheckPhokha = false);  // The switch of the MCtruth check
		declareProperty("m_CheckPythia",   m_CheckPythia = false);  // The switch of the MCtruth check
		declareProperty("m_JpsiTagCtrl",   m_JpsiTagCtrl = false);  // The switch of the Jpsi tag check
		declareProperty("m_ExcluChannel",  m_ExcluChannel = false); // The switch of tagging exclusive channel with specific m_ExcluID: false) no, true) yes.

		declareProperty("m_Ecm",           m_Ecm = 2.2324);  // The c.m. energy
		declareProperty("m_ExcluID",       m_ExcluID = 0);   // The switch of the [0] pi+pi-, [1] pi+pi-pi0, [2] pi+pi-2pi0, [3] 2(pi+pi-), [4] pi+pi-3pi0, [5] 2(pi+pi-)pi0,
		// [6] 2(pi+pi-pi0), [7] 3(pi+pi-), [8] pi+pi-4pi0, [9] 2(pi+pi-)3pi0, [10] K+K-, [11] K+K-pi+pi-
		declareProperty("m_Eratio",        m_Eratio = 0.65); // energy limit of second shower
		declareProperty("m_eedeg",         m_eedeg = 10.0);  // angle limit between two shower

		// Set Vr, Vz and costheta cut value //
		declareProperty("m_vr0cut",        m_vr0cut = 0.5);
		declareProperty("m_vz0cut",        m_vz0cut = 5.0);
		declareProperty("m_cosThetacut",   m_cosThetacut = 0.93);
		declareProperty("m_MomRatiocut",   m_MomRatiocut = 0.94);

		declareProperty("m_GammaBarrelEth",m_GammaBarrelEth = 0.025);
		declareProperty("m_GammaEndCapEth",m_GammaEndCapEth = 0.050);

		declareProperty("m_GamConEcut",    m_GamConEcut = 0.1);
		declareProperty("m_GamConAngcut",  m_GamConAngcut = 15);

		declareProperty("m_EOPcut",        m_EOPcut = 0.8);   // Energy over momentum cut for charged tracks and 3-prong events.
		declareProperty("m_PIDRatioCut",   m_PIDRatioCut = 0.25);  // PID ratio cut for charged tracks and 3-prong events.
		declareProperty("m_EBhabhaRatio",  m_EBhabhaRatio = 0.65); // momentum cut for E/p larger than 0.8 tracks
		declareProperty("m_dEdxChipCut",   m_dEdxChipCut = 10); // cut to remove non-collision tracks

		declareProperty("m_NGoodCut",      m_NGoodCut = 1);  // reserve the 1-prong events

		declareProperty("m_EVISratio",     m_EVISratio = 0.4);     // Visible energy cut
		declareProperty("m_IsogamECut",    m_IsogamECut = 0.1);    // cut of isolated photon energy
		declareProperty("m_IsogamAngCut",  m_IsogamAngCut = 20);   // cut of angle between isolated photon and charged track
		declareProperty("m_NIsogamCut",    m_NIsogamCut = 2);      // cut of number of isolated photon when Ngood = 2

		declareProperty("m_d2ThetaCut",    m_d2ThetaCut = 10.0);
		declareProperty("m_d2PhiCut",      m_d2PhiCut = 15.0);
		declareProperty("m_d3ThetaCut",    m_d3ThetaCut = 10.0);
		declareProperty("m_d3PhiCut",      m_d3PhiCut = 15.0);
		declareProperty("m_NlargeeopCut",  m_NlargeeopCut = 1);
		declareProperty("m_NlargeprobECut",m_NlargeprobECut = 1);

		declareProperty("m_NGood0EvtECut", m_NGood0EvtECut = 1.0); // Event energy cut for 0-prong events, the ratio w.r.t. Ebeam rather than Ecm.
		declareProperty("m_NGood0NPi0Cut", m_NGood0NPi0Cut = 2);    // The requirement of number of good pi0, applied for 0-prong events
		declareProperty("m_EvtERatio",     m_EvtERatio = 1.0);     // Event energy cut for 1-prong events, the ratio w.r.t. Ebeam rather than Ecm.
		declareProperty("m_Pi0LowMass",    m_Pi0LowMass = 0.110);  // The low limit of pi0 mass, applied for 1-prong events
		declareProperty("m_Pi0HigMass",    m_Pi0HigMass = 0.150);  // The high limit of pi0 mass, applied for 1-prong events
		declareProperty("m_NGoodPi0Cut",   m_NGoodPi0Cut = 1);     // The requirement of number of good pi0, applied for 1-prong events
		declareProperty("m_BalanceCut",    m_BalanceCut = 0.750);  // The high limit of event Balance, applied for 1-prong events

		declareProperty("m_MCTruthPi0Width",    m_MCTruthPi0Width = 2.5e-6);  // The width requirement of pi0 in MCtruth level
		declareProperty("m_MCTruthEtaWidth",    m_MCTruthEtaWidth = 5.0e-7);  // The width requirement of eta in MCtruth level

	}

///////////////////////////////////////////////////////////////////////////////////////////////
StatusCode HadronSel::initialize(){
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;

	// MC truth
	NTuplePtr nt0(ntupleSvc(), "FILE1/MCtruth");
	if ( nt0 ) m_tuple0 = nt0;
	else {
		m_tuple0 = ntupleSvc()->book ("FILE1/MCtruth", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple0 ){

			status = m_tuple0->addItem("NchMC",             m_NchMC);
			status = m_tuple0->addItem("NtyMC",             m_NtyMC);
			status = m_tuple0->addItem("JpsiTagMC",         m_JpsiTagMC);
			status = m_tuple0->addItem("ExcluTagMC",        m_ExcluTagMC); // Tag of Exclusive channels

			status = m_tuple0->addItem("indexmc",           m_idxmc,0,100);
			status = m_tuple0->addItem("pdgid",             m_idxmc,m_pdgidmc);
			status = m_tuple0->addItem("motheridx",         m_idxmc,m_motheridxmc);

			status = m_tuple0->addItem("npipMC",            m_npipMC,0,100);
			status = m_tuple0->addIndexedItem("pippxMC",    m_npipMC,m_pippxMC);
			status = m_tuple0->addIndexedItem("pippyMC",    m_npipMC,m_pippyMC);
			status = m_tuple0->addIndexedItem("pippzMC",    m_npipMC,m_pippzMC);
			status = m_tuple0->addIndexedItem("pipEnMC",    m_npipMC,m_pipEnMC);

			status = m_tuple0->addItem("npimMC",            m_npimMC,0,100);
			status = m_tuple0->addIndexedItem("pimpxMC",    m_npimMC,m_pimpxMC);
			status = m_tuple0->addIndexedItem("pimpyMC",    m_npimMC,m_pimpyMC);
			status = m_tuple0->addIndexedItem("pimpzMC",    m_npimMC,m_pimpzMC);
			status = m_tuple0->addIndexedItem("pimEnMC",    m_npimMC,m_pimEnMC);

			status = m_tuple0->addItem("nkpMC",             m_nkpMC,0,100);
			status = m_tuple0->addIndexedItem("kppxMC",     m_nkpMC,m_kppxMC);
			status = m_tuple0->addIndexedItem("kppyMC",     m_nkpMC,m_kppyMC);
			status = m_tuple0->addIndexedItem("kppzMC",     m_nkpMC,m_kppzMC);
			status = m_tuple0->addIndexedItem("kpEnMC",     m_nkpMC,m_kpEnMC);

			status = m_tuple0->addItem("nkmMC",             m_nkmMC,0,100);
			status = m_tuple0->addIndexedItem("kmpxMC",     m_nkmMC,m_kmpxMC);
			status = m_tuple0->addIndexedItem("kmpyMC",     m_nkmMC,m_kmpyMC);
			status = m_tuple0->addIndexedItem("kmpzMC",     m_nkmMC,m_kmpzMC);
			status = m_tuple0->addIndexedItem("kmEnMC",     m_nkmMC,m_kmEnMC);

			status = m_tuple0->addItem("nppMC",             m_nppMC,0,100);
			status = m_tuple0->addIndexedItem("pppxMC",     m_nppMC,m_pppxMC);
			status = m_tuple0->addIndexedItem("pppyMC",     m_nppMC,m_pppyMC);
			status = m_tuple0->addIndexedItem("pppzMC",     m_nppMC,m_pppzMC);
			status = m_tuple0->addIndexedItem("ppEnMC",     m_nppMC,m_ppEnMC);

			status = m_tuple0->addItem("npmMC",             m_npmMC,0,100);
			status = m_tuple0->addIndexedItem("pmpxMC",     m_npmMC,m_pmpxMC);
			status = m_tuple0->addIndexedItem("pmpyMC",     m_npmMC,m_pmpyMC);
			status = m_tuple0->addIndexedItem("pmpzMC",     m_npmMC,m_pmpzMC);
			status = m_tuple0->addIndexedItem("pmEnMC",     m_npmMC,m_pmEnMC);

			status = m_tuple0->addItem("ngamMC",            m_ngamMC);
			status = m_tuple0->addItem("nn0MC",             m_nn0MC);
			status = m_tuple0->addItem("nantin0MC",         m_nantin0MC);
			status = m_tuple0->addItem("nklongMC",          m_nklongMC);
			status = m_tuple0->addItem("ngaminiMC",         m_ngaminiMC);

			status = m_tuple0->addItem("npi0MC",            m_npi0MC,0,100);
			status = m_tuple0->addIndexedItem("pi0pxMC",    m_npi0MC,m_pi0pxMC);
			status = m_tuple0->addIndexedItem("pi0pyMC",    m_npi0MC,m_pi0pyMC);
			status = m_tuple0->addIndexedItem("pi0pzMC",    m_npi0MC,m_pi0pzMC);
			status = m_tuple0->addIndexedItem("pi0EnMC",    m_npi0MC,m_pi0EnMC);
			status = m_tuple0->addIndexedItem("pi0MaMC",    m_npi0MC,m_pi0MaMC);
			status = m_tuple0->addIndexedItem("pi0IdMC",    m_npi0MC,m_pi0IdMC);

			status = m_tuple0->addItem("netaMC",            m_netaMC,0,100);
			status = m_tuple0->addIndexedItem("etapxMC",    m_netaMC,m_etapxMC);
			status = m_tuple0->addIndexedItem("etapyMC",    m_netaMC,m_etapyMC);
			status = m_tuple0->addIndexedItem("etapzMC",    m_netaMC,m_etapzMC);
			status = m_tuple0->addIndexedItem("etaEnMC",    m_netaMC,m_etaEnMC);
			status = m_tuple0->addIndexedItem("etaMaMC",    m_netaMC,m_etaMaMC);
			status = m_tuple0->addIndexedItem("etaIdMC",    m_netaMC,m_etaIdMC);

			status = m_tuple0->addItem("nkstMC",            m_nkstMC,0,100);
			status = m_tuple0->addIndexedItem("kstpxMC",    m_nkstMC,m_kstpxMC);
			status = m_tuple0->addIndexedItem("kstpyMC",    m_nkstMC,m_kstpyMC);
			status = m_tuple0->addIndexedItem("kstpzMC",    m_nkstMC,m_kstpzMC);
			status = m_tuple0->addIndexedItem("kstEnMC",    m_nkstMC,m_kstEnMC);
			status = m_tuple0->addIndexedItem("kstMaMC",    m_nkstMC,m_kstMaMC);
			status = m_tuple0->addIndexedItem("kstIdMC",    m_nkstMC,m_kstIdMC);

			status = m_tuple0->addItem("nkstar0MC",         m_nkstar0MC,0,100);
			status = m_tuple0->addIndexedItem("kstar0pxMC", m_nkstar0MC,m_kstar0pxMC);
			status = m_tuple0->addIndexedItem("kstar0pyMC", m_nkstar0MC,m_kstar0pyMC);
			status = m_tuple0->addIndexedItem("kstar0pzMC", m_nkstar0MC,m_kstar0pzMC);
			status = m_tuple0->addIndexedItem("kstar0EnMC", m_nkstar0MC,m_kstar0EnMC);
			status = m_tuple0->addIndexedItem("kstar0MaMC", m_nkstar0MC,m_kstar0MaMC);
			status = m_tuple0->addIndexedItem("kstar0IdMC", m_nkstar0MC,m_kstar0IdMC);

			status = m_tuple0->addItem("nkstar0barMC",         m_nkstar0barMC,0,100);
			status = m_tuple0->addIndexedItem("kstar0barpxMC", m_nkstar0barMC,m_kstar0barpxMC);
			status = m_tuple0->addIndexedItem("kstar0barpyMC", m_nkstar0barMC,m_kstar0barpyMC);
			status = m_tuple0->addIndexedItem("kstar0barpzMC", m_nkstar0barMC,m_kstar0barpzMC);
			status = m_tuple0->addIndexedItem("kstar0barEnMC", m_nkstar0barMC,m_kstar0barEnMC);
			status = m_tuple0->addIndexedItem("kstar0barMaMC", m_nkstar0barMC,m_kstar0barMaMC);
			status = m_tuple0->addIndexedItem("kstar0barIdMC", m_nkstar0barMC,m_kstar0barIdMC);

			status = m_tuple0->addItem("nkstarpMC",         m_nkstarpMC,0,100);
			status = m_tuple0->addIndexedItem("kstarppxMC", m_nkstarpMC,m_kstarppxMC);
			status = m_tuple0->addIndexedItem("kstarppyMC", m_nkstarpMC,m_kstarppyMC);
			status = m_tuple0->addIndexedItem("kstarppzMC", m_nkstarpMC,m_kstarppzMC);
			status = m_tuple0->addIndexedItem("kstarpEnMC", m_nkstarpMC,m_kstarpEnMC);
			status = m_tuple0->addIndexedItem("kstarpMaMC", m_nkstarpMC,m_kstarpMaMC);
			status = m_tuple0->addIndexedItem("kstarpIdMC", m_nkstarpMC,m_kstarpIdMC);

			status = m_tuple0->addItem("nkstarmMC",         m_nkstarmMC,0,100); //kstarmbar
			status = m_tuple0->addIndexedItem("kstarmpxMC", m_nkstarmMC,m_kstarmpxMC);
			status = m_tuple0->addIndexedItem("kstarmpyMC", m_nkstarmMC,m_kstarmpyMC);
			status = m_tuple0->addIndexedItem("kstarmpzMC", m_nkstarmMC,m_kstarmpzMC);
			status = m_tuple0->addIndexedItem("kstarmEnMC", m_nkstarmMC,m_kstarmEnMC);
			status = m_tuple0->addIndexedItem("kstarmMaMC", m_nkstarmMC,m_kstarmMaMC);
			status = m_tuple0->addIndexedItem("kstarmIdMC", m_nkstarmMC,m_kstarmIdMC);

			status = m_tuple0->addItem("nChgMC",            m_nChgMC,0,100);
			status = m_tuple0->addIndexedItem("p3MC",       m_nChgMC,m_p3MC);
			status = m_tuple0->addIndexedItem("costhMC",    m_nChgMC,m_costhMC);
			status = m_tuple0->addIndexedItem("p3MCRest",   m_nChgMC,m_p3MCRest);
			status = m_tuple0->addIndexedItem("costhMCRest",m_nChgMC,m_costhMCRest);

			status = m_tuple0->addItem("nNeuMC",            m_nNeuMC,0,10000);

			status = m_tuple0->addItem("EcmsMC",            m_EcmsMC);
			status = m_tuple0->addItem("ElabMC",            m_ElabMC);
			status = m_tuple0->addItem("EisrMC",            m_EisrMC);
			status = m_tuple0->addItem("EeffMC",            m_EeffMC);
			status = m_tuple0->addItem("p3ElabMC",          m_p3ElabMC);

			status = m_tuple0->addItem("nisrMC",            m_nisrMC,0,20);
			status = m_tuple0->addIndexedItem("cosisrMC",   m_nisrMC,m_cosisrMC);
			status = m_tuple0->addIndexedItem("pzisrMC",    m_nisrMC,m_pzisrMC);
			status = m_tuple0->addIndexedItem("p3isrMC",    m_nisrMC,m_p3isrMC);
		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple0) << endmsg;
			return StatusCode::FAILURE;
		}
	}

	// detector level information
	NTuplePtr nt1(ntupleSvc(), "FILE1/Rinfo");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book("FILE1/Rinfo", CLID_ColumnWiseTuple, "track N-Tuple example");
		if ( m_tuple1 )    {

			status = m_tuple1->addItem("Nch",               m_Nch);
			status = m_tuple1->addItem("Nty",               m_Nty);
			status = m_tuple1->addItem("Eisr",              m_Eisr);
			status = m_tuple1->addItem("Ecms",              m_Ecms);
			status = m_tuple1->addItem("Eeff",              m_Eeff);
			status = m_tuple1->addItem("JpsiTag",           m_JpsiTag);
			status = m_tuple1->addItem("ExcluTag",          m_ExcluTag); // Double counting tag

			status = m_tuple1->addItem("indexmc",           m_idx,0,100);
			status = m_tuple1->addItem("pdgid",             m_idx,m_pdgid);
			status = m_tuple1->addItem("motheridx",         m_idx,m_motheridx);
			status = m_tuple1->addIndexedItem("TruthP4",	m_idx, 4, m_TruthP4);

			status = m_tuple1->addItem("runno",             m_runno);
			status = m_tuple1->addItem("event",             m_event);
			status = m_tuple1->addItem("ebeam",             m_ebeam);

			status = m_tuple1->addItem("dthe",              m_dthe);  // delta theta between two showers
			status = m_tuple1->addItem("eraw1",             m_eraw1); // energy of first shower
			status = m_tuple1->addItem("eraw2",             m_eraw2); // energy of second shower

			status = m_tuple1->addItem("ngam",              m_ngam,0,10000);
			status = m_tuple1->addIndexedItem("GamThe",     m_ngam,m_GamThe);
			status = m_tuple1->addIndexedItem("GamPhi",     m_ngam,m_GamPhi);
			status = m_tuple1->addIndexedItem("GamEne",     m_ngam,m_GamEne);

			status = m_tuple1->addItem("NIsogam",           m_NIsogam,0,1000);
			status = m_tuple1->addIndexedItem("IsoGamThe",  m_NIsogam,m_IsoGamThe);
			status = m_tuple1->addIndexedItem("IsoGamPhi",  m_NIsogam,m_IsoGamPhi);
			status = m_tuple1->addIndexedItem("IsoGamEne",  m_NIsogam,m_IsoGamEne);

			status = m_tuple1->addItem("Ngamcon",           m_Ngamcon); // Number of Gamma conversion
			status = m_tuple1->addItem("Egamcon",           m_Egamcon); // Energy of Gamma conversion

			status = m_tuple1->addItem("nInipi0",           m_nInipi0,0,100);
			status = m_tuple1->addIndexedItem("pi0px",      m_nInipi0,m_pi0px);
			status = m_tuple1->addIndexedItem("pi0py",      m_nInipi0,m_pi0py);
			status = m_tuple1->addIndexedItem("pi0pz",      m_nInipi0,m_pi0pz);
			status = m_tuple1->addIndexedItem("pi0En",      m_nInipi0,m_pi0En);

			status = m_tuple1->addItem("nIniGood",          m_IninGood);
			status = m_tuple1->addItem("nGood",             m_nGood,0,100);
			status = m_tuple1->addIndexedItem("EMCE",       m_nGood,m_EMCE); // Energy deposited in EMC of each good charged track
			status = m_tuple1->addIndexedItem("TrkP",       m_nGood,m_TrkP);
			status = m_tuple1->addIndexedItem("dEdx",       m_nGood,m_dEdx);
			status = m_tuple1->addIndexedItem("dEdxchip",   m_nGood,m_dEdxchip);
			status = m_tuple1->addIndexedItem("costheta",   m_nGood,m_costheta);
			status = m_tuple1->addIndexedItem("phi",        m_nGood,m_phi);
			status = m_tuple1->addIndexedItem("charge",     m_nGood,m_charge);
			status = m_tuple1->addIndexedItem("vr0",        m_nGood,m_vr0);
			status = m_tuple1->addIndexedItem("vz0",        m_nGood,m_vz0);
			status = m_tuple1->addIndexedItem("px_mdc",     m_nGood,m_px_mdc);
			status = m_tuple1->addIndexedItem("py_mdc",     m_nGood,m_py_mdc);
			status = m_tuple1->addIndexedItem("pz_mdc",     m_nGood,m_pz_mdc);
			status = m_tuple1->addIndexedItem("E_mdc",      m_nGood,m_E_mdc);
			status = m_tuple1->addIndexedItem("mass_pid",   m_nGood,m_mass_pid); // Particle identified of the good charged track
			status = m_tuple1->addIndexedItem("EOPRatio",   m_nGood,m_EOPRatio); // E(EMC energy) over P of the good charged track

			status = m_tuple1->addItem("N_large_eop",       m_N_large_eop);
			status = m_tuple1->addItem("N_large_probE",     m_N_large_probE);

			status = m_tuple1->addItem("EEmcNeu",           m_EEmcNeu);//EMC energy of all good photons
			status = m_tuple1->addItem("EEmcChg",           m_EEmcChg);//EMC energy of total charged tracks
			status = m_tuple1->addItem("TotVisE",           m_TotVisE);//Total visable energy of an event
			status = m_tuple1->addItem("TotEvtE",           m_TotEvtE);//Total estimated energy of an event, in which the charged tracks are regared as pion
			status = m_tuple1->addItem("Balance",           m_Balance);//The Balance of an event
			status = m_tuple1->addItem("avz",               m_avz);

			status = m_tuple1->addItem("d2ang",             m_d2ang); // angle between the two charged tracks
			status = m_tuple1->addItem("d2the",             m_d2the);
			status = m_tuple1->addItem("d2phi",             m_d2phi);

			status = m_tuple1->addItem("d3ang",             m_d3ang);
			status = m_tuple1->addItem("d3the",             m_d3the);
			status = m_tuple1->addItem("d3phi",             m_d3phi);

			status = m_tuple1->addItem("ETrkAftPID",        m_ETrkAftPID);  // Calculated energy of charged tacks with PID
			status = m_tuple1->addItem("npplus",            m_npplus);
			status = m_tuple1->addItem("npminus",           m_npminus);
			status = m_tuple1->addItem("nkplus",            m_nkplus);
			status = m_tuple1->addItem("nkminus",           m_nkminus);
			status = m_tuple1->addItem("npiplus",           m_npiplus);
			status = m_tuple1->addItem("npiminus",          m_npiminus);

			status = m_tuple1->addItem("nInikst",           m_nInikst,0,100);
			status = m_tuple1->addIndexedItem("kstpx",      m_nInikst,m_kstpx);
			status = m_tuple1->addIndexedItem("kstpy",      m_nInikst,m_kstpy);
			status = m_tuple1->addIndexedItem("kstpz",      m_nInikst,m_kstpz);
			status = m_tuple1->addIndexedItem("kstEn",      m_nInikst,m_kstEn);

			status = m_tuple1->addItem("nInikstar0",        m_nInikstar0,0,100);
			status = m_tuple1->addIndexedItem("kstar0px",   m_nInikstar0,m_kstar0px);
			status = m_tuple1->addIndexedItem("kstar0py",   m_nInikstar0,m_kstar0py);
			status = m_tuple1->addIndexedItem("kstar0pz",   m_nInikstar0,m_kstar0pz);
			status = m_tuple1->addIndexedItem("kstar0En",   m_nInikstar0,m_kstar0En);

			status = m_tuple1->addItem("nInikstar0bar",     m_nInikstar0bar,0,100);
			status = m_tuple1->addIndexedItem("kstar0barpx",m_nInikstar0bar,m_kstar0barpx);
			status = m_tuple1->addIndexedItem("kstar0barpy",m_nInikstar0bar,m_kstar0barpy);
			status = m_tuple1->addIndexedItem("kstar0barpz",m_nInikstar0bar,m_kstar0barpz);
			status = m_tuple1->addIndexedItem("kstar0barEn",m_nInikstar0bar,m_kstar0barEn);

			status = m_tuple1->addItem("nInikstarp",        m_nInikstarp,0,100);
			status = m_tuple1->addIndexedItem("kstarppx",   m_nInikstarp,m_kstarppx);
			status = m_tuple1->addIndexedItem("kstarppy",   m_nInikstarp,m_kstarppy);
			status = m_tuple1->addIndexedItem("kstarppz",   m_nInikstarp,m_kstarppz);
			status = m_tuple1->addIndexedItem("kstarpEn",   m_nInikstarp,m_kstarpEn);

			status = m_tuple1->addItem("nInikstarm",        m_nInikstarm,0,100);
			status = m_tuple1->addIndexedItem("kstarmpx",   m_nInikstarm,m_kstarmpx);
			status = m_tuple1->addIndexedItem("kstarmpy",   m_nInikstarm,m_kstarmpy);
			status = m_tuple1->addIndexedItem("kstarmpz",   m_nInikstarm,m_kstarmpz);
			status = m_tuple1->addIndexedItem("kstarmEn",   m_nInikstarm,m_kstarmEn);

			status = m_tuple1->addItem(			"npi0",			m_npi0, 0, 100000);
			status = m_tuple1->addIndexedItem(	"chi2pi0",		m_npi0,			m_chi2pi0);
			status = m_tuple1->addIndexedItem(	"pxpi0",		m_npi0,			m_pxpi0);
			status = m_tuple1->addIndexedItem(	"pypi0",		m_npi0,			m_pypi0);
			status = m_tuple1->addIndexedItem(	"pzpi0",		m_npi0,			m_pzpi0);
			status = m_tuple1->addIndexedItem(	"Epi0",			m_npi0,			m_Epi0);
			status = m_tuple1->addIndexedItem(	"idgam1",		m_npi0,			m_idgam1);
			status = m_tuple1->addIndexedItem(	"idgam2",		m_npi0,			m_idgam2);

			status = m_tuple1->addItem(			"nkap",     	m_nkap,			0,		100);
			status = m_tuple1->addIndexedItem(	"idkap",		m_nkap,			m_idkap);
			status = m_tuple1->addIndexedItem(	"mkap",			m_nkap,			m_mkap);
			status = m_tuple1->addIndexedItem(	"pxkap",		m_nkap,			m_pxkap);
			status = m_tuple1->addIndexedItem(	"pykap",		m_nkap,			m_pykap);
			status = m_tuple1->addIndexedItem(	"pzkap",		m_nkap,			m_pzkap);
			status = m_tuple1->addIndexedItem(	"Ekap",			m_nkap,			m_Ekap);
			status = m_tuple1->addIndexedItem(	"KapProbK",		m_nkap,			m_KapProbK);
                        status = m_tuple1->addIndexedItem(      "KapProbPi",            m_nkap,                 m_KapProbPi);
                        status = m_tuple1->addIndexedItem(      "KapProbP",             m_nkap,                 m_KapProbP);

			status = m_tuple1->addItem(			"nkam",     	m_nkam,			0,		100);
			status = m_tuple1->addIndexedItem(	"idkam",		m_nkam,			m_idkam);
			status = m_tuple1->addIndexedItem(	"mkam",			m_nkam,			m_mkam);
			status = m_tuple1->addIndexedItem(	"pxkam",		m_nkam,			m_pxkam);
			status = m_tuple1->addIndexedItem(	"pykam",		m_nkam,			m_pykam);
			status = m_tuple1->addIndexedItem(	"pzkam",		m_nkam,			m_pzkam);
			status = m_tuple1->addIndexedItem(	"Ekam",			m_nkam,			m_Ekam);
                        status = m_tuple1->addIndexedItem(      "KamProbK",             m_nkam,                 m_KamProbK);
                        status = m_tuple1->addIndexedItem(      "KamProbPi",            m_nkam,                 m_KamProbPi);
                        status = m_tuple1->addIndexedItem(      "KamProbP",             m_nkam,                 m_KamProbP);

			status = m_tuple1->addItem(			"npip",     	m_npip,			0,		100);
			status = m_tuple1->addIndexedItem(	"idpip",		m_npip,			m_idpip);
			status = m_tuple1->addIndexedItem(	"mpip",			m_npip,			m_mpip);
			status = m_tuple1->addIndexedItem(	"pxpip",		m_npip,			m_pxpip);
			status = m_tuple1->addIndexedItem(	"pypip",		m_npip,			m_pypip);
			status = m_tuple1->addIndexedItem(	"pzpip",		m_npip,			m_pzpip);
			status = m_tuple1->addIndexedItem(	"Epip",			m_npip,			m_Epip);
                        status = m_tuple1->addIndexedItem(      "PipProbK",             m_npip,                 m_PipProbK);
                        status = m_tuple1->addIndexedItem(      "PipProbPi",            m_npip,                 m_PipProbPi);
                        status = m_tuple1->addIndexedItem(      "PipProbP",             m_npip,                 m_PipProbP);

			status = m_tuple1->addItem(			"npim",     	m_npim,			0,		100);
			status = m_tuple1->addIndexedItem(	"idpim",		m_npim,			m_idpim);
			status = m_tuple1->addIndexedItem(	"mpim",			m_npim,			m_mpim);
			status = m_tuple1->addIndexedItem(	"pxpim",		m_npim,			m_pxpim);
			status = m_tuple1->addIndexedItem(	"pypim",		m_npim,			m_pypim);
			status = m_tuple1->addIndexedItem(	"pzpim",		m_npim,			m_pzpim);
			status = m_tuple1->addIndexedItem(	"Epim",			m_npim,			m_Epim);
                        status = m_tuple1->addIndexedItem(      "PimProbK",             m_npim,                 m_PimProbK);
                        status = m_tuple1->addIndexedItem(      "PimProbPi",            m_npim,                 m_PimProbPi);
                        status = m_tuple1->addIndexedItem(      "PimProbP",             m_npim,                 m_PimProbP);

			status = m_tuple1->addItem(			"nKS0",			m_nKS0, 0, 100);
			status = m_tuple1->addIndexedItem(	"pxKS0",		m_nKS0,			m_pxKS0);
			status = m_tuple1->addIndexedItem(	"pyKS0",		m_nKS0,			m_pyKS0);
			status = m_tuple1->addIndexedItem(	"pzKS0",		m_nKS0,			m_pzKS0);
			status = m_tuple1->addIndexedItem(	"EKS0",			m_nKS0,			m_EKS0);
			status = m_tuple1->addIndexedItem(	"mKS0",			m_nKS0,			m_mKS0);
			status = m_tuple1->addIndexedItem(	"chi2KS0IP",	m_nKS0,			m_chi2KS0IP);
			status = m_tuple1->addIndexedItem(	"chi2KS0SP",	m_nKS0,			m_chi2KS0SP);
			status = m_tuple1->addIndexedItem(	"LifeKS0",		m_nKS0,			m_LifeKS0);
			status = m_tuple1->addIndexedItem(	"LenKS0",		m_nKS0,			m_LenKS0);
			status = m_tuple1->addIndexedItem(	"LenErrKS0",	m_nKS0,			m_LenErrKS0);
			status = m_tuple1->addIndexedItem(	"IdTrp",		m_nKS0,			m_IdTrp);
			status = m_tuple1->addIndexedItem(	"pxTrp",		m_nKS0,			m_pxTrp);
			status = m_tuple1->addIndexedItem(	"pyTrp",		m_nKS0,			m_pyTrp);
			status = m_tuple1->addIndexedItem(	"pzTrp",		m_nKS0,			m_pzTrp);
			status = m_tuple1->addIndexedItem(	"ETrp",			m_nKS0,			m_ETrp);
			status = m_tuple1->addIndexedItem(	"IdTrm",		m_nKS0,			m_IdTrm);
			status = m_tuple1->addIndexedItem(	"pxTrm",		m_nKS0,			m_pxTrm);
			status = m_tuple1->addIndexedItem(	"pyTrm",		m_nKS0,			m_pyTrm);
			status = m_tuple1->addIndexedItem(	"pzTrm",		m_nKS0,			m_pzTrm);
			status = m_tuple1->addIndexedItem(	"ETrm",			m_nKS0,			m_ETrm);

		}
		else    { 
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}

	//////////////////////////////////--------end of book--------//////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;

}

//***********************************************************************
StatusCode HadronSel::execute() {

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	int runNumber = eventHeader->runNumber();  
	int event=eventHeader->eventNumber();

	int runNo = fabs(runNumber);

	if(runNo >= 47467 && runNo <= 47493) m_Ecm = 3.49;
	else if(runNo >= 51657 && runNo <= 51893) m_Ecm = 3.5080;
	else if(runNo >= 51584 && runNo <= 51656) m_Ecm = 3.5097;
	else if(runNo >= 51894 && runNo <= 52090) m_Ecm = 3.5104;
	else if(runNo >= 52298 && runNo <= 52332) m_Ecm = 3.5146;

	unsigned int Nch, Nty;
	int jpsitag=0;
	int exclutag=-1;

	Nch = eventHeader->flag1();
	Nty = eventHeader->flag2();

	m_Nch = Nch; m_NchMC = Nch;
	m_Nty = Nty; m_NtyMC = Nty;
	m_runno = runNumber;
	m_event = event;

	NcutIni++;
	RscanDQ rdq = RscanDQ(runNo);
	double Ebeam = rdq.getEbeam();
	if(m_InputEcm==false) { Ebeam = rdq.getEbeam(); }
	if(m_InputEcm==true) { Ebeam = m_Ecm/2.0; }
	m_ebeam = Ebeam;

	int Status = rdq.getStatus();
	if (m_RealData&&(Status==-1))  return StatusCode::SUCCESS;

	//  if (m_CheckHybrid&&Nty!=11) return StatusCode::SUCCESS; // topology analysis only for Luarlw component in Hybrid sample

	int Number = rdq.getNumber();
	//  if (Number == -13)  return StatusCode::SUCCESS;

	NcutTot++;
	//cout<<endl;
	//cout<<endl;
	//cout<<"Ebeam = "<<Ebeam<<", Status = "<<Status<<", Number = "<<Number<<", NcutTot = "<<NcutTot<<", Nch = "<<Nch<<", Nty = " <<Nty<< ", runNo = " <<runNo<< endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// ##### for MC truth information #####
	////

	if(eventHeader->runNumber()<0)
	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		if (!mcParticleCol)
		{
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}

		// record for topology information
		int m_numParticle = 1;
		int nprmary = 0;
		Event::McParticleCol::iterator iter_mc1 = mcParticleCol->begin();
		for(; iter_mc1 != mcParticleCol->end(); iter_mc1++)
		{
			if(!(*iter_mc1)->decayFromGenerator()) continue;
			if((*iter_mc1)->primaryParticle()) nprmary++;
		}

		Event::McParticleCol::iterator iter_mc2 = mcParticleCol->begin();
		if(nprmary==1)
		{
			m_numParticle = 0;
			for(; iter_mc2 != mcParticleCol->end(); iter_mc2++)
			{
				if(!(*iter_mc2)->decayFromGenerator()) continue;
				if((*iter_mc2)->primaryParticle())
				{
					m_pdgidmc[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridxmc[m_numParticle] = 0;
					m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridx[m_numParticle] = 0;
				}
				else
				{
					m_pdgidmc[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridxmc[m_numParticle] = ((*iter_mc2)->mother()).trackIndex();
					m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex();
				}
				m_numParticle += 1;
			}
			m_idxmc = m_numParticle;
			m_idx = m_numParticle;
		}
		if(nprmary>1)
		{
			m_numParticle = 1;
			for(; iter_mc2 != mcParticleCol->end(); iter_mc2++)
			{
				if(!(*iter_mc2)->decayFromGenerator()) continue;
				if((*iter_mc2)->primaryParticle())
				{
					m_pdgidmc[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridxmc[m_numParticle] = 0;
					m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridx[m_numParticle] = 0;
				}
				else
				{
					m_pdgidmc[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridxmc[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1;
					m_pdgid[m_numParticle] = (*iter_mc2)->particleProperty();
					m_motheridx[m_numParticle] = ((*iter_mc2)->mother()).trackIndex() + 1;
				}
				m_numParticle += 1;
				m_pdgidmc[0] = 11111;
				m_motheridxmc[0] = 0;
				m_pdgid[0] = 11111;
				m_motheridx[0] = 0;
			}
			m_idxmc = m_numParticle;
			m_idx = m_numParticle;
		}

		int nmc_pp=0, nmc_pip=0, nmc_kp=0;
		int nmc_pm=0, nmc_pim=0, nmc_km=0;

		int nmc_ep_ini=0, nmc_em_ini=0, nmc_mp_ini=0, nmc_mm_ini=0;
		int nmc_ep_pi0=0, nmc_em_pi0=0;
		int nmc_ep_oth=0, nmc_em_oth=0, nmc_mp_oth=0, nmc_mm_oth=0;

		int nmc_gamma=0, nmc_n0=0, nmc_antin0=0, nmc_klong=0;
		int nmc_pi0=0, nmc_eta=0, nmc_kst=0, nmc_kstar0=0, nmc_kstar0bar=0, nmc_kstarp=0, nmc_kstarm=0;

		int mmc_gamma_pi0=0, mmc_gamma_ini=0, mmc_gamma_oth=0;

		int nmc_nChg=0, nmc_nNeu=0, nmc_nisr=0;

		int iflagchg=0, iflagneu=0, iflagisr=0;

		int tempid = 0;
		HepLorentzVector mc_p4elab, mc_p4eisr;
		double mc_ecms=0.0, mc_elab=0.0, mc_eisr=0.0, mc_eeff=0.0, mc_p3elab=0.0;

		Vp4 pGamMC;
		pGamMC.clear();
		Vint iGamMC;
		iGamMC.clear();

		Event::McParticleCol::iterator iter_mc0 = mcParticleCol->begin();
		for(; iter_mc0 != mcParticleCol->end(); iter_mc0++) {
			if(!(*iter_mc0)->decayFromGenerator()) continue;

			HepLorentzVector mctrue_track0 = (*iter_mc0)->initialFourMomentum();

			if(m_CheckLUARLW==true)
			{
				if((*iter_mc0)->primaryParticle()){ mc_p4elab += mctrue_track0; }
				if((*iter_mc0)->primaryParticle()&&(*iter_mc0)->particleProperty()==-22)
				{ 
					mc_p4eisr += mctrue_track0;
					m_cosisrMC[iflagisr] = mctrue_track0.cosTheta();
					m_pzisrMC[iflagisr] = mctrue_track0.pz();
					m_p3isrMC[iflagisr] = mctrue_track0.rho();
					iflagisr++;
				}
			}

			if(m_CheckHybrid==true)
			{
				if((*iter_mc0)->primaryParticle()){ mc_p4elab += mctrue_track0; }
				if((*iter_mc0)->primaryParticle()&&((*iter_mc0)->particleProperty()==-22 || (*iter_mc0)->particleProperty()==22))
				{ 
					mc_p4eisr += mctrue_track0;
					m_cosisrMC[iflagisr] = mctrue_track0.cosTheta();
					m_p3isrMC[iflagisr] = mctrue_track0.rho();
					iflagisr++;
				}
			}

			if(m_CheckPythia==true)
			{
				if((*iter_mc0)->primaryParticle()&&(*iter_mc0)->particleProperty()==11&&((*iter_mc0)->mother()).particleProperty()==11) continue;
				if((*iter_mc0)->primaryParticle()){ mc_p4elab += mctrue_track0; }
				if((*iter_mc0)->primaryParticle()&&(*iter_mc0)->particleProperty()==22&&((*iter_mc0)->mother()).particleProperty()==22)
				{
					mc_p4eisr += mctrue_track0;
					m_cosisrMC[iflagisr] = mctrue_track0.cosTheta();
					m_p3isrMC[iflagisr] = mctrue_track0.rho();
					iflagisr++;
				}
			}
		}

		nmc_nisr = iflagisr;

		int ParNum = 0;
		double Ecms = 0.0;
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		for(; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if(!(*iter_mc)->decayFromGenerator()) continue;
			if(m_CheckPythia==true&&(*iter_mc)->primaryParticle()&&(*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()==11) continue;

			HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
			Ecms += mctrue_track.e();

			HepLorentzVector p4Truth = (*iter_mc)->initialFourMomentum();
			for(int ll = 0; ll < 4; ll++) m_TruthP4[ParNum][ll] = p4Truth[ll];
			ParNum = ParNum + 1;
			//cout << " event = " << event << ", primary = " << (*iter_mc)->primaryParticle() << ", iter_mc = " << tempid << ", PDGID = " << (*iter_mc)->particleProperty() << ", Mother PDGID = " << ((*iter_mc)->mother()).particleProperty() << ", ene = " << mctrue_track.e() << ", mass = " << mctrue_track.m() << endl;

			//      if((*iter_mc)->primaryParticle())
			//      {
			////        cout << " primary = " << setw(2) << (*iter_mc)->primaryParticle() << " it = " << setw(2) << tempid << ", ID = " << setw(6) << (*iter_mc)->particleProperty() << ", MID = " << setw(6) << ((*iter_mc)->mother()).particleProperty() << setiosflags(ios::fixed)<<setprecision(10) << ", ene = " << setw(12) << mctrue_track.e() << ", mass = " << setw(13) << mctrue_track.m() << ", cos = " << setw(13) << mctrue_track.cosTheta() << ", p3MC = " << setw(12) << mctrue_track.rho() << ", Ecms = " << setw(12) << Ecms << endl;
			//        cout << " primary = " << setw(2) << (*iter_mc)->primaryParticle() << " it = " << setw(2) << tempid << ", ID = " << setw(6) << (*iter_mc)->particleProperty() << ", MID = " << setw(6) << ((*iter_mc)->mother()).particleProperty() << setiosflags(ios::fixed)<<setprecision(10) << ", ene = " << setw(12) << mctrue_track.e() << ", mass = " << setw(13) << mctrue_track.m() << ", px = " << setw(13) << mctrue_track.px() << ", py = " << setw(13) << mctrue_track.py() << ", pz = " << setw(13) << mctrue_track.pz() << endl;
			//      }

			if((*iter_mc)->particleProperty()==22&&((*iter_mc)->mother()).particleProperty()==22) { mmc_gamma_ini++; }
			if((*iter_mc)->particleProperty()==22&&((*iter_mc)->mother()).particleProperty()==111){ mmc_gamma_pi0++; }
			if((*iter_mc)->particleProperty()==22&&((*iter_mc)->mother()).particleProperty()!=22&&((*iter_mc)->mother()).particleProperty()!=111) { mmc_gamma_oth++; }

			if((*iter_mc)->particleProperty()==22)
			{
				iGamMC.push_back(tempid);
				pGamMC.push_back(mctrue_track);

				nmc_nNeu++;
				nmc_gamma++;
			}

			if((*iter_mc)->particleProperty()==211)
			{
				nmc_nChg++;
				m_pippxMC[nmc_pip] = mctrue_track.px();
				m_pippyMC[nmc_pip] = mctrue_track.py();
				m_pippzMC[nmc_pip] = mctrue_track.pz();
				m_pipEnMC[nmc_pip] = mctrue_track.e();
				nmc_pip++;
			}

			if((*iter_mc)->particleProperty()==-211)
			{
				nmc_nChg++;
				m_pimpxMC[nmc_pim] = mctrue_track.px();
				m_pimpyMC[nmc_pim] = mctrue_track.py();
				m_pimpzMC[nmc_pim] = mctrue_track.pz();
				m_pimEnMC[nmc_pim] = mctrue_track.e();
				nmc_pim++;
			}

			if((*iter_mc)->particleProperty()==321)
			{
				nmc_nChg++;
				m_kppxMC[nmc_kp] = mctrue_track.px();
				m_kppyMC[nmc_kp] = mctrue_track.py();
				m_kppzMC[nmc_kp] = mctrue_track.pz();
				m_kpEnMC[nmc_kp] = mctrue_track.e();
				nmc_kp++;
			}

			if((*iter_mc)->particleProperty()==-321)
			{
				nmc_nChg++;
				m_kmpxMC[nmc_km] = mctrue_track.px();
				m_kmpyMC[nmc_km] = mctrue_track.py();
				m_kmpzMC[nmc_km] = mctrue_track.pz();
				m_kmEnMC[nmc_km] = mctrue_track.e();
				nmc_km++;
			}

			if((*iter_mc)->particleProperty()==2212)
			{
				nmc_nChg++;
				m_pppxMC[nmc_pp] = mctrue_track.px();
				m_pppyMC[nmc_pp] = mctrue_track.py();
				m_pppzMC[nmc_pp] = mctrue_track.pz();
				m_ppEnMC[nmc_pp] = mctrue_track.e();
				nmc_pp++;
			}

			if((*iter_mc)->particleProperty()==-2212)
			{
				nmc_nChg++;
				m_pmpxMC[nmc_pm] = mctrue_track.px();
				m_pmpyMC[nmc_pm] = mctrue_track.py();
				m_pmpzMC[nmc_pm] = mctrue_track.pz();
				m_pmEnMC[nmc_pm] = mctrue_track.e();
				nmc_pm++;
			}

			if((*iter_mc)->particleProperty()==-11&&((*iter_mc)->mother()).particleProperty()==-11){ nmc_nChg++; nmc_ep_ini++; }
			if((*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()==11){ nmc_nChg++; nmc_em_ini++; }

			if((*iter_mc)->particleProperty()==-11&&((*iter_mc)->mother()).particleProperty()==111){ nmc_nChg++; nmc_ep_pi0++; }
			if((*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()==111){ nmc_nChg++; nmc_em_pi0++; }

			if((*iter_mc)->particleProperty()==-11&&((*iter_mc)->mother()).particleProperty()!=-11&&((*iter_mc)->mother()).particleProperty()!=111){ nmc_nChg++; nmc_ep_oth++; }
			if((*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()!=11&&((*iter_mc)->mother()).particleProperty()!=111){ nmc_nChg++; nmc_em_oth++; }

			if((*iter_mc)->particleProperty()==-13&&((*iter_mc)->mother()).particleProperty()==-13){ nmc_nChg++; nmc_mp_ini++; }
			if((*iter_mc)->particleProperty()==13&&((*iter_mc)->mother()).particleProperty()==13){ nmc_nChg++; nmc_mm_ini++; }

			if((*iter_mc)->particleProperty()==-13&&((*iter_mc)->mother()).particleProperty()!=-13){ nmc_nChg++; nmc_mp_oth++; }
			if((*iter_mc)->particleProperty()==13&&((*iter_mc)->mother()).particleProperty()!=13){ nmc_nChg++; nmc_mm_oth++; }

			if((*iter_mc)->particleProperty()==2112){ nmc_nNeu++; nmc_n0++; }
			if((*iter_mc)->particleProperty()==-2112){ nmc_nNeu++; nmc_antin0++; }

			if((*iter_mc)->particleProperty()==111)
			{
				m_pi0pxMC[nmc_pi0] = mctrue_track.px();
				m_pi0pyMC[nmc_pi0] = mctrue_track.py();
				m_pi0pzMC[nmc_pi0] = mctrue_track.pz();
				m_pi0EnMC[nmc_pi0] = mctrue_track.e();
				m_pi0MaMC[nmc_pi0] = mctrue_track.m();

				m_pi0px[nmc_pi0] = mctrue_track.px();
				m_pi0py[nmc_pi0] = mctrue_track.py();
				m_pi0pz[nmc_pi0] = mctrue_track.pz();
				m_pi0En[nmc_pi0] = mctrue_track.e();

				if((*iter_mc)->primaryParticle()) m_pi0IdMC[nmc_pi0] = 1;
				if(!(*iter_mc)->primaryParticle()) m_pi0IdMC[nmc_pi0] = 0;
				nmc_pi0++;
			}

			if((*iter_mc)->particleProperty()==221)
			{
				m_etapxMC[nmc_eta] = mctrue_track.px();
				m_etapyMC[nmc_eta] = mctrue_track.py();
				m_etapzMC[nmc_eta] = mctrue_track.pz();
				m_etaEnMC[nmc_eta] = mctrue_track.e();
				m_etaMaMC[nmc_eta] = mctrue_track.m();

				if((*iter_mc)->primaryParticle()) m_etaIdMC[nmc_eta] = 1;
				if(!(*iter_mc)->primaryParticle()) m_etaIdMC[nmc_eta] = 0;
				nmc_eta++;
			}

			if((*iter_mc)->particleProperty()==310)
			{
				m_kstpxMC[nmc_kst] = mctrue_track.px();
				m_kstpyMC[nmc_kst] = mctrue_track.py();
				m_kstpzMC[nmc_kst] = mctrue_track.pz();
				m_kstEnMC[nmc_kst] = mctrue_track.e();
				m_kstMaMC[nmc_kst] = mctrue_track.m();

				m_kstpx[nmc_kst] = mctrue_track.px();
				m_kstpy[nmc_kst] = mctrue_track.py();
				m_kstpz[nmc_kst] = mctrue_track.pz();
				m_kstEn[nmc_kst] = mctrue_track.e();

				if((*iter_mc)->primaryParticle()) m_kstIdMC[nmc_kst] = 1;
				if(!(*iter_mc)->primaryParticle()) m_kstIdMC[nmc_kst] = 0;
				nmc_kst++;
			}

			if((*iter_mc)->particleProperty()==313)
			{
				m_kstar0pxMC[nmc_kstar0] = mctrue_track.px();
				m_kstar0pyMC[nmc_kstar0] = mctrue_track.py();
				m_kstar0pzMC[nmc_kstar0] = mctrue_track.pz();
				m_kstar0EnMC[nmc_kstar0] = mctrue_track.e();
				m_kstar0MaMC[nmc_kstar0] = mctrue_track.m();

				m_kstar0px[nmc_kstar0] = mctrue_track.px();
				m_kstar0py[nmc_kstar0] = mctrue_track.py();
				m_kstar0pz[nmc_kstar0] = mctrue_track.pz();
				m_kstar0En[nmc_kstar0] = mctrue_track.e();

				if((*iter_mc)->primaryParticle()) m_kstar0IdMC[nmc_kstar0] = 1;
				if(!(*iter_mc)->primaryParticle()) m_kstar0IdMC[nmc_kstar0] = 0;
				nmc_kstar0++;
			}

			if((*iter_mc)->particleProperty()==-313)
			{
				m_kstar0barpxMC[nmc_kstar0bar] = mctrue_track.px();
				m_kstar0barpyMC[nmc_kstar0bar] = mctrue_track.py();
				m_kstar0barpzMC[nmc_kstar0bar] = mctrue_track.pz();
				m_kstar0barEnMC[nmc_kstar0bar] = mctrue_track.e();
				m_kstar0barMaMC[nmc_kstar0bar] = mctrue_track.m();

				m_kstar0barpx[nmc_kstar0bar] = mctrue_track.px();
				m_kstar0barpy[nmc_kstar0bar] = mctrue_track.py();
				m_kstar0barpz[nmc_kstar0bar] = mctrue_track.pz();
				m_kstar0barEn[nmc_kstar0bar] = mctrue_track.e();

				if((*iter_mc)->primaryParticle()) m_kstar0barIdMC[nmc_kstar0bar] = 1;
				if(!(*iter_mc)->primaryParticle()) m_kstar0barIdMC[nmc_kstar0bar] = 0;
				nmc_kstar0bar++;
			}

			if((*iter_mc)->particleProperty()==323)
			{
				m_kstarppxMC[nmc_kstarp] = mctrue_track.px();
				m_kstarppyMC[nmc_kstarp] = mctrue_track.py();
				m_kstarppzMC[nmc_kstarp] = mctrue_track.pz();
				m_kstarpEnMC[nmc_kstarp] = mctrue_track.e();
				m_kstarpMaMC[nmc_kstarp] = mctrue_track.m();

				m_kstarppx[nmc_kstarp] = mctrue_track.px();
				m_kstarppy[nmc_kstarp] = mctrue_track.py();
				m_kstarppz[nmc_kstarp] = mctrue_track.pz();
				m_kstarpEn[nmc_kstarp] = mctrue_track.e();

				if((*iter_mc)->primaryParticle()) m_kstarpIdMC[nmc_kstarp] = 1;
				if(!(*iter_mc)->primaryParticle()) m_kstarpIdMC[nmc_kstarp] = 0;
				nmc_kstarp++;
			}

			if((*iter_mc)->particleProperty()==-323)
			{
				m_kstarmpxMC[nmc_kstarm] = mctrue_track.px();
				m_kstarmpyMC[nmc_kstarm] = mctrue_track.py();
				m_kstarmpzMC[nmc_kstarm] = mctrue_track.pz();
				m_kstarmEnMC[nmc_kstarm] = mctrue_track.e();
				m_kstarmMaMC[nmc_kstarm] = mctrue_track.m();

				m_kstarmpx[nmc_kstarm] = mctrue_track.px();
				m_kstarmpy[nmc_kstarm] = mctrue_track.py();
				m_kstarmpz[nmc_kstarm] = mctrue_track.pz();
				m_kstarmEn[nmc_kstarm] = mctrue_track.e();

				if((*iter_mc)->primaryParticle()) m_kstarmIdMC[nmc_kstarm] = 1;
				if(!(*iter_mc)->primaryParticle()) m_kstarmIdMC[nmc_kstarm] = 0;
				nmc_kstarm++;
			}

			if((*iter_mc)->particleProperty()==130){ nmc_klong++; }

			if((*iter_mc)->particleProperty()==443&&((*iter_mc)->mother()).particleProperty()==443){ jpsitag=1; }

			if(nmc_nChg==(iflagchg+1))
			{
				m_p3MC[iflagchg] = mctrue_track.rho();
				m_costhMC[iflagchg] = mctrue_track.cosTheta();

				// further boost to rest frame of virtual photon
				Hep3Vector virtualgam = -(mc_p4elab-mc_p4eisr).boostVector();
				mctrue_track = mctrue_track.boost(virtualgam);
				// store the results
				m_p3MCRest[iflagchg] = mctrue_track.rho();
				m_costhMCRest[iflagchg] = mctrue_track.cosTheta();

				iflagchg = nmc_nChg;
			}

			tempid++; 
		}

		m_JpsiTagMC = jpsitag;
		m_JpsiTag = jpsitag;

		//    cout << setiosflags(ios::fixed)<<setprecision(10) << "Total: ene = " << setw(12) << mc_p4elab.e() << ", mass = " << setw(13) << mc_p4elab.m() << ", px = " << setw(13) << mc_p4elab.px() << ", py = " << setw(13) << mc_p4elab.py() << ", pz = " << setw(13) << mc_p4elab.pz() << endl;

		mc_ecms = mc_p4elab.m();
		mc_elab = mc_p4elab.e();
		mc_eisr = mc_p4eisr.e();
		mc_eeff = (mc_p4elab-mc_p4eisr).m();
		mc_p3elab = mc_p4elab.rho();

		m_npipMC = nmc_pip;
		m_npimMC = nmc_pim;
		m_nkpMC = nmc_kp;
		m_nkmMC = nmc_km;
		m_nppMC = nmc_pp;
		m_npmMC = nmc_pm;
		m_ngamMC = nmc_gamma;
		m_ngaminiMC = mmc_gamma_ini;
		m_nn0MC = nmc_n0;
		m_nantin0MC = nmc_antin0;

		m_npi0MC = nmc_pi0;
		m_nInipi0 = nmc_pi0;
		m_netaMC = nmc_eta;
		m_nkstMC = nmc_kst;
		m_nInikst = nmc_kst;
		m_nklongMC = nmc_klong;

		m_nkstar0MC = nmc_kstar0;
		m_nkstar0barMC = nmc_kstar0bar;
		m_nkstarpMC = nmc_kstarp;
		m_nkstarmMC = nmc_kstarm;
		m_nInikstar0 = nmc_kstar0;
		m_nInikstar0bar = nmc_kstar0bar;
		m_nInikstarp = nmc_kstarp;
		m_nInikstarm = nmc_kstarm;

		m_nChgMC = nmc_nChg;
		m_IninGood = nmc_nChg;
		m_nNeuMC = nmc_nNeu;

		m_EcmsMC = mc_ecms; m_Ecms = mc_ecms;
		m_ElabMC = mc_elab;
		m_EisrMC = mc_eisr; m_Eisr = mc_eisr;
		m_nisrMC = nmc_nisr;
		m_EeffMC = mc_eeff; m_Eeff = mc_eeff;
		m_p3ElabMC = mc_p3elab;
		//cout<<"  ecms = "<<mc_ecms<<", eisr = "<<mc_eisr<<", eeff = "<<mc_eeff<<endl;
		//cout<<  setiosflags(ios::fixed)<<setprecision(12) <<"  labene = "<<m_ElabMC<<", cmsene = "<<m_EcmsMC<<", nisrMC = "<<m_nisrMC<<", eisr = "<<mc_eisr<<", eeff = "<<mc_eeff<<endl;

		if(m_CheckLUARLW==true)
		{
			if(nmc_pp==0&&nmc_pm==0&&nmc_n0==0&&nmc_antin0==0&&nmc_ep_ini==0&&nmc_em_ini==0&&nmc_mp_ini==0&&nmc_mm_ini==0&&nmc_ep_oth==0&&nmc_em_oth==0&&nmc_mp_oth==0&&nmc_mm_oth==0&&mmc_gamma_oth==0&&mmc_gamma_ini==0)
			{
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 0;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==1&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 1;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==2&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 2;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 3;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==3&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 4;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==1&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 5;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==2&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 6;
				if(nmc_pip==3&&nmc_pim==3&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 7;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==4&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 8;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==3&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 9;
				if(nmc_pip==0&&nmc_pim==0&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==1&&nmc_km==1) exclutag = 10;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==1&&nmc_km==1) exclutag = 11;
			}
		}

		if(m_CheckHybrid==true)
		{
			if(Nty!=11&&nmc_pp==0&&nmc_pm==0&&nmc_n0==0&&nmc_antin0==0&&nmc_ep_ini==0&&nmc_em_ini==0&&nmc_mp_ini==0&&nmc_mm_ini==0&&nmc_ep_oth==0&&nmc_em_oth==0&&nmc_mp_oth==0&&nmc_mm_oth==0&&mmc_gamma_oth==0)
			{
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=2) exclutag = 0;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==1&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=2) exclutag = 1;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==2&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=2) exclutag = 2;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=2) exclutag = 3;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==3&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=1) exclutag = 4;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==1&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=1) exclutag = 5;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==2&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=1) exclutag = 6;
				if(nmc_pip==3&&nmc_pim==3&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=1) exclutag = 7;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==4&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=1) exclutag = 8;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==3&&nmc_klong==0&&nmc_kp==0&&nmc_km==0&&mmc_gamma_ini<=1) exclutag = 9;
				if(nmc_pip==0&&nmc_pim==0&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==1&&nmc_km==1&&mmc_gamma_ini<=1) exclutag = 10;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==1&&nmc_km==1&&mmc_gamma_ini<=1) exclutag = 11;
			}

			if(Nty==11&&nmc_pp==0&&nmc_pm==0&&nmc_n0==0&&nmc_antin0==0&&nmc_ep_ini==0&&nmc_em_ini==0&&nmc_mp_ini==0&&nmc_mm_ini==0&&nmc_ep_oth==0&&nmc_em_oth==0&&nmc_mp_oth==0&&nmc_mm_oth==0&&mmc_gamma_oth==0&&mmc_gamma_ini==0)
			{
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 0;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==1&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 1;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==2&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 2;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 3;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==3&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 4;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==1&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 5;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==2&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 6;
				if(nmc_pip==3&&nmc_pim==3&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 7;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==4&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 8;
				if(nmc_pip==2&&nmc_pim==2&&nmc_pi0==3&&nmc_klong==0&&nmc_kp==0&&nmc_km==0) exclutag = 9;
				if(nmc_pip==0&&nmc_pim==0&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==1&&nmc_km==1) exclutag = 10;
				if(nmc_pip==1&&nmc_pim==1&&nmc_pi0==0&&nmc_klong==0&&nmc_kp==1&&nmc_km==1) exclutag = 11;
			}
		}

		m_ExcluTagMC = exclutag;
		m_ExcluTag = exclutag;

	}

	if(m_JpsiTagCtrl==true&&jpsitag!=1) return StatusCode::SUCCESS;
	if(m_ExcluChannel==true&&exclutag!=m_ExcluID) return StatusCode::SUCCESS; // tag specific exclusive channel with ID = exclutag

	VetoDoubleCount++;

	if(m_QQbarMC==true || m_QEDBKGMC==true || m_CheckLUARLW==true || m_CheckHybrid==true) m_tuple0->write(); // write the MCtruth for the Signal and background MC

	// finish the MC truth block
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
	if(!evtRecEvent)
	{
		log << MSG::FATAL << "Could not find EvtRecEvent" << endreq;
		return StatusCode::FAILURE;
	}

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
	if(!evtRecTrkCol)
	{
		log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
		return StatusCode::FAILURE;
	}

	Hep3Vector xorigin(0,0,0);

	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid())
	{
		double* dbv = vtxsvc->PrimaryVertex();
		double*  vv = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	//// ##### the 1st step -->  reject Bhabha ####
	//// #####      EMC information only       #### 

	double Efst = 0.0, Esec = 0.0;
	int ifst = -1, isec = -1;
	for(int i=0; i<evtRecEvent->totalTracks(); i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin()+i;
		if(!(*itTrk)->isEmcShowerValid()) continue;

		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		double eraw = emcTrk->energy();

		if(eraw>Efst)
		{
			isec = ifst;
			Esec = Efst;

			ifst = i;
			Efst = eraw;
		}
		else if(eraw>Esec)
		{
			isec = i;
			Esec = eraw;
		}
	}

	if(ifst!=-1&&isec!=-1)
	{
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+ifst))->emcShower();
		RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+isec))->emcShower();

		double eraw1 = g1Trk->energy();
		double the1  = g1Trk->theta();
		double phi1  = g1Trk->phi();
		double eraw2 = g2Trk->energy();
		double the2  = g2Trk->theta();
		double phi2  = g2Trk->phi();

		HepLorentzVector ptrk1, ptrk2;
		ptrk1.setPx(eraw1*sin(the1)*cos(phi1));
		ptrk1.setPy(eraw1*sin(the1)*sin(phi1));
		ptrk1.setPz(eraw1*cos(the1));
		ptrk1.setE(eraw1);
		ptrk2.setPx(eraw2*sin(the2)*cos(phi2));
		ptrk2.setPy(eraw2*sin(the2)*sin(phi2));
		ptrk2.setPz(eraw2*cos(the2));
		ptrk2.setE(eraw2);

		double dthe = fabs((the1+the2)*180.0/CLHEP::pi-180);

		//--> Veto Bhabha events
		if(dthe<m_eedeg && eraw2>m_Eratio*Ebeam) return StatusCode::SUCCESS;

		m_dthe  = dthe;
		m_eraw1 = eraw1;
		m_eraw2 = eraw2;
	}

	NcutVetobb++;
	//cout<<"NcutVetobb = "<<NcutVetobb<<endl; 
	//////////////////////////////////////////////////////////////////////////////////////
	// Selection of the good photons
	Vp4 pGam;
	pGam.clear();
	Vint iGam;
	iGam.clear();

	int nGam = 0;
	for(int i=evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++) // loop among all the neutral tracks
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin()+i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();

		double emct = emcTrk->time();
		double fcos = cos(emcTrk->theta());
		double eraw = emcTrk->energy();

		double theta = emcTrk->theta();
		double   phi = emcTrk->phi();

		//-> angle(tracks), minimun deposited energy, time
		if(emct<0||emct>14) continue;

		if(fabs(fcos)<0.8)
		{
			if(eraw<m_GammaBarrelEth) continue;
		}
		else if(fabs(fcos)>0.86&&fabs(fcos)<0.92)
		{
			if(eraw<m_GammaEndCapEth) continue;
		}
		else
		{
			continue;
		}

		iGam.push_back(i);
	}

	nGam = iGam.size();
	m_ngam = nGam;

	double ShwEneu = 0.0;
	double Balance = 0.0;

	HepLorentzVector gamptrk;
	HepLorentzVector Totgamptrk;

	for(int i=0; i<nGam; i++)
	{
		EvtRecTrackIterator  itTrk = evtRecTrkCol->begin()+iGam[i];
		RecEmcShower *emcTrk = (*itTrk)->emcShower();

		double eraw  = emcTrk->energy();
		double theta = emcTrk->theta();
		double phi   = emcTrk->phi();

		ShwEneu += eraw;  ///// sum of energy of good photons

		gamptrk.setPx(eraw*sin(theta)*cos(phi));
		gamptrk.setPy(eraw*sin(theta)*sin(phi));
		gamptrk.setPz(eraw*cos(theta)); 
		gamptrk.setE(eraw);

		Totgamptrk += gamptrk;

		m_GamThe[i] = theta*180.0/CLHEP::pi;
		m_GamPhi[i] = phi*180.0/CLHEP::pi;
		m_GamEne[i] = eraw;
		pGam.push_back(gamptrk);
	}

	Hep3Vector Totgamp = Totgamptrk.getV();
	double pneutrk = Totgamp.mag();

	if(ShwEneu!=0.0) Balance = pneutrk/ShwEneu;

	m_Balance = Balance;

	// cout << "nGam = " << nGam << endl; 
	// end of selecting good photons
	///////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////
	// selection of Isolated photons
	Vint iIsoGam;  // isolated photon
	iIsoGam.clear();

	int nIsoGam = 0;
	for(int i=evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;

		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		double eraw = emcTrk->energy();

		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

		double dthe = 200.; 
		double dphi = 200.;
		double dang = 200.;

		for(int j=0; j<evtRecEvent->totalCharged(); j++)
		{
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;

			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();

			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();

			double angd = extpos.angle(emcpos);
			double thed = extpos.theta() - emcpos.theta();
			double phid = extpos.deltaPhi(emcpos);

			thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+CLHEP::pi, CLHEP::twopi) - CLHEP::pi;
			phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+CLHEP::pi, CLHEP::twopi) - CLHEP::pi;

			if(angd < dang)
			{
				dang = angd;
				dthe = thed;
				dphi = phid;
			}
		}

		if(dang>=200) continue;
		double emct = emcTrk->time();
		double fcos = cos(emcTrk->theta());

		dthe = fabs(dthe*180/(CLHEP::pi));
		dphi = fabs(dphi*180/(CLHEP::pi));
		dang = fabs(dang*180/(CLHEP::pi));

		//-> angle(tracks), minimun deposited energy, time
		if(dang<m_IsogamAngCut) continue;
		if(eraw<m_IsogamECut ) continue;
		if(emct<0||emct>14) continue;

		iIsoGam.push_back(i);
	}

	nIsoGam = iIsoGam.size();
	m_NIsogam = nIsoGam;

	for(int i=0; i<nIsoGam; i++)
	{
		EvtRecTrackIterator  itTrk=evtRecTrkCol->begin()+iIsoGam[i];
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		double eraw  = emcTrk->energy();
		double theta = emcTrk->theta();
		double phi   = emcTrk->phi();

		m_IsoGamThe[i] = theta*180/(CLHEP::pi);
		m_IsoGamPhi[i] = phi*180/(CLHEP::pi);
		m_IsoGamEne[i] = eraw;
	}

	// cout << "nIso = " << nIsoGam << endl; 
	// end of selecting Isolated photons
	//////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////
	//// ##### the 2nd step -->  select good track  #####
	////       
	Vint iGoodp,iGoodm,iGood;
	iGood.clear(); iGoodp.clear(); iGoodm.clear();

	Vp4 Vp4ep, Vp4em;
	Vp4ep.clear(); Vp4em.clear(); 
	Vint Viep, Viem;
	Viep.clear(); Viem.clear();

	////////////////////////////////////////////////////////////////////////////////////
	// loop among the total charged tracks in each event

	int nGood = 0;
	for(int i=0; i<evtRecEvent->totalCharged(); i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin()+i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;

		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]); // Vr in xy plane
		double  Rvz0=vecipa[3]; // Vz

		double cost = cos(mdcTrk->theta());
		double    p = mdcTrk->p();
		double   pt = mdcTrk->pxy();
		double chi2 = mdcTrk->chi2(); // used for test
		int    ndof = mdcTrk->ndof(); // used for test

		//--> cuts for track selection
		if(fabs(Rvz0) >= m_vz0cut) continue;
		if(fabs(Rvxy0)>= m_vr0cut) continue;
		if(fabs(cost) >= m_cosThetacut) continue; 
		if(p > m_MomRatiocut*Ebeam) continue;

		/////////////////////////////////////////////////////////////////////////////
		// veto particle whose mass is higher than proton
		// use chi_{dE/dx}(proton) information

		if((*itTrk)->isMdcDedxValid())
		{
			RecMdcDedx* dedxTrk = (*itTrk)->mdcDedx();
			double dEdxChiProton = dedxTrk->chiP();

			//--> cut on dE/dx with proton hypothesis
			if(dEdxChiProton>m_dEdxChipCut) continue;
		}

		if((*itTrk)->isEmcShowerValid())
		{
			RecEmcShower *emcTrk = (*itTrk)->emcShower();
			double eraw = emcTrk->energy();

			//--> veto e+/e- track of Bhabha
			if(p>m_EBhabhaRatio*Ebeam && eraw/p>m_EOPcut) continue;

			Hep3Vector p3trk = mdcTrk->p3();
			double BhabhaE = sqrt(p3trk.mag2()+xmass[0]*xmass[0]); // regarding the track as the electron.
			HepLorentzVector ptrk(p3trk, BhabhaE);

			// record e+/e- tracks for gamma conversion
			if(eraw/p>m_EOPcut)
			{
				if(mdcTrk->charge()==1)
				{
					Vp4em.push_back(ptrk);
					Viem.push_back(i);
				}
				else if(mdcTrk->charge()==-1)
				{
					Vp4ep.push_back(ptrk);
					Viep.push_back(i);
				} 
			}
		}

		iGood.push_back(i); // iGood[i], where i=0; i<evtRecEvent->totalCharged(); i++

		if(mdcTrk->charge()==1)
		{
			iGoodp.push_back(i);
		}
		else if(mdcTrk->charge()==-1)
		{
			iGoodm.push_back(i);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////
	// veto e+/e- from gamma conversion
	Vint iCand; iCand.clear(); // candidate after excluding the gamma conversion events
	Vint iGamConEpEm; iGamConEpEm.clear(); // track ID of the electron and positron of gamma conversion 
	double Egamcon = 0;
	int Ngamcon = 0;

	for(int i=0; i<Vp4ep.size(); i++)
	{
		for(int j=0; j<Vp4em.size(); j++)
		{
			double massee = (Vp4ep[i]+Vp4em[j]).m();
			double opangl = Vp4ep[i].angle(Vp4em[j].vect())*180.0/CLHEP::pi;

			// gamma conversion events
			if(massee<m_GamConEcut && opangl<m_GamConAngcut)
			{
				double energyee = (Vp4ep[i]+Vp4em[j]).e();
				Egamcon += energyee;
				Ngamcon++;

				iGamConEpEm.push_back(Viep[i]);
				iGamConEpEm.push_back(Viem[j]);
			}
		}
	}

	m_Ngamcon = Ngamcon;
	m_Egamcon = Egamcon;

	for(int i=0; i<iGood.size(); i++)
	{
		int tag = 0;
		for(int j=0; j<iGamConEpEm.size(); j++)
		{
			if(iGood[i]==iGamConEpEm[j]) tag = 1;
		}

		//--> Excluding the gamma conversion events
		if(tag==0) iCand.push_back(iGood[i]); // iCand[i], where i = iGood[0], iGood[1], ..., iGood[nGood]
	}

	// end of veto e+/e- from gamma conversion
	// end of the 2nd step -->  select good track 
	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////

	nGood = iCand.size();
	m_nGood = nGood;

	// cout << "nGood = " << nGood << endl; 
	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	// The calculating of the N_large_eop and N_large_probE

	int N_large_eop = 0;
	int N_large_probE = 0;
	double ShwEtrk = 0.0;   
	double EneEtrk = 0.0;   
	double TotVisE = 0.0;   
	double TotEvtE = 0.0;   

	double p_fst = 0, p_sec = 0;
	Hep3Vector p3_fst, p3_sec;

	ParticleID *pid = ParticleID::instance();

	if(nGood>0)
	{
		for(int i=0; i<nGood; i++)
		{
			EvtRecTrackIterator itTrk = evtRecTrkCol->begin()+iCand[i];
			RecMdcTrack *mdcTrk =(*itTrk)->mdcTrack();

			// accumulate the energy of an event
			double pTrk = 0.0;
			double TrkEInMDC = 0.0;
			Hep3Vector p3Trk;

			pTrk = mdcTrk->p();
			p3Trk = mdcTrk->p3();
			TrkEInMDC = sqrt(p3Trk.mag2()+xmass[2]*xmass[2]); // regarding the track as pions.

			EneEtrk += TrkEInMDC;

			// Indexing the tracks with EOP larger than 0.8 (used for nGood==3)
			double TrkEInEMC = 0.0;
			double eop = 0.0;
			if((*itTrk)->isEmcShowerValid())
			{
				RecEmcShower *emcTrk = (*itTrk)->emcShower();
				TrkEInEMC = emcTrk->energy();
				eop = TrkEInEMC/pTrk;

				if(eop>m_EOPcut)
				{
					N_large_eop++;
				}

				ShwEtrk += TrkEInEMC;
			}

			// using PID to index the tracks with PID ratio larger than 0.25 (used for nGood==3)
			pid->init();
			pid->setMethod(pid->methodProbability());
			pid->setChiMinCut(4);
			pid->setRecTrack(*itTrk);
			pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());
			pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton() | pid->onlyElectron() | pid->onlyMuon());
			pid->calculate();

			double pid_ratio = 0.0;
			if(pid->IsPidInfoValid())
			{
				pid_ratio = pid->probElectron()/(pid->probElectron()+pid->probPion()+pid->probKaon()+pid->probProton());
				if(pid_ratio>m_PIDRatioCut)
				{
					N_large_probE++;
				}
			}

			// reserve two tracks with highest momentum
			if(pTrk>p_fst)
			{
				p_sec = p_fst;
				p3_sec = p3_fst;
				p_fst = pTrk;
				p3_fst = p3Trk;
			}
			else if(pTrk>p_sec)
			{
				p_sec = pTrk;
				p3_sec = p3Trk;
			}

			// dE/dx info for good charged tracks
			double dEdx = -1.0;
			double dEdxchip = -100.0;
			if((*itTrk)->isMdcDedxValid())
			{
				RecMdcDedx* dedxTrk = (*itTrk)->mdcDedx();
				dEdx = dedxTrk->probPH();
				dEdxchip = dedxTrk->chiP();
			}

			// record the infos of good charged tracks
			double Rvxy0 = 0.0;
			double Rvz0 = 0.0;
			double trkphi = 0.0;
			HepVector a = mdcTrk->helix();
			HepSymMatrix Ea = mdcTrk->err();
			HepPoint3D point0(0.,0.,0.);
			HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
			VFHelix helixip(point0,a,Ea);
			helixip.pivot(IP);
			HepVector vecipa = helixip.a();
			Rvxy0 = vecipa[0];
			Rvz0 = vecipa[3];
			trkphi = mdcTrk->phi()*180.0/CLHEP::pi;

			m_vr0[i]  = Rvxy0;
			m_vz0[i]  = Rvz0;
			m_TrkP[i] = pTrk;
			m_EMCE[i] = TrkEInEMC;
			m_dEdx[i] = dEdx;
			m_dEdxchip[i] = dEdxchip;
			m_phi[i]  = trkphi;

			m_charge[i]   = mdcTrk->charge();
			m_costheta[i] = cos(mdcTrk->theta());
		}
	}

	m_N_large_eop = N_large_eop;
	m_N_large_probE = N_large_probE;

	TotVisE = ShwEtrk + ShwEneu;
	TotEvtE = ShwEtrk + EneEtrk + ShwEneu;

	m_EEmcNeu = ShwEneu;
	m_EEmcChg = ShwEtrk;
	m_TotVisE = TotVisE;
	m_TotEvtE = TotEvtE;

	// end of recording the info of good charged tracks
	//////////////////////////////////////////////////////////////////////////////////////

	//// Reconstruction of pi0
	int npi0 = 0;
	Vp4 p4pi0A; p4pi0A.clear();
	Vint gam1ID; gam1ID.clear();
	Vp4 p4gam1; p4gam1.clear();
	Vint gam2ID; gam2ID.clear();
	Vp4 p4gam2; p4gam2.clear();
	Vdb chi2pi0; chi2pi0.clear();

	if(nGam>=2) {
		KalmanKinematicFit* kmfitpi0 = KalmanKinematicFit::instance();
		for(int i=0; i<nGam-1; i++)
		{
			for(int j=i+1; j<nGam; j++)
			{
				RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
				RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
				HepLorentzVector emcp1 = pGam[i];
				HepLorentzVector emcp2 = pGam[j];
				double mgg = (emcp1+emcp2).m();
				if(mgg < m_Pi0LowMass || mgg > m_Pi0HigMass) continue;
				kmfitpi0->init();
				//kmfitpi0->setIterNumber(5);
				kmfitpi0->setChisqCut(2500);
				kmfitpi0->AddTrack(0,0.0,g1Trk);
				kmfitpi0->AddTrack(1,0.0,g2Trk);
				kmfitpi0->AddResonance(0,mpion0,0,1);
				if(!kmfitpi0->Fit(0)) continue;
				bool oksq = kmfitpi0->Fit();
				if(oksq)
				{
					double chisqpi0 = kmfitpi0->chisq();
					HepLorentzVector pGam1 = kmfitpi0->pfit(0);
					HepLorentzVector pGam2 = kmfitpi0->pfit(1);
					HepLorentzVector ppi0 = kmfitpi0->pfit(0) + kmfitpi0->pfit(1);
					if(chisqpi0 <=2500)
					{
						p4pi0A.push_back(ppi0);
						chi2pi0.push_back(chisqpi0);
						gam1ID.push_back(i);
						p4gam1.push_back(pGam1);
						gam2ID.push_back(j);
						p4gam2.push_back(pGam2);
						npi0++;
					}
				}
			}
		}
	}

	m_npi0 = npi0;

	for(int i = 0; i < npi0; i++) {

		m_chi2pi0[i] = chi2pi0[i];
		m_pxpi0[i] = p4pi0A[i][0];
		m_pypi0[i] = p4pi0A[i][1];
		m_pzpi0[i] = p4pi0A[i][2];
		m_Epi0[i] = p4pi0A[i][3];
		m_idgam1[i] = gam1ID[i];
		m_idgam2[i] = gam2ID[i];
	}

	// cout << "npi0 = " << npi0 << endl; 
	// end of selecting of pi0
	//////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////
	// the 3nd step -->  select events with at least 1 tracks

	//--> Cut on number of Good tracks
	if(nGood<m_NGoodCut) return StatusCode::SUCCESS;
	NcutNGoodgeZero++;

	//////////////////////////////////////////////////////////////////////////////////////
	// information on 0-prong events
	if(nGood==0) // Note that here nGood = iCand.size();
	{
		NcutNGoodeqZero++;

		// cut on event energy (total estimated energy of the event)
		if(TotEvtE<m_NGood0EvtECut*Ebeam) return StatusCode::SUCCESS;
		NcutNGoodeqZero_EvtECut++;

		//--> requirment on number of good pi0
		if(npi0<m_NGood0NPi0Cut) return StatusCode::SUCCESS;
		NcutNGoodeqZero_NPi0Cut++;

		//--> requirment on number of isolated photons of the event
		if(nIsoGam<m_NIsogamCut) return StatusCode::SUCCESS;
		NcutNGoodeqZero_NIsogamCut++;
	}

	//////////////////////////////////////////////////////////////////////////////////////
	// information on 1-prong events
	if(nGood==1) // Note that here nGood = iCand.size();
	{
		NcutNGoodeqOne++;

		// cut on event energy (total estimated energy of the event)
		if(TotEvtE<m_EvtERatio*Ebeam) return StatusCode::SUCCESS;
		NcutNGoodeqOne_EvtECut++;

		//--> requirment on number of good pi0
		if(npi0<m_NGoodPi0Cut) return StatusCode::SUCCESS;
		NcutNGoodeqOne_NPi0Cut++;

		//--> requirment on Balance of the event
		if(Balance>m_BalanceCut || Balance==0) return StatusCode::SUCCESS;
		NcutNGoodeqOne_BalanceCut++;
	}

	// information on 2-prong events
	if(nGood==2) // Note that here nGood = iCand.size();
	{
		NcutNGoodeqTwo++;
		RecMdcTrack *mdcTrk1 = (*(evtRecTrkCol->begin()+iCand[0]))->mdcTrack();
		RecMdcTrack *mdcTrk2 = (*(evtRecTrkCol->begin()+iCand[1]))->mdcTrack();

		Hep3Vector mtrk1 = mdcTrk1->p3();
		Hep3Vector mtrk2 = mdcTrk2->p3();

		double d2ang = mtrk1.angle(mtrk2)*180.0/CLHEP::pi;
		double the21 = mtrk1.theta()*180.0/CLHEP::pi;
		double the22 = mtrk2.theta()*180.0/CLHEP::pi;
		double phi21 = mtrk1.phi()*180.0/CLHEP::pi;
		double phi22 = mtrk2.phi()*180.0/CLHEP::pi;
		double d2the = fabs(the21+the22-180.0);
		double d2phi = fabs(fabs(phi21-phi22)-180.0);

		//--> angle requirements: the tracks are not back-to-back
		//|||||| Note that the requirements should be satisfied simultaneously.
		if(d2the<m_d2ThetaCut && d2phi<m_d2PhiCut) return StatusCode::SUCCESS; 
		NcutNGoodeqTwo_AngCut++;

		//--> requirment on number of isolated photon
		if(nIsoGam<m_NIsogamCut) return StatusCode::SUCCESS;
		NcutNGoodeqTwo_Vetobb++;

		m_d2ang = d2ang;
		m_d2the = d2the;
		m_d2phi = d2phi;
	}

	// information on 3-prong events
	if(nGood==3)
	{
		NcutNGoodeqThree++;
		if(p_fst>0&&p_sec>0)
		{
			double d3ang = p3_fst.angle(p3_sec)*180.0/CLHEP::pi;
			double the31 = p3_fst.theta()*180.0/CLHEP::pi;
			double the32 = p3_sec.theta()*180.0/CLHEP::pi;
			double phi31 = p3_fst.phi()*180.0/CLHEP::pi;
			double phi32 = p3_sec.phi()*180.0/CLHEP::pi;
			double d3the = fabs(the31+the32-180.0);
			double d3phi = fabs(fabs(phi31-phi32)-180.0);

			//--> cuts on the angles between tracks
			//|||||| Note that the requirements should be satisfied simultaneously.
			if(d3the<m_d3ThetaCut && d3phi<m_d3PhiCut) return StatusCode::SUCCESS;      
			NcutNGoodeqThree_AngCut++;

			//--> cut on the EOP of the tracks
			if(N_large_eop > m_NlargeeopCut) return StatusCode::SUCCESS;
			NcutNGoodeqThree_EOPCut++;

			//--> cut on the probability of electron
			if(N_large_probE > m_NlargeprobECut) return StatusCode::SUCCESS;
			NcutNGoodeqThree_Vetobb++;

			m_d3ang = d3ang;
			m_d3the = d3the;
			m_d3phi = d3phi;
		}
	}

	// information on 4,5,6,...-prong events

	if(nGood>=4){ NcutNGoodgtThree++; }

	NcutFinal++;
	// end of counting the final good charged tracks
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////
	//// ##### the 3rd step --> resonance reconstruction  #####

	///////////////////////////////////////////////////////////////////////////
	// particle identification

	Vp4 ppip, ppim, pkap, pkam, pprop, pprom;
	ppip.clear(); ppim.clear();	pkap.clear(); pkam.clear(); pprop.clear(); pprom.clear();
	Vint kapID; kapID.clear();
	Vint kamID; kamID.clear();
	Vint pipID; pipID.clear();
	Vint pimID; pimID.clear();
	Vint propID; propID.clear();
	Vint promID; promID.clear();
	Vdb ProbK; ProbK.clear();
	Vdb ProbPi; ProbPi.clear();
	Vdb ProbP; ProbP.clear();

	double Etrk_pid = 0.0;
	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin()+iCand[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());      
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
		pid->calculate();

		if(!(pid->IsPidInfoValid())) continue;
		double  prob_pi = pid->probPion();
		double  prob_K  = pid->probKaon();
		double  prob_p  = pid->probProton();
		ProbK.push_back(prob_K);
		ProbPi.push_back(prob_pi);
		ProbP.push_back(prob_p);

		// pi+/-
		if (prob_pi > prob_K && prob_pi > prob_p && prob_pi > 0.001) {

			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
			RecMdcKalTrack::setPidType (RecMdcKalTrack::pion);

			HepLorentzVector ptrk;
			ptrk.setPx(mdcKalTrk->px());
			ptrk.setPy(mdcKalTrk->py());
			ptrk.setPz(mdcKalTrk->pz());
			double p3 = ptrk.mag();
			ptrk.setE(sqrt(p3*p3+mpi*mpi)); 

			Etrk_pid += sqrt(p3*p3+mpi*mpi);

			if(mdcKalTrk->charge() > 0) ppip.push_back(ptrk);
			else ppim.push_back(ptrk);

			if(mdcKalTrk->charge() > 0) pipID.push_back(i);
			else pimID.push_back(i);
		}

		// K+/-
		if (prob_K > prob_pi && prob_K > prob_p && prob_K > 0.001) {

			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
			RecMdcKalTrack::setPidType (RecMdcKalTrack::kaon);

			HepLorentzVector ptrk;
			ptrk.setPx(mdcKalTrk->px());
			ptrk.setPy(mdcKalTrk->py());
			ptrk.setPz(mdcKalTrk->pz());
			double p3 = ptrk.mag();
			ptrk.setE(sqrt(p3*p3+mk*mk));

			Etrk_pid += sqrt(p3*p3+mk*mk);

			if (mdcKalTrk->charge() > 0) pkap.push_back(ptrk);
			else pkam.push_back(ptrk);

			if (mdcKalTrk->charge() > 0) kapID.push_back(i);
			else kamID.push_back(i);
		}

		if (prob_p > prob_K && prob_p > prob_pi && prob_p > 0.001) {

			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
			RecMdcKalTrack::setPidType (RecMdcKalTrack::proton);

			HepLorentzVector ptrk;
			ptrk.setPx(mdcKalTrk->px());
			ptrk.setPy(mdcKalTrk->py());
			ptrk.setPz(mdcKalTrk->pz());
			double p3 = ptrk.mag();
			ptrk.setE(sqrt(p3*p3+mp*mp));

			Etrk_pid += sqrt(p3*p3+mp*mp);

			if (mdcKalTrk->charge() > 0) pprop.push_back(ptrk);
			else pprom.push_back(ptrk);

			if (mdcKalTrk->charge() > 0) propID.push_back(i);
			else promID.push_back(i);
		}
	}

	m_ETrkAftPID = Etrk_pid;

	int npip = ppip.size();
	int npim = ppim.size();
	int nprop= pprop.size();
	int nprom= pprom.size();
	int nkap = pkap.size();
	int nkam = pkam.size();

	m_npip = npip;
	m_npim = npim;
	m_nkap = nkap;
	m_nkam = nkam;

	// Saving the information of each tracks

	for (int i = 0; i < nkap; i++) {
		int id = kapID[i];
		RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + id;
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		HepLorentzVector p4 = mdcKalTrk->p4(mk);

		m_idkap[i] = id;
		m_pxkap[i] = mdcKalTrk->px();
		m_pykap[i] = mdcKalTrk->py();
		m_pzkap[i] = mdcKalTrk->pz();
		m_Ekap[i] = p4[3];
		m_mkap[i] = mk;

		m_KapProbK[i] = ProbK[id];
		m_KapProbPi[i] = ProbPi[id];
		m_KapProbP[i] = ProbP[id];
	}

	for (int i = 0; i < nkam; i++) {
		int id = kamID[i];
		RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + id;
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		HepLorentzVector p4 = mdcKalTrk->p4(mk);

		m_idkam[i] = id;
		m_pxkam[i] = mdcKalTrk->px();
		m_pykam[i] = mdcKalTrk->py();
		m_pzkam[i] = mdcKalTrk->pz();
		m_Ekam[i] = p4[3];
		m_mkam[i] = mk;

		m_KamProbK[i] = ProbK[id];
		m_KamProbPi[i] = ProbPi[id];
		m_KamProbP[i] = ProbP[id];
	}

	for (int i = 0; i < npip; i++) {
		int id = pipID[i];
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + id;
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		HepLorentzVector p4 = mdcKalTrk->p4(mpi);

		m_idpip[i] = id;
		m_pxpip[i] = mdcKalTrk->px();
		m_pypip[i] = mdcKalTrk->py();
		m_pzpip[i] = mdcKalTrk->pz();
		m_Epip[i] = p4[3];
		m_mpip[i] = mpi;

		m_PipProbK[i] = ProbK[id];
		m_PipProbPi[i] = ProbPi[id];
		m_PipProbP[i] = ProbP[id];
	}

	for (int i = 0; i < npim; i++) {
		int id = pimID[i];
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + id;
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		HepLorentzVector p4 = mdcKalTrk->p4(mpi);

		m_idpim[i] = id;
		m_pxpim[i] = mdcKalTrk->px();
		m_pypim[i] = mdcKalTrk->py();
		m_pzpim[i] = mdcKalTrk->pz();
		m_Epim[i] = p4[3];
		m_mpim[i] = mpi;

		m_PimProbK[i] = ProbK[id];
		m_PimProbPi[i] = ProbPi[id];
		m_PimProbP[i] = ProbP[id];
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	double aveVz=0.0;

	if(nGood>=1) {
		for(int i=0; i<nGood; i++) { //charge must be +-1 
			RecMdcTrack *mdcTrk = (*(evtRecTrkCol->begin()+iCand[i]))->mdcTrack();

			HepVector a = mdcTrk->helix();
			HepSymMatrix Ea = mdcTrk->err();
			HepPoint3D point0(0.,0.,0.);
			HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
			VFHelix helixip(point0,a,Ea);
			helixip.pivot(IP);
			HepVector vecipa = helixip.a();
			double  Rvz0=vecipa[3];
			aveVz += Rvz0;
		}
	}

	m_avz = double(aveVz)/nGood;
	//cout<<"m_avz = "<<m_avz<<endl;

	//////////////////////////////////////////////////////////////////////////////////////////
	//// Reconstruction of KS0 -> pi+ pi- with vertex fit

	//nGood loop for Without PID and All tracks are considered as pion
	Vint IdTrp; IdTrp.clear();
	Vint IdTrm; IdTrm.clear();
	Vp4 p4Trp; p4Trp.clear();
	Vp4 p4Trm; p4Trm.clear();
	for(int i = 0; i < evtRecEvent->totalCharged(); i++) {

		EvtRecTrackIterator itTrk = evtRecTrkCol->begin()+i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;

		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]); // Vr in xy plane
		double  Rvz0=vecipa[3]; // Vz

		double cost = cos(mdcTrk->theta());
		double    p = mdcTrk->p();
		double   pt = mdcTrk->pxy();
		double chi2 = mdcTrk->chi2(); // used for test
		int    ndof = mdcTrk->ndof(); // used for test

		//--> cuts for track selection
		if(fabs(Rvz0) >= 20) continue;
		if(fabs(cost) >= m_cosThetacut) continue; 

		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		HepLorentzVector p4 = mdcKalTrk->p4(mpi);

		if(mdcTrk->charge()==1) {
			IdTrp.push_back(i);
			p4Trp.push_back(p4);			
		}
		else if(mdcTrk->charge()==-1) {
			IdTrm.push_back(i);
			p4Trm.push_back(p4);	
		}
	}
	int NumTrp = IdTrp.size();
	int NumTrm = IdTrm.size();
	//cout << NumTrp << "	" << NumTrm << endl;

	int NumKS0 = 0;
	Vint IDTrp; IDTrp.clear();
	Vint IDTrm; IDTrm.clear();
	Vp4 P4Trp; P4Trp.clear();
	Vp4 P4Trm; P4Trm.clear();
	Vp4 P4KS0; P4KS0.clear();
	Vdb Chi2KS0IP; Chi2KS0IP.clear();
	Vdb Chi2KS0SP; Chi2KS0SP.clear();
	Vdb LIFEKS0; LIFEKS0.clear();
	Vdb LENKS0; LENKS0.clear();
	Vdb LENERRKS0; LENERRKS0.clear();
	if(NumTrp > 0 && NumTrm > 0) {
		for(int ia = 0; ia < NumTrp; ia++) {

			RecMdcKalTrack *pip1Trk = (*(evtRecTrkCol->begin()+IdTrp[ia]))->mdcKalTrack();
			pip1Trk->setPidType(RecMdcKalTrack::pion);
			WTrackParameter wvpip1Trk = WTrackParameter(mpi, pip1Trk->getZHelix(), pip1Trk->getZError());

			for(int ib = 0; ib < NumTrm; ib++) {

				RecMdcKalTrack *pim1Trk = (*(evtRecTrkCol->begin()+IdTrm[ib]))->mdcKalTrack();
				pim1Trk->setPidType(RecMdcKalTrack::pion);
				WTrackParameter wvpim1Trk = WTrackParameter(mpi, pim1Trk->getZHelix(), pim1Trk->getZError());

				Hep3Vector ip(0,0,0);
				HepSymMatrix ipEx(3,0);
				IVertexDbSvc* vtxsvc;
				Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
				if (vtxsvc->isVertexValid()) {
					double* dbv = vtxsvc->PrimaryVertex();
					double* vv  = vtxsvc->SigmaPrimaryVertex();
					ip.setX(dbv[0]);
					ip.setY(dbv[1]);
					ip.setZ(dbv[2]);
					ipEx[0][0] = vv[0] * vv[0];
					ipEx[1][1] = vv[1] * vv[1];
					ipEx[2][2] = vv[2] * vv[2];
				}
				else return StatusCode::SUCCESS; // if cannot load vertex information, will go to another event

				VertexParameter bs;
				bs.setVx(ip);
				bs.setEvx(ipEx);
				///////////////////////////////////////////////////////////////////
				//  set a common vertex with huge error
				///////////////////////////////////////////////////////////////////
				HepPoint3D    vx(0., 0., 0.);
				HepSymMatrix  Evx(3, 0);
				double bx = 1E+6;
				double by = 1E+6;
				double bz = 1E+6;
				Evx[0][0] = bx * bx;
				Evx[1][1] = by * by;
				Evx[2][2] = bz * bz;
				VertexParameter vxpar;
				vxpar.setVx(vx);
				vxpar.setEvx(Evx);
				///////////////////////////////////////////////////////////////////
				//  do vertex fit
				///////////////////////////////////////////////////////////////////
				VertexFit *vtxfit = VertexFit::instance();
				vtxfit->init();
				vtxfit->AddTrack(0, wvpip1Trk);
				vtxfit->AddTrack(1, wvpim1Trk);
				vtxfit->AddVertex(0, vxpar, 0, 1);
				if (!(vtxfit->Fit(0))) continue;
				vtxfit->Swim(0);
				vtxfit->BuildVirtualParticle(0);
				double vtx_chisq = vtxfit->chisq(0);
				///////////////////////////////////////////////////////////////////
				//  do second vertex fit
				///////////////////////////////////////////////////////////////////
				SecondVertexFit *svtxfit = SecondVertexFit::instance();
				svtxfit->init();
				svtxfit->setPrimaryVertex(bs);
				svtxfit->AddTrack(0, vtxfit->wVirtualTrack(0));
				svtxfit->setVpar(vtxfit->vpar(0));
				if (!svtxfit->Fit()) continue;
				double svtx_chisq = svtxfit->chisq();
				wvpip1Trk = vtxfit->wtrk(0);
				wvpim1Trk = vtxfit->wtrk(1);
				HepLorentzVector pKs = svtxfit->p4par();
				HepVector vtxKs = svtxfit->vpar().Vx();
				double ctau = svtxfit->ctau();
				double len = svtxfit->decayLength();
				double lenerr = svtxfit->decayLengthError();
				HepLorentzVector p4pip = vtxfit->pfit(0);
				HepLorentzVector p4pim = vtxfit->pfit(1);
				NumKS0 = NumKS0 + 1;

				IDTrp.push_back(ia);
				IDTrm.push_back(ib);
				P4Trp.push_back(p4pip);
				P4Trm.push_back(p4pim);
				P4KS0.push_back(pKs);
				Chi2KS0IP.push_back(vtx_chisq);
				Chi2KS0SP.push_back(svtx_chisq);
				LIFEKS0.push_back(ctau);
				LENKS0.push_back(len);
				LENERRKS0.push_back(lenerr);
			}
		}
	}
	m_nKS0 = NumKS0;
	//cout << "NumKS0 = " << NumKS0 << endl;

	for(int i = 0; i < NumKS0; i++) {
		m_chi2KS0IP[i] = Chi2KS0IP[i];
		m_chi2KS0SP[i] = Chi2KS0SP[i];
		m_LifeKS0[i] = LIFEKS0[i];
		m_LenKS0[i] = LENKS0[i];
		m_LenErrKS0[i] = LENERRKS0[i];

		m_pxKS0[i] = P4KS0[i][0];
		m_pyKS0[i] = P4KS0[i][1];
		m_pzKS0[i] = P4KS0[i][2];
		m_EKS0[i] = P4KS0[i][3];
		m_mKS0[i] = sqrt(abs(P4KS0[i][3]*P4KS0[i][3]-P4KS0[i][0]*P4KS0[i][0]-P4KS0[i][1]*P4KS0[i][1]-P4KS0[i][2]*P4KS0[i][2]));

		m_IdTrp[i] = IDTrp[i];
		m_pxTrp[i] = P4Trp[i][0];
		m_pyTrp[i] = P4Trp[i][1];
		m_pzTrp[i] = P4Trp[i][2];
		m_ETrp[i] = P4Trp[i][3];

		m_IdTrm[i] = IDTrm[i];
		m_pxTrm[i] = P4Trm[i][0];
		m_pyTrm[i] = P4Trm[i][1];
		m_pzTrm[i] = P4Trm[i][2];
		m_ETrm[i] = P4Trm[i][3];
	}

	m_tuple1->write();
	return StatusCode::SUCCESS;

} 
//end of execute()
/////////////////////////////////////////////////////////////////////////////////////////


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode HadronSel::finalize(){
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;

	cout << "initial events   :    " << NcutIni << endl;
	cout << "total events     :    " << NcutTot << endl;
	cout << "double count tag :    " << VetoDoubleCount << endl;
	cout << "veto Bhabha      :    " << NcutVetobb << endl;
	cout << "nGood>=0         :    " << NcutNGoodgeZero << endl;
	cout << "crude nGood = 0  :    " << NcutNGoodeqZero << endl;
	cout << "crude nGood = 1  :    " << NcutNGoodeqOne << endl;
	cout << "crude nGood = 2  :    " << NcutNGoodeqTwo << endl;
	cout << "crude nGood = 3  :    " << NcutNGoodeqThree << endl;
	cout << "crude nGood > 3  :    " << NcutNGoodgtThree << endl;
	cout << "clean nGood = 0  :    " << NcutNGoodeqZero_NIsogamCut << endl;
	cout << "clean nGood = 1  :    " << NcutNGoodeqOne_BalanceCut << endl;
	cout << "clean nGood = 2  :    " << NcutNGoodeqTwo_Vetobb << endl;
	cout << "clean nGood = 3  :    " << NcutNGoodeqThree_Vetobb << endl;
	cout << "clean nGood > 3  :    " << NcutNGoodgtThree << endl;
	cout << "Final events     :    " << NcutNGoodeqZero_NIsogamCut <<" + "<< NcutNGoodeqOne_BalanceCut <<" + "<< NcutNGoodeqTwo_Vetobb <<" + "<< NcutNGoodeqThree_Vetobb <<" + "<< NcutNGoodgtThree <<" = "<< NcutNGoodeqZero_NIsogamCut+NcutNGoodeqOne_BalanceCut+NcutNGoodeqTwo_Vetobb+NcutNGoodeqThree_Vetobb+NcutNGoodgtThree << "::" << NcutFinal << endl;
	cout << endl;
	cout << "cut flow nGood=0 :    " << NcutNGoodeqZero << " --> " << NcutNGoodeqZero_EvtECut << " --> " << NcutNGoodeqZero_NPi0Cut << " --> " << NcutNGoodeqZero_NIsogamCut << endl;
	cout << "cut flow nGood=1 :    " << NcutNGoodeqOne << " --> " << NcutNGoodeqOne_EvtECut << " --> " << NcutNGoodeqOne_NPi0Cut << " --> " << NcutNGoodeqOne_BalanceCut << endl;
	cout << "cut flow nGood=2 :    " << NcutNGoodeqTwo << " --> " << NcutNGoodeqTwo_AngCut << " --> " << NcutNGoodeqTwo_Vetobb << endl;
	cout << "cut flow nGood=3 :    " << NcutNGoodeqThree << " --> " << NcutNGoodeqThree_AngCut << " --> " << NcutNGoodeqThree_EOPCut << " --> " << NcutNGoodeqThree_Vetobb << endl;
	cout << endl;


	return StatusCode::SUCCESS;
}
