#include "RevisedSampleAlg/JpsiToKsKlEtaAlg.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "DatabaseSvc/IDatabaseSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "DstEvent/TofHitStatus.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"

#include "McTruth/McParticle.h"
#include "EvTimeEvent/RecEsTime.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#include "DTagTool/DTagTool.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "CLHEP/Geometry/Point3D.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
#include "GammaConv/GammaConv.h"
#include "MeasuredEcmsSvc/IMeasuredEcmsSvc.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>
#include <set>
const double mpion=0.13957018; //pi+-
const double mpi0=0.1349766;
const double meta=0.547853;
const double mep=0.95778;
const double mk0=0.497611; //ks&kl
const double melectron=0.000510999;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
int no_evt1=0;
JpsiToKsKlEtaAlg::JpsiToKsKlEtaAlg(const std::string& name, ISvcLocator* pSvcLocator):Algorithm(name,pSvcLocator){
	declareProperty("Ecms",            m_Ecms=4.260);
	declareProperty("Debug",           m_debug=0);
}
JpsiToKsKlEtaAlg::~JpsiToKsKlEtaAlg(){
}
StatusCode JpsiToKsKlEtaAlg::initialize(){
	cout<<"test initialize () --- "<<endl; 
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToKsKlEtaAlg::initialize()"<<endreq;
	//add your code here
	StatusCode status;
	NTuplePtr nt1(ntupleSvc(), "FILE1/tree");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/tree", CLID_ColumnWiseTuple, "N-Tuple example");
		if ( m_tuple1 )    {
			status = m_tuple1->addItem ("run",   m_runNo);
			status = m_tuple1->addItem ("evt",   m_evtNo);
			//status = m_tuple1->addItem ("no_pip",    m_no_pip,0,20);
			//status = m_tuple1->addIndexedItem ("p4_pip", m_no_pip, 4, m_p4_pip);
			//status = m_tuple1->addIndexedItem ("pip_id", m_no_pip, m_pip_id);
			//status = m_tuple1->addItem ("no_pim",    m_no_pim,0,20);
			//status = m_tuple1->addIndexedItem ("p4_pim", m_no_pim, 4, m_p4_pim);
			//status = m_tuple1->addIndexedItem ("pim_id", m_no_pim, m_pim_id);
			status = m_tuple1->addItem ("no_gam",  m_no_gam,0,50);
			status = m_tuple1->addIndexedItem ("gam_par", m_no_gam, 6,  m_gam_par);
			//status = m_tuple1->addItem ("no_pi0", m_no_pi0, 0,50);
			//status = m_tuple1->addIndexedItem ("pi0_par",  m_no_pi0, 8, m_pi0_par);
			status = m_tuple1->addItem ("no_eta", m_no_eta, 0,50);
			status = m_tuple1->addIndexedItem ("eta_par",  m_no_eta, 8, m_eta_par);
			//status = m_tuple1->addItem ("no_chrp",   m_no_chrp,0,20);
			//status = m_tuple1->addIndexedItem ("chrp_p3",  m_no_chrp, 3, m_chrp_p3);
			//status = m_tuple1->addIndexedItem ("chrp_id",  m_no_chrp, m_chrp_id);
			//status = m_tuple1->addIndexedItem ("chrp_Eemc",  m_no_chrp, m_chrp_Eemc);
			//status = m_tuple1->addIndexedItem ("chrp_prob",  m_no_chrp, 25,m_chrp_prob);
			//status = m_tuple1->addItem ("no_chrm",   m_no_chrm,0,20);
			//status = m_tuple1->addIndexedItem ("chrm_p3",  m_no_chrm, 3, m_chrm_p3);
			//status = m_tuple1->addIndexedItem ("chrm_id",  m_no_chrm, m_chrm_id);
			//status = m_tuple1->addIndexedItem ("chrm_Eemc",  m_no_chrm, m_chrm_Eemc);
			//status = m_tuple1->addIndexedItem ("chrm_prob",  m_no_chrm, 25,m_chrm_prob);
			//status = m_tuple1->addItem ("no_etaTOpipipi0", m_no_etaTOpipipi0,0,50);
			//status = m_tuple1->addIndexedItem ("M_etaTOpipipi0", m_no_etaTOpipipi0, m_M_etaTOpipipi0);
			status = m_tuple1->addIndexedItem ("pipFeta_id", m_no_etaTOpipipi0, m_pipFeta_id);
			status = m_tuple1->addIndexedItem ("pimFeta_id", m_no_etaTOpipipi0, m_pimFeta_id);
			status = m_tuple1->addIndexedItem ("pi0Feta_id", m_no_etaTOpipipi0, m_pi0Feta_id);
			status = m_tuple1->addIndexedItem("chis_etaTOpipipi0",m_no_etaTOpipipi0,m_chis_etaTOpipipi0);
			status = m_tuple1->addIndexedItem ("p4_pipFeta_kf",m_no_etaTOpipipi0,4,m_p4_pipFeta_kf);
			status = m_tuple1->addIndexedItem ("p4_pimFeta_kf",m_no_etaTOpipipi0,4,m_p4_pimFeta_kf);
			status = m_tuple1->addIndexedItem ("p4_gam1Feta_kf",m_no_etaTOpipipi0,4,m_p4_gam1Feta_kf);
			status = m_tuple1->addIndexedItem ("p4_gam2Feta_kf",m_no_etaTOpipipi0,4,m_p4_gam2Feta_kf);
			status = m_tuple1->addItem ("no_etaTO3pi0", m_no_etaTO3pi0,0,50);
			status = m_tuple1->addIndexedItem ("M_etaTO3pi0", m_no_etaTO3pi0, m_M_etaTO3pi0);
			status = m_tuple1->addIndexedItem ("pi01Feta_id", m_no_etaTO3pi0, m_pi01Feta_id);
			status = m_tuple1->addIndexedItem ("pi02Feta_id", m_no_etaTO3pi0, m_pi02Feta_id);
			status = m_tuple1->addIndexedItem ("pi03Feta_id", m_no_etaTO3pi0, m_pi03Feta_id);
			status = m_tuple1->addIndexedItem ("chis_etaTO3pi0", m_no_etaTO3pi0, m_chis_etaTO3pi0);
			status = m_tuple1->addIndexedItem("p4_gam1Fpi01Feta_kf",m_no_etaTO3pi0,4,m_p4_gam1Fpi01Feta_kf);
			status = m_tuple1->addIndexedItem("p4_gam2Fpi01Feta_kf",m_no_etaTO3pi0,4,m_p4_gam2Fpi01Feta_kf);
			status = m_tuple1->addIndexedItem("p4_gam1Fpi02Feta_kf",m_no_etaTO3pi0,4,m_p4_gam1Fpi02Feta_kf);
			status = m_tuple1->addIndexedItem("p4_gam2Fpi02Feta_kf",m_no_etaTO3pi0,4,m_p4_gam2Fpi02Feta_kf);
			status = m_tuple1->addIndexedItem("p4_gam1Fpi03Feta_kf",m_no_etaTO3pi0,4,m_p4_gam1Fpi03Feta_kf);
			status = m_tuple1->addIndexedItem("p4_gam2Fpi03Feta_kf",m_no_etaTO3pi0,4,m_p4_gam2Fpi03Feta_kf);
			//status = m_tuple1->addItem ("no_phiTOkk", m_no_phiTOkk,0,100);
			//status = m_tuple1->addItem("no_ksTOpipi", m_no_ksTOpipi,0,100);
			//status = m_tuple1->addIndexedItem("chis_phiTOkk",m_no_phiTOkk,m_chis_phiTOkk);
			//status = m_tuple1->addIndexedItem("chis_ksTOpipi", m_no_ksTOpipi, m_chis_ksTOpipi);
			//status = m_tuple1->addIndexedItem ("p4_kapFphi_kf", m_no_phiTOkk,4, m_p4_kapFphi_kf);
			//status = m_tuple1->addIndexedItem("p4_pipFks_kf", m_no_ksTOpipi,4, m_p4_pipFks_kf);
			//status = m_tuple1->addIndexedItem ("p4_kamFphi_kf", m_no_phiTOkk,4, m_p4_kamFphi_kf);
            //status = m_tuple1->addIndexedItem("p4_pimFks_kf", m_no_ksTOpipi,4, m_p4_pimFks_kf);			
			//status = m_tuple1->addIndexedItem ("M_phiTOkk", m_no_phiTOkk, m_M_phiTOkk);
			//status = m_tuple1->addIndexedItem("M_ksTOpipi", m_no_ksTOpipi, m_M_ksTOpipi);
			//status = m_tuple1->addIndexedItem ("rconv", m_no_phiTOkk, m_rconv);
		    //status = m_tuple1->addIndexedItem("rconv", m_no_ksTOpipi, m_rconv);
			//status = m_tuple1->addIndexedItem ("kapFphi_id", m_no_phiTOkk, m_kapFphi_id);
			//status = m_tuple1->addIndexedItem("pipFks_id", m_no_ksTOpipi, m_pipFks_id);
			//status = m_tuple1->addIndexedItem ("kamFphi_id", m_no_phiTOkk, m_kamFphi_id);
            //status = m_tuple1->addIndexedItem("pimFks_id", m_no_ksTOpipi, m_pimFks_id);
			//status = m_tuple1->addIndexedItem ("cos_ee", m_no_phiTOkk, m_cos_ee);
    		//status = m_tuple1->addIndexedItem("cos_ee", m_no_ksTOpipi, m_cos_ee);

			status = m_tuple1->addItem("indexmc", m_idxmc, 0, 100);
			status = m_tuple1->addIndexedItem("pdgid" , m_idxmc, m_pdgid);
			status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);

			//pip, pim from ks0
			status = m_tuple1->addItem(         "npip",         m_npip,         0,      100);
            status = m_tuple1->addIndexedItem(  "idpip",        m_npip,         m_idpip);
            status = m_tuple1->addIndexedItem(  "mpip",         m_npip,         m_mpip);
            status = m_tuple1->addIndexedItem(  "pxpip",        m_npip,         m_pxpip);
            status = m_tuple1->addIndexedItem(  "pypip",        m_npip,         m_pypip);
            status = m_tuple1->addIndexedItem(  "pzpip",        m_npip,         m_pzpip);
            status = m_tuple1->addIndexedItem(  "Epip",         m_npip,         m_Epip);
                        status = m_tuple1->addIndexedItem(      "PipProbK",             m_npip,                 m_PipProbK);
                        status = m_tuple1->addIndexedItem(      "PipProbPi",            m_npip,                 m_PipProbPi);
                        status = m_tuple1->addIndexedItem(      "PipProbP",             m_npip,                 m_PipProbP);

            status = m_tuple1->addItem(         "npim",         m_npim,         0,      100);
            status = m_tuple1->addIndexedItem(  "idpim",        m_npim,         m_idpim);
            status = m_tuple1->addIndexedItem(  "mpim",         m_npim,         m_mpim);
            status = m_tuple1->addIndexedItem(  "pxpim",        m_npim,         m_pxpim);
            status = m_tuple1->addIndexedItem(  "pypim",        m_npim,         m_pypim);
            status = m_tuple1->addIndexedItem(  "pzpim",        m_npim,         m_pzpim);
            status = m_tuple1->addIndexedItem(  "Epim",         m_npim,         m_Epim);
                        status = m_tuple1->addIndexedItem(      "PimProbK",             m_npim,                 m_PimProbK);
                        status = m_tuple1->addIndexedItem(      "PimProbPi",            m_npim,                 m_PimProbPi);
                        status = m_tuple1->addIndexedItem(      "PimProbP",             m_npim,                 m_PimProbP);

            status = m_tuple1->addItem(         "nKS0",         m_nKS0, 0, 100);
            status = m_tuple1->addIndexedItem(  "pxKS0",        m_nKS0,         m_pxKS0);
            status = m_tuple1->addIndexedItem(  "pyKS0",        m_nKS0,         m_pyKS0);
            status = m_tuple1->addIndexedItem(  "pzKS0",        m_nKS0,         m_pzKS0);
            status = m_tuple1->addIndexedItem(  "EKS0",         m_nKS0,         m_EKS0);
            status = m_tuple1->addIndexedItem(  "mKS0",         m_nKS0,         m_mKS0);
            status = m_tuple1->addIndexedItem(  "chi2KS0IP",    m_nKS0,         m_chi2KS0IP);
            status = m_tuple1->addIndexedItem(  "chi2KS0SP",    m_nKS0,         m_chi2KS0SP);
            status = m_tuple1->addIndexedItem(  "LifeKS0",      m_nKS0,         m_LifeKS0);
            status = m_tuple1->addIndexedItem(  "LenKS0",       m_nKS0,         m_LenKS0);
            status = m_tuple1->addIndexedItem(  "LenErrKS0",    m_nKS0,         m_LenErrKS0);
            status = m_tuple1->addIndexedItem(  "IdTrp",        m_nKS0,         m_IdTrp);
            status = m_tuple1->addIndexedItem(  "pxTrp",        m_nKS0,         m_pxTrp);
            status = m_tuple1->addIndexedItem(  "pyTrp",        m_nKS0,         m_pyTrp);
            status = m_tuple1->addIndexedItem(  "pzTrp",        m_nKS0,         m_pzTrp);
            status = m_tuple1->addIndexedItem(  "ETrp",         m_nKS0,         m_ETrp);
            status = m_tuple1->addIndexedItem(  "IdTrm",        m_nKS0,         m_IdTrm);
            status = m_tuple1->addIndexedItem(  "pxTrm",        m_nKS0,         m_pxTrm);
            status = m_tuple1->addIndexedItem(  "pyTrm",        m_nKS0,         m_pyTrm);
            status = m_tuple1->addIndexedItem(  "pzTrm",        m_nKS0,         m_pzTrm);
            status = m_tuple1->addIndexedItem(  "ETrm",         m_nKS0,         m_ETrm);

		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}
	//--------end of book--------

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;

}
//******************************************************
StatusCode JpsiToKsKlEtaAlg::beginRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToKsKlEtaAlg::beginRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;
}

StatusCode JpsiToKsKlEtaAlg::execute(){
	if(m_debug){
		no_evt1++;
		cout<<"no_evt1="<<no_evt1<<endl;
	}  
	if(m_debug) cout<<"begin execute()--"<<endl;
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToKsKlEtaAlg::execute()"<<endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	IMeasuredEcmsSvc* ecmsSvc;
	StatusCode status=serviceLocator()->service("MeasuredEcmsSvc", ecmsSvc, true);
	if(!status.isSuccess()) return status;
	int runNo = eventHeader->runNumber();
	int eventNo = eventHeader->eventNumber();

	if(m_debug) cout<<"runNo="<<runNo<<" eventNo="<<eventNo<<endl;

	//MC truth
	SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
	if(mcParticleCol){
		Event::McParticleCol::iterator iter_mc=mcParticleCol->begin();
		m_idxmc=8;
		const int incPid=443; //443 is the PDG code of J/psi
		bool incPdcy(false);
		int rootIndex(-1);
		for (;iter_mc!=mcParticleCol->end();iter_mc++)
		{
			if((*iter_mc)->primaryParticle()) continue;
			if(!(*iter_mc)->decayFromGenerator()) continue;
			if((*iter_mc)->particleProperty()==incPid)
			{
				incPdcy=true;
				rootIndex=(*iter_mc)->trackIndex();
			}
			if(!incPdcy) continue; 
			m_pdgid[m_idxmc]=(*iter_mc)->particleProperty();
			if((*iter_mc)->particleProperty()==incPid||((*iter_mc)->mother()).particleProperty()==incPid) 
				m_motheridx[m_idxmc]=((*iter_mc)->mother()).trackIndex()-rootIndex;
			else 
				m_motheridx[m_idxmc]=((*iter_mc)->mother()).trackIndex()-rootIndex-1;
			m_idxmc++;
		}
	}

	//select good tracks
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	log << MSG::DEBUG <<"ncharg, nneu, tottks = "
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() <<endreq;

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

	Vint pipID; pipID.clear();
	Vint pimID; pimID.clear();
	//Vint kapID; kapID.clear();
	//Vint kamID; kamID.clear();
	Vint chrpID; chrpID.clear();
	Vint chrmID; chrmID.clear();
	for(int i=0; i<evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		bool bool_tp=goodTrk(itTrk); //select good trks:vr\vz\cos-theta
		if(m_debug) cout<<i<<" bool_tp="<<bool_tp<<endl;
		if(!bool_tp) continue;

		//distinguish chr+ chr-
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		int chr = mdcTrk->charge();
		if(m_debug) cout<<i<<" chr="<<chr<<endl;
		if(chr>0)  chrpID.push_back(i);
		else  chrmID.push_back(i);

		//PID
		ISimplePIDSvc*  m_simplePIDSvc;
		Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);
		m_simplePIDSvc->preparePID(*itTrk);
		bool pion_pid= m_simplePIDSvc->ispion();
		double prob[25];

		PID(itTrk,prob);
		if (prob[2] > prob[3] && prob[2] > prob[4] && prob[3] > 0.001) bool pion_pid = 1;

		//bool kaon_pid= m_simplePIDSvc->iskaon();
		if(pion_pid) {
			if(chr>0) pipID.push_back(i);
			else pimID.push_back(i); 
		} 
		//else if(kaon_pid){
		//	if(chr>0) kapID.push_back(i);
		//	else kamID.push_back(i);
		//}
	}
	if(m_debug){
		cout<<"totaltrk: no_chrp="<<chrpID.size()<<"; no_chrm="<<chrmID.size()<<endl;
		cout<<"   pions: no_pip="<<pipID.size()<<"; no_pim="<<pimID.size()<<endl;
		//cout<<"   kaons: no_kap="<<kapID.size()<<"; no_kam="<<kamID.size()<<endl;
	}  

	int no_chrp=chrpID.size();
	int no_chrm=chrmID.size();
	if(no_chrp>20||no_chrm>20) return StatusCode::SUCCESS;
	int no_pip=pipID.size();
	int no_pim=pimID.size();
	//int no_kap=kapID.size();
	//int no_kam=kamID.size();

	m_runNo = runNo;
	m_evtNo= eventNo;

	m_no_chrp=no_chrp;  m_no_chrm=no_chrm;
	m_no_pip=no_pip;  m_no_pim=no_pim;
	//m_no_kap=no_kap;  m_no_kam=no_kam;

	//***save information for K+K-
	/*
	for(int i=0;i<no_kap;i++){
		int kap_id=kapID[i];
		if(m_debug) cout<<"kap_id="<<kap_id<<endl;
		RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
		EvtRecTrackIterator itTrk_kap=evtRecTrkCol->begin() + kap_id;
		RecMdcKalTrack* mdcKalTrk_kap = (*itTrk_kap)->mdcKalTrack();
		HepLorentzVector p4_kap_tp=mdcKalTrk_kap->p4(mkaon);
		for(int bb=0;bb<4;bb++) { m_p4_kap[i][bb]=p4_kap_tp[bb];}
		m_kap_id[i]=kap_id;
	}   
	for(int i=0;i<no_kam;i++){
		int kam_id=kamID[i];
		if(m_debug) cout<<"kam_id="<<kam_id<<endl;
		RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
		EvtRecTrackIterator itTrk_kam=evtRecTrkCol->begin() + kam_id;
		RecMdcKalTrack* mdcKalTrk_kam = (*itTrk_kam)->mdcKalTrack();
		HepLorentzVector p4_kam_tp=mdcKalTrk_kam->p4(mkaon);
		for(int bb=0;bb<4;bb++) {m_p4_kam[i][bb]=p4_kam_tp[bb];}
		m_kam_id[i]=kam_id;
	}
	*/

	//*****save information for pi+pi-
	/*
	for(int i=0;i<no_pip;i++){
		int pip_id=pipID[i];
		if(m_debug) cout<<"pip_id="<<pip_id<<endl;
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		EvtRecTrackIterator itTrk_pip=evtRecTrkCol->begin() + pip_id;
		RecMdcKalTrack* mdcKalTrk_pip = (*itTrk_pip)->mdcKalTrack();
		HepLorentzVector p4_pip_tp=mdcKalTrk_pip->p4(mpion);
		for(int bb=0;bb<4;bb++)  {m_p4_pip[i][bb]=p4_pip_tp[bb];}
		m_pip_id[i]=pip_id;
	}
	for(int i=0;i<no_pim;i++){
		int pim_id=pimID[i];
		if(m_debug) cout<<"pim_id="<<pim_id<<endl;
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		EvtRecTrackIterator itTrk_pim=evtRecTrkCol->begin() + pim_id;
		RecMdcKalTrack* mdcKalTrk_pim = (*itTrk_pim)->mdcKalTrack();
		HepLorentzVector p4_pim_tp=mdcKalTrk_pim->p4(mpion);
		for(int bb=0;bb<4;bb++)  {m_p4_pim[i][bb]=p4_pim_tp[bb];}
		m_pim_id[i]=pim_id;
	}
	*/

	//good showers
	Vint goodGam; goodGam.clear();
	int no_gam=0;
	for(int i=evtRecEvent->totalCharged();i<evtRecEvent->totalTracks();i++){
		//if(m_debug) cout<<"i="<<i<<" no_gam="<<no_gam<<endl;
		EvtRecTrackIterator itTrk_shower=evtRecTrkCol->begin() + i;
		RecEmcShower* shower = (*itTrk_shower)->emcShower();
		double dang_min=-100;
		if(good_gam(shower,dang_min)) {
			if(m_debug) cout<<" dang_min="<<dang_min<<endl;
			goodGam.push_back(i);
			HepLorentzVector p4_shw = getP4(shower);
			for(int ll=0;ll<4;ll++) m_gam_par[no_gam][ll]=p4_shw[ll];
			m_gam_par[no_gam][4]=i;
			m_gam_par[no_gam][5]=dang_min;
			no_gam++;
		}
		if(no_gam>48) break;
	}
	if(m_debug) cout<<"no_gam="<<no_gam<<endl;
	if(no_gam>48) return StatusCode::SUCCESS;
	m_no_gam=no_gam;

	//*****eta->gamma gamma
	if(m_debug) cout<<"#####eta----"<<endl;
	if(m_debug) cout<<"no_gam="<<no_gam<<endl;
	ss=-1;
	if(m_debug) cout<<"m_no_pi0="<<m_no_pi0<<endl;
	for(int gg=0;gg<no_gam;gg++){
		if(m_debug) cout<<"gg="<<gg<<endl;
		EvtRecTrackIterator itTrk_shr1=evtRecTrkCol->begin() + goodGam[gg];

		double px_tp=m_gam_par[gg][0];
		RecEmcShower* shr1 = (*itTrk_shr1)->emcShower();
		HepLorentzVector p4_shr1 = getP4(shr1);
		double px_shr1=p4_shr1.px();
		double dang_min_tp=m_gam_par[gg][5];
		if(m_debug) cout<<"px_tp="<<px_tp<<" px_shr1="<<px_shr1<<"; dang_min_tp="<<dang_min_tp<<endl;
		if(dang_min_tp>-10&&dang_min_tp<10) continue;
		for(int hh=gg+1;hh<no_gam;hh++){
			if(m_debug) cout<<"hh="<<hh<<endl;
			EvtRecTrackIterator itTrk_shr2=evtRecTrkCol->begin() + goodGam[hh];
			double px_tp2=m_gam_par[hh][0];
			RecEmcShower* shr2 = (*itTrk_shr2)->emcShower();
			HepLorentzVector p4_shr2 = getP4(shr2);
			double px_shr2=p4_shr2.px();
			double dang_min_tp2=m_gam_par[hh][5];
			if(m_debug) cout<<"px_tp2="<<px_tp2<<" px_shr2="<<px_shr2<<"; dang_min_tp2="<<dang_min_tp2<<endl;
			if(dang_min_tp2>-10&&dang_min_tp2<10) continue;

			HepLorentzVector p4_eta=p4_shr1+p4_shr2;
			double eta_m=p4_eta.m();
			if(m_debug) cout<<"eta_m="<<eta_m<<endl;
			if(eta_m<0.48||eta_m>0.59) continue;
			else ss++;
			if(m_debug) cout<<"ss="<<ss<<endl;
			if(ss>48) return StatusCode::SUCCESS;

			HepLorentzVector p4_eta_tp,p4_eta_1c;
			double eta_chis=-100;
			if(m_debug) cout<<"before saveeta: eta_chis="<<eta_chis<<endl;
			saveeta(shr1,shr2,eta_chis,p4_eta_tp,p4_eta_1c);
			if(m_debug) cout<<"after saveeta: eta_chis="<<eta_chis<<endl;

			for(int ll=0;ll<4;ll++) m_eta_par[ss][ll]=p4_eta_1c[ll];
			m_eta_par[ss][4]=eta_m;
			m_eta_par[ss][5]=eta_chis;
			m_eta_par[ss][6]=goodGam[gg];
			m_eta_par[ss][7]=goodGam[hh];
		}
	}
	int no_eta=ss+1;
	m_no_eta=no_eta;

	bool etaTOgg=false;
	if(no_eta>0) etaTOgg=true;

	//Ks->pi+pi-
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
    //cout << NumTrp << "   " << NumTrm << endl;

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

}//end of execute()
//**************************************************************
StatusCode JpsiToKsKlEtaAlg::endRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToKsKlEtaAlg::endRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;

}
StatusCode JpsiToKsKlEtaAlg::finalize(){
	if(m_debug)  cout<<"finalize: no_evt1="<<no_evt1<<endl;
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToKsKlEtaAlg::finalize()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;
}

//add your code here,for chrer member-functions
bool JpsiToKsKlEtaAlg::goodTrk(EvtRecTrackIterator itTrk) {
	double VrCut      = 1.0;
	double VzCut      = 10.0;
	double CosThetaCut= 0.93;
	if ( !(*itTrk)->isMdcTrackValid() ) return false;
	RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
	if ( mdcTrk->charge()==0 ) return false;

	Hep3Vector xorigin(0,0,0);
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
		double* dbv = vtxsvc->PrimaryVertex();
		double*  vv = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}
	HepVector a = mdcTrk->helix();
	HepSymMatrix Ea = mdcTrk->err();
	HepPoint3D point0(0.,0.,0.);
	HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
	VFHelix helixip3(point0,a,Ea);
	helixip3.pivot(IP);
	HepVector  vecipa = helixip3.a();

	double dr=fabs(vecipa[0]);
	double dz=fabs(vecipa[3]);
	double costheta=cos(mdcTrk->theta());
	if (  dr>= VrCut ) return false;
	if (  dz>= VzCut ) return false;
	if ( fabs(costheta) >= CosThetaCut ) return false;

	return true;
}

bool JpsiToKsKlEtaAlg::good_gam(RecEmcShower *emcTrk,double& dang_min){

	SmartDataPtr<EvtRecEvent> recEvt(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> recTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);

	double  m_minEnergy = 0.025;
	bool  m_useBarrelEndcap   = true;
	double  m_maxCosThetaBarrel = 0.8;
	double  m_minCosThetaEndcap = 0.86;
	double  m_maxCosThetaEndcap = 0.92;
	double  m_minEndcapEnergy   = 0.050;
	bool  m_applyTimeCut = true;
	double  m_minTime      = 0.;
	double  m_maxTime      = 14.;
	double  m_applyDangCut = false;
	double  m_minDang      = 10.0;

	double eraw = emcTrk->energy();
	double phi =  emcTrk->phi();
	double the =  emcTrk->theta();
	HepLorentzVector shP4( eraw * sin(the) * cos(phi),
			eraw * sin(the) * sin(phi),
			eraw * cos(the),
			eraw );
	double cosThetaSh = shP4.vect().cosTheta();

	/// Minimum energy
	if (shP4.e() <= m_minEnergy) return( false );
	/// Barrel/Endcap
	bool inBarrelEndcap = false;
	if(fabs(cosThetaSh) < m_maxCosThetaBarrel) inBarrelEndcap = true;
	if((fabs(cosThetaSh) > m_minCosThetaEndcap)
			&&(fabs(cosThetaSh) < m_maxCosThetaEndcap)
			&&(shP4.e() > m_minEndcapEnergy)) inBarrelEndcap = true;
	if(m_useBarrelEndcap&&!inBarrelEndcap) return( false );


	/// Time, only apply timing cuts if "recEvt->totalCharged() > 0"
	if ( m_applyTimeCut ) {
		double time = emcTrk->time();
		if ( recEvt->totalCharged() > 0 ) {
			if ( time < m_minTime || time > m_maxTime ) return false;
		}
		else {
			RecEmcShower* firstG = (*(recTrkCol->begin()))->emcShower();
			double deltaTime = fabs(time - firstG->time());
			if ( deltaTime > 10 ) return false;
		}
	}
	/// Dang
	Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
	double dang = 200.;
	if (recEvt->totalCharged() > 0)
	{
		for (int j = 0; j < recEvt->totalCharged(); j++) {
			EvtRecTrackIterator jtTrk = recTrkCol->begin() + j;
			if ( !(*jtTrk)->isExtTrackValid() ) continue;
			RecExtTrack* extTrk = (*jtTrk)->extTrack();
			if ( extTrk->emcVolumeNumber() == -1 ) continue;
			Hep3Vector extpos = extTrk->emcPosition();
			double  angd1 = extpos.angle(emcpos);
			if ( angd1 < dang ) dang = angd1;
		}
		if(dang<200) dang = dang * 180 / (CLHEP::pi);
		else dang=-200;
		if (m_applyDangCut&&(fabs(dang) <= m_minDang)) return( false );
	}  // End of "recEvt->totalCharged() > 0" IF

	dang_min=dang;

	return( true );
}
HepLorentzVector JpsiToKsKlEtaAlg::getP4(RecEmcShower* gTrk){

	double eraw = gTrk->energy();
	double phi =  gTrk->phi();
	double the =  gTrk->theta();

	return HepLorentzVector( eraw * sin(the) * cos(phi),
			eraw * sin(the) * sin(phi),
			eraw * cos(the),
			eraw );
}
void JpsiToKsKlEtaAlg::savepi0(RecEmcShower *shr1,RecEmcShower *shr2,double& pi0_chis,HepLorentzVector& p4_pi0,HepLorentzVector& p4_pi0_1c){
	double xmpi0=0.1349766;
	pi0_chis=-100;
	bool pi0_sta=0;
	HepLorentzVector g1P4 = getP4(shr1);
	HepLorentzVector g2P4 = getP4(shr2);
	p4_pi0 = g1P4 + g2P4;

	KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
	kmfit->init();
	kmfit->setIterNumber(5);
	kmfit->AddTrack(0, 0.0, shr1);
	kmfit->AddTrack(1, 0.0, shr2);
	kmfit->AddResonance(0, xmpi0, 0, 1);

	bool oksq =kmfit->Fit(0);
	pi0_sta=oksq;
	if(m_debug)  cout<<"oksq = "<<oksq<<endl;
	kmfit->BuildVirtualParticle(0);
	pi0_chis = kmfit->chisq(0);
	p4_pi0_1c=kmfit->pfit(0)+kmfit->pfit(1);

	if(m_debug) {
		cout<<"pi0_chis ="<<pi0_chis<<endl;
		cout<<"pi0_m_1c = "<<p4_pi0_1c.m()<<endl;
	}
}
void JpsiToKsKlEtaAlg::saveeta(RecEmcShower *shr1,RecEmcShower *shr2,double& eta_chis,HepLorentzVector& p4_eta,HepLorentzVector& p4_eta_1c){
	double xmeta=0.547853;
	eta_chis=-100;
	bool eta_sta=0;
	HepLorentzVector g1P4 = getP4(shr1);
	HepLorentzVector g2P4 = getP4(shr2);
	p4_eta = g1P4 + g2P4;

	KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
	kmfit->init();
	kmfit->setIterNumber(5);
	kmfit->AddTrack(0, 0.0, shr1);
	kmfit->AddTrack(1, 0.0, shr2);
	kmfit->AddResonance(0, xmeta, 0, 1);

	bool oksq =kmfit->Fit(0);
	eta_sta=oksq;
	if(m_debug)  cout<<"oksq = "<<oksq<<endl;
	kmfit->BuildVirtualParticle(0);
	eta_chis = kmfit->chisq(0);
	if(m_debug) cout<<"eta_chis = "<<eta_chis<<endl;
	p4_eta_1c=kmfit->pfit(0)+kmfit->pfit(1);

	if(m_debug) {
		cout<<"eta_chis ="<<eta_chis<<endl;
		cout<<"eta_m_1c = "<<p4_eta_1c.m()<<endl;
	}
}
void JpsiToKsKlEtaAlg::PID(EvtRecTrackIterator itTrk,double prob[25]){
	ParticleID *pid = ParticleID::instance();
	pid->init();
	pid->setMethod(pid->methodProbability());
	pid->setChiMinCut(4);
	pid->setRecTrack(*itTrk);
	pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE());
	pid->identify(pid->onlyElectron() | pid->onlyMuon() | pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
	pid->calculate();
	if( (pid->IsPidInfoValid()) ){
		prob[0] = pid->probElectron();
		prob[1] = pid->probMuon();
		prob[2] = pid->probPion();
		prob[3] = pid->probKaon();
		prob[4] = pid->probProton();
		if(m_debug) cout<<"prob. of e/mu/pi/K/P : "<<prob[0]<<"; "<<prob[1]<<"; "<<prob[2]<<"; "<<prob[3]<<"; "<<prob[4]<<endl;
		prob[5]= pid->chiDedx(0);
		prob[6]= pid->chiTof1(0);
		prob[7]= pid->chiTof2(0);
		prob[8]= pid->chiTofE(0);

		prob[9]= pid->chiDedx(1);
		prob[10]= pid->chiTof1(1);
		prob[11]= pid->chiTof2(1);
		prob[12]= pid->chiTofE(1);

		prob[13]= pid->chiDedx(2);
		prob[14]= pid->chiTof1(2);
		prob[15]= pid->chiTof2(2);
		prob[16]= pid->chiTofE(2);

		prob[17]= pid->chiDedx(3);
		prob[18]= pid->chiTof1(3);
		prob[19]= pid->chiTof2(3);
		prob[20]= pid->chiTofE(3);

		prob[21]= pid->chiDedx(4);
		prob[22]= pid->chiTof1(4);
		prob[23]= pid->chiTof2(4);
		prob[24]= pid->chiTofE(4);
	}
	else {
		for(int nn=0;nn<25;nn++) prob[nn]=1000;
	}
	if(m_debug) cout<<"in PID(): prob[0]="<<prob[0]<<"; prob[24]="<<prob[24]<<endl;
}
