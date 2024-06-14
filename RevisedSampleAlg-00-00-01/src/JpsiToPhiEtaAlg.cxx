#include "SampleAlg/JpsiToPhiEtaAlg.h"

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
const double mkaon=0.493677; //k+-
const double mpi0=0.1349766;
const double meta=0.547853;
const double mep=0.95778;
const double mphi=1.019455;
const double melectron=0.000510999;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
int no_evt1=0;
JpsiToPhiEtaAlg::JpsiToPhiEtaAlg(const std::string& name, ISvcLocator* pSvcLocator):Algorithm(name,pSvcLocator){
	declareProperty("Ecms",            m_Ecms=4.260);
	declareProperty("Debug",           m_debug=0);
}
JpsiToPhiEtaAlg::~JpsiToPhiEtaAlg(){
}
StatusCode JpsiToPhiEtaAlg::initialize(){
	cout<<"test initialize () --- "<<endl; 
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToPhiEtaAlg::initialize()"<<endreq;
	//add your code here
	StatusCode status;
	NTuplePtr nt1(ntupleSvc(), "FILE1/tree");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/tree", CLID_ColumnWiseTuple, "N-Tuple example");
		if ( m_tuple1 )    {
			status = m_tuple1->addItem ("run",   m_runNo);
			status = m_tuple1->addItem ("evt",   m_evtNo);
			status = m_tuple1->addItem ("no_pip",    m_no_pip,0,20);
			status = m_tuple1->addIndexedItem ("p4_pip", m_no_pip, 4, m_p4_pip);
			status = m_tuple1->addIndexedItem ("pip_id", m_no_pip, m_pip_id);
			status = m_tuple1->addItem ("no_pim",    m_no_pim,0,20);
			status = m_tuple1->addIndexedItem ("p4_pim", m_no_pim, 4, m_p4_pim);
			status = m_tuple1->addIndexedItem ("pim_id", m_no_pim, m_pim_id);
			status = m_tuple1->addItem ("no_kap",    m_no_kap,0,20);
			status = m_tuple1->addIndexedItem ("p4_kap", m_no_kap, 4, m_p4_kap);
			status = m_tuple1->addIndexedItem ("kap_id", m_no_kap, m_kap_id);
			status = m_tuple1->addItem ("no_kam",    m_no_kam,0,20);
			status = m_tuple1->addIndexedItem ("p4_kam", m_no_kam, 4, m_p4_kam);
			status = m_tuple1->addIndexedItem ("kam_id", m_no_kam, m_kam_id);
			status = m_tuple1->addItem ("no_gam",  m_no_gam,0,50);
			status = m_tuple1->addIndexedItem ("gam_par", m_no_gam, 6,  m_gam_par);
			status = m_tuple1->addItem ("no_pi0", m_no_pi0, 0,50);
			status = m_tuple1->addIndexedItem ("pi0_par",  m_no_pi0, 8, m_pi0_par);
			status = m_tuple1->addItem ("no_eta", m_no_eta, 0,50);
			status = m_tuple1->addIndexedItem ("eta_par",  m_no_eta, 8, m_eta_par);
			status = m_tuple1->addItem ("no_chrp",   m_no_chrp,0,20);
			status = m_tuple1->addIndexedItem ("chrp_p3",  m_no_chrp, 3, m_chrp_p3);
			status = m_tuple1->addIndexedItem ("chrp_id",  m_no_chrp, m_chrp_id);
			status = m_tuple1->addIndexedItem ("chrp_Eemc",  m_no_chrp, m_chrp_Eemc);
			status = m_tuple1->addIndexedItem ("chrp_prob",  m_no_chrp, 25,m_chrp_prob);
			status = m_tuple1->addItem ("no_chrm",   m_no_chrm,0,20);
			status = m_tuple1->addIndexedItem ("chrm_p3",  m_no_chrm, 3, m_chrm_p3);
			status = m_tuple1->addIndexedItem ("chrm_id",  m_no_chrm, m_chrm_id);
			status = m_tuple1->addIndexedItem ("chrm_Eemc",  m_no_chrm, m_chrm_Eemc);
			status = m_tuple1->addIndexedItem ("chrm_prob",  m_no_chrm, 25,m_chrm_prob);
			status = m_tuple1->addItem ("no_etaTOpipipi0", m_no_etaTOpipipi0,0,50);
			status = m_tuple1->addIndexedItem ("M_etaTOpipipi0", m_no_etaTOpipipi0, m_M_etaTOpipipi0);
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
			status = m_tuple1->addItem ("no_phiTOkk", m_no_phiTOkk,0,100);
			status = m_tuple1->addIndexedItem("chis_phiTOkk",m_no_phiTOkk,m_chis_phiTOkk);
			status = m_tuple1->addIndexedItem ("p4_kapFphi_kf", m_no_phiTOkk,4, m_p4_kapFphi_kf);
			status = m_tuple1->addIndexedItem ("p4_kamFphi_kf", m_no_phiTOkk,4, m_p4_kamFphi_kf);
			status = m_tuple1->addIndexedItem ("M_phiTOkk", m_no_phiTOkk, m_M_phiTOkk);
			status = m_tuple1->addIndexedItem ("rconv", m_no_phiTOkk, m_rconv);
			status = m_tuple1->addIndexedItem ("kapFphi_id", m_no_phiTOkk, m_kapFphi_id);
			status = m_tuple1->addIndexedItem ("kamFphi_id", m_no_phiTOkk, m_kamFphi_id);
			status = m_tuple1->addIndexedItem ("cos_ee", m_no_phiTOkk, m_cos_ee);

			status = m_tuple1->addItem("indexmc", m_idxmc, 0, 100);
			status = m_tuple1->addIndexedItem("pdgid" , m_idxmc, m_pdgid);
			status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);

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
StatusCode JpsiToPhiEtaAlg::beginRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToPhiEtaAlg::beginRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;
}

StatusCode JpsiToPhiEtaAlg::execute(){
	if(m_debug){
		no_evt1++;
		cout<<"no_evt1="<<no_evt1<<endl;
	}  
	if(m_debug) cout<<"begin execute()--"<<endl;
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToPhiEtaAlg::execute()"<<endreq;
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
		m_idxmc=0;
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
	Vint kapID; kapID.clear();
	Vint kamID; kamID.clear();
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
		bool kaon_pid= m_simplePIDSvc->iskaon();
		if(pion_pid) {
			if(chr>0) pipID.push_back(i);
			else pimID.push_back(i); 
		} 
		else if(kaon_pid){
			if(chr>0) kapID.push_back(i);
			else kamID.push_back(i);
		}
	}
	if(m_debug){
		cout<<"totaltrk: no_chrp="<<chrpID.size()<<"; no_chrm="<<chrmID.size()<<endl;
		cout<<"   pions: no_pip="<<pipID.size()<<"; no_pim="<<pimID.size()<<endl;
		cout<<"   kaons: no_kap="<<kapID.size()<<"; no_kam="<<kamID.size()<<endl;
	}  

	int no_chrp=chrpID.size();
	int no_chrm=chrmID.size();
	if(no_chrp>20||no_chrm>20) return StatusCode::SUCCESS;
	int no_pip=pipID.size();
	int no_pim=pimID.size();
	int no_kap=kapID.size();
	int no_kam=kamID.size();

	m_runNo = runNo;
	m_evtNo= eventNo;

	m_no_chrp=no_chrp;  m_no_chrm=no_chrm;
	m_no_pip=no_pip;  m_no_pim=no_pim;
	m_no_kap=no_kap;  m_no_kam=no_kam;

	//***save information for K+K-
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

	//*****save information for pi+pi-
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

	//*****pi0
	if(m_debug) cout<<"#####pi0----"<<endl;
	int ss=-1;
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

			HepLorentzVector p4_pi0=p4_shr1+p4_shr2;
			double pi0_m=p4_pi0.m();
			if(m_debug) cout<<"pi0_m="<<pi0_m<<endl;
			if(pi0_m<0.105||pi0_m>0.160) continue;
			else ss++;
			if(ss>48) return StatusCode::SUCCESS;

			HepLorentzVector p4_pi0_tp,p4_pi0_1c;
			double pi0_chis=-100;
			if(m_debug) cout<<"before savepi0: pi0_chis="<<pi0_chis<<endl;
			savepi0(shr1,shr2,pi0_chis,p4_pi0_tp,p4_pi0_1c);
			if(m_debug) cout<<"after savepi0: pi0_chis="<<pi0_chis<<endl;

			for(int ll=0;ll<4;ll++) m_pi0_par[ss][ll]=p4_pi0_1c[ll];
			m_pi0_par[ss][4]=pi0_m;
			m_pi0_par[ss][5]=pi0_chis;
			m_pi0_par[ss][6]=goodGam[gg];
			m_pi0_par[ss][7]=goodGam[hh];
		}
	}
	int no_pi0=ss+1;
	m_no_pi0=no_pi0;
	if(m_debug) cout<<"ss="<<ss<<endl;

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

	//all chr+
	for(int tt=0;tt<no_chrp;tt++){
		EvtRecTrackIterator itTrk_tt = evtRecTrkCol->begin() + chrpID[tt];
		RecMdcTrack* mdcTrk = (*itTrk_tt)->mdcTrack();
		Hep3Vector p3_chrp=mdcTrk->p3();
		for(int cc=0;cc<3;cc++) m_chrp_p3[tt][cc]=p3_chrp[cc];
		m_chrp_id[tt]=chrpID[tt];
		//deposit energy in EMC
		if((*itTrk_tt)->isEmcShowerValid()){
			RecEmcShower *emcTrk_tt = (*itTrk_tt)->emcShower();
			double Eemc_tt=emcTrk_tt->energy();
			m_chrp_Eemc[tt]=Eemc_tt;
		}
		else m_chrp_Eemc[tt]=-100;
		//PID for chr+
		double prob_p[25];
		for(int aa=0;aa<25;aa++) prob_p[aa]=-100;
		PID(itTrk_tt,prob_p);
		for(int ll=0;ll<25;ll++)  m_chrp_prob[tt][ll]=prob_p[ll];
	}
	//all chr-
	for(int tt=0;tt<no_chrm;tt++){
		EvtRecTrackIterator itTrk_tt = evtRecTrkCol->begin() + chrmID[tt];
		RecMdcTrack* mdcTrk = (*itTrk_tt)->mdcTrack();
		Hep3Vector p3_chrm=mdcTrk->p3();
		for(int cc=0;cc<3;cc++) m_chrm_p3[tt][cc]=p3_chrm[cc];
		m_chrm_id[tt]=chrmID[tt];
		//deposit energy in EMC
		if((*itTrk_tt)->isEmcShowerValid()){
			RecEmcShower *emcTrk_tt = (*itTrk_tt)->emcShower();
			double Eemc_tt=emcTrk_tt->energy();
			m_chrm_Eemc[tt]=Eemc_tt;
		}
		else m_chrm_Eemc[tt]=-100;
		//PID for chr-
		double prob_m[25];
		for(int aa=0;aa<25;aa++) prob_m[aa]=-100;
		PID(itTrk_tt,prob_m);
		for(int ll=0;ll<25;ll++)  m_chrm_prob[tt][ll]=prob_m[ll];
	}


	//eta->pi+pi-pi0
	ss=-1;
	for(int i=0;i<no_pip;i++){
		int pip_id=pipID[i];
		if(m_debug) cout<<"pip_id="<<pip_id<<endl;
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		EvtRecTrackIterator itTrk_pip=evtRecTrkCol->begin() + pip_id;
		RecMdcKalTrack* mdcKalTrk_pip = (*itTrk_pip)->mdcKalTrack();
		HepLorentzVector p4_pip_tp=mdcKalTrk_pip->p4(mpion);
		WTrackParameter wvpipFeta;
		wvpipFeta=WTrackParameter(mpion, mdcKalTrk_pip->getZHelix(), mdcKalTrk_pip->getZError());

		for(int j=0;j<no_pim;j++){
			int pim_id=pimID[j];
			if(m_debug) cout<<"pim_id="<<pim_id<<endl;
			RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
			EvtRecTrackIterator itTrk_pim=evtRecTrkCol->begin() + pim_id;
			RecMdcKalTrack* mdcKalTrk_pim = (*itTrk_pim)->mdcKalTrack();
			HepLorentzVector p4_pim_tp=mdcKalTrk_pim->p4(mpion);
			WTrackParameter wvpimFeta;
			wvpimFeta=WTrackParameter(mpion, mdcKalTrk_pim->getZHelix(), mdcKalTrk_pim->getZError());

			for(int k=0;k<no_pi0;k++){
				HepLorentzVector p4_pi0_tp;
				for(int ll=0;ll<4;ll++) p4_pi0_tp[ll]=m_pi0_par[k][ll];

				HepLorentzVector p4_etaTOpipipi0=p4_pip_tp+p4_pim_tp+p4_pi0_tp;
				double m_etaTOpipipi0_tp=p4_etaTOpipipi0.m();
				if(m_debug) cout<<"m_etaTOpipipi0_tp="<<m_etaTOpipipi0_tp<<endl;
				if(m_etaTOpipipi0_tp>0.48&&m_etaTOpipipi0_tp<0.59) {
					ss++;
					if(ss>48) return StatusCode::SUCCESS;
					m_M_etaTOpipipi0[ss]=m_etaTOpipipi0_tp;  
					m_pipFeta_id[ss]=i;
					m_pimFeta_id[ss]=j;				      
					m_pi0Feta_id[ss]=k;
					if(m_debug) cout<<"ss="<<ss<<endl;
				}
				else continue;

				int gam1Fpi0_id=m_pi0_par[k][6];
				int gam2Fpi0_id=m_pi0_par[k][7];
				EvtRecTrackIterator itTrk_shower1=evtRecTrkCol->begin() + gam1Fpi0_id;
				RecEmcShower* shower1 = (*itTrk_shower1)->emcShower();
				EvtRecTrackIterator itTrk_shower2=evtRecTrkCol->begin() + gam2Fpi0_id;
				RecEmcShower* shower2 = (*itTrk_shower2)->emcShower();

				//kinematic fit for eta->pi+pi-pi0
				KalmanKinematicFit * kmfit_eta = KalmanKinematicFit::instance();
				kmfit_eta->init();
				kmfit_eta->setChisqCut(2500);
				kmfit_eta->AddTrack(0, wvpipFeta);
				kmfit_eta->AddTrack(1, wvpimFeta);
				kmfit_eta->AddTrack(2, 0.0, shower1);
				kmfit_eta->AddTrack(3, 0.0, shower2);
				kmfit_eta->AddResonance(0, mpi0, 2,3);
				kmfit_eta->AddResonance(1, meta, 0,1,2,3);
				kmfit_eta->Fit();
				double chis_eta = kmfit_eta->chisq();
				if(m_debug) cout<<"chis_eta="<<chis_eta<<endl;
				m_chis_etaTOpipipi0[ss]=chis_eta;

				HepLorentzVector p4_pipFeta_kf=kmfit_eta->pfit(0);
				HepLorentzVector p4_pimFeta_kf=kmfit_eta->pfit(1);
				HepLorentzVector p4_gam1Feta_kf=kmfit_eta->pfit(2);
				HepLorentzVector p4_gam2Feta_kf=kmfit_eta->pfit(3);
				for(int xx=0;xx<4;xx++) { 
					m_p4_pipFeta_kf[ss][xx]=p4_pipFeta_kf[xx];
					m_p4_pimFeta_kf[ss][xx]=p4_pimFeta_kf[xx];
					m_p4_gam1Feta_kf[ss][xx]=p4_gam1Feta_kf[xx];
					m_p4_gam2Feta_kf[ss][xx]=p4_gam2Feta_kf[xx];
				}
				HepLorentzVector p4_etaTOpipipi0_kf=p4_pipFeta_kf+p4_pimFeta_kf+p4_gam1Feta_kf+p4_gam2Feta_kf;
				double m_eta_1c=p4_etaTOpipipi0_kf.m();
				if(m_debug) cout<<"eta->pipipi0: m_eta_1c="<<m_eta_1c<<endl;
			}
		}
	}
	m_no_etaTOpipipi0=ss+1;
	bool etaTOpipipi0=false;
	if(m_no_etaTOpipipi0>0) etaTOpipipi0=true;

	//eta->3pi0
	int mm=-1;
	for(int i=0;i<no_pi0;i++){
		for(int j=i+1;j<no_pi0;j++){
			for(int k=j+1;k<no_pi0;k++){
				HepLorentzVector p4_pi0_1,p4_pi0_2,p4_pi0_3;
				for(int ll=0;ll<4;ll++){
					p4_pi0_1[ll]=m_pi0_par[i][ll];
					p4_pi0_2[ll]=m_pi0_par[j][ll];
					p4_pi0_3[ll]=m_pi0_par[k][ll];
				}	
				//ensure not the same gammas
				int gam1_Fpi01=m_pi0_par[i][6]; 
				int gam2_Fpi01=m_pi0_par[i][7]; 
				if(m_debug) cout<<"gam1_Fpi01="<<gam1_Fpi01<<" gam2_Fpi01="<<gam2_Fpi01<<endl;
				int gam1_Fpi02=m_pi0_par[j][6]; 
				int gam2_Fpi02=m_pi0_par[j][7];
				if(m_debug) cout<<"gam1_Fpi02="<<gam1_Fpi02<<" gam2_Fpi02="<<gam2_Fpi02<<endl;
				int gam1_Fpi03=m_pi0_par[k][6]; 
				int gam2_Fpi03=m_pi0_par[k][7];
				if(m_debug) cout<<"gam1_Fpi03="<<gam1_Fpi03<<" gam2_Fpi03="<<gam2_Fpi03<<endl;

				if((gam1_Fpi01==gam1_Fpi02)||(gam1_Fpi01==gam2_Fpi02)||(gam1_Fpi01==gam1_Fpi03)||(gam1_Fpi01==gam2_Fpi03)) continue;
				if((gam2_Fpi01==gam1_Fpi02)||(gam2_Fpi01==gam2_Fpi02)||(gam2_Fpi01==gam1_Fpi03)||(gam2_Fpi01==gam2_Fpi03)) continue;
				if((gam1_Fpi02==gam1_Fpi01)||(gam1_Fpi02==gam2_Fpi01)||(gam1_Fpi02==gam1_Fpi03)||(gam1_Fpi02==gam2_Fpi03)) continue;
				if((gam2_Fpi02==gam1_Fpi01)||(gam2_Fpi02==gam2_Fpi01)||(gam2_Fpi02==gam1_Fpi03)||(gam2_Fpi02==gam2_Fpi03)) continue;
				if((gam1_Fpi03==gam1_Fpi01)||(gam1_Fpi03==gam2_Fpi01)||(gam1_Fpi03==gam1_Fpi02)||(gam1_Fpi03==gam2_Fpi02)) continue;
				if((gam2_Fpi03==gam1_Fpi01)||(gam2_Fpi03==gam2_Fpi01)||(gam2_Fpi03==gam1_Fpi02)||(gam2_Fpi03==gam2_Fpi02)) continue;

				HepLorentzVector p4_etaTO3pi0=p4_pi0_1+p4_pi0_2+p4_pi0_3;
				double m_etaTO3pi0_tp=p4_etaTO3pi0.m();
				if(m_debug) cout<<"m_etaTO3pi0_tp="<<m_etaTO3pi0_tp<<endl;
				if(m_etaTO3pi0_tp>0.48&&m_etaTO3pi0_tp<0.59) {
					mm++;
					if(mm>48) return StatusCode::SUCCESS;
					m_M_etaTO3pi0[mm]=m_etaTO3pi0_tp;
					m_pi01Feta_id[mm]=i;
					m_pi02Feta_id[mm]=j;
					m_pi03Feta_id[mm]=k;
					if(m_debug) cout<<"mm="<<mm<<endl;
				}
				else continue;

				//kinematic fit for eta->3pi0
				EvtRecTrackIterator itTrk_shower1_Fpi01=evtRecTrkCol->begin() + gam1_Fpi01;
				RecEmcShower* shower1_Fpi01 = (*itTrk_shower1_Fpi01)->emcShower();
				EvtRecTrackIterator itTrk_shower2_Fpi01=evtRecTrkCol->begin() + gam2_Fpi01;
				RecEmcShower* shower2_Fpi01 = (*itTrk_shower2_Fpi01)->emcShower();

				EvtRecTrackIterator itTrk_shower1_Fpi02=evtRecTrkCol->begin() + gam1_Fpi02;
				RecEmcShower* shower1_Fpi02 = (*itTrk_shower1_Fpi02)->emcShower();
				EvtRecTrackIterator itTrk_shower2_Fpi02=evtRecTrkCol->begin() + gam2_Fpi02;
				RecEmcShower* shower2_Fpi02 = (*itTrk_shower2_Fpi02)->emcShower();

				EvtRecTrackIterator itTrk_shower1_Fpi03=evtRecTrkCol->begin() + gam1_Fpi03;
				RecEmcShower* shower1_Fpi03 = (*itTrk_shower1_Fpi03)->emcShower();
				EvtRecTrackIterator itTrk_shower2_Fpi03=evtRecTrkCol->begin() + gam2_Fpi03;
				RecEmcShower* shower2_Fpi03 = (*itTrk_shower2_Fpi03)->emcShower();

				KalmanKinematicFit * kmfit_etaTO3pi0 = KalmanKinematicFit::instance();
				kmfit_etaTO3pi0->init();
				kmfit_etaTO3pi0->setChisqCut(2500);
				kmfit_etaTO3pi0->AddTrack(0, 0.0, shower1_Fpi01);
				kmfit_etaTO3pi0->AddTrack(1, 0.0, shower2_Fpi01);
				kmfit_etaTO3pi0->AddTrack(2, 0.0, shower1_Fpi02);
				kmfit_etaTO3pi0->AddTrack(3, 0.0, shower2_Fpi02);
				kmfit_etaTO3pi0->AddTrack(4, 0.0, shower1_Fpi03);
				kmfit_etaTO3pi0->AddTrack(5, 0.0, shower2_Fpi03);
				kmfit_etaTO3pi0->AddResonance(0, mpi0, 0,1);
				kmfit_etaTO3pi0->AddResonance(1, mpi0, 2,3);
				kmfit_etaTO3pi0->AddResonance(2, mpi0, 4,5);
				kmfit_etaTO3pi0->AddResonance(3, meta, 0,1,2,3,4,5);
				kmfit_etaTO3pi0->Fit();
				double chis_eta = kmfit_etaTO3pi0->chisq();
				if(m_debug) cout<<"chis_eta="<<chis_eta<<endl;
				m_chis_etaTO3pi0[mm]=chis_eta;

				HepLorentzVector p4_gam1Fpi01Feta_kf=kmfit_etaTO3pi0->pfit(0);
				HepLorentzVector p4_gam2Fpi01Feta_kf=kmfit_etaTO3pi0->pfit(1);
				HepLorentzVector p4_gam1Fpi02Feta_kf=kmfit_etaTO3pi0->pfit(2);
				HepLorentzVector p4_gam2Fpi02Feta_kf=kmfit_etaTO3pi0->pfit(3);
				HepLorentzVector p4_gam1Fpi03Feta_kf=kmfit_etaTO3pi0->pfit(4);
				HepLorentzVector p4_gam2Fpi03Feta_kf=kmfit_etaTO3pi0->pfit(5);

				for(int xx=0;xx<4;xx++) {
					m_p4_gam1Fpi01Feta_kf[mm][xx]=p4_gam1Fpi01Feta_kf[xx];
					m_p4_gam2Fpi01Feta_kf[mm][xx]=p4_gam2Fpi01Feta_kf[xx];
					m_p4_gam1Fpi02Feta_kf[mm][xx]=p4_gam1Fpi02Feta_kf[xx];
					m_p4_gam2Fpi02Feta_kf[mm][xx]=p4_gam2Fpi02Feta_kf[xx];
					m_p4_gam1Fpi03Feta_kf[mm][xx]=p4_gam1Fpi03Feta_kf[xx];
					m_p4_gam2Fpi03Feta_kf[mm][xx]=p4_gam2Fpi03Feta_kf[xx];
				}
				HepLorentzVector p4_etaTO3pi0_kf=p4_gam1Fpi01Feta_kf+p4_gam2Fpi01Feta_kf+p4_gam1Fpi02Feta_kf+p4_gam2Fpi02Feta_kf+p4_gam1Fpi03Feta_kf+p4_gam2Fpi03Feta_kf;
				double m_eta_1c=p4_etaTO3pi0_kf.m();
				if(m_debug) cout<<"eta->3pi0: m_eta_1c="<<m_eta_1c<<endl;
			}
		}
	}
	m_no_etaTO3pi0=mm+1;
	bool etaTO3pi0=false;
	if(m_no_etaTO3pi0>0) etaTO3pi0=true;

	//phi->K+K-
	ss=-1;
	for(int i=0;i<no_kap;i++){
		int kap_id=kapID[i];
		if(m_debug) cout<<"kap_id="<<kap_id<<endl;
		RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
		EvtRecTrackIterator itTrk_kap=evtRecTrkCol->begin() + kap_id;
		RecMdcKalTrack* mdcKalTrk_kap = (*itTrk_kap)->mdcKalTrack();
		HepLorentzVector p4_kap_tp=mdcKalTrk_kap->p4(mkaon);
		WTrackParameter wvkapFphi;
		wvkapFphi=WTrackParameter(mkaon, mdcKalTrk_kap->getZHelixK(), mdcKalTrk_kap->getZErrorK());

		for(int j=0;j<no_kam;j++){
			int kam_id=kamID[j];
			if(m_debug) cout<<"kam_id="<<kam_id<<endl;
			RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
			EvtRecTrackIterator itTrk_kam=evtRecTrkCol->begin() + kam_id;
			RecMdcKalTrack* mdcKalTrk_kam = (*itTrk_kam)->mdcKalTrack();

			HepLorentzVector p4_kam_tp=mdcKalTrk_kam->p4(mkaon);
			WTrackParameter wvkamFphi;
			wvkamFphi=WTrackParameter(mkaon, mdcKalTrk_kam->getZHelixK(), mdcKalTrk_kam->getZErrorK());

			HepLorentzVector p4_phiTOkk=p4_kap_tp+p4_kam_tp;
			double m_phiTOkk_tp=p4_phiTOkk.m();
			if(m_debug) cout<<"m_phiTOkk_tp="<<m_phiTOkk_tp<<endl;

			//Rxy of K+K-
			RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
			EvtRecTrackIterator itTrk_ep=evtRecTrkCol->begin() + kap_id;
			RecMdcKalTrack* mdcKalTrk_ep = (*itTrk_ep)->mdcKalTrack();
			RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
			EvtRecTrackIterator itTrk_em=evtRecTrkCol->begin() + kam_id;
			RecMdcKalTrack* mdcKalTrk_em = (*itTrk_em)->mdcKalTrack();

			Hep3Vector xorigin(0,0,0);
			IVertexDbSvc*  vtxsvc;
			Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
			if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
			}
			HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
			GammaConv gconv = GammaConv(mdcKalTrk_ep->helix(), mdcKalTrk_em->helix(),IP);
			double rconv = gconv.getRXY1();
			if(m_debug) cout<<"rconv="<<rconv<<endl;

			HepLorentzVector p4_ep=mdcKalTrk_ep->p4(melectron);
			HepLorentzVector p4_em=mdcKalTrk_em->p4(melectron);

			Hep3Vector p3_ep=p4_ep.v();
			Hep3Vector p3_em=p4_em.v();
			double ang_ee=p3_ep.angle(p3_em);
			double cos_ee=cos(ang_ee);

			//end Rxy
			if(m_phiTOkk_tp>0.98&&m_phiTOkk_tp<1.06) {
				ss++;
				if(ss>98) return StatusCode::SUCCESS;
				m_rconv[ss]=rconv;
				m_M_phiTOkk[ss]=m_phiTOkk_tp;  
				m_kapFphi_id[ss]=i;
				m_kamFphi_id[ss]=j;		
				m_cos_ee[ss]=cos_ee;
				if(m_debug) cout<<"ss="<<ss<<endl;
			}
			else continue;
			//kinematic fit for phi->K+K-
			KalmanKinematicFit * kmfit_phi = KalmanKinematicFit::instance();
			kmfit_phi->init();
			kmfit_phi->setChisqCut(2500);
			kmfit_phi->AddTrack(0, wvkapFphi);
			kmfit_phi->AddTrack(1, wvkamFphi);
			kmfit_phi->AddResonance(0, mphi, 0,1);
			kmfit_phi->Fit();
			double chis_phi = kmfit_phi->chisq();
			if(m_debug) cout<<"chis_phi="<<chis_phi<<endl;
			m_chis_phiTOkk[ss]=chis_phi;

			HepLorentzVector p4_kapFphi_kf=kmfit_phi->pfit(0);
			HepLorentzVector p4_kamFphi_kf=kmfit_phi->pfit(1);

			if(m_debug) cout<<"p4_kapFphi_kf.px()="<<p4_kapFphi_kf.px()<<endl;

			for(int xx=0;xx<4;xx++) {
				m_p4_kapFphi_kf[ss][xx]=p4_kapFphi_kf[xx];
				m_p4_kamFphi_kf[ss][xx]=p4_kamFphi_kf[xx];
			}
			HepLorentzVector p4_phiTOkk_kf=p4_kapFphi_kf+p4_kamFphi_kf;
			double m_phi_1c=p4_phiTOkk_kf.m();
			if(m_debug) cout<<"m_phi_1c="<<m_phi_1c<<endl;

		}
	}
	m_no_phiTOkk=ss+1;
	bool phiTOkk=false;
	if(m_no_phiTOkk>0) phiTOkk=true;

	if(m_debug) {
		cout<<"etaTOgg:"<<etaTOgg<<endl;
		cout<<"etaTOpipipi0:"<<etaTOpipipi0<<endl;
		cout<<"etaTO3pi0:"<<etaTO3pi0<<endl;
		cout<<"phiTOkk:"<<phiTOkk<<endl;
	} 

	if(phiTOkk||etaTOgg||etaTOpipipi0||etaTO3pi0)  m_tuple1->write();

	return StatusCode::SUCCESS;
}//end of execute()
//**************************************************************
StatusCode JpsiToPhiEtaAlg::endRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToPhiEtaAlg::endRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;

}
StatusCode JpsiToPhiEtaAlg::finalize(){
	if(m_debug)  cout<<"finalize: no_evt1="<<no_evt1<<endl;
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"JpsiToPhiEtaAlg::finalize()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;
}

//add your code here,for chrer member-functions
bool JpsiToPhiEtaAlg::goodTrk(EvtRecTrackIterator itTrk) {
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

bool JpsiToPhiEtaAlg::good_gam(RecEmcShower *emcTrk,double& dang_min){

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
HepLorentzVector JpsiToPhiEtaAlg::getP4(RecEmcShower* gTrk){

	double eraw = gTrk->energy();
	double phi =  gTrk->phi();
	double the =  gTrk->theta();

	return HepLorentzVector( eraw * sin(the) * cos(phi),
			eraw * sin(the) * sin(phi),
			eraw * cos(the),
			eraw );
}
void JpsiToPhiEtaAlg::savepi0(RecEmcShower *shr1,RecEmcShower *shr2,double& pi0_chis,HepLorentzVector& p4_pi0,HepLorentzVector& p4_pi0_1c){
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
void JpsiToPhiEtaAlg::saveeta(RecEmcShower *shr1,RecEmcShower *shr2,double& eta_chis,HepLorentzVector& p4_eta,HepLorentzVector& p4_eta_1c){
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
void JpsiToPhiEtaAlg::PID(EvtRecTrackIterator itTrk,double prob[25]){
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
