#ifndef JpsiToPhiEtaAlg_Header
#define JpsiToPhiEtaAlg_Header

#include "GaudiKernel/Algorithm.h"
//you can add oher necessary header files
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/AlgFactory.h"
#include "VertexFit/IVertexDbSvc.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

class JpsiToPhiEtaAlg:public Algorithm {
  public:
    JpsiToPhiEtaAlg(const std::string& name, ISvcLocator* pSvcLocator);
    ~JpsiToPhiEtaAlg();
    StatusCode initialize();
    StatusCode beginRun();   
    StatusCode execute();
    StatusCode endRun();
    StatusCode finalize();
    bool goodTrk(EvtRecTrackIterator itTrk);
    bool good_gam(RecEmcShower *emcTrk, double& dang_min);
    HepLorentzVector getP4(RecEmcShower* gTrk);
    void savepi0(RecEmcShower *shr1,RecEmcShower *shr2,double& pi0_chis,HepLorentzVector& p4_pi0,HepLorentzVector& p4_pi0_1c);
    void saveeta(RecEmcShower *shr1,RecEmcShower *shr2,double& eta_chis,HepLorentzVector& p4_eta,HepLorentzVector& p4_eta_1c);
    void PID(EvtRecTrackIterator itTrk,double prob[25]);

  private:
    bool	m_debug;
    double      m_Ecms;

    NTuple::Tuple* m_tuple1;                    
    NTuple::Item<int>  m_runNo;
    NTuple::Item<int>  m_evtNo;

    NTuple::Item<int>  m_no_pip;
    NTuple::Matrix<double> m_p4_pip;
    NTuple::Array<int> m_pip_id;
    NTuple::Item<int>  m_no_pim;
    NTuple::Matrix<double> m_p4_pim;
    NTuple::Array<int> m_pim_id;
    NTuple::Item<int>  m_no_kap;
    NTuple::Matrix<double> m_p4_kap;
    NTuple::Array<int> m_kap_id;
    NTuple::Item<int>  m_no_kam;
    NTuple::Matrix<double> m_p4_kam;
    NTuple::Array<int> m_kam_id;
    NTuple::Item<int>  m_no_gam;
    NTuple::Matrix<double> m_gam_par;
    NTuple::Item<int>  m_no_pi0;
    NTuple::Matrix<double> m_pi0_par;
    NTuple::Item<int>  m_no_eta;
    NTuple::Matrix<double> m_eta_par;
    NTuple::Item<int>  m_no_chrp;
    NTuple::Matrix<double> m_chrp_p3;
    NTuple::Array<int> m_chrp_id;
    NTuple::Array<double> m_chrp_Eemc;
    NTuple::Matrix<double> m_chrp_prob;
    NTuple::Item<int>  m_no_chrm;
    NTuple::Matrix<double> m_chrm_p3;
    NTuple::Array<int> m_chrm_id;
    NTuple::Array<double> m_chrm_Eemc;
    NTuple::Matrix<double> m_chrm_prob;
    NTuple::Item<int>  m_no_etaTOpipipi0;
    NTuple::Array<double> m_chis_etaTOpipipi0;
    NTuple::Matrix<double> m_p4_pipFeta_kf;
    NTuple::Matrix<double> m_p4_pimFeta_kf;
    NTuple::Matrix<double> m_p4_gam1Feta_kf;
    NTuple::Matrix<double> m_p4_gam2Feta_kf;
    NTuple::Array<double> m_M_etaTOpipipi0;
    NTuple::Array<int> m_pipFeta_id;
    NTuple::Array<int> m_pimFeta_id;
    NTuple::Array<int> m_pi0Feta_id;
    NTuple::Item<int>  m_no_etaTO3pi0;
    NTuple::Array<double> m_chis_etaTO3pi0;
    NTuple::Matrix<double> m_p4_gam1Fpi01Feta_kf;
    NTuple::Matrix<double> m_p4_gam2Fpi01Feta_kf;
    NTuple::Matrix<double> m_p4_gam1Fpi02Feta_kf;
    NTuple::Matrix<double> m_p4_gam2Fpi02Feta_kf;
    NTuple::Matrix<double> m_p4_gam1Fpi03Feta_kf;
    NTuple::Matrix<double> m_p4_gam2Fpi03Feta_kf;
    NTuple::Array<double> m_M_etaTO3pi0;
    NTuple::Array<int> m_pi01Feta_id;
    NTuple::Array<int> m_pi02Feta_id;
    NTuple::Array<int> m_pi03Feta_id;
    NTuple::Item<int>  m_no_phiTOkk;
    NTuple::Array<double> m_chis_phiTOkk;
    NTuple::Matrix<double> m_p4_kapFphi_kf;
    NTuple::Matrix<double> m_p4_kamFphi_kf;
    NTuple::Array<double> m_M_phiTOkk;
    NTuple::Array<double> m_rconv;
    NTuple::Array<int> m_kapFphi_id;
    NTuple::Array<int> m_kamFphi_id;
    NTuple::Array<double> m_cos_ee;

    //MC info. 
    NTuple::Item<int>  m_idxmc;
    NTuple::Array<int>  m_pdgid;
    NTuple::Array<int>  m_motheridx;
};
#endif//JpsiToPhiEtaAlg_Header
