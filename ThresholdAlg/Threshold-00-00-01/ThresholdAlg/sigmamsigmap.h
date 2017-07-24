#ifndef Physics_Analysis_sigmamsigmap
#define Physics_Analysis_sigmamsigmap

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

#include "TH1.h"
#include "TH2.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "VertexFit/WTrackParameter.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"


class sigmamsigmap: public Algorithm {

		public:
				sigmamsigmap(const std::string& name, ISvcLocator* pSvcLocator);
				StatusCode initialize();
				StatusCode execute();
				StatusCode finalize();

		private:
				double pnbar(const double mnbar, Hep3Vector &nbarpos, HepLorentzVector &pi0);
				double m_energyThreshold;
				double m_energyThreshold1;
				double m_energyThreshold2;
				double m_gammaAngleCut;
				double cms;

				//


				// define NTuple here

				NTuple::Tuple* m_tuple0;
				NTuple::Item<double> m_mcpnbar;
				NTuple::Item<double> m_mcpn;
				NTuple::Item<double> m_mcppip;
				NTuple::Item<double> m_mcppim;
				NTuple::Item<double> m_mcpsigma;
				NTuple::Item<double> m_mcpsigma_b;
				NTuple::Item<double> m_mcpxsigma;
				NTuple::Item<double> m_mcpysigma;
				NTuple::Item<double> m_mcpzsigma;
				NTuple::Item<double> m_mcesigma;

				NTuple::Item<double> m_mcpp;
				NTuple::Item<double> m_mcppbar;
				NTuple::Item<double> m_mcpxpbar;
				NTuple::Item<double> m_mcpypbar;
				NTuple::Item<double> m_mcpzpbar;
				NTuple::Item<double> m_mcepbar;
				NTuple::Item<double> m_mcpposx;
				NTuple::Item<double> m_mcpposy;
				NTuple::Item<double> m_mcpposz;

				NTuple::Item<double> m_mcppi0;
				NTuple::Item<double> m_mcpgamma;

				NTuple::Item<double> m_mccossigma;

				NTuple::Item<long>   m_nmc_p;
				NTuple::Item<long>   m_nmc_pi0;
				NTuple::Item<long>   m_nmc_pip;
				NTuple::Item<long>   m_nmc_pim;
				NTuple::Item<long>   m_nmc_nbar;
				NTuple::Item<long>   m_nmc_n;
				NTuple::Item<long>   m_nmc_sigma;
				NTuple::Item<long>   m_nmc_gamma;
				NTuple::Item<long>   m_nmc_p_loc;
				NTuple::Item<long>   m_nmc_pi0_loc;
				NTuple::Item<long>   m_nmc_gamma_loc;
				NTuple::Item<long>	 m_runNo;
				NTuple::Item<long>	 m_event;

				NTuple::Tuple* m_tuple1;
				NTuple::Item<int> 	 m_idxmc;
				NTuple::Array<int> 	 m_pdgid;
				NTuple::Array<int> 	 m_motheridx;
				NTuple::Item<double> m_mcpgamma_loc;
				NTuple::Item<double> m_mcenbar_loc;
				NTuple::Item<double> m_mcpxnbar_loc;
				NTuple::Item<double> m_mcpynbar_loc;
				NTuple::Item<double> m_mcpznbar_loc;
				NTuple::Item<double> m_mcen_loc;
				NTuple::Item<double> m_mcepip_loc;
				NTuple::Item<double> m_mcepim_loc;
				NTuple::Item<double> m_mcpxsigma_loc;
				NTuple::Item<double> m_mcpysigma_loc;
				NTuple::Item<double> m_mcpzsigma_loc;
				NTuple::Item<double> m_mcesigma_loc;
				NTuple::Item<double> m_mcpxpbar_loc;
				NTuple::Item<double> m_mcpypbar_loc;
				NTuple::Item<double> m_mcpzpbar_loc;
				NTuple::Item<double> m_mcepbar_loc;

				NTuple::Item<double> m_mcpposx_loc;
				NTuple::Item<double> m_mcpposy_loc;
				NTuple::Item<double> m_mcpposz_loc;
				NTuple::Item<double> m_mcpimposx_loc;
				NTuple::Item<double> m_mcpimposy_loc;
				NTuple::Item<double> m_mcpimposz_loc;

				NTuple::Item<double> m_apRvz0;
				NTuple::Item<double> m_apRvxy0;
				NTuple::Item<double> m_pRvz0;
				NTuple::Item<double> m_pRvxy0;
				NTuple::Item<double> m_pimRvz0;
				NTuple::Item<double> m_pimRvxy0;
				NTuple::Item<double> m_pipRvz0;
				NTuple::Item<double> m_pipRvxy0;
				NTuple::Item<double> m_pid_ppip;
				NTuple::Item<double> m_pid_pip_px;
				NTuple::Item<double> m_pid_pip_py;
				NTuple::Item<double> m_pid_pip_pz;
				NTuple::Item<double> m_pip_ep;
				NTuple::Item<double> m_pid_ppim;
				NTuple::Item<double> m_pid_pim_px;
				NTuple::Item<double> m_pid_pim_py;
				NTuple::Item<double> m_pid_pim_pz;
				NTuple::Item<double> m_pim_ep;
				NTuple::Item<double> m_pid_pm;

				NTuple::Item<double> m_nbar_energy;
				NTuple::Item<double> m_nbar_hits;
				NTuple::Item<double> m_nbar_cos;
				NTuple::Item<double> m_nbar_eseed;
				NTuple::Item<double> m_nbar_shape;
				NTuple::Item<double> m_nbar_latmom;
				NTuple::Item<double> m_nbar_secmom;
				NTuple::Item<double> m_nbar_hit;
				NTuple::Item<double> m_nbar_match;
				NTuple::Item<double> m_nbar_costheta;
				NTuple::Item<double> m_nbar_phi;
				NTuple::Item<double> m_nbar_time;

				NTuple::Item<double> m_chisq;
				NTuple::Item<double> m_Sigmac_m;
				NTuple::Item<double> m_Sigmac_p;
				NTuple::Item<double> m_Sigmac_cos;
				NTuple::Item<double> m_aSigmac_m;
				NTuple::Item<double> m_aSigmac_p;
				NTuple::Item<double> m_aSigamc_cos;

				
				NTuple::Item<double> m_pxpim;
				NTuple::Item<double> m_pypim;
				NTuple::Item<double> m_pzpim;
				NTuple::Item<double> m_epim;
				NTuple::Item<double> m_pxpip;
				NTuple::Item<double> m_pypip;
				NTuple::Item<double> m_pzpip;
				NTuple::Item<double> m_epip;
				NTuple::Item<double> m_pxnbar;
				NTuple::Item<double> m_pynbar;
				NTuple::Item<double> m_pznbar;
				NTuple::Item<double> m_enbar;
				NTuple::Item<double> m_pxn;
				NTuple::Item<double> m_pyn;
				NTuple::Item<double> m_pzn;
				NTuple::Item<double> m_en;

				NTuple::Item<double> m_chisqv2;
				NTuple::Item<double> m_Sigmacv2_m;
				NTuple::Item<double> m_Sigmacv2_p;
				NTuple::Item<double> m_Sigmacv2_cos;
				NTuple::Item<double> m_aSigmacv2_m;
				NTuple::Item<double> m_aSigmacv2_p;
				NTuple::Item<double> m_aSigmacv2_cos;

				
				NTuple::Item<double> m_pxpimv2;
				NTuple::Item<double> m_pypimv2;
				NTuple::Item<double> m_pzpimv2;
				NTuple::Item<double> m_epimv2;
				NTuple::Item<double> m_pxpipv2;
				NTuple::Item<double> m_pypipv2;
				NTuple::Item<double> m_pzpipv2;
				NTuple::Item<double> m_epipv2;
				NTuple::Item<double> m_pxnbarv2;
				NTuple::Item<double> m_pynbarv2;
				NTuple::Item<double> m_pznbarv2;
				NTuple::Item<double> m_enbarv2;
				NTuple::Item<double> m_pxnv2;
				NTuple::Item<double> m_pynv2;
				NTuple::Item<double> m_pznv2;
				NTuple::Item<double> m_env2;

				NTuple::Item<double> m_mpipi;
				NTuple::Item<double> m_mpipirecoil;

				TH1D* hmdc_nGood;
				TH1D* hmdc_nCharge;
				TH1D* hnpip;
				TH1D* hnpim;
				TH1D* hecm_egam;
				TH1D* hchisqpi0;
				TH1D* hppi0;
};

#endif
