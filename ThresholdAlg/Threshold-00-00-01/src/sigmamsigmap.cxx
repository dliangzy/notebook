//author Liu Liang
//version 0.1
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/IHistogramSvc.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "DstEvent/TofHitStatus.h"

#include "McTruth/McParticle.h"
#include "Identifier/TofID.h"
#include "ParticleID/ParticleID.h"
#include "TrigEvent/TrigEvent.h"
#include "TrigEvent/TrigData.h"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"

#include "ThresholdAlg/sigmamsigmap.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#ifndef ENABLE_BACKWARDS_COMPATIBITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>

const double mpi = 0.13957;
const double mkaon = 0.493677;
const double mproton = 0.938272;
const double mpion0 = 0.1395766;
const double mpion = 0.139570;
const double mEta = 0.547853;
const double n_mass = 0.939565;
const double nbar_mass = 0.939565;
const double Sigmam_mass = 1.197436;
const double Sigma_mass = 1.18937;
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
const double velc = 299.792458;			// tof path unit in mm

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

static int Ncut_none=0, 	Ncut_ncharge=0, 
		   Ncut_pid=0, 		Ncut_nshower=0, 
		   Ncut_npi0=0;
static int Ncut_pppi0pi0=0,	Ncut_ppi0=0,
		   Ncut_pbarpi0=0,	Ncut_nnbarpipi,
		   Ncut_mc;
static int Ncut_kin=0,		Ncut_nbar=0,
		   Ncut_asigma=0, 	Ncut_kin_flag1=0,
		   Ncut_kin_flag2=0;
static int Ncut_pid_flag1=0,	Ncut_pid_flag2=0,
		   Ncut_pid_flag3=0;

//------------------------------------------------------------------------------------------------
sigmamsigmap::sigmamsigmap(const std::string& name, ISvcLocator* pSvcLocator):
	Algorithm(name, pSvcLocator){
		declareProperty("energyTreshold",m_energyThreshold=0.025);
		declareProperty("energyTreshold1",m_energyThreshold2=0.025);
		declareProperty("energyTreshold2",m_energyThreshold2=0.050);
		declareProperty("gammaAngleCut",m_gammaAngleCut=20.0);
		declareProperty("cms",cms=3.65);
	}

//-------------------------------------------------------------------------------

StatusCode sigmamsigmap::initialize(){
	MsgStream log(msgSvc(), name());

	log << MSG::INFO << "in initialize()" << endmsg;

	StatusCode status;

	NTuplePtr nt0(ntupleSvc(), "FILE1/MCtruth");
	if ( nt0 ) m_tuple0 = nt0;
	else	{
		m_tuple0 = ntupleSvc()->book ("FILE1/MCtruth", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if( m_tuple0 )	{
			status  = m_tuple0->addItem ("nmc_pip", 	m_nmc_pip);
			status  = m_tuple0->addItem ("nmc_nbar", 	m_nmc_nbar);
			status  = m_tuple0->addItem ("nmc_pim",		m_nmc_pim);
			status  = m_tuple0->addItem ("nmc_n",		m_nmc_n);
			status  = m_tuple0->addItem ("nmc_sigma",	m_nmc_sigma);
			status  = m_tuple0->addItem ("nmc_pi0",		m_nmc_pi0);
			status  = m_tuple0->addItem ("nmc_p",		m_nmc_p);

			status  = m_tuple0->addItem ("mcpnbar", 	m_mcpnbar);
			status  = m_tuple0->addItem ("mcppip",		m_mcppip);
			status  = m_tuple0->addItem ("mcpn",		m_mcpn);
			status  = m_tuple0->addItem ("mcppim",		m_mcppim);
			status  = m_tuple0->addItem ("mcppi0",		m_mcppi0);
			status  = m_tuple0->addItem ("mcpp",		m_mcpp);
			status  = m_tuple0->addItem ("mcpsigma",	m_mcpsigma);
			status  = m_tuple0->addItem ("mccossogma",	m_mccossigma);
			status  = m_tuple0->addItem ("nmc_gamma",	m_nmc_gamma);
			status  = m_tuple0->addItem ("mcgamma", 	m_mcpgamma);
		}
		else	{
			log << MSG::ERROR << "Connot book N-tuple:" << long(m_tuple0) << endmsg;
			return StatusCode::FAILURE;
		}
	}

/*
	NTuplePtr nt1(ntupleSvc(), "FILE1/shower");
	if ( nt1 ) m_tuple1 = nt1;
	else	{
		m_tuple1 = ntupleSvc()->book ("FILE1/shower", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if( m_tuple1 )	{
			status  = m_tuple1->addItem ();
			status  = m_tuple1->addItem ();
			status  = m_tuple1->addItem ();
			status  = m_tuple1->addItem ();
			status  = m_tuple1->addItem ();
			status  = m_tuple1->addItem ();
			status  = m_tuple1->addItem ();
			status  = m_tuple1->addItem ();
		}
		else	{
			log << MSG::ERROR << "Connot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}
*/

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;
}



//========================================================================================




StatusCode sigmamsigmap::execute()	{
	MsgStream log(msgSvc(),name());
	log << MSG::INFO << "in execute()" << endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	int runNo=eventHeader->runNumber();
	int event=eventHeader->eventNumber();
	m_runNo = eventHeader->runNumber();
	m_event = eventHeader->eventNumber();
	log << MSG::DEBUG << "run, evtnum = "
		<< runNo << " , "
		<< event << endreq;
	Ncut_none++;
	
	int RadTrig = 0;
	if ( runNo>0 )	{
		SmartDataPtr<TrigData> trigData(eventSvc(),EventModel::Trig::TrigData);
		if( !trigData )	{
			cout << "Could not find Trigger Data for physics analysis" << endl;
			return StatusCode::FAILURE;
		}
		RadTrig = trigData->getTrigChannel(9);
	}

//-----------------------------------------------------------------------------
	if ( runNo<0 )	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
		if ( mcParticleCol )	{
			int m_numParticle = 0;
			bool psipDecay = false;
			int rootIndex = -1;
			int nmc_eta = 0, nmc_pi0 = 0, nmc_gampi0 = 0, nmc_gameta = 0, nmc_antip = 0, nmc_a0 = 0;
			Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
			for ( ; iter_mc != mcParticleCol->end(); iter_mc++)	{
				//don't understand
				//psipDecay will always be true.
				if ((*iter_mc)->primaryParticle()) continue;
				if (!(*iter_mc)->decayFromGenerator()) continue;
				if ((*iter_mc)->particleProperty()==100443)	{
					psipDecay = true;
					rootIndex = (*iter_mc)->trackIndex();
				}
				if ( !psipDecay ) continue;
				int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
				int pdgid = (*iter_mc)->particleProperty();
				m_pdgid[m_numParticle] = pdgid;
				m_motheridx[m_numParticle] = mcidx;
				m_numParticle += 1;
			}
			m_idxmc = m_numParticle;
		}
	}

//----------------------------------------------------------
	HepLorentzVector mc_p4pip, mc_p4nbar, mc_p4n, mc_p4sigma, mc_p4pim, mc_p4sigmabar;
	if ( runNo < 0 )	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		if ( mcParticleCol )	{
			bool vphoDecay = false;
			int rootIndex = -1;
			int nmc_pip   = 0,		nmc_pim = 0,
				nmc_nbar  = 0,		nmc_n   = 0,
				nmc_sigma = 0,		nmc_pi0 = 0,
				nmc_gamma,			nmc_p	= 0;
			Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
			for ( ; iter_mc != mcParticleCol->end(); iter_mc++)	{
				if (!(*iter_mc)->decayFromGenerator()) continue;
				vphoDecay = true;
				HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
				if((*iter_mc)->particleProperty()==-3112)	{
					mc_p4sigmabar = mctrue_track;
				}
				if ((*iter_mc)->particleProperty()==-2112 && ((*iter_mc)->mother()).particleProperty() == -3112)	{
					mc_p4nbar      = mctrue_track;
					m_mcpnbar      = mctrue_track.rho();
					m_mcenbar_loc  = mctrue_track.e();
					m_mcpxnbar_loc = mctrue_track.px();			//px
					m_mcpynbar_loc = mctrue_track.px();			//px?
					m_mcpznbar_loc = mctrue_track.px();			//px?
					nmc_nbar++;
				}

				if ((*iter_mc)->particleProperty()==2112 && ((*iter_mc)->mother()).particleProperty() == 3112)	{
					mc_p4n		= mctrue_track;
					m_mcpn		= mctrue_track.rho();
					m_mcen_loc	= mctrue_track.e();
					nmc_n++;
				}
				
				if ((*iter_mc)->particleProperty()==211 && ((*iter_mc)->mother()).particleProperty() == -3112)	{
					mc_p4pip		= mctrue_track;
					m_mcppip		= mctrue_track.rho();
					m_mcepip_loc	= mctrue_track.e();
					nmc_pip++;
				}

				if ((*iter_mc)->particleProperty()==-211 && ((*iter_mc)->mother()).particleProperty() == 3112)	{
					mc_p4pim = mctrue_track;
					m_mcppim = mctrue_track.rho();
					m_mcepim_loc = mctrue_track.e();
					m_mcpimposx_loc = (*iter_mc)->initialPosition().x();
					m_mcpimposy_loc = (*iter_mc)->initialPosition().y();
					m_mcpimposz_loc = (*iter_mc)->initialPosition().z();
					nmc_pim++;
				}
				if ((*iter_mc)->particleProperty()==-3112)	{
					mc_p4sigma		= mctrue_track;
					m_mcpsigma		= mctrue_track.rho();
					m_mccossigma	= cos(mctrue_track.theta());
					nmc_sigma++;
				}
				if((*iter_mc)->particleProperty()==111 && ( ((*iter_mc)->mother()).particleProperty()== -3112 || ((*iter_mc)->mother()).particleProperty() == 3112 ))	{
					nmc_pi0++;
					m_mcppi0 = mctrue_track.rho();
				}
				
				if ((*iter_mc)->particleProperty()==2112 && ((*iter_mc)->mother()).particleProperty() == 3112)	{
					nmc_p++;
					m_mcpp = mctrue_track.rho();
				}//same with nmc_n
				if ((*iter_mc)->particleProperty()==22)	{
					m_mcpgamma = mctrue_track.rho();
					m_mcpgamma_loc = mctrue_track.rho();
					nmc_gamma++;
				}

			}
			m_nmc_nbar		= nmc_nbar;
			m_nmc_n			= nmc_n;
			m_nmc_pip		= nmc_pip;
			m_nmc_pim		= nmc_pim;
			m_nmc_sigma		= nmc_sigma;
			m_nmc_pi0		= nmc_pi0;
			m_nmc_p			= nmc_p;
			m_nmc_pi0_loc	= nmc_pi0;
			m_nmc_p_loc		= nmc_p;
			m_nmc_gamma		= nmc_gamma;
			m_nmc_gamma_loc = nmc_gamma;
			m_tuple0->write();

			if(nmc_pi0 == 2) { Ncut_pppi0pi0++; }
			if(nmc_pi0 == 1 && nmc_p == 1 ) { Ncut_ppi0++; }
			if(nmc_pi0 == 1 && nmc_p == 0 ) { Ncut_pbarpi0++; }
			if(nmc_pi0 == 0 ) {Ncut_nnbarpipi++; }
			Ncut_mc++;
		}
	}

/*
	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	Hep3Vector xorigin(0,0,0);
	IVertexDbSvc* vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid())	{
		double* dbv = vtxsvc->PrimaryVertex();
		double* vv  = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}

	Vint iGood;
	iGood.clear();
	int nCharge = 0;
	for ( int i = 0; i < evtRecEvent->totalCharged(); i++ )	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *mdcTrk =(*itTrk)->mdcTrack();
		double xv	= xorigin.x();
		double yv	= xorigin.y();
		double phi0	= mdcTrk->helix(1);
		double x0	= mdcTrk->x();
		double y0	= mdcTrk->y();
		double z0	= mdcTrk->z();
		double Rxy	= (x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosnt ruction
		HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
		VFHelix helixip(piont0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double Rvxy0 = fabs(vecipa[0]);
		double Rvz0 = vecipa[0];
		double Rvphi0 = vecipa[0];
		double cost = cos(mdcTrk->theta());
		if(fabs(cost) >= 0.93) continue;
		if(fabs(Rvz0) >= 10.0) continue;
		if(fabs(Rvxy0)>= 2.0 ) continue;
		iGood.push_back(i);
		nCharge += mdcTrk->charge();
	}

	//
	//========================================================

	//
	int nGood = iGood.size();
	m_nGood = nGood;
	m_nCharge = nCharge;
	hmdc_nGood->Fill(nGood);
	hmdc_nCharge->Fill(nCharge);

	if(nGood != 2) return SUCCESS;
	Ncut_ncharge++;
	//======================pid========================

	Vint ipip, ipim, ipp, ipm, iKp, iKm;
	ipip.clear();
	ipim.clear();
	ipp.clear();
	ipm.clear();
	iKp.clear();
	iKm.clear();

	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nGood; i++ )	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2 | pid->useTofE());		//pid->useTofE();
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isMdcDedxValid()) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		RecMdcTrack* dedxTrk = (*itTrk)->mdcDedx();
		SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
		SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
		double prob_pi = pid->probPion();
		double prob_K  = pid->probKaon();
		double prob_p  = pid->probProton();
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		if(mdcKalTrk->charge()>0)	{
			if( prob_pi > prob_K && prob_pi > prob_p )	{
				ipip.push_back(iGood[i]);
				RecMdcKalTrack::setPidType (RecMdcKalTrack:pion);
				HepLorentzVector ptrk_pip;
				ptrk_pip.setPx(mdcKalTrk->px());
				ptrk_pip.setPy(mdcKalTrk->py());
				ptrk_pip.setPz(mdcKalTrk->pz());
				double p3 = ptrk_pip.mag();
				ptrk_pip.setE(sqrt(p3*p3+mpion*mpion));

				m_pid_ppip = ptrk_pip.rho();
				m_pid_pip_px = ptrk_pip.px();
				m_pid_pip_py = ptrk_pip.py();
				m_pid_pip_pz = ptrk_pip.pz();
				double pip_ep=0;
				if((*itTrk)->isEmcShowerValid())	{
					RecEmcShower *emcTrk = (*itTrk)->emcShower();
					pip_ep = emcTrk->energy()/ptrk_pip.rho();
				}
				m_pip_ep = pip_ep;

				HepVector a = mdcTrk->helix();
				HepSymMatrix Ea = mdcTrk->err();
				HepPoint3D point0(0.,0.,0.);	//The initial point for MDC recosntruction
				HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
				VFHelix helixip(point0, a, Ea);
				helixip.pivot(IP);
				HepVector vecipa = helixip.a();
				m_pipRvxy0 = fabs(vecipa[0]);	//the nearest distance to IP in xy plane
				m_pipRvz0 = vecipa[3];		//the nearest distance to IP in z direction
			}
			if( prob_p > prob_K && prob_p > prob_pi )	{
				ipp.push_back(iGood[i]);
				RecMdcKalTrack::setPidType (RecMdcKalTrack:proton);
				HepLorentzVector ptrk_pp;
				ptrk_pp.setPx(mdcKalTrk->px());
				ptrk_pp.setPy(mdcKalTrk->py());
				ptrk_pp.setPz(mdcKalTrk->pz());
				double p3 = ptrk_pp.mag();
				ptrk_pp.setE(sqrt(p3*p3+mproton*mproton));
				m_pip_pp = ptrk_pp.rho();
				HepVector a = mdcTrk->helix();
				HepSymMatrix Ea = mdcTrk->err();
				HepPoint3D point0(0.,0.,0.);	//The initial point for MDC recosntruction
				HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
				VFHelix helixip(point0, a, Ea);
				helixip.pivot(IP);
				HepVector vecipa = helixip.a();
				m_pRvxy0 = fabs(vecipa[0]);	//the nearest distance to IP in xy plane
				m_pRvz0 = vecipa[3];		//the nearest distance to IP in z direction
			}

			if( prob_K > prob_p && prob_K > prob_pi )	{
				iKp.push_back(iGood[i]);
			}
		}
		

		if(mdcKalTrk->charge()<0 )	{
			if( prob_pi > prob_K && prob_pi > prob_p )	{
				ipim.push_back(iGood[i]);
				RecMdcKalTrack::setPidType (RecMdcKalTrack:pion);
				HepLorentzVector ptrk_pim;
				ptrk_pim.setPx(mdcKalTrk->px());
				ptrk_pim.setPy(mdcKalTrk->py());
				ptrk_pim.setPz(mdcKalTrk->pz());
				double p3 = ptrk_pim.mag();
				ptrk_pim.setE(sqrt(p3*p3+mpion*mpion));

				m_pid_ppim = ptrk_pim.rho();
				m_pid_pim_px = ptrk_pim.px();
				m_pid_pim_py = ptrk_pim.py();
				m_pid_pim_pz = ptrk_pim.pz();
				double pim_ep=0;
				if((*itTrk)->isEmcShowerValid())	{
					RecEmcShower *emcTrk = (*itTrk)->emcShower();
					pim_ep = emcTrk->energy()/ptrk_pim.rho();
				}
				m_pim_ep = pim_ep;

				HepVector a = mdcTrk->helix();
				HepSymMatrix Ea = mdcTrk->err();
				HepPoint3D point0(0.,0.,0.);	//The initial point for MDC recosntruction
				HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
				VFHelix helixip(point0, a, Ea);
				helixip.pivot(IP);
				HepVector vecipa = helixip.a();
				m_pimRvxy0 = fabs(vecipa[0]);	//the nearest distance to IP in xy plane
				m_pimRvz0 = vecipa[3];		//the nearest distance to IP in z direction
			}
			if( prob_p > prob_K && prob_p > prob_pi )	{
				ipm.push_back(iGood[i]);
				RecMdcKalTrack::setPidType (RecMdcKalTrack:proton);
				HepLorentzVector ptrk_pm;
				ptrk_pm.setPx(mdcKalTrk->px());
				ptrk_pm.setPy(mdcKalTrk->py());
				ptrk_pm.setPz(mdcKalTrk->pz());
				double p3 = ptrk_pm.mag();
				ptrk_pm.setE(sqrt(p3*p3+mproton*mproton));
				m_pip_pm = ptrk_pm.rho();
				HepVector a = mdcTrk->helix();
				HepSymMatrix Ea = mdcTrk->err();
				HepPoint3D point0(0.,0.,0.);	//The initial point for MDC recosntruction
				HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
				VFHelix helixip(point0, a, Ea);
				helixip.pivot(IP);
				HepVector vecipa = helixip.a();
				m_apRvxy0 = fabs(vecipa[0]);	//the nearest distance to IP in xy plane
				m_apRvz0 = vecipa[3];		//the nearest distance to IP in z direction
			}
			if( prob_K>prob_p && prob_K > prob_pi)	{
				iKm.push_back(iGood[i]);
			}
		}
	}

	int npip = ipip.size();
	int npim = ipip.size();
	int npp  = ipp.size();
	int npm  = ipm.size();
	m_npip = npip;
	m_npim = npim;
	m_npp = ipp.size();
	m_npm = ipm.size();
	m_nKp = iKp.size();
	m_nKm = iKm.size();
	hnpip->Fill(npip);
	hnpim->Fill(npim);

	if ( npim*npip != 1 ) return SUCCESS;
	Ncut_pid++;
	//==============================================================
	Vint ishower;
	double etot = 0;
	ishower.clear();
	int goodgam = 0;

	for(int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++)	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
		if ( !(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
		double dang = 200;
		for ( int j = 0; j < evtRecEvent->totalCharged(); j++)	{
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();
			double angd = extpos.angle(emcpos);
			if(angd < dang)	{
				dang = angd;
			}
		}
		if( dang != 200)	{
			dang = dang*180/ (CLHEP::pi);
		}
		double ctht = cos(emcTrk->theta());
		if((emcTrk->module() == 1) && (fabs(ctht) > 0.8 )) continue;
		if((emcTrk->module() == 0 || emcTrk->module() == 2) && (fabs(ctht) > 0.92 || fabs(ctht) < 0.86)) continue;
		if(fabs(dang) < m_gammaAngleCut) continue;

		//??
		if (emcTrk->module() == 1){
			if(emcTrk->energy()<m_energyThreshold1) continue;
		}
		if (emcTrk->module() == 0 || emcTrk->module() == 2){
			if(emcTrk->energy()<m_energyThreshold2) continue;
		}
		//??
		etot+=emcTrk->energy();
		ishower.push_back(i);
	}
	int nshower = ishower.size();
	if(nshower<1) return SUCCESS;
	Ncut_nshower++;


	double energy = 0;
	int index = -1;
	for(int k = 0; k<nshower; k++)	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin+ishower[k];
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		double eraw = emcTrk->energy();
		if ( eraw > 2.0 ) continue;
		if ( eraw > energy )	{
			energy = eraw;
			index = k;
		}
	}
	if ( index == -1 ) return SUCCESS;
	RecEmcShower *nbarTrk = (*(evtRecTrkCol->begin() + ishower[index]))->emcShower();
	Hep3Vector nbarpos(nbarTrk->x(),nbarTrk->y(), nbarTrk->z());
	Hep3Vector nbar_pos_40d;
	double nbar_ene_40d = 0;
	m_nbar_energy = nbarTrk->energy();
	m_nbar_hits   = nbarTrk->numHits();
	m_nbar_secmom = nbarTrk->secondMoment();
	m_nbar_cos    = cos(nbarTrk->theta());
	for( int i = 0; i < nshower; i++)	{
		if(i == index) continue;
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ishower[i];
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
		double angle = emcpos.angle(nbarpos)*180.0/3.1415926;
		if(angle<40)	{
			nbar_ene_40d += emcTrk->energy();
			nbar_pos_40d += emcTrk->energy()*emcpos;
		}
	}
	nbar_pos_40d = (nbar_pos_40d + nbarTrk->energy()*nbarpos)/(nbar_ene_40d + nbarTrk->energy());
	if (nbarTrk->energy() < 0.2 || nbarTrk->numHits() < 20 || nbarTrk->secondMoment() < 20) return SUCCESS;
	Ncut_nbar++;
	//=======================================================
	HepLorentzVector ecms;
	ecms.setPy(0); 
	ecms.setPz(0);
	ecms.setE(cms);
	ecms.setPx(cms*sin(0.022*0.5));

	RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack();
	RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();
	
	WTrackParameter wvpimTrk = WTrackParameter(mpion, pimTrk->getZHelix(), pimTrk->getZError());
	HepLorentzVector p4pim = wvpimTrk.p();
	WTrackParameter wvpipTrk = WTrackParameter(mpion, pipTrk->getZHelix(), pipTrk->getZError());
	HepLorentzVector p4pip = wvpipTrk.p();

	HepLorentzVector p4pipi = p4pip + p4pim;
	HepLorentzVector p4pipirecoil = ecms - p4pip - p4pim;

	m_mpipi = p4pipi.m();
	m_mpipirecoil = p4pipirecoil.m();

	double chi2 = 999.;
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	kmfit->init();
	kmfit->AddTrack(0,wvpimTrk);
	kmfit->AddTrack(1,wvpipTrk);
	kmfit->AddMissTrack(2, nbar_mass, nbarTrk);
	kmfit->AddMissTrack(3, n_mass);
	kmfit->AddResonance(0, Sigmam_mass, 1, 2);
	kmfit->AddFourMomentum(1, ecms);
	bool oksq = kmfit->Fit();
	if(oksq)	{
		chi2 = kmfit->chisq();
		HepLorentzVector p4pim = kmfit->pfit(0);
		HepLorentzVector p4pip = kmfit->pfit(1);
		HepLorentzVector p4nbar = kmfit->pfit(2);
		HepLorentzVector p4n = kmfit->pfit(3);

		HepLorentzVector p4sigmabar = p4pip + p4nbar;
		HepLorentzVector p4sigma = ecms - p4sigmabar;
		m_Sigmac_m = p4sigma.m();
		m_Sigmac_m = p4sigma.rho();
		m_Sigmac_cos = cos(p4sigma.theta());
		
		m_aSigmac_m = p4sigmabar.m();
		m_aSigmac_m = p4sigmabar.rho();
		m_aSigmac_cos = cos(p4sigmabar.theta());

		m_pxpim = p4pim.px();		m_pypim = p4pim.py();
		m_pzpim = p4pim.pz();		m_epim  = p4pim.e();

		*/




			

/*

	if(chi2 == 999 && chi2v2 == 999) return SUCCESS;
	Ncut_kin++;
*/
//	m_tuple1->write();
	return StatusCode::SUCCESS;

}

//---------------------------------------------------------------------------------
StatusCode sigmamsigmap::finalize()	{
	cout << "Selection criteria:	" << "Survived"    << endl;
	cout << "Total number:			" << Ncut_none     << endl;
	cout << "cut mc:				" << Ncut_mc       << endl;
	cout << "Charged trk:			" << Ncut_ncharge  << endl;
	cout << "Pid sel:				" << Ncut_pid      << endl;
//	cout << "" <<   << endl;
//	cout << "" <<   << endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}



