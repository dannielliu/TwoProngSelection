#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"



#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"




#include "McTruth/McParticle.h"
#include "McTruth/DecayMode.h"
#include "McTruth/MdcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/MucMcHit.h"
#include "McTruth/McEvent.h"


#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "TwoProngAlg/TwoProng.h"
#include "VertexFit/Helix.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "ParticleID/ParticleID.h"
#include "VertexFit/SecondVertexFit.h"

#include <vector>
//const double twopi = 6.2831853;
//const double pi = 3.1415927;
//const double mpi0 = 0.1349766;
const double me = 0.000511;
const double mu = 0.105658;
const double mpi = 0.13957;
const double mka = 0.493677;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
//const double velc = 29.9792458;  tof_path unit in cm.
const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Hep3Vector> Vp3;
//const HepLorentzVector p_cms(0.034067,0.0,0.0,3.097);
//const Hep3Vector u_cms = -p_cms.boostVector();
static int counter[10]={0,0,0,0,0,0,0,0,0,0};
/////////////////////////////////////////////////////////////////////////////

TwoProng::TwoProng(const std::string& name, ISvcLocator* pSvcLocator) :
		Algorithm(name, pSvcLocator) {

				//Declare the properties 
				declareProperty("runmc",m_runmc = false);
				declareProperty("mcmatch",m_mcmatch = true);  
				declareProperty("ecms",m_ecms = 3.686);
				declareProperty("espread",m_espread=0.0013);
				declareProperty("beamangle",m_beamangle = 0.022); 

				declareProperty("mpi0",mpi0 = 0.1349766); 
				declareProperty("meta",meta = 0.547853); 
				declareProperty("metap",metap = 0.95766);
				declareProperty("mrho",mrho = 0.77549);
                //charged tracks
				declareProperty("Vr0cut", m_vr0cut=1.0);
				declareProperty("Vz0cut", m_vz0cut=10.0);
				//declareProperty("Vr0cut", m_vr0cut=10.0);
				//declareProperty("Vz0cut", m_vz0cut=30.0);
				declareProperty("Coscut", m_coscut=0.93);
                                declareProperty("MomentumCut", m_p0cut=2.0);               
                //gamma  
				declareProperty("UseBarrelEmc",    m_useBarrel        = true);
				declareProperty("UseEndcapEmc",    m_useEndcap        = true);
				declareProperty("BarrelEnergyThreshold", m_energyBThreshold=0.025);
				declareProperty("EndcapEnergyThreshold", m_energyEThreshold=0.025);
				declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
				declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
				declareProperty("GammaTrkCut", m_gammaTrkCut=0.0);
				declareProperty("GammaTimeWinHCut", m_gammathCut=14.0);
				declareProperty("GammaTimeWinLCut", m_gammatlCut=0.0);
				declareProperty("GammaNumLCut", m_gammanumlCut=2);
				declareProperty("GammaNumHCut", m_gammanumhCut=8);

	     //final 6c
				declareProperty("chisq6ccut", m_chisq6ccut=200.0);

		}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode TwoProng::initialize(){
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;

	StatusCode status;
	//   static const bool CREATEIFNOTTHERE(true);
	//   StatusCode EventTagSvcStatus = service("EventTagSvc", m_EventTagSvc, CREATEIFNOTTHERE);
	//   if (!EventTagSvcStatus.isSuccess() || 0 ==m_EventTagSvc ) {
	//     log << MSG::ERROR << " Could not initialize Decay code service" << endreq;
	//     return EventTagSvcStatus;
	//   }
	NTuplePtr nt0(ntupleSvc(), "GETAPETAP/mctruth");
	if ( nt0 ) m_tuple0 = nt0;
	else {
		m_tuple0 = ntupleSvc()->book ("GETAPETAP/mctruth", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple0 )    {
			status = m_tuple0->addItem ("truthppip", 4,  m_truthppip);
			status = m_tuple0->addItem ("truthppim", 4,  m_truthppim);
			status = m_tuple0->addItem ("mcgpx", 3,  m_mcgpx);
			status = m_tuple0->addItem ("mcgpy" , 3,  m_mcgpy);
			status = m_tuple0->addItem ("mcgpz" , 3,  m_mcgpz);
			status = m_tuple0->addItem ("mcge" , 3,  m_mcge);
			status = m_tuple0->addItem ("dang_min", m_dang_min);
			status = m_tuple0->addItem ("dang_e", m_dang_e);
		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple0) << endmsg;
			return StatusCode::FAILURE;
		}
	}        





	NTuplePtr nt1(ntupleSvc(), "GETAPETAP/vxyz");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("GETAPETAP/vxyz", CLID_ColumnWiseTuple, "kk N-Tuple example");
		if ( m_tuple1 )    {
			status = m_tuple1->addItem ("vx0",   m_vx0);
			status = m_tuple1->addItem ("vy0",   m_vy0);
			status = m_tuple1->addItem ("vz0",   m_vz0);
			status = m_tuple1->addItem ("vr0",   m_vr0);
			status = m_tuple1->addItem ("rvxy0",  m_rvxy0);
			status = m_tuple1->addItem ("rvz0",   m_rvz0);
			status = m_tuple1->addItem ("rvphi0", m_rvphi0);
			//status = m_tuple1->addItem ("costheta",m_costheta);
			status = m_tuple1->addItem ("phi", m_phi_mdc);
			status = m_tuple1->addItem ("p", m_pt);
			status = m_tuple1->addItem ("rvx0", m_rvx0);
			status = m_tuple1->addItem ("rvy0", m_rvy0);
			status = m_tuple1->addItem ("phi0", m_phi0);
		} 
		else    { 
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}



	NTuplePtr nt4(ntupleSvc(), "GETAPETAP/TwoProng");
	if ( nt4 ) m_tuple4 = nt4;
	else {
		m_tuple4 = ntupleSvc()->book ("GETAPETAP/TwoProng", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple4 )    {
			status = m_tuple4->addItem ("run",  m_run); 
			status = m_tuple4->addItem ("rec",  m_rec); 
			status = m_tuple4->addItem ("evttag",  m_evttag); 
			status = m_tuple4->addItem("indexmc",          m_idxmc, 0, 100);
			status = m_tuple4->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
			status = m_tuple4->addIndexedItem("motheridx", m_idxmc, m_motheridx);    
                        status = m_tuple4->addItem ("ngch" ,   m_ngch);
			status = m_tuple4->addItem ("ncharg",   m_ncharg);
			status = m_tuple4->addItem ("nneu",    m_nneu);

			status = m_tuple4->addItem ("ntof1", m_ntof1,0,5);
			status = m_tuple4->addIndexedItem("toflayer1", m_ntof1, m_toflayer1);
			status = m_tuple4->addIndexedItem("tof1", m_ntof1, m_tof1);
			status = m_tuple4->addIndexedItem("t01", m_ntof1, m_t01);
			status = m_tuple4->addItem ("ntof2", m_ntof2,0,5);
			status = m_tuple4->addIndexedItem("toflayer2", m_ntof2, m_toflayer2);
			status = m_tuple4->addIndexedItem("tof2", m_ntof2, m_tof2);
			status = m_tuple4->addIndexedItem("t02", m_ntof2, m_t02);
			status = m_tuple4->addItem ("mustatus", m_mustatus);
			status = m_tuple4->addItem ("emcstatus", m_emcstatus);
			status = m_tuple4->addItem ("emctrk1", m_emctrk1);
			status = m_tuple4->addItem ("emctrk2", m_emctrk2);
			status = m_tuple4->addItem ("epratio1", m_epratio1);
			status = m_tuple4->addItem ("epratio2", m_epratio2);

			status = m_tuple4->addItem ("kappx" , m_kappx);
			status = m_tuple4->addItem ("kappy" , m_kappy);
			status = m_tuple4->addItem ("kappz" , m_kappz);
			status = m_tuple4->addItem ("kape" ,  m_kape) ;

			status = m_tuple4->addItem ("kampx" , m_kampx);
			status = m_tuple4->addItem ("kampy" , m_kampy);
			status = m_tuple4->addItem ("kampz" , m_kampz);
			status = m_tuple4->addItem ("kame" ,  m_kame) ;

			status = m_tuple4->addItem ("costheta",  m_costheta) ;

		}
		else    { 
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
			return StatusCode::FAILURE;
		}
	}
	//
	//--------end of book--------
	//

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode TwoProng::execute() {
	const double beamEnergy = m_ecms/2.;
	const HepLorentzVector p_cms(m_ecms*sin(m_beamangle*0.5),0.0,0.0,m_ecms);
	//const HepLorentzVector p_cms(0.040547,0,0,3.68632);
	const Hep3Vector u_cms = -p_cms.boostVector(); 



	StatusCode sc=StatusCode::SUCCESS;

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	counter[0]++;


	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	if (!eventHeader) 
	{
		log << MSG::FATAL << "Could not find Event Header" << endreq;
		return StatusCode::SUCCESS;
	}

	m_run = eventHeader->runNumber();
	m_rec = eventHeader->eventNumber();
	m_evttag=eventHeader->eventTag();
	//if( m_run >= 9755 && m_run<= 9766 ) return StatusCode::SUCCESS;
	//if( m_run >= 34298 && m_run<= 34313 ) return StatusCode::SUCCESS;
	//if( m_run >= 34461 && m_run<= 34477 ) return StatusCode::SUCCESS;
	//if( m_run >= 35086 && m_run<= 35098 ) return StatusCode::SUCCESS;


	Hep3Vector vgammamc_ext, vpi0_gammamc_high_ext, vpi0_gammamc_low_ext;
	if (eventHeader->runNumber()<0)
	{
		//MC information
		Vp4 mcpGam,mcpGam1,mcpGam2, mcppip,mcppim;

		mcpGam.clear();
		mcpGam1.clear();
		mcpGam2.clear();
		mcppip.clear();
		mcppim.clear();
		Vint mcIEta,mcIGam;
		Vp3 vmcGam;
		mcIEta.clear();
		mcIGam.clear();
		vmcGam.clear();
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");

		int m_numParticle = 0;
		if (!mcParticleCol)
		{
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}

		bool psipDecay = false;
		int rootIndex = -1;

		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			HepLorentzVector  mctrue_track = (*iter_mc)->initialFourMomentum();
			HepLorentzVector mctrack_iniposition = (*iter_mc)->initialPosition();
			HepLorentzVector mctrack_finposition = (*iter_mc)->finalPosition();

			if ((*iter_mc)->primaryParticle()) continue;
			if (!(*iter_mc)->decayFromGenerator()) continue;
			//if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
			int imom=443;
			if(m_ecms<3.2)imom=443;
			else if (m_ecms>3.6)imom=100443;
			if ((*iter_mc)->particleProperty()==imom)  
				//if ((*iter_mc)->particleProperty()==100443) 
			{
				psipDecay = true;
				rootIndex = (*iter_mc)->trackIndex();
				//    m_evttag=(m_EventTagSvc->getChainCode(*iter_mc))<<8;

			}
			if (!psipDecay) continue;
			int mcidxi=(*iter_mc)->trackIndex() - rootIndex;	
			int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
			int pdgid = (*iter_mc)->particleProperty();
			m_pdgid[m_numParticle] = pdgid;
			m_motheridx[m_numParticle] = mcidx;
			m_numParticle += 1; 
			if(m_mcmatch){
				if( (*iter_mc)->particleProperty() == 321){
					mcppip.push_back(mctrue_track);
					m_truthppip[0] = mctrue_track.px();
					m_truthppip[1] = mctrue_track.py();
					m_truthppip[2] = mctrue_track.pz();
					m_truthppip[3] = mctrue_track.e();
				}

				if( (*iter_mc)->particleProperty() == -321){
					mcppim.push_back(mctrue_track);
					m_truthppim[0] = mctrue_track.px();
					m_truthppim[1] = mctrue_track.py();
					m_truthppim[2] = mctrue_track.pz();
					m_truthppim[3] = mctrue_track.e();
				}

				if( (*iter_mc)->particleProperty() == 22 &&
						((*iter_mc)->mother()).particleProperty()==100443){//443 for J/psi;100443 for psi(2S)
							m_mcgpx[2]  = mctrue_track.px();
							m_mcgpy[2]  = mctrue_track.py();
							m_mcgpz[2]  = mctrue_track.pz();
							m_mcge[2]  = mctrue_track.e();
							vgammamc_ext=(*iter_mc)->finalPosition().vect();
							vgammamc_ext[0] = mctrue_track.px();
							vgammamc_ext[1] = mctrue_track.py();
							vgammamc_ext[2] = mctrue_track.pz();
						}
				if( (*iter_mc)->particleProperty() == 22 &&
						((*iter_mc)->mother()).particleProperty()==111){//111 for pi0; 221 for eta

					mcIGam.push_back(mcidx);
					mcpGam.push_back(mctrue_track);
					vmcGam.push_back((*iter_mc)->finalPosition().vect());
					//					cout<<"(*iter_mc)->finalPosition().vect() ="<<(*iter_mc)->finalPosition().vect()<<endl;
				}
				if( (*iter_mc)->particleProperty() == 111 &&
						((((*iter_mc)->mother()).particleProperty()==10441)
						 ||(((*iter_mc)->mother()).particleProperty()==20443)||(((*iter_mc)->mother()).particleProperty()==445))
				  ){
					mcIEta.push_back(mcidxi);
				}
			} 
		}
		//cout<<"mcIGam.size() ="<<mcIGam.size()<<endl;
		//cout<<"mcIEta.size() ="<<mcIEta.size()<<endl;


		m_idxmc = m_numParticle;

		if(m_mcmatch&&mcIGam.size()==2&&mcIEta.size()==1){

			mcpGam1.push_back(mcpGam[0]);
			mcpGam2.push_back(mcpGam[1]);
			if(mcpGam1.size()==1&&mcpGam2.size()==1){

				//   std::cout<<":::"<<mcIEta.size()<<", "<<mcIGam.size()<<", "<<mcpGam1.size()<<", "<<mcpGam2.size()<<std::endl;
				double mce1=mcpGam1[0].e();
				double mce2=mcpGam2[0].e();
				int mcee=(mce1>=mce2)?1:2;

				if (mcee==1){

					m_mcgpx[0]  = mcpGam1[0].px();
					m_mcgpy[0]  = mcpGam1[0].py();
					m_mcgpz[0]  = mcpGam1[0].pz();
					m_mcge[0]  = mcpGam1[0].e();
					vpi0_gammamc_high_ext = vmcGam[0];      
					vpi0_gammamc_low_ext = vmcGam[1];
					m_mcgpx[1]  = mcpGam2[0].px();
					m_mcgpy[1]  = mcpGam2[0].py();
					m_mcgpz[1]  = mcpGam2[0].pz();
					m_mcge[1]  = mcpGam2[0].e();
				}
				else{
					m_mcgpx[0]  = mcpGam2[0].px();
					m_mcgpy[0]  = mcpGam2[0].py();
					m_mcgpz[0]  = mcpGam2[0].pz();
					m_mcge[0]  = mcpGam2[0].e();
					vpi0_gammamc_high_ext = vmcGam[1];
					vpi0_gammamc_low_ext = vmcGam[0];
					m_mcgpx[1]  = mcpGam1[0].px();
					m_mcgpy[1]  = mcpGam1[0].py();
					m_mcgpz[1]  = mcpGam1[0].pz();
					m_mcge[1]  = mcpGam1[0].e();
				}
			}

		}

		//		m_tuple0->write();               
	}

	//cout<<"vpi0_gammamc_high_ext ="<<vpi0_gammamc_high_ext<<endl;
	//cout<<"vpi0_gammamc_low_ext ="<<vpi0_gammamc_low_ext<<endl;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	if (!evtRecEvent) 
	{
		log << MSG::FATAL << "Could not find EvtRecEvent" << endreq;
		return StatusCode::SUCCESS;
	}
	log << MSG::INFO <<"ncharg, nneu, tottks = " 
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() <<endreq;
	//  if(evtRecEvent->totalNeutral()>30)return sc;

	m_ncharg  = evtRecEvent->totalCharged();

	m_nneu = evtRecEvent->totalNeutral();



        HepPoint3D vx(0., 0., 0.);
        HepSymMatrix Evx(3, 0);

        IVertexDbSvc*  vtxsvc;
        Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
        if(vtxsvc->isVertexValid()){
        	double* dbv = vtxsvc->PrimaryVertex();
        	double*  vv = vtxsvc->SigmaPrimaryVertex();
        	//	if (m_reader.isRunNumberValid( m_run)) {
        	//   HepVector dbv = m_reader.PrimaryVertex( m_run);
        	//    HepVector vv = m_reader.SigmaPrimaryVertex( m_run);
        	vx.setX(dbv[0]);
        	vx.setY(dbv[1]);
        	vx.setZ(dbv[2]);
        	Evx[0][0]=vv[0]*vv[0];
        	Evx[0][1]=vv[0]*vv[1];
        	Evx[1][1]=vv[1]*vv[1];
        	Evx[1][2]=vv[1]*vv[2];
        	Evx[2][2]=vv[2]*vv[2];
        }

        VertexParameter vx_db;
        vx_db.setVx(vx);
        vx_db.setEvx(Evx);




	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
	if (!evtRecTrkCol) 
	{
		log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
		return StatusCode::SUCCESS;
	}

	Vint iGood, ikap, ikam,itrkp, itrkm, ipip, ipim;
	iGood.clear();
	ikap.clear();
	ikam.clear();
	itrkp.clear();
	itrkm.clear();
	ipip.clear();
	ipim.clear();

	Vp4 pkap, pkam;
	pkap.clear();
	pkam.clear();

	int nCharge = 0;
	for(int i = 0; i < evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		//if(!(*itTrk)->isMdcKalTrackValid()) continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		double pch=mdcTrk->p();
		double x0=mdcTrk->x();
		double y0=mdcTrk->y();
		double z0=mdcTrk->z();
		double costheta = mdcTrk->theta();
		double phi_mdc = mdcTrk->phi();
		double phi0=mdcTrk->helix(1);
		double xv=vx.x();
		double yv=vx.y();
		double zv=vx.z();
		double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
		double m_vx0 = x0;
		double m_vy0 = y0;
		double m_vz0 = z0;
		double m_vr0 = Rxy;
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xv,yv,zv); 
		VFHelix helixip(point0,a,Ea); 
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=vecipa[1];
		m_rvxy0=Rvxy0;
		m_rvz0=Rvz0;
		//		m_rvphi0=Rvphi0;
		//m_costheta = costheta;
		//		m_phi_mdc = phi_mdc;
		m_pt = pch;
		m_rvx0 = Rvxy0 * cos(Rvphi0);
		m_rvy0 = Rvxy0 * sin(Rvphi0);
		m_phi0 = Rvphi0;	
		//		m_tuple1->write();
		if(fabs(Rvz0) >= m_vz0cut) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;
		//		if(pch >= m_p0cut) continue;
		double cost = cos(mdcTrk->theta());
		if(fabs(cost) >= m_coscut ) continue;

		iGood.push_back(i);//   iGood.push_back((*itTrk)->trackId());
		nCharge += mdcTrk->charge();

	}


	//
	// Finish Good Charged Track Selection
	//
	int nGood = iGood.size();
	m_ngch=nGood;
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;

	//if( !((nGood ==3) ||(nGood ==6 && nCharge ==0))){
	if( nGood != 2 || nCharge != 0){
		return StatusCode::SUCCESS;
	}
	counter[1]++;

	int ntof1=0;
	int ntof2=0;
	int emcstatus=0x0;
	int mustatus=0x0;
	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		if (!(*itTrk)->isMdcTrackValid()) continue;
		if (!(*itTrk)->isTofTrackValid()) continue;
		SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
		SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();

		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack

		if(mdcTrk->charge() >0) {
			ikap.push_back(iGood[i]);
			if ((*itTrk)->isEmcShowerValid()){
				emcstatus = emcstatus | 0x1;
				double ep_ratio =  (*itTrk)->emcShower()->energy()/mdcTrk->p();
				double emc_p = (*itTrk)->emcShower()->energy();
				m_epratio1 = ep_ratio;
				m_emctrk1 = (*itTrk)->emcShower()->energy();
			}
			for (;iter_tof!=tofTrkCol.end();iter_tof++){
				TofHitStatus *status = new TofHitStatus;
				status->setStatus((*iter_tof)->status());
				if(!(status->is_counter())) continue;
				//if(status->layer()!=0) continue; // for endcap
				double tof = (*iter_tof)->tof();
				m_toflayer1[ntof1] = status->layer();
				m_tof1[ntof1] = tof;
				m_t01[ntof1] = (*iter_tof)->t0();
				ntof1++;
				delete status;
			}
			if ( (*itTrk)->isMucTrackValid()) mustatus = mustatus | 0x1;;
		} else {
			ikam.push_back(iGood[i]);
			if ((*itTrk)->isEmcShowerValid()){
				emcstatus = emcstatus | 0x2;
				double ep_ratio =  (*itTrk)->emcShower()->energy()/mdcTrk->p();
				double emc_p = (*itTrk)->emcShower()->energy();
				m_epratio2 = ep_ratio;
				m_emctrk2 = (*itTrk)->emcShower()->energy();
			}
			for (;iter_tof!=tofTrkCol.end();iter_tof++){
				TofHitStatus *status = new TofHitStatus;
				status->setStatus((*iter_tof)->status());
				if(!(status->is_counter())) continue;
				//if(status->layer()!=0) continue; // for endcap
				double tof = (*iter_tof)->tof();
				m_toflayer2[ntof2] = status->layer();
				m_tof2[ntof2] = tof;
				m_t02[ntof1] = (*iter_tof)->t0();
				ntof2++;
				delete status;
			}
			if ( (*itTrk)->isMucTrackValid()) mustatus = mustatus | 0x2;;

		}
	}
	m_ntof1 = ntof1;
	m_ntof2 = ntof2;
	//std::cout<<"ntof "<<ntof<<std::endl;
	m_emcstatus = emcstatus;
	m_mustatus = mustatus;

	int nkap = ikap.size();
	int nkam = ikam.size();

	//log << MSG::INFO << "PID ka+, ka- = " << nkap << " , " << nkam << endreq;
	if(nkap != 1) return sc;
	if(nkam != 1) return sc;

	counter[2]++;


	RecMdcKalTrack *kamTrk = (*(evtRecTrkCol->begin()+ikam[0]))->mdcKalTrack();
	WTrackParameter wvkamTrk = WTrackParameter(mka, kamTrk->getZHelix(), kamTrk->getZError());
	HepLorentzVector p4_kam = wvkamTrk.p();
	m_kampx = p4_kam.px();
	m_kampy = p4_kam.py();
	m_kampz = p4_kam.pz();
	m_kame  = p4_kam.e();

	RecMdcKalTrack *kapTrk = (*(evtRecTrkCol->begin()+ikap[0]))->mdcKalTrack();
	WTrackParameter wvkapTrk = WTrackParameter(mka, kapTrk->getZHelix(), kapTrk->getZError());
	HepLorentzVector p4_kap = wvkapTrk.p();
	m_kappx = p4_kap.px();
	m_kappy = p4_kap.py();
	m_kappz = p4_kap.pz();
	m_kape  = p4_kap.e();

	m_costheta = p4_kam.vect().cosTheta(p4_kap.vect());

	m_tuple4->write();

	return StatusCode::SUCCESS;
	}


	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	StatusCode TwoProng::finalize() {
		cout<<"ee to K K : TwoProng  "<<endl;
		cout<<"total number:         "<<counter[0]<<endl;
		cout<<"good track:           "<<counter[1]<<endl;
		cout<<"2 track:              "<<counter[2]<<endl;
		MsgStream log(msgSvc(), name());
		log << MSG::INFO << "in finalize()" << endmsg;
		return StatusCode::SUCCESS;
	}

