#ifndef Physics_Analysis_TwoProng_H
#define Physics_Analysis_TwoProng_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"
//#include "EventTag/EventTagSvc.h"
class TwoProng : public Algorithm {

public:
  TwoProng(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

private:



  //   ReadBeamParFromDb m_reader;
  bool m_runmc;
  bool m_mcmatch;
  double m_ecms;
  double m_espread;
  double m_beamangle;
  double mpi0;
  double mrho;
  double meta;
  double metap;
  double m_vr0cut;
  double m_vz0cut;
  double m_coscut;
  double m_p0cut;

  bool m_useBarrel;
  bool m_useEndcap;
  double m_energyBThreshold;
  double m_energyEThreshold;
  double m_gammaPhiCut;
  double m_gammaThetaCut;
  double m_gammaTrkCut;
  double m_gammathCut;
  double m_gammatlCut;
  int  m_gammanumlCut;
  int  m_gammanumhCut;


  double m_chisq6ccut;

  NTuple::Tuple*  m_tuple0;      // truth
  NTuple::Array<double>  m_truthppip;
  NTuple::Array<double>  m_truthppim;
  NTuple::Item<double>  m_dang_min;
  NTuple::Item<double>  m_dang_e;
  
  NTuple::Tuple*  m_tuple1;      // charged track vertex
  NTuple::Item<double>  m_vx0;
  NTuple::Item<double>  m_vy0;
  NTuple::Item<double>  m_vz0;
  NTuple::Item<double>  m_vr0;
  NTuple::Item<double>  m_rvxy0;
  NTuple::Item<double>  m_rvz0;
  NTuple::Item<double>  m_rvphi0;
  //NTuple::Item<double>  m_costheta;
  NTuple::Item<double>  m_phi_mdc;
  NTuple::Item<double>  m_pt;
  NTuple::Item<double>  m_rvx0;
  NTuple::Item<double>  m_rvy0;
  NTuple::Item<double>  m_phi0;
 

  NTuple::Tuple*  m_tuple4;  

  NTuple::Array<double> m_mcgpx;
  NTuple::Array<double> m_mcgpy; 
  NTuple::Array<double> m_mcgpz;
  NTuple::Array<double> m_mcge;

  NTuple::Item<long>  m_ncharg;
  NTuple::Item<long>  m_nneu;
  NTuple::Item<long>  m_ngch; 
  NTuple::Item<long>  m_run; 
  NTuple::Item<long>  m_rec;
  NTuple::Item<long>  m_evttag; 
  NTuple::Item<long>   m_idxmc;
  NTuple::Array<int>  m_pdgid;
  NTuple::Array<int>  m_motheridx;
  NTuple::Item<int>   m_ISRtag;

  NTuple::Item<long>    m_ntof1;
  NTuple::Array<int>   m_toflayer1;
  NTuple::Array<double> m_tof1;
  NTuple::Array<double> m_t01;
  NTuple::Item<long>    m_ntof2;
  NTuple::Array<int>   m_toflayer2;
  NTuple::Array<double> m_tof2;
  NTuple::Array<double> m_t02;
  NTuple::Item<long>   m_mustatus;
  NTuple::Item<long>   m_emcstatus;
  NTuple::Item<double>   m_emctrk1;
  NTuple::Item<double>   m_emctrk2;
  NTuple::Item<double>   m_epratio1;
  NTuple::Item<double>   m_epratio2;

  NTuple::Item<double> m_kappx;
  NTuple::Item<double> m_kappy; 
  NTuple::Item<double> m_kappz;
  NTuple::Item<double> m_kape;

  NTuple::Item<double> m_kampx;
  NTuple::Item<double> m_kampy; 
  NTuple::Item<double> m_kampz;
  NTuple::Item<double> m_kame;

  NTuple::Item<double> m_costheta;
//NTuple::Array<double> m_pippx;
//NTuple::Array<double> m_pippy; 
//NTuple::Array<double> m_pippz;
//NTuple::Array<double> m_pipe;

//NTuple::Array<double> m_pimpx;
//NTuple::Array<double> m_pimpy; 
//NTuple::Array<double> m_pimpz;
//NTuple::Array<double> m_pime;



  //  IEventTagSvc* m_EventTagSvc;


};

#endif 
