#include <utility>
#include <limits>
#include <TLorentzVector.h>
#include <MEUtils/CFourVector.hpp>
#include <MElCalc/MElCalc.hpp>
#include "KinFit5Pi.hpp"

const std::vector<std::string> KinFit5Pi::_seqTracks = {"pi-_0", "pi-_1", "pi+_0", "pi+_1"};
const std::vector<std::string> KinFit5Pi::_seqPhotons = {"g0", "g1"};
const std::set<std::string> KinFit5Pi::_permPhotons = {"g0", "g1"};
const std::set<std::string> KinFit5Pi::_permTracks = {"pi-_0", "pi-_1", "pi+_0", "pi+_1"};
const std::set<std::string> KinFit5Pi::_permAll = {"pi-_0", "pi-_1", "pi+_0", "pi+_1", "g0", "g1"};
const std::vector<std::set<std::string>> KinFit5Pi::_perm3Pi =
  {{"pi-_0", "pi+_0", "g0", "g1"},
   {"pi-_0", "pi+_1", "g0", "g1"},
   {"pi-_1", "pi+_0", "g0", "g1"},
   {"pi-_1", "pi+_1", "g0", "g1"}};
const std::vector<std::set<std::string>> KinFit5Pi::_permPiPlPiMi =
  {{"pi-_0", "pi+_0"},
   {"pi-_0", "pi+_1"},
   {"pi-_1", "pi+_0"},
   {"pi-_1", "pi+_1"}};

KinFit5Pi::KinFit5Pi(TTree* tree) :
  KFCmd::TrPh(tree) {
}

KinFit5Pi::~KinFit5Pi() {}

bool KinFit5Pi::cutTracks() {
  _trackIndices.clear();
  for (int i = 0; i < nt; i++) {
    bool point =
      (std::fabs(tz[i]) < 12.0) &&
      (std::fabs(trho[i]) < 1.0);
    bool dedx =
      (tdedx[i] > 0.0) &&
      (tdedx[i] < 15000.0);
    bool ptot =
      (tptot[i] > 5.0) &&
      (tptot[i] < 1000.0);
    if (point && dedx && ptot) _trackIndices.push_back(i);
  }
  if (_trackIndices.size() == 4) {
    int totalCharge = 0;
    for (int i = 0; i < 4; ++i) totalCharge += tcharge[_trackIndices[i]];
    return (totalCharge == 0);
  }
  return false;
}

bool KinFit5Pi::cutPhotons() {
  _photonIndices.clear();
  for (int i = 0; i < nph; i++)
    if ((phen[i] > 20.0) && (phen[i] < 1000.0))
      _photonIndices.push_back(i);
  if (_photonIndices.size() >= 2) return true;
  return false;
}

bool KinFit5Pi::cut() {
  if (!cutTracks()) return false;
  if (!cutPhotons()) return false;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(_trackIndices.begin(), _trackIndices.end(),
	    [&charges](int i, int j) { return charges[i] < charges[j]; });
  return true;
}

void KinFit5Pi::fill_fit5pi(KFCmd::Hypo4ChPions2Photons* hypo, std::size_t iph, std::size_t jph) {
  kf_err_5pi = hypo->getErrorCode();
  kf_chi2_5pi = hypo->getChiSquare();
  kf_chi2ph_5pi = hypo->getChiSquare(_permPhotons);
  kf_chi2t_5pi = hypo->getChiSquare(_permTracks);
  std::copy(_trackIndices.begin(),
	    _trackIndices.end(), tind);
  phind[0] = _photonIndices[iph];
  phind[1] = _photonIndices[jph];
  TLorentzVector in_tm_recoil(0, 0, 0, 2 * emeas);
  in_tm_recoil -= hypo->getInitialMomentum(_permTracks);
  TLorentzVector kf_tm_recoil(0, 0, 0, 2 * emeas);
  kf_tm_recoil -= hypo->getFinalMomentum(_permTracks);
  in_tprecoil = in_tm_recoil.P();
  in_terecoil = in_tm_recoil.E();
  in_tm2recoil = in_tm_recoil.M2();
  kf_tprecoil = kf_tm_recoil.P();
  kf_terecoil = kf_tm_recoil.E();
  kf_tm2recoil = kf_tm_recoil.M2();
  in_mgg = hypo->getInitialMomentum(_permPhotons).M();
  kf_mgg = hypo->getFinalMomentum(_permPhotons).M();
  TVector3 ivertex;
  TVector3 fvertex;
  TLorentzVector timomentum;
  TLorentzVector tfmomentum;
  TLorentzVector phimomentum;
  TLorentzVector phfmomentum;
  ivertex = hypo->getInitialVertex("vtx0");
  in_vx = ivertex.X();
  in_vy = ivertex.Y();
  in_vz = ivertex.Z();
  fvertex = hypo->getFinalVertex("vtx0");
  kf_vx = fvertex.X();
  kf_vy = fvertex.Y();
  kf_vz = fvertex.Z();
  for (int i = 0; i < 2; ++i) {
    phimomentum = hypo->getInitialMomentum(_seqPhotons[i]);
    in_phth[i] = phimomentum.Theta();
    in_phphi[i] = phimomentum.Phi();
    in_phen[i] = phimomentum.E();
    phfmomentum = hypo->getFinalMomentum(_seqPhotons[i]);
    kf_phth[i] = phfmomentum.Theta();
    kf_phphi[i] = phfmomentum.Phi();
    kf_phen[i] = phfmomentum.E();
  }
  for (int i = 0; i < 4; ++i) {
    in_m3pi[i] = hypo->getInitialMomentum(_perm3Pi[i]).M();
    kf_m3pi[i] = hypo->getFinalMomentum(_perm3Pi[i]).M();
    in_mpiplpimi[i] = hypo->getInitialMomentum(_permPiPlPiMi[i]).M();
    kf_mpiplpimi[i] = hypo->getFinalMomentum(_permPiPlPiMi[i]).M();
    timomentum = hypo->getInitialMomentum(_seqTracks[i]);
    in_tth[i] = timomentum.Theta();
    in_tphi[i] = timomentum.Phi();
    in_tz[i] = tz[tind[i]];
    in_ten[i] = ten[tind[i]];
    in_tetot[i] = timomentum.E();
    in_tptot[i] = timomentum.P();
    in_tcharge[i] = tcharge[tind[i]];
    tfmomentum = hypo->getFinalMomentum(_seqTracks[i]);
    kf_tth[i] = tfmomentum.Theta();
    kf_tphi[i] = tfmomentum.Phi();
    kf_tetot[i] = tfmomentum.E();
    kf_tptot[i] = tfmomentum.P();
  }
  CFourVector cfv_piMi0(hypo->getFinalMomentum("pi-_0"), 1.e-3);
  CFourVector cfv_piMi1(hypo->getFinalMomentum("pi-_1"), 1.e-3);
  CFourVector cfv_piPl0(hypo->getFinalMomentum("pi+_0"), 1.e-3);
  CFourVector cfv_piPl1(hypo->getFinalMomentum("pi+_1"), 1.e-3);
  CFourVector cfv_pi0(hypo->getFinalMomentum(_permPhotons), 1.e-3);
  mel2_omegapipi = MElCalc::getOmega2PiMEl2(std::make_pair(cfv_piMi0, cfv_piMi1),
					    std::make_pair(cfv_piPl0, cfv_piPl1), cfv_pi0);
  std::vector<CFourVector> t_momenta;
  t_momenta.reserve(5);
  t_momenta.push_back(cfv_piMi0);
  t_momenta.push_back(cfv_piMi1);
  t_momenta.push_back(cfv_piPl0);
  t_momenta.push_back(cfv_piPl1);
  t_momenta.push_back(cfv_pi0);
  ph_space_wt = phase_space_weight(t_momenta);

  TLorentzVector in_total_dmomentum = hypo->getInitialMomentum(_permAll);
  in_total_de = in_total_dmomentum.E() - 2 * emeas;
  in_total_px = in_total_dmomentum.X();
  in_total_py = in_total_dmomentum.Y();
  in_total_pz = in_total_dmomentum.Z();

  TLorentzVector kf_total_dmomentum = hypo->getFinalMomentum(_permAll);
  kf_total_de = kf_total_dmomentum.E() - 2 * emeas;
  kf_total_px = kf_total_dmomentum.X();
  kf_total_py = kf_total_dmomentum.Y();
  kf_total_pz = kf_total_dmomentum.Z();
}

bool KinFit5Pi::fit5pi(KFCmd::Hypo4ChPions2Photons* hypo) {
  if (!hypo->fillTrack("pi-_0", _trackIndices[0], *this)) return false;
  if (!hypo->fillTrack("pi-_1", _trackIndices[1], *this)) return false;
  if (!hypo->fillTrack("pi+_0", _trackIndices[2], *this)) return false;
  if (!hypo->fillTrack("pi+_1", _trackIndices[3], *this)) return false;
  bool result = false;
  for (std::size_t iph = 0; iph + 1 < _photonIndices.size(); ++iph) {
    if (!hypo->fillPhoton("g0", _photonIndices[iph], *this)) continue;
    for (std::size_t jph = iph + 1; jph < _photonIndices.size(); ++jph) {
        if (!hypo->fillPhoton("g1", _photonIndices[jph], *this)) continue;
	hypo->optimize();
	if (hypo->getErrorCode() != 0) continue;
	if (hypo->getChiSquare() > kf_chi2_5pi) continue;
	fill_fit5pi(hypo, iph, jph);
	result = true;
    }
  }
  return result;
}

void KinFit5Pi::fill_fit5pi_mpi0(KFCmd::Hypo4ChPions2Photons* hypo) {
  kf_err_5pi_mpi0 = hypo->getErrorCode();
  kf_chi2_5pi_mpi0 = hypo->getChiSquare();
  kf_chi2ph_5pi_mpi0 = hypo->getChiSquare(_permPhotons);
  kf_chi2t_5pi_mpi0 = hypo->getChiSquare(_permTracks);
  for (int i = 0; i < 4; ++i) {
    kf_mpi0_m3pi[i] = hypo->getFinalMomentum(_perm3Pi[i]).M();
    TLorentzVector tfmomentum = hypo->getFinalMomentum(_seqTracks[i]);
    kf_mpi0_tth[i] = tfmomentum.Theta();
    kf_mpi0_tphi[i] = tfmomentum.Phi();
  }
  kf_mpi0_mgg_dbl = hypo->getFinalMomentum(_permPhotons).M();
  kf_mpi0_mgg2_dbl = hypo->getFinalMomentum(_permPhotons).M2();
  kf_mpi0_mgg = hypo->getFinalMomentum(_permPhotons).M();
}

void KinFit5Pi::fit5pi_mpi0(KFCmd::Hypo4ChPions2Photons* hypo) {
  if (!hypo->fillTrack("pi-_0", _trackIndices[0], *this)) return;
  if (!hypo->fillTrack("pi-_1", _trackIndices[1], *this)) return;
  if (!hypo->fillTrack("pi+_0", _trackIndices[2], *this)) return;
  if (!hypo->fillTrack("pi+_1", _trackIndices[3], *this)) return;
  if (!hypo->fillPhoton("g0", phind[0], *this)) return;
  if (!hypo->fillPhoton("g1", phind[1], *this)) return;
  hypo->optimize();
  if (hypo->getErrorCode() != 0) return;
  fill_fit5pi_mpi0(hypo);
}

void KinFit5Pi::fill_fit4pi(KFCmd::Hypo4ChPions* hypo) {
  kf_err_4pi = hypo->getErrorCode();
  kf_chi2_4pi = hypo->getChiSquare();
}

void KinFit5Pi::fit4pi(KFCmd::Hypo4ChPions* hypo) {
  if (!hypo->fillTrack("pi-_0", _trackIndices[0], *this)) return;
  if (!hypo->fillTrack("pi-_1", _trackIndices[1], *this)) return;
  if (!hypo->fillTrack("pi+_0", _trackIndices[2], *this)) return;
  if (!hypo->fillTrack("pi+_1", _trackIndices[3], *this)) return;
  hypo->optimize();
  if (hypo->getErrorCode() != 0) return;
  fill_fit4pi(hypo);
}

void KinFit5Pi::Loop(const std::string& outpath, double mfield) {
  _time.Start();
  setStatus(0);
  if (fChain == 0) return;
  TFile* outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutptuBranches(out_tree);
  fChain->GetEntry(0);
  KFCmd::Hypo4ChPions2Photons hypo5pi(2 * emeas, mfield);
  KFCmd::Hypo4ChPions2Photons hypo5pi_mpi0(2 * emeas, mfield);
  hypo5pi_mpi0.enableConstraint("m-pi0-constraint");
  KFCmd::Hypo4ChPions hypo4pi(2 * emeas, mfield);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  _timePerEntry.Start();
  int npassed = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (!cut()) continue;
    clean();
    hypo5pi.setBeamXY(xbeam, ybeam);
    hypo5pi.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
    hypo5pi.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);

    // hypo5pi.disableVertexConstraintX("pi-_0");
    // hypo5pi.disableVertexConstraintX("pi-_1");
    // hypo5pi.disableVertexConstraintX("pi+_0");
    // hypo5pi.disableVertexConstraintX("pi+_1");

    // hypo5pi.disableVertexConstraintY("pi-_0");
    // hypo5pi.disableVertexConstraintY("pi-_1");
    // hypo5pi.disableVertexConstraintY("pi+_0");
    // hypo5pi.disableVertexConstraintY("pi+_1");
    
    
    // hypo5pi.fixVertexComponent("vtx0", z0, KFBase::VERTEX_Z);
    // hypo5pi.disableVertexConstraintXYZ("pi-_0");
    // hypo5pi.disableVertexConstraintXYZ("pi-_1");
    // hypo5pi.disableVertexConstraintXYZ("pi+_0");
    // hypo5pi.disableVertexConstraintXYZ("pi+_1");
        
    hypo5pi_mpi0.setBeamXY(xbeam, ybeam);
    hypo5pi_mpi0.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
    hypo5pi_mpi0.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);

    // hypo5pi_mpi0.disableVertexConstraintX("pi-_0");
    // hypo5pi_mpi0.disableVertexConstraintX("pi-_1");
    // hypo5pi_mpi0.disableVertexConstraintX("pi+_0");
    // hypo5pi_mpi0.disableVertexConstraintX("pi+_1");

    // hypo5pi_mpi0.disableVertexConstraintY("pi-_0");
    // hypo5pi_mpi0.disableVertexConstraintY("pi-_1");
    // hypo5pi_mpi0.disableVertexConstraintY("pi+_0");
    // hypo5pi_mpi0.disableVertexConstraintY("pi+_1");
    
    // hypo5pi_mpi0.fixVertexComponent("vtx0", z0, KFBase::VERTEX_Z);
    // hypo5pi_mpi0.disableVertexConstraintXYZ("pi-_0");
    // hypo5pi_mpi0.disableVertexConstraintXYZ("pi-_1");
    // hypo5pi_mpi0.disableVertexConstraintXYZ("pi+_0");
    // hypo5pi_mpi0.disableVertexConstraintXYZ("pi+_1");

    hypo4pi.setBeamXY(xbeam, ybeam);
    hypo4pi.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
    hypo4pi.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);

    // hypo4pi.disableVertexConstraintX("pi-_0");
    // hypo4pi.disableVertexConstraintX("pi-_1");
    // hypo4pi.disableVertexConstraintX("pi+_0");
    // hypo4pi.disableVertexConstraintX("pi+_1");

    // hypo4pi.disableVertexConstraintY("pi-_0");
    // hypo4pi.disableVertexConstraintY("pi-_1");
    // hypo4pi.disableVertexConstraintY("pi+_0");
    // hypo4pi.disableVertexConstraintY("pi+_1");
    
    // hypo4pi.fixVertexComponent("vtx0", z0, KFBase::VERTEX_Z);
    // hypo4pi.disableVertexConstraintXYZ("pi-_0");
    // hypo4pi.disableVertexConstraintXYZ("pi-_1");
    // hypo4pi.disableVertexConstraintXYZ("pi+_0");
    // hypo4pi.disableVertexConstraintXYZ("pi+_1");
    
    if (!fit5pi(&hypo5pi)) continue;
    npassed++;
    fit5pi_mpi0(&hypo5pi_mpi0);
    fit4pi(&hypo4pi);
    out_tree->Fill();
    printStatus(jentry, nentries);
  }
  printStatus(nentries, nentries);
  out_tree->Write();
  outfl->Close();
  delete outfl;
  printSummary(npassed, nentries);
  
}

void KinFit5Pi::clean() {
  in_total_de = std::numeric_limits<float>::infinity();
  in_total_px = std::numeric_limits<float>::infinity();
  in_total_py = std::numeric_limits<float>::infinity();
  in_total_pz = std::numeric_limits<float>::infinity();
  
  kf_total_de = std::numeric_limits<float>::infinity();
  kf_total_px = std::numeric_limits<float>::infinity();
  kf_total_py = std::numeric_limits<float>::infinity();
  kf_total_pz = std::numeric_limits<float>::infinity();

  ph_space_wt = 0;
  mel2_omegapipi = 0;
  std::fill(tind, tind + 4, 0);
  std::fill(phind, phind + 2, 0);
  in_tprecoil = 0;
  in_terecoil = 0;
  in_tm2recoil = 0;
  kf_tprecoil = 0;
  kf_terecoil = 0;
  kf_tm2recoil = 0;
  kf_chi2_5pi = std::numeric_limits<float>::infinity();
  kf_chi2ph_5pi = std::numeric_limits<float>::infinity();
  kf_chi2t_5pi = std::numeric_limits<float>::infinity();
  kf_err_5pi = 1;
  kf_chi2_5pi_mpi0 = std::numeric_limits<float>::infinity();
  kf_chi2ph_5pi_mpi0 = std::numeric_limits<float>::infinity();
  kf_chi2t_5pi_mpi0 = std::numeric_limits<float>::infinity();
  kf_err_5pi_mpi0 = 1;
  kf_chi2_4pi = std::numeric_limits<float>::infinity();
  kf_err_4pi = 1;
  in_mgg = 0;
  kf_mgg = 0;
  kf_mpi0_mgg = 0;
  kf_mpi0_mgg_dbl = 0;
  kf_mpi0_mgg2_dbl = 0;
  std::fill(in_m3pi, in_m3pi + 4, 0);
  std::fill(kf_m3pi, kf_m3pi + 4, 0);
  std::fill(kf_mpi0_m3pi, kf_mpi0_m3pi + 4, 0);
  std::fill(in_mpiplpimi, in_mpiplpimi + 4, 0);
  std::fill(kf_mpiplpimi, kf_mpiplpimi + 4, 0);
  in_vx = 0;
  in_vy = 0;
  in_vz = 0;
  std::fill(in_tth, in_tth + 4, 0);
  std::fill(in_tphi, in_tphi + 4, 0);
  std::fill(in_tz, in_tz + 4, 0);
  std::fill(in_ten, in_ten + 4, 0);
  std::fill(in_tptot, in_tptot + 4, 0);
  std::fill(in_tetot, in_tetot + 4, 0);
  std::fill(in_tcharge, in_tcharge + 4, 0);
  kf_vx = 0;
  kf_vy = 0;
  kf_vz = 0;
  std::fill(kf_tth, kf_tth + 4, 0);
  std::fill(kf_mpi0_tth, kf_mpi0_tth + 4, 0);
  
  std::fill(kf_tphi, kf_tphi + 4, 0);
  std::fill(kf_mpi0_tphi, kf_mpi0_tphi + 4, 0);
  
  std::fill(kf_tetot, kf_tetot + 4, 0);
  std::fill(kf_tptot, kf_tptot + 4, 0);
  std::fill(in_phth, in_phth + 2, 0);
  std::fill(in_phphi, in_phphi + 2, 0);
  std::fill(in_phen, in_phen + 2, 0);
  std::fill(kf_phth, kf_phth + 2, 0);
  std::fill(kf_phphi, kf_phphi + 2, 0);
  std::fill(kf_phen, kf_phen + 2, 0);
}

void KinFit5Pi::printSummary(int npassed, int nentries) {
  _time.Stop();
  std::cout << "_________________________________________" << std::endl;
  std::cout << "TOTAL NUMBER OF EVENTS :\t" << nentries << std::endl;
  std::cout << "NUMBER OF PASSED EVENTS :\t" << npassed << std::endl;
  std::cout << "TOTAL CPU TIME :\t" << _time.CpuTime() << std::endl;
  std::cout << "TOTAL REAL TIME :\t" << _time.RealTime() << std::endl;
  std::cout << "_________________________________________" << std::endl;
}

void KinFit5Pi::printStatus(int entry_number, int nentries) {
  int status =  (100.0 * entry_number) / nentries;
  if (getStatus() == status) return;
  setStatus(status);
  _timePerEntry.Stop();
  std::cout << " [STATUS : " <<
    std::setfill('0') << std::setw(2) <<
    getStatus()  << "%]" <<
    std::setw(0) <<
    "\tCPU TIME: " <<
    std::setprecision(3) <<
    std::fixed <<
    _timePerEntry.CpuTime() <<
    "\tREAL TIME: " <<
    _timePerEntry.RealTime() << std::endl;
  _timePerEntry.Start();
}

void KinFit5Pi::setStatus(int status) {
  _status = status;
}

int KinFit5Pi::getStatus() const {
  return _status;
}

void KinFit5Pi::setupOutptuBranches(TTree* tree) {
  tree->Branch("in_total_de", &in_total_de, "in_total_de/F");
  tree->Branch("in_total_px", &in_total_px, "in_total_px/F");
  tree->Branch("in_total_py", &in_total_py, "in_total_py/F");
  tree->Branch("in_total_pz", &in_total_pz, "in_total_pz/F");
  
  tree->Branch("kf_total_de", &kf_total_de, "kf_total_de/F");
  tree->Branch("kf_total_px", &kf_total_px, "kf_total_px/F");
  tree->Branch("kf_total_py", &kf_total_py, "kf_total_py/F");
  tree->Branch("kf_total_pz", &kf_total_pz, "kf_total_pz/F");
  
  tree->Branch("in_vx", &in_vx, "in_vx/F");
  tree->Branch("in_vy", &in_vy, "in_vy/F");
  tree->Branch("in_vz", &in_vz, "in_vz/F");
  tree->Branch("in_tth", in_tth, "in_tth[4]/F");
  tree->Branch("in_tphi", in_tphi, "in_tphi[4]/F");
  tree->Branch("in_tz", in_tz, "in_tz[4]/F");
  tree->Branch("in_ten", in_ten, "in_ten[4]/F");
  tree->Branch("in_tptot", in_tptot, "in_tptot[4]/F");
  tree->Branch("in_tetot", in_tetot, "in_tetot[4]/F");
  tree->Branch("in_tcharge", in_tcharge, "in_tcharge[4]/I");
  tree->Branch("kf_vx", &kf_vx, "kf_vx/F");
  tree->Branch("kf_vy", &kf_vy, "kf_vy/F");
  tree->Branch("kf_vz", &kf_vz, "kf_vz/F");
  tree->Branch("kf_tth", kf_tth, "kf_tth[4]/F");
  tree->Branch("kf_mpi0_tth", kf_mpi0_tth, "kf_mpi0_tth[4]/F");
 
  tree->Branch("kf_tetot", kf_tetot, "kf_tetot[4]/F");
  tree->Branch("kf_tphi", kf_tphi, "kf_tphi[4]/F");
  tree->Branch("kf_mpi0_tphi", kf_mpi0_tphi, "kf_mpi0_tphi[4]/F");
  
  tree->Branch("kf_tptot", kf_tptot, "kf_tptot[4]/F");
  tree->Branch("in_phth", in_phth, "in_phth[2]/F");
  tree->Branch("in_phphi", in_phphi, "in_phphi[2]/F");
  tree->Branch("in_phen", in_phen, "in_phen[2]/F");
  tree->Branch("kf_phth", kf_phth, "kf_phth[2]/F");
  tree->Branch("kf_phphi", kf_phphi, "kf_phphi[2]/F");
  tree->Branch("kf_phen", kf_phen, "kf_phen[2]/F");
 
  tree->Branch("ph_space_wt", &ph_space_wt, "ph_space_wt/F");
  tree->Branch("mel2_omegapipi", &mel2_omegapipi, "mel2_omegapipi/F");
  tree->Branch("tind", tind, "tind[4]/I");
  tree->Branch("phind", phind, "phind[2]/I");
  tree->Branch("in_tprecoil", &in_tprecoil, "in_tprecoil/F");
  tree->Branch("in_terecoil", &in_terecoil, "in_terecoil/F");
  tree->Branch("in_tm2recoil", &in_tm2recoil, "in_tm2recoil/F");
  tree->Branch("kf_tprecoil", &kf_tprecoil, "kf_tprecoil/F");
  tree->Branch("kf_terecoil", &kf_terecoil, "kf_terecoil/F");
  tree->Branch("kf_tm2recoil", &kf_tm2recoil, "kf_tm2recoil/F");
  tree->Branch("kf_chi2_5pi", &kf_chi2_5pi, "kf_chi2_5pi/F");
  tree->Branch("kf_chi2ph_5pi", &kf_chi2ph_5pi, "kf_chi2ph_5pi/F");
  tree->Branch("kf_chi2t_5pi", &kf_chi2t_5pi, "kf_chi2t_5pi/F");
  tree->Branch("kf_err_5pi", &kf_err_5pi, "kf_err_5pi/I");
  tree->Branch("kf_chi2_5pi_mpi0", &kf_chi2_5pi_mpi0, "kf_chi2_5pi_mpi0/F");
  tree->Branch("kf_chi2ph_5pi_mpi0", &kf_chi2ph_5pi_mpi0, "kf_chi2ph_5pi_mpi0/F");
  tree->Branch("kf_chi2t_5pi_mpi0", &kf_chi2t_5pi_mpi0, "kf_chi2t_5pi_mpi0/F");
  tree->Branch("kf_err_5pi_mpi0", &kf_err_5pi_mpi0, "kf_err_5pi_mpi0/I");
  tree->Branch("kf_chi2_4pi", &kf_chi2_4pi, "kf_chi2_4pi/F");
  tree->Branch("kf_err_4pi", &kf_err_4pi, "kf_err_4pi/I");
  tree->Branch("in_mgg", &in_mgg, "in_mgg/F");
  tree->Branch("kf_mgg", &kf_mgg, "kf_mgg/F");

  tree->Branch("kf_mpi0_mgg", &kf_mpi0_mgg, "kf_mpi0_mgg/F");
  tree->Branch("kf_mpi0_mgg_dbl", &kf_mpi0_mgg_dbl, "kf_mpi0_mgg_dbl/D");
  tree->Branch("kf_mpi0_mgg2_dbl", &kf_mpi0_mgg2_dbl, "kf_mpi0_mgg2_dbl/D");
  
  tree->Branch("in_m3pi", in_m3pi, "in_m3pi[4]/F");
  tree->Branch("kf_m3pi", kf_m3pi, "kf_m3pi[4]/F");
  tree->Branch("kf_mpi0_m3pi", kf_mpi0_m3pi, "kf_mpi0_m3pi[4]/F");
  tree->Branch("in_mpiplpimi", in_mpiplpimi, "in_mpiplpimi[4]/F");
  tree->Branch("kf_mpiplpimi", kf_mpiplpimi, "kf_mpiplpimi[4]/F");
 
  tree->Branch("ebeam", &ebeam, "ebeam/F");
  tree->Branch("emeas", &emeas, "emeas/F");
  tree->Branch("demeas", &demeas, "demeas/F");
  tree->Branch("emeas0", &emeas0, "emeas0/F");
  tree->Branch("demeas0", &demeas0, "demeas0/F");
  tree->Branch("xbeam", &xbeam, "xbeam/F");
  tree->Branch("ybeam", &ybeam, "ybeam/F");
  tree->Branch("runnum", &runnum, "runnum/I");
  tree->Branch("finalstate_id", &finalstate_id, "finalstate_id/I");
  tree->Branch("evnum", &evnum, "evnum/I");
  tree->Branch("trigbits",  &trigbits,  "trigbits/I");
  tree->Branch("trigmchs",  &trigmchs,  "trigmchs/I");
  tree->Branch("trigtime",  &trigtime,  "trigtime/F");
  tree->Branch("time",  &time,  "time/F");
  tree->Branch("dcfittime", &dcfittime, "dcfittime/F");
  tree->Branch("anttime",   &anttime,   "anttime/F");
  tree->Branch("mutime", &mutime, "mutime/F");
  tree->Branch("is_coll",   &is_coll,   "is_coll/I");
  tree->Branch("is_bhabha", &is_bhabha, "is_bhabha/I");
  tree->Branch("nt_total",  &nt_total,  "nt_total/I");
  tree->Branch("ecaltot",   &ecaltot,   "ecaltot/F");
  tree->Branch("ecalneu",   &ecalneu,   "ecalneu/F");
  tree->Branch("z0", &z0, "z0/F");
  tree->Branch("psumch", &psumch, "psumch/F");
  tree->Branch("psumnu", &psumnu, "psumnu/F");
  tree->Branch("lumoff", &lumoff, "lumoff/F");
  tree->Branch("lumofferr", &lumofferr, "lumofferr/F");
  tree->Branch("nv_total", &nv_total,   "nv_total/I");
  tree->Branch("nv",   &nv,  "nv/I");
  tree->Branch("vtrk",  vtrk,   "vtrk[nv]/I");
  tree->Branch("vind",  vind,   "vind[nv][10]/I");
  tree->Branch("vchi",  vchi,   "vchi[nv]/F");
  tree->Branch("vxyz",  vxyz,   "vxyz[nv][3]/F");
  tree->Branch("nt", &nt, "nt/I");
  tree->Branch("it",  it, "it[2]/I");
  tree->Branch("tnhit",  tnhit, "tnhit[nt]/I");
  tree->Branch("tlength", tlength,  "tlength[nt]/F");
  tree->Branch("tphi",   tphi,  "tphi[nt]/F");
  tree->Branch("tth", tth,   "tth[nt]/F");
  tree->Branch("tptot",  tptot, "tptot[nt]/F");
  tree->Branch("tphiv",  tphiv, "tphiv[nt]/F");
  tree->Branch("tthv",   tthv,  "tthv[nt]/F");
  tree->Branch("tptotv", tptotv, "tptotv[nt]/F");
  tree->Branch("trho",   trho,  "trho[nt]/F");
  tree->Branch("tdedx",  tdedx, "tdedx[nt]/F");
  tree->Branch("tz",  tz, "tz[nt]/F");
  tree->Branch("tt0",    tt0,   "tt0[nt]/F");
  tree->Branch("tant",   tant,  "tant[nt]/F");
  tree->Branch("tchi2r", tchi2r,    "tchi2r[nt]/F");
  tree->Branch("tchi2z", tchi2z,    "tchi2z[nt]/F");
  tree->Branch("tchi2ndf",   tchi2ndf,  "tchi2ndf[nt]/F");
  tree->Branch("tcharge",    tcharge,   "tcharge[nt]/I");
  tree->Branch("ten",    ten,   "ten[nt]/F");
  tree->Branch("tfc",    tfc,   "tfc[nt]/F");
  tree->Branch("tenlxe", tenlxe,    "tenlxe[nt]/F");
  tree->Branch("tlengthlxe",tlengthlxe,"tlengthlxe[nt]/F");
  tree->Branch("tenslxe_layers",tenslxe_layers,  "tenslxe_layers[nt][14]/F");
  tree->Branch("tencsi", tencsi,    "tencsi[nt]/F");
  tree->Branch("tenbgo", tenbgo,    "tenbgo[nt]/F");
  tree->Branch("tclth",  tclth, "tclth[nt]/F");
  tree->Branch("tclphi", tclphi,    "tclphi[nt]/F");
  tree->Branch("terr",   terr,  "terr[nt][3][3]/F");
  tree->Branch("terr0",  terr0, "terr0[nt][5][5]/F");
  tree->Branch("tindlxe",    tindlxe,   "tindlxe[nt]/I");
  tree->Branch("tzcc",   tzcc,  "tzcc[nt][2]/F");
  tree->Branch("txyzatcl",   txyzatcl,  "txyzatcl[nt][3]/F");
  tree->Branch("txyzatlxe",  txyzatlxe, "txyzatlxe[nt][3]/F");
  tree->Branch("tenconv",    tenconv,   "tenconv[nt]/I");
  tree->Branch("nks_total", &nks_total, "nks_total/I");
  tree->Branch("nks",   &nks,   "nks/I");
  tree->Branch("ksvind", ksvind,    "ksvind[nks][2]/I");
  tree->Branch("kstype", kstype,    "kstype[nks]/I");
  tree->Branch("ksfstatus", ksfstatus,    "ksfstatus[nks]/I");
  tree->Branch("ksvchi", ksvchi,    "ksvchi[nks]/F");
  tree->Branch("ksvxyz", ksvxyz,    "ksvxyz[nks][3]/F");
  tree->Branch("ksminv", ksminv,    "ksminv[nks]/F");
  tree->Branch("ksalign",    ksalign,   "ksalign[nks]/F");
  tree->Branch("kstlen", kstlen,    "kstlen[nks]/F");
  tree->Branch("ksdpsi", ksdpsi,    "ksdpsi[nks]/F");
  tree->Branch("kslen",  kslen, "kslen[nks]/F");
  tree->Branch("ksz0",   ksz0,  "ksz0[nks]/F");
  tree->Branch("ksphi",  ksphi, "ksphi[nks]/F");
  tree->Branch("ksth",   ksth,  "ksth[nks]/F");
  tree->Branch("ksptot", ksptot,    "ksptot[nks]/F");
  tree->Branch("kspiphi",    kspiphi,   "kspiphi[nks][2]/F");
  tree->Branch("kspith", kspith,    "kspith[nks][2]/F");
  tree->Branch("kspipt", kspipt,    "kspipt[nks][2]/F");
  tree->Branch("kserr",  kserr, "kserr[nks][3][3]/F");
  // LXe treacks
  tree->Branch("ntlxe_total", &ntlxe_total, "ntlxe_total/I");
  tree->Branch("ntlxe",   &ntlxe,   "ntlxe/I");
  tree->Branch("ntlxelayers",   ntlxelayers,   "ntlxelayers[ntlxe]/I");
  tree->Branch("tlxenhit",    tlxenhit, "tlxenhit[ntlxe]/I");
  tree->Branch("tlxelength",    tlxelength,  "tlxelength[ntlxe]/F");
  tree->Branch("tlxededx",    tlxededx,  "tlxededx[ntlxe]/F");
  tree->Branch("tlxeir",    tlxeir,  "tlxeir[ntlxe]/F");
  tree->Branch("tlxeitheta",    tlxeitheta,  "tlxeitheta[ntlxe]/F");
  tree->Branch("tlxeiphi",    tlxeiphi,  "tlxeiphi[ntlxe]/F");
  tree->Branch("tlxevtheta",    tlxevtheta,  "tlxevtheta[ntlxe]/F");
  tree->Branch("tlxevphi",    tlxevphi,  "tlxevphi[ntlxe]/F");
  tree->Branch("tlxechi2",    tlxechi2,  "tlxechi2[ntlxe]/F");
  tree->Branch("tlxesen",    tlxesen,  "tlxesen[ntlxe]/F");
  tree->Branch("tlxesen_layers",tlxesen_layers,  "tlxesen_layers[ntlxe][14]/F");
  tree->Branch("nph_total",  &nph_total,"nph_total/I");
  tree->Branch("nph",    &nph,  "nph/I");
  tree->Branch("phen",   phen,  "phen[nph]/F");
  tree->Branch("phth",   phth,  "phth[nph]/F");
  tree->Branch("phphi",  phphi, "phphi[nph]/F");
  tree->Branch("phrho",  phrho, "phrho[nph]/F");
  tree->Branch("phen0",   phen0,  "phen0[nph]/F");
  tree->Branch("phth0",   phth0,  "phth0[nph]/F");
  tree->Branch("phphi0",  phphi0, "phphi0[nph]/F");
  tree->Branch("phlxe",  phlxe, "phlxe[nph]/F");
  tree->Branch("phslxe_layers",  phslxe_layers, "phslxe_layers[nph][14]/F");
  tree->Branch("pherr",  pherr, "pherr[nph][3]/F");
  tree->Branch("phcsi",  phcsi, "phcsi[nph]/F");
  tree->Branch("phbgo",  phbgo, "phbgo[nph]/F");
  tree->Branch("phflag", phflag,    "phflag[nph]/I");
  tree->Branch("phconv", phconv,    "phconv[nph]/I");
  tree->Branch("phfc",   phfc,  "phfc[nph]/I");
  tree->Branch("nzcs_total", &nzcs_total, "nzcs_total/I");
  tree->Branch("nzcs",   &nzcs,   "nzcs/I");
  tree->Branch("zcsch",   zcsch,  "zcsch[nzcs]/I");
  tree->Branch("zcsstat", zcsstat,    "zcsstat[nzcs]/I");
  tree->Branch("zcsamp",  zcsamp, "zcsamp[nzcs]/F");
  tree->Branch("zcstime", zcstime,    "zcstime[nzcs]/F");
  tree->Branch("zcsphi",  zcsphi, "zcsphi[nzcs]/F");
  tree->Branch("nzcc_total", &nzcc_total, "nzcc_total/I");
  tree->Branch("nzcc",   &nzcc,   "nzcc/I");
  tree->Branch("zccl",   zccl,    "zccl[nzcc]/I");
  tree->Branch("zccns",  zccns,   "zccns[nzcc]/I");
  tree->Branch("zccamp", zccamp,  "zccamp[nzcc]/F");
  tree->Branch("zcct",   zcct,    "zcct[nzcc]/I");
  tree->Branch("zccz",   zccz,    "zccz[nzcc]/F");
  tree->Branch("zccvalid",   zccvalid,    "zccvalid[nzcc]/I");
  tree->Branch("nant",   &nant,   "nant/I");
  tree->Branch("antch",  antch, "antch[nant]/I");
  tree->Branch("antt0",  antt0, "antt0[nant]/F");
  tree->Branch("antt1",  antt1, "antt1[nant]/F");
  tree->Branch("anta0",  anta0, "anta0[nant]/F");
  tree->Branch("anta1",  anta1, "anta1[nant]/F");
  tree->Branch("antst",  antst, "antst[nant]/I");
  tree->Branch("nmu",   &nmu,   "nmu/I");
  tree->Branch("much",  much, "much[nmu]/I");
  tree->Branch("mut0",  mut0, "mut0[nmu]/F");
  tree->Branch("mut1",  mut1, "mut1[nmu]/F");
  tree->Branch("mut2",  mut2, "mut2[nmu]/F");
  tree->Branch("mut3",  mut3, "mut3[nmu]/F");
  tree->Branch("mua0",  mua0, "mua0[nmu]/F");
  tree->Branch("mua1",  mua1, "mua1[nmu]/F");
  tree->Branch("mua2",  mua2, "mua2[nmu]/F");
  tree->Branch("mua3",  mua3, "mua3[nmu]/F");
  tree->Branch("must",  must, "must[nmu]/I");
  tree->Branch("nsim",&nsim,"nsim/I");
  tree->Branch("simtype",simtype,"simtype[nsim]/I");
  tree->Branch("simorig",simorig,"simorig[nsim]/I");
  tree->Branch("simmom",simmom,"simmom[nsim]/F");
  tree->Branch("simphi",simphi,"simphi[nsim]/F");
  tree->Branch("simtheta",simtheta,"simtheta[nsim]/F");
  tree->Branch("simvtx",simvtx,"simvtx[nsim]/F");
  tree->Branch("simvty",simvty,"simvty[nsim]/F");
  tree->Branch("simvtz",simvtz,"simvtz[nsim]/F");
  tree->Branch("ncorr", &ncorr, "ncorr/I");
  tree->Branch("idcorr", idcorr, "idcorr[ncorr]/I");
  tree->Branch("bitcorr", bitcorr, "bitcorr[ncorr]/I");
  tree->Branch("nbadbank", &nbadbank, "nbadbank/I");
  tree->Branch("nbadbankg", &nbadbankg, "nbadbankg/I");
  tree->Branch("nbadbanks", nbadbanks, "nbadbanks[nbadbankg]/I");
  tree->Branch("nlostbanks", &nlostbanks, "nlostbanks/I");
  tree->Branch("ncorruptedbanks", &ncorruptedbanks, "ncorruptedbanks/I");
}
