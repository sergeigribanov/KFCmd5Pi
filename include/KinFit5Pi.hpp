#ifndef __KINFIT5PI_HPP__
#define __KINFIT5PI_HPP__
#include <string>
#include <vector>
#include <set>
#include <TStopwatch.h>
#include <KFCmd/TrPh.hpp>
#include <KFCmd/Hypo4ChPions2Photons.hpp>
#include <KFCmd/Hypo4ChPions.hpp>

class KinFit5Pi : public KFCmd::TrPh {
public:
  explicit KinFit5Pi(TTree* = 0);
  virtual ~KinFit5Pi();
  virtual void Loop(const std::string&, double) override final;
private:
  bool cutTracks();
  bool cutPhotons();
  void setupOutptuBranches(TTree*);
  bool cut();
  void clean();
  bool fit5pi(KFCmd::Hypo4ChPions2Photons*);
  void fit5pi_mpi0(KFCmd::Hypo4ChPions2Photons*);
  void fit4pi(KFCmd::Hypo4ChPions*);
  void fill_fit5pi(KFCmd::Hypo4ChPions2Photons*, std::size_t, std::size_t);
  void fill_fit5pi_mpi0(KFCmd::Hypo4ChPions2Photons*);
  void fill_fit4pi(KFCmd::Hypo4ChPions*);
  void printSummary(int, int);
  void printStatus(int, int);
  void setStatus(int);
  int getStatus() const;
  static const std::vector<std::string> _seqTracks;
  static const std::vector<std::string> _seqPhotons;
  static const std::set<std::string> _permPhotons;
  static const std::set<std::string> _permTracks;
  static const std::set<std::string> _permAll;
  static const std::vector<std::set<std::string>> _perm3Pi;
  static const std::vector<std::set<std::string>> _permPiPlPiMi;
  int _status;
  TStopwatch _time;
  TStopwatch _timePerEntry;
  std::vector<std::size_t> _trackIndices;
  std::vector<std::size_t> _photonIndices;
  float ph_space_wt;
  float mel2_omegapipi;
  int tind[4];
  int phind[2];
  float in_tprecoil;
  float in_terecoil;
  float in_tm2recoil;
  float kf_tprecoil;
  float kf_terecoil;
  float kf_tm2recoil;
  float kf_chi2_5pi;
  float kf_chi2ph_5pi;
  float kf_chi2t_5pi;
  int kf_err_5pi;
  float kf_chi2_5pi_mpi0;
  float kf_chi2ph_5pi_mpi0;
  float kf_chi2t_5pi_mpi0;
  int kf_err_5pi_mpi0;
  float kf_chi2_4pi;
  int kf_err_4pi;
  float in_mgg;
  float kf_mgg;
  float in_m3pi[4];
  float kf_m3pi[4];
  float kf_mpi0_m3pi[4];
  float in_mpiplpimi[4];
  float kf_mpiplpimi[4];

  float in_vx;
  float in_vy;
  float in_vz;
  float kf_vx;
  float kf_vy;
  float kf_vz;
  
  float in_tth[4];
  float in_tphi[4];
  float in_tz[4];
  float in_ten[4];
  float in_tptot[4];
  float in_tetot[4];
  int in_tcharge[4];
  float kf_tth[4];
  float kf_tphi[4];
  float kf_tetot[4];
  float kf_tptot[4];

  float kf_mpi0_tth[4];
  float kf_mpi0_tphi[4];

  float in_phth[2];
  float in_phphi[2];
  float in_phen[2];
  float kf_phth[2];
  float kf_phphi[2];
  float kf_phen[2];

  float in_total_de;
  float in_total_px;
  float in_total_py;
  float in_total_pz;

  float kf_total_de;
  float kf_total_px;
  float kf_total_py;
  float kf_total_pz;

  float kf_mpi0_mgg;
  double kf_mpi0_mgg_dbl;
  double kf_mpi0_mgg2_dbl;

};

#endif
