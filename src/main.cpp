#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <TFile.h>
#include <TTree.h>
#include "KinFit5Pi.hpp"
namespace po = boost::program_options;

typedef struct {
  std::string ifname;
  std::string ofname;
  double mfield;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()(
      "help",
      "CMD-3 kinematic fits")(
      "ifname", po::value<std::string>(&(opts->ifname)),
      "The path to input .root file.")
    ("ofname", po::value<std::string>(&(opts->ofname)),
     "The path to output .root file.")
    ("mfield", po::value<double>(&(opts->mfield)),
     "Magnetic field");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

template <class T>
T* find_object(TFile* fl, const std::string& object_name) {
  T* object = dynamic_cast<T*>(fl->Get(object_name.c_str()));
  if (!object) {
    std::cerr << "[!] Object with name \"" << object_name << "\" of class "
              << T::Class_Name() << " is not found in file " << fl->GetName()
              << std::endl;
    exit(1);
  }
  return object;
}

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options:");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
   if (vmap.count("ifname") &&
      vmap.count("ofname") &&
      vmap.count("mfield")) {
     auto fl = TFile::Open(opts.ifname.c_str(), "read");
     auto tree = find_object<TTree>(fl, "tr_ph");
     KinFit5Pi fitter(tree);
     fitter.Loop(opts.ofname, opts.mfield);
     fl->Close();
   } else {
     help(desc);
   }
  return 0;
}
