#ifndef DPGAnalysis_SiStripTools_APVShotFinder_h
#define DPGAnalysis_SiStripTools_APVShotFinder_h

#include <vector>
#include "DPGAnalysis/SiStripTools/interface/APVShot.h"

class SiStripDigi;


namespace edm {
  template <class T> class DetSet;
  template <class T> class DetSetVector;
}

class APVShotFinder {

 public:
  APVShotFinder(const bool zs=true);
  APVShotFinder(const edm::DetSet<SiStripDigi>& digis, const bool zs=true);
  APVShotFinder(const edm::DetSetVector<SiStripDigi>& digicoll, const bool zs=true);

  void computeShots(const edm::DetSet<SiStripDigi>& digis);
  void computeShots(const edm::DetSetVector<SiStripDigi>& digicoll);

  const std::vector<APVShot>& getShots() const;

 private:

  void addShots(const edm::DetSet<SiStripDigi>& digis);

  bool _zs;
  std::vector<APVShot> _shots;

};

#endif // DPGAnalysis_SiStripTools_APVShotFinder_h
