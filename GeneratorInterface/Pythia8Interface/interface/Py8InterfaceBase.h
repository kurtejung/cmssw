#ifndef gen_Py8InterfaceBase_h
#define gen_Py8InterfaceBase_h

#include <vector>
#include <string>

#include "GeneratorInterface/Core/interface/ParameterCollector.h"

#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC2.h>

namespace gen {

   class Py8InterfaceBase {

      public:
         
	 Py8InterfaceBase( edm::ParameterSet const& ps );
	 ~Py8InterfaceBase() {}
	 
         virtual bool generatePartonsAndHadronize() = 0;
         bool decay() { return true; } // NOT used - let's call it "design imperfection"
         bool readSettings( int ); // common func
         virtual bool initializeForInternalPartons() = 0;
         bool declareStableParticles( const std::vector<int>& ); // common func
         bool declareSpecialSettings( const std::vector<std::string>& ); // common func
         virtual void finalizeEvent() = 0; 
         virtual void statistics();
         virtual const char* classname() const = 0;
	       
      protected:
         
	 std::auto_ptr<Pythia8::Pythia> fMasterGen;
	 std::auto_ptr<Pythia8::Pythia> fDecayer;
	 HepMC::Pythia8ToHepMC          toHepMC;
	 ParameterCollector	        fParameters;
	 
	 unsigned int                   pythiaPylistVerbosity;
         bool                           pythiaHepMCVerbosity;
	 unsigned int                   maxEventsToPrint;

   };

}

#endif
