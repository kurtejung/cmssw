// -*- C++ -*-
//
//

// class template HadronizerFilter<HAD> provides an EDFilter which uses
// the hadronizer type HAD to read in external partons and hadronize them, 
// and decay the resulting particles, in the CMS framework.

#ifndef gen_HadronizerFilter_h
#define gen_HadronizerFilter_h

#include <memory>

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

// #include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

// LHE Run
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"

// LHE Event
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

namespace edm
{
  template <class HAD, class DEC> class HadronizerFilter : public EDFilter
  {
  public:
    typedef HAD Hadronizer;
    typedef DEC Decayer;

    // The given ParameterSet will be passed to the contained
    // Hadronizer object.
    explicit HadronizerFilter(ParameterSet const& ps);

    virtual ~HadronizerFilter();

    virtual bool filter(Event& e, EventSetup const& es);
    virtual void endJob();
    virtual bool beginRun(Run &, EventSetup const&);
    virtual bool endRun(Run &, EventSetup const&);
    virtual bool beginLuminosityBlock(LuminosityBlock &, EventSetup const&);
    virtual bool endLuminosityBlock(LuminosityBlock &, EventSetup const&);
    virtual void respondToOpenInputFile(FileBlock const& fb);
    virtual void respondToCloseInputFile(FileBlock const& fb);
    virtual void respondToOpenOutputFiles(FileBlock const& fb);
    virtual void respondToCloseOutputFiles(FileBlock const& fb);

  private:
    Hadronizer hadronizer_;
    // gen::ExternalDecayDriver* decayer_;
    Decayer* decayer_;
    bool fromSource_;

  };

  //------------------------------------------------------------------------
  //
  // Implementation

  template <class HAD, class DEC>
  HadronizerFilter<HAD,DEC>::HadronizerFilter(ParameterSet const& ps) :
    EDFilter(),
    hadronizer_(ps),
    decayer_(0),
    fromSource_(true)
  {
    // TODO:
    // Put the list of types produced by the filters here.
    // The current design calls for:
    //   * LHEGeneratorInfo
    //   * LHEEvent
    //   * HepMCProduct
    // But I can not find the LHEGeneratorInfo class; it might need to
    // be invented.

    if ( ps.exists("ExternalDecays") )
    {
       //decayer_ = new gen::ExternalDecayDriver(ps.getParameter<ParameterSet>("ExternalDecays"));
       ParameterSet ps1 = ps.getParameter<ParameterSet>("ExternalDecays");
       decayer_ = new Decayer(ps1);
    }

    produces<edm::HepMCProduct>();
    produces<GenEventInfoProduct>();
    produces<GenLumiInfoProduct, edm::InLumi>();
    produces<GenRunInfoProduct, edm::InRun>();
  }

  template <class HAD, class DEC>
  HadronizerFilter<HAD,DEC>::~HadronizerFilter()
  { if (decayer_) delete decayer_; }

  template <class HAD, class DEC>
  bool
  HadronizerFilter<HAD, DEC>::filter(Event& ev, EventSetup const& /* es */)
  {
    hadronizer_.setEDMEvent(ev);

    // get LHE stuff and pass to hadronizer!
    //
    edm::Handle<LHEEventProduct> product;
    if ( fromSource_ ) 
      ev.getByLabel("source", product);
    else
      ev.getByLabel("externalLHEProducer", product);

    lhef::LHEEvent *lheEvent =
		new lhef::LHEEvent(hadronizer_.getLHERunInfo(), *product);
    hadronizer_.setLHEEvent( lheEvent );
    
    // hadronizer_.generatePartons();
    if ( !hadronizer_.hadronize() ) return false ;

    //  this is "fake" stuff
    // in principle, decays are done as part of full event generation,
    // except for particles that are marked as to be kept stable
    // but we currently keep in it the design, because we might want
    // to use such feature for other applications
    //
    if ( !hadronizer_.decay() ) return false;
    
    std::auto_ptr<HepMC::GenEvent> event (hadronizer_.getGenEvent());
    if( !event.get() ) return false; 

    // The external decay driver is being added to the system,
    // it should be called here
    //
    if ( decayer_ ) 
    {
      event.reset( decayer_->decay( event.get() ) );
    }

    if ( !event.get() ) return false;

    // check and perform if there're any unstable particles after 
    // running external decay packges
    //
    hadronizer_.resetEvent( event.release() );
    if ( !hadronizer_.residualDecay() ) return false;

    hadronizer_.finalizeEvent();

    event.reset( hadronizer_.getGenEvent() );
    if ( !event.get() ) return false;

    event->set_event_number( ev.id().event() );

    std::auto_ptr<GenEventInfoProduct> genEventInfo(hadronizer_.getGenEventInfo());
    if (!genEventInfo.get())
    { 
      // create GenEventInfoProduct from HepMC event in case hadronizer didn't provide one
      genEventInfo.reset(new GenEventInfoProduct(event.get()));
    }
    ev.put(genEventInfo);

    std::auto_ptr<HepMCProduct> bare_product(new HepMCProduct());
    bare_product->addHepMCData( event.release() );
    ev.put(bare_product);

    return true;
  }

  template <class HAD, class DEC>
  void
  HadronizerFilter<HAD,DEC>::endJob()
  { }

  template <class HAD, class DEC>
  bool
  HadronizerFilter<HAD,DEC>::beginRun(Run& run, EventSetup const& es)
  {
        
    // this is run-specific
    
    // get LHE stuff and pass to hadronizer!

    std::vector< edm::Handle<LHERunInfoProduct> > productV;

    run.getManyByType(productV); 

    if ( productV.size() > 1 ) 
      throw edm::Exception(errors::EventCorruption) 
        << "More than one LHERunInfoProduct present";
    else
      {
        if ( productV[0].provenance()->moduleLabel() == "externalLHEProducer" ) 
          fromSource_ = false;
      }

    hadronizer_.setLHERunInfo( new lhef::LHERunInfo(*productV[0]) ) ;
    lhef::LHERunInfo* lheRunInfo = hadronizer_.getLHERunInfo().get();
    lheRunInfo->initLumi();
    
    return true;
  
  }

  template <class HAD, class DEC>
  bool
  HadronizerFilter<HAD,DEC>::endRun(Run& r, EventSetup const&)
  {
    // Retrieve the LHE run info summary and transfer determined
    // cross-section into the generator run info

    const lhef::LHERunInfo* lheRunInfo = hadronizer_.getLHERunInfo().get();
    lhef::LHERunInfo::XSec xsec = lheRunInfo->xsec();

    GenRunInfoProduct& genRunInfo = hadronizer_.getGenRunInfo();
    genRunInfo.setInternalXSec( GenRunInfoProduct::XSec(xsec.value(), xsec.error()) );

    // If relevant, record the integrated luminosity for this run
    // here.  To do so, we would need a standard function to invoke on
    // the contained hadronizer that would report the integrated
    // luminosity.

    hadronizer_.statistics();
    if ( decayer_ ) decayer_->statistics();
    lheRunInfo->statistics();

    std::auto_ptr<GenRunInfoProduct> griproduct( new GenRunInfoProduct(genRunInfo) );
    r.put(griproduct);
    
    return true;
  }

  template <class HAD, class DEC>
  bool
  HadronizerFilter<HAD,DEC>::beginLuminosityBlock(LuminosityBlock &, EventSetup const& es)
  {
    
    lhef::LHERunInfo* lheRunInfo = hadronizer_.getLHERunInfo().get();
    lheRunInfo->initLumi();
    
    if ( !hadronizer_.readSettings(1) )
       throw edm::Exception(errors::Configuration) 
	 << "Failed to read settings for the hadronizer "
	 << hadronizer_.classname() << " \n";

    if ( decayer_ )
    {
       decayer_->init(es);
       if ( !hadronizer_.declareStableParticles( decayer_->operatesOnParticles() ) )
          throw edm::Exception(errors::Configuration)
            << "Failed to declare stable particles in hadronizer "
            << hadronizer_.classname()
	    << " for internal parton generation\n";
       if ( !hadronizer_.declareSpecialSettings( decayer_->specialSettings() ) )
          throw edm::Exception(errors::Configuration)
            << "Failed to declare special settings in hadronizer "
            << hadronizer_.classname()
            << "\n";
    }

    if (! hadronizer_.initializeForExternalPartons())
      throw edm::Exception(errors::Configuration) 
	<< "Failed to initialize hadronizer "
	<< hadronizer_.classname()
	<< " for external parton generation\n";

    return true;

  }

  template <class HAD, class DEC>
  bool
  HadronizerFilter<HAD,DEC>::endLuminosityBlock(LuminosityBlock &lumi, EventSetup const&)
  {

    const lhef::LHERunInfo* lheRunInfo = hadronizer_.getLHERunInfo().get();


    std::vector<lhef::LHERunInfo::Process> LHELumiProcess = lheRunInfo->getLumiProcesses();
    std::vector<GenLumiInfoProduct::ProcessInfo> GenLumiProcess;
    for(unsigned int i=0; i < LHELumiProcess.size(); i++){
      lhef::LHERunInfo::Process thisProcess=LHELumiProcess[i];

      GenLumiInfoProduct::ProcessInfo temp;      
      temp.setProcess(thisProcess.process());
      temp.setLheXSec(thisProcess.getLHEXSec().value(),thisProcess.getLHEXSec().error());
      temp.setNPassPos(thisProcess.nPassPos());
      temp.setNPassNeg(thisProcess.nPassNeg());
      temp.setNTotalPos(thisProcess.nTotalPos());
      temp.setNTotalNeg(thisProcess.nTotalNeg());
      temp.setTried(thisProcess.tried().n(), thisProcess.tried().sum(), thisProcess.tried().sum2());
      temp.setSelected(thisProcess.selected().n(), thisProcess.selected().sum(), thisProcess.selected().sum2());
      temp.setKilled(thisProcess.killed().n(), thisProcess.killed().sum(), thisProcess.killed().sum2());
      temp.setAccepted(thisProcess.accepted().n(), thisProcess.accepted().sum(), thisProcess.accepted().sum2());
      temp.setAcceptedBr(thisProcess.acceptedBr().n(), thisProcess.acceptedBr().sum(), thisProcess.acceptedBr().sum2());
      GenLumiProcess.push_back(temp);
    }
    std::auto_ptr<GenLumiInfoProduct> genLumiInfo(new GenLumiInfoProduct());
    genLumiInfo->setHEPIDWTUP(lheRunInfo->getHEPRUP()->IDWTUP);
    genLumiInfo->setProcessInfo( GenLumiProcess );
    lumi.put(genLumiInfo);    
    
    return true;
  }

  template <class HAD, class DEC>
  void
  HadronizerFilter<HAD,DEC>::respondToOpenInputFile(FileBlock const& fb)
  { }

  template <class HAD, class DEC>
  void
  HadronizerFilter<HAD,DEC>::respondToCloseInputFile(FileBlock const& fb)
  { }

  template <class HAD, class DEC>
  void
  HadronizerFilter<HAD,DEC>::respondToOpenOutputFiles(FileBlock const& fb)
  { }

  template <class HAD, class DEC>
  void
  HadronizerFilter<HAD,DEC>::respondToCloseOutputFiles(FileBlock const& fb)
  { }

}

#endif // gen_HadronizerFilter_h
