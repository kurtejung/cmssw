import FWCore.ParameterSet.Config as cms

# make heavyIonPatCandidates
from PhysicsTools.PatAlgos.producersHeavyIons.heavyIonPatCandidates_cff import *

# make selectedLayer1Objects
from PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff import *

#get the specific jet containers for HI
from PhysicsTools.PatAlgos.selectionLayer1.hiJetSelector_cff import *

patHeavyIonDefaultSequence = cms.Sequence(
    heavyIonPatCandidates  * 
    selectedhiPatJets        *
    selectedPatMuons       *
    selectedPatPhotons
)
