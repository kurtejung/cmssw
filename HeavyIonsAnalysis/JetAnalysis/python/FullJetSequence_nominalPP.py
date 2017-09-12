import FWCore.ParameterSet.Config as cms

### PP RECO does not include R=3 or R=5 jets.
### re-RECO is only possible for PF, RECO is missing calotowers
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
ak5PFJets.doAreaFastjet = True
ak3PFJets = ak5PFJets.clone(rParam = 0.3)
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
ak3GenJets = ak5GenJets.clone(rParam = 0.3)
ak4GenJets = ak5GenJets.clone(rParam = 0.4)

#SoftDrop PF jets
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
akSoftDrop4PFJets = cms.EDProducer(
    "SoftDropJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(0.4),
    useOnlyCharged = cms.bool(False),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
<<<<<<< HEAD
=======
akSoftDrop5PFJets = akSoftDrop4PFJets.clone(rParam = cms.double(0.5), R0 = cms.double(0.5))
akSoftDrop4PFz01bm1Jets = akSoftDrop4PFJets.clone(beta = cms.double(-1))
akSoftDrop4PFz01b1Jets = akSoftDrop4PFJets.clone(beta = cms.double(1))
akSoftDrop4PFz005bm1Jets = akSoftDrop4PFJets.clone(beta = cms.double(-1), zcut=cms.double(0.05))
akSoftDrop4PFz005bm2Jets = akSoftDrop4PFJets.clone(beta = cms.double(-2), zcut=cms.double(0.05))
>>>>>>> fully-working PbPb configs and included soft-drop variations in pp

from HeavyIonsAnalysis.JetAnalysis.akSoftDrop4GenJets_cfi import akSoftDrop4GenJets

#Filter PF jets
akFilter4PFJets = cms.EDProducer(
    "FastjetJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    useFiltering = cms.bool(True),
    nFilt = cms.int32(4),
    rFilt = cms.double(0.15),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

from RecoJets.Configuration.GenJetParticles_cff import *
from RecoHI.HiJetAlgos.HiGenJets_cff import *
from HeavyIonsAnalysis.JetAnalysis.makePartons_cff import myPartons

from HeavyIonsAnalysis.JetAnalysis.jets.ak3PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4CaloJetSequence_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_pp_mc_cff import *
<<<<<<< HEAD
=======
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_z01_bm1_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_z01_b1_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_z005_bm1_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_z005_bm2_pp_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop5PFJetSequence_pp_mc_cff import *
>>>>>>> fully-working PbPb configs and included soft-drop variations in pp

highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("generalTracks"),
                                cut = cms.string('quality("highPurity")')
)

from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
ak3GenNjettiness = Njettiness.clone(
                    src = cms.InputTag("ak3GenJets"),
                    R0  = cms.double( 0.3 )
)

ak4GenNjettiness = Njettiness.clone(
                    src = cms.InputTag("ak4GenJets"),
                    R0  = cms.double( 0.4 )
)


# Other radii jets and calo jets need to be reconstructed
jetSequences = cms.Sequence(
    myPartons +
    genParticlesForJets +
    #ak3GenJets +
    ak4GenJets +
    ak3GenNjettiness + 
    ak4GenNjettiness + 
    #ak5GenJets +
    ak3PFJets +
    #ak5PFJets +
    akSoftDrop4PFJets +
    #akSoftDrop5PFJets +
    #akSoftDrop4PFz01bm1Jets +
    #akSoftDrop4PFz01b1Jets +
    #akSoftDrop4PFz005bm1Jets +
    #akSoftDrop4PFz005bm2Jets + 
    #akFilter4PFJets +
    #akFilter5PFJets +
    akSoftDrop4GenJets +
    #akSoftDrop5GenJets +
    highPurityTracks +
#    ak3PFJetSequence +
#    ak4PFJetSequence +
#    ak5PFJetSequence +
    #ak4CaloJetSequence +
    akSoftDrop4PFJetSequence 
    #akSoftDrop4PFz01bm1JetSequence +
    #akSoftDrop4PFz01b1JetSequence +
    #akSoftDrop4PFz005bm1JetSequence +
    #akSoftDrop4PFz005bm2JetSequence
#    akSoftDrop5PFJetSequence
)
