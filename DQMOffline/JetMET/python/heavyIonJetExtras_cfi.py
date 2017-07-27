from RecoHI.HiCentralityAlgos.CentralityBin_cfi import *
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
from RecoHI.HiJetAlgos.hiFJRhoProducer import hiFJRhoProducer
from RecoHI.HiJetAlgos.hiFJGridEmptyAreaCalculator_cff import hiFJGridEmptyAreaCalculator

centralityBin.Centrality = cms.InputTag("hiCentrality")
centralityBin.centralityVariable = cms.string("HFtowers")
centralityBin.nonDefaultGlauberModel = cms.string("Hydjet_Drum")

kt4PFJets.src = cms.InputTag('particleFlowTmp')
kt4PFJets.doAreaFastjet = True
kt4PFJets.jetPtMin      = cms.double(0.0)
kt4PFJets.GhostArea     = cms.double(0.005)

from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

akCs4PFJets = cms.EDProducer(
    "CSJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm  = cms.string("AntiKt"),
    rParam        = cms.double(0.4),
    etaMap    = cms.InputTag('hiFJRhoProducer','mapEtaEdges'),
    rho       = cms.InputTag('hiFJRhoProducer','mapToRho'),
    rhom      = cms.InputTag('hiFJRhoProducer','mapToRhoM'),
    csAlpha   = cms.double(1.),
    writeCompound = cms.bool(True),
    jetCollInstanceName = cms.string("pfParticlesCs")
)
akCs4PFJets.src           = cms.InputTag('particleFlowTmp')
akCs4PFJets.doAreaFastjet = cms.bool(True)
akCs4PFJets.jetPtMin      = cms.double(0.0)
akCs4PFJets.useExplicitGhosts = cms.bool(True)
akCs4PFJets.GhostArea     = cms.double(0.005)

akCs3PFJets = akCs4PFJets.clone(rParam=cms.double(0.3))
akCs5PFJets = akCs4PFJets.clone(rParam=cms.double(0.5))

hiFJRhoProducer.etaRanges = cms.vdouble(-5.0, -3.0, -2.1, -1.3, 1.3, 2.1, 3.0, 5.0)
hiFJRhoProducer.jetSource = cms.InputTag("kt4PFJets")

hiFJGridEmptyAreaCalculator.jetSource = cms.InputTag("kt4PFJets")

jetDQMAnalyzerPrequel = cms.Sequence(	centralityBin
					* kt4PFJets
					* hiFJRhoProducer
                                        * hiFJGridEmptyAreaCalculator
					* akCs3PFJets * akCs4PFJets * akCs5PFJets
)
