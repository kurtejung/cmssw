import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.jets.HiReRecoJets_HI_cff import PFTowers, kt4PFJets, hiFJRhoProducer, hiFJGridEmptyAreaCalculator, akPu2CaloJets, akPu2PFJets, akCs2PFJets, akPu3CaloJets, akPu3PFJets, akCs3PFJets, akPu4CaloJets, akPu4PFJets, akCs4PFJets, akCsSoftDrop4PFJets, akCsSoftDropZ05B154PFJets, ak4PFJets

#jet analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.akPu2CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu2PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs2PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs3PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs4PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_PbPb_mc_cff import *

from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDropZ05B154PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCsFilter4PFJetSequence_PbPb_mc_cff import *
#from HeavyIonsAnalysis.JetAnalysis.jets.akCsFilter5PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDrop4PFJetSequence_PbPb_mc_cff import *
#from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDrop5PFJetSequence_PbPb_mc_cff import *

highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("hiGeneralTracks"),
                                cut = cms.string('quality("highPurity")'))

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
offlinePrimaryVertices.TrackLabel = 'highPurityTracks'

#the following lines are in the wrong python config
#they should be in a python config handling reco, not analyzers. To be fixed
akPu3PFJets.minimumTowersFraction = cms.double(0.5)
akPu4PFJets.minimumTowersFraction = cms.double(0.5)

akPu3CaloJets.minimumTowersFraction = cms.double(0.)
akPu4CaloJets.minimumTowersFraction = cms.double(0.)


jetSequences = cms.Sequence(
    PFTowers + 
    kt4PFJets +
    hiFJRhoProducer +
    hiFJGridEmptyAreaCalculator + 

#    akPu4PFJetsNoLimits +

    #akPu2CaloJets +
    #akPu2PFJets +
    #akCs2PFJets +

#    akPu3CaloJets +
#    akPu3PFJets +
    #akVs3CaloJets +
    #akVs3PFJets +
#    akCs3PFJets +

    ak4PFJets +
    akPu4CaloJets +
    akPu4PFJets +
    akCs4PFJets +

    #akCsSoftDropZ05B154PFJets +
    
    #akPu5CaloJets +
    #akPu5PFJets +
    #akVs5CaloJets +
    #akVs5PFJets +
    #akCs5PFJets +

    #akCsFilter4PFJets +
    #akCsFilter5PFJets +
    akCsSoftDrop4PFJets +
    #akCsSoftDrop5PFJets +

    highPurityTracks +
    offlinePrimaryVertices +

    #akPu2CaloJetSequence +
    #akPu2PFJetSequence +
    #akCs2PFJetSequence +

    #akCsSoftDropZ05B154PFJetSequence
    #akPu3CaloJetSequence +
    #akVs3CaloJetSequence +
    #akVs3PFJetSequence +
    #akPu3PFJetSequence +
    #akCs3PFJetSequence +

    #akPu4CaloJetSequence +
    #akVs4CaloJetSequence +
    #akVs4PFJetSequence +
    #akPu4PFJetSequence +
    akCs4PFJetSequence +

    #akPu5CaloJetSequence +
    #akVs5CaloJetSequence +
    #akVs5PFJetSequence +
    #akPu5PFJetSequence +
    #akCs5PFJetSequence +

    #akCsFilter4PFJetSequence +
    #akCsFilter5PFJetSequence +
    akCsSoftDrop4PFJetSequence #+
    #akCsSoftDrop5PFJetSequence
)
