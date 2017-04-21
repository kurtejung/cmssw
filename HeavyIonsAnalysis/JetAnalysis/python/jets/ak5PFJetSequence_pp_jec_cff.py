

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

ak5PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("ak5PFJets"),
    matched = cms.InputTag("ak5GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.5
    )

ak5PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("ak5GenJets"),
    matched = cms.InputTag("ak5GenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.5
    )

ak5PFparton = patJetPartonMatch.clone(src = cms.InputTag("ak5PFJets")
                                                        )

ak5PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("ak5PFJets"),
    payload = "AK5PF_offline"
    )

ak5PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('ak5CaloJets'))

#ak5PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5GenJets'))

ak5PFbTagger = bTaggers("ak5PF",0.5)

#create objects locally since they dont load properly otherwise
#ak5PFmatch = ak5PFbTagger.match
ak5PFparton = patJetPartonMatch.clone(src = cms.InputTag("ak5PFJets"), matched = cms.InputTag("genParticles"))
ak5PFPatJetFlavourAssociationLegacy = ak5PFbTagger.PatJetFlavourAssociationLegacy
ak5PFPatJetPartons = ak5PFbTagger.PatJetPartons
ak5PFJetTracksAssociatorAtVertex = ak5PFbTagger.JetTracksAssociatorAtVertex
#ak5PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
ak5PFSimpleSecondaryVertexHighEffBJetTags = ak5PFbTagger.SimpleSecondaryVertexHighEffBJetTags
ak5PFSimpleSecondaryVertexHighPurBJetTags = ak5PFbTagger.SimpleSecondaryVertexHighPurBJetTags
ak5PFCombinedSecondaryVertexBJetTags = ak5PFbTagger.CombinedSecondaryVertexBJetTags
ak5PFCombinedSecondaryVertexV2BJetTags = ak5PFbTagger.CombinedSecondaryVertexV2BJetTags
ak5PFJetBProbabilityBJetTags = ak5PFbTagger.JetBProbabilityBJetTags
ak5PFSoftPFMuonByPtBJetTags = ak5PFbTagger.SoftPFMuonByPtBJetTags
ak5PFSoftPFMuonByIP3dBJetTags = ak5PFbTagger.SoftPFMuonByIP3dBJetTags
ak5PFTrackCountingHighEffBJetTags = ak5PFbTagger.TrackCountingHighEffBJetTags
ak5PFTrackCountingHighPurBJetTags = ak5PFbTagger.TrackCountingHighPurBJetTags
ak5PFPatJetPartonAssociationLegacy = ak5PFbTagger.PatJetPartonAssociationLegacy

ak5PFImpactParameterTagInfos = ak5PFbTagger.ImpactParameterTagInfos
ak5PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
ak5PFJetProbabilityBJetTags = ak5PFbTagger.JetProbabilityBJetTags

ak5PFSecondaryVertexTagInfos = ak5PFbTagger.SecondaryVertexTagInfos
ak5PFSimpleSecondaryVertexHighEffBJetTags = ak5PFbTagger.SimpleSecondaryVertexHighEffBJetTags
ak5PFSimpleSecondaryVertexHighPurBJetTags = ak5PFbTagger.SimpleSecondaryVertexHighPurBJetTags
ak5PFCombinedSecondaryVertexBJetTags = ak5PFbTagger.CombinedSecondaryVertexBJetTags
ak5PFCombinedSecondaryVertexV2BJetTags = ak5PFbTagger.CombinedSecondaryVertexV2BJetTags

ak5PFSecondaryVertexNegativeTagInfos = ak5PFbTagger.SecondaryVertexNegativeTagInfos
ak5PFNegativeSimpleSecondaryVertexHighEffBJetTags = ak5PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
ak5PFNegativeSimpleSecondaryVertexHighPurBJetTags = ak5PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
ak5PFNegativeCombinedSecondaryVertexBJetTags = ak5PFbTagger.NegativeCombinedSecondaryVertexBJetTags
ak5PFPositiveCombinedSecondaryVertexBJetTags = ak5PFbTagger.PositiveCombinedSecondaryVertexBJetTags
ak5PFNegativeCombinedSecondaryVertexV2BJetTags = ak5PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
ak5PFPositiveCombinedSecondaryVertexV2BJetTags = ak5PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

ak5PFSoftPFMuonsTagInfos = ak5PFbTagger.SoftPFMuonsTagInfos
ak5PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
ak5PFSoftPFMuonBJetTags = ak5PFbTagger.SoftPFMuonBJetTags
ak5PFSoftPFMuonByIP3dBJetTags = ak5PFbTagger.SoftPFMuonByIP3dBJetTags
ak5PFSoftPFMuonByPtBJetTags = ak5PFbTagger.SoftPFMuonByPtBJetTags
ak5PFNegativeSoftPFMuonByPtBJetTags = ak5PFbTagger.NegativeSoftPFMuonByPtBJetTags
ak5PFPositiveSoftPFMuonByPtBJetTags = ak5PFbTagger.PositiveSoftPFMuonByPtBJetTags
ak5PFPatJetFlavourIdLegacy = cms.Sequence(ak5PFPatJetPartonAssociationLegacy*ak5PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#ak5PFPatJetFlavourAssociation = ak5PFbTagger.PatJetFlavourAssociation
#ak5PFPatJetFlavourId = cms.Sequence(ak5PFPatJetPartons*ak5PFPatJetFlavourAssociation)

ak5PFJetBtaggingIP       = cms.Sequence(ak5PFImpactParameterTagInfos *
            (ak5PFTrackCountingHighEffBJetTags +
             ak5PFTrackCountingHighPurBJetTags +
             ak5PFJetProbabilityBJetTags +
             ak5PFJetBProbabilityBJetTags 
            )
            )

ak5PFJetBtaggingSV = cms.Sequence(ak5PFImpactParameterTagInfos
            *
            ak5PFSecondaryVertexTagInfos
            * (ak5PFSimpleSecondaryVertexHighEffBJetTags+
                ak5PFSimpleSecondaryVertexHighPurBJetTags+
                ak5PFCombinedSecondaryVertexBJetTags+
                ak5PFCombinedSecondaryVertexV2BJetTags
              )
            )

ak5PFJetBtaggingNegSV = cms.Sequence(ak5PFImpactParameterTagInfos
            *
            ak5PFSecondaryVertexNegativeTagInfos
            * (ak5PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                ak5PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                ak5PFNegativeCombinedSecondaryVertexBJetTags+
                ak5PFPositiveCombinedSecondaryVertexBJetTags+
                ak5PFNegativeCombinedSecondaryVertexV2BJetTags+
                ak5PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

ak5PFJetBtaggingMu = cms.Sequence(ak5PFSoftPFMuonsTagInfos * (ak5PFSoftPFMuonBJetTags
                +
                ak5PFSoftPFMuonByIP3dBJetTags
                +
                ak5PFSoftPFMuonByPtBJetTags
                +
                ak5PFNegativeSoftPFMuonByPtBJetTags
                +
                ak5PFPositiveSoftPFMuonByPtBJetTags
              )
            )

ak5PFJetBtagging = cms.Sequence(ak5PFJetBtaggingIP
            *ak5PFJetBtaggingSV
            *ak5PFJetBtaggingNegSV
#            *ak5PFJetBtaggingMu
            )

ak5PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("ak5PFJets"),
        genJetMatch          = cms.InputTag("ak5PFmatch"),
        genPartonMatch       = cms.InputTag("ak5PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5PFcorr")),
        JetPartonMapSource   = cms.InputTag("ak5PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("ak5PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("ak5PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("ak5PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("ak5PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("ak5PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("ak5PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("ak5PFJetBProbabilityBJetTags"),
            cms.InputTag("ak5PFJetProbabilityBJetTags"),
            #cms.InputTag("ak5PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("ak5PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("ak5PFTrackCountingHighEffBJetTags"),
            cms.InputTag("ak5PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("ak5PFJetID"),
        addBTagInfo = True,
        addTagInfos = True,
        addDiscriminators = True,
        addAssociatedTracks = True,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = True,
        addGenPartonMatch = True,
        addGenJetMatch = True,
        embedGenJetMatch = True,
        embedGenPartonMatch = True,
        # embedCaloTowers = False,
        # embedPFCandidates = True
        )

ak5PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("ak5PFJets"),
           	    R0  = cms.double( 0.5)
)
ak5PFpatJetsWithBtagging.userData.userFloats.src += ['ak5PFNjettiness:tau1','ak5PFNjettiness:tau2','ak5PFNjettiness:tau3']

ak5PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("ak5PFpatJetsWithBtagging"),
                                                             genjetTag = 'ak5GenJets',
                                                             rParam = 0.5,
                                                             matchJets = cms.untracked.bool(False),
                                                             matchTag = 'patJetsWithBtagging',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlow'),
                                                             trackTag = cms.InputTag("generalTracks"),
                                                             fillGenJets = True,
                                                             isMC = True,
							     doSubEvent = False,
                                                             useHepMC = cms.untracked.bool(True),
							     genParticles = cms.untracked.InputTag("genParticles"),
							     eventInfoTag = cms.InputTag("generator"),
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(True),
                                                             bTagJetName = cms.untracked.string("ak5PF"),
                                                             jetName = cms.untracked.string("ak5PF"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(False),
							     doSubJets = cms.untracked.bool(False),
                                                             doGenSubJets = cms.untracked.bool(False),     
                                                             subjetGenTag = cms.untracked.InputTag("ak5GenJets"),
                                                             doGenTaus = False
                                                             )

ak5PFJetSequence_mc = cms.Sequence(
                                                  #ak5PFclean
                                                  #*
                                                  ak5PFmatch
                                                  #*
                                                  #ak5PFmatchGroomed
                                                  *
                                                  ak5PFparton
                                                  *
                                                  ak5PFcorr
                                                  *
                                                  #ak5PFJetID
                                                  #*
                                                  ak5PFPatJetFlavourIdLegacy
                                                  #*
			                          #ak5PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  ak5PFJetTracksAssociatorAtVertex
                                                  *
                                                  ak5PFJetBtagging
                                                  *
                                                  ak5PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  ak5PFpatJetsWithBtagging
                                                  *
                                                  ak5PFJetAnalyzer
                                                  )

ak5PFJetSequence_data = cms.Sequence(ak5PFcorr
                                                    *
                                                    #ak5PFJetID
                                                    #*
                                                    ak5PFJetTracksAssociatorAtVertex
                                                    *
                                                    ak5PFJetBtagging
                                                    *
                                                    ak5PFNjettiness 
                                                    *
                                                    ak5PFpatJetsWithBtagging
                                                    *
                                                    ak5PFJetAnalyzer
                                                    )

ak5PFJetSequence_jec = cms.Sequence(ak5PFJetSequence_mc)
ak5PFJetSequence_mb = cms.Sequence(ak5PFJetSequence_mc)

ak5PFJetSequence = cms.Sequence(ak5PFJetSequence_jec)
ak5PFJetAnalyzer.genPtMin = cms.untracked.double(1)
ak5PFJetAnalyzer.jetPtMin = cms.double(1)
