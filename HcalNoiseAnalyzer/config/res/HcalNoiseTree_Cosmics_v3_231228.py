# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECO --data -s RAW2DIGI,RECO --scenario cosmics --filein file:5C1B1DE5-9B38-E211-A048-001D09F24FBA.root --fileout DummyOutput.root --conditions GR_R_72_V1::All --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# in Gobinda's script
#process.load("Configuration.Geometry.GeometryIdeal_cff") 
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi") 
#process.load("RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi") 
#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi") 
#process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff") 

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(10000)
)

#process.metFilter = cms.EDFilter("METFilter",
#                                 CaloMETsrc = cms.InputTag("caloMet"),
#                                 MinMET     = cms.double(7.0)
#)

process.options = cms.untracked.PSet(
     wantSummary = cms.untracked.bool(True)
)
process.MessageLogger = cms.Service("MessageLogger",
 cout = cms.untracked.PSet( 
   default = cms.untracked.PSet( ## kill all messages in the log 
   limit = cms.untracked.int32(0) 
  ), 
  FwkJob = cms.untracked.PSet( ## but FwkJob category - those unlimitted 
   limit = cms.untracked.int32(-1) 
  ),
  FwkReport = cms.untracked.PSet(
   reportEvery = cms.untracked.int32(100), ## print event record number
   limit = cms.untracked.int32(-1) 
  ),
  FwkSummary = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
  #  reportEvery = cms.untracked.int32(100),
  #  limit = cms.untracked.int32(10000000)
  )
 ),
 categories = cms.untracked.vstring('FwkJob','FwkReport','FwkSummary'), 
 destinations = cms.untracked.vstring('cout')
)



# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    #fileNames = cms.untracked.vstring('/store/data//Commissioning2014/Cosmics/RAW/v3/000/225/125/00000/888ADCC7-352D-E411-B401-02163E00A091.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/MinimumBias/RAW/v3/000/224/512/00000/3CC0EA67-1727-E411-B11A-02163E008EFD.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v1/000/221/107/00000/521DB9A3-1EC1-E311-A9E4-02163E00BA2A.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/227/391/00000/90DAFD68-0750-E411-8E0A-02163E008BE3.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/MinimumBias/RAW/v3/000/227/391/00000/A2CD1AE3-0650-E411-A39E-02163E008CFE.root')
    #    fileNames = cms.untracked.vstring('/store/data/Commissioning2014/HcalHPDNoise/RAW/v3/000/227/489/00000/0CF75C7A-5855-E411-92A1-02163E00A129.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/22D7FAD0-6490-E411-891E-02163E0104D6.root')
    fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/265CEF40-5E90-E411-A5F4-02163E011C1F.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/26672730-5E90-E411-8C85-02163E011BE3.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/62B946C7-6490-E411-B25E-02163E01193D.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/6451164E-4790-E411-BB28-02163E011945.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/682BE4C5-6490-E411-8E8B-02163E011C45.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/6A76831A-6590-E411-B89D-02163E011BDE.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/8A0391CE-6490-E411-B4EF-02163E00FC3C.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/8E75352E-5E90-E411-9AD8-02163E00FB18.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/AA8D8C70-5E90-E411-9CAB-02163E00FC3C.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/B81762C5-6490-E411-A1CA-02163E0119E8.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/C21A4ECC-6490-E411-97E6-02163E011BE3.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/CAFE674F-5E90-E411-8A5F-02163E00FB9F.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/E06FFF2B-5E90-E411-8F72-02163E0104D6.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/E2F36DC2-6490-E411-84C7-02163E01192A.root')
    #fileNames = cms.untracked.vstring('/store/data/Commissioning2014/Cosmics/RAW/v3/000/231/228/00000/FC8BBCC1-6490-E411-BBB7-02163E00FDB9.root')
                            )


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('RECO nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

#process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    outputCommands = process.RECOSIMEventContent.outputCommands,
#    fileName = cms.untracked.string('DummyOutput.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('')
#    )
#)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#needed for 720 version
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_R_72_V1::All', '')
#needed for 703 version
#process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_R_70_V2::All', '')

# Hcal noise analyzers
process.HBHENoiseFilterResultProducer = cms.EDProducer(
   'HBHENoiseFilterResultProducer',
   noiselabel = cms.InputTag('hcalnoise'),
   minRatio = cms.double(-999),
   maxRatio = cms.double(999),
   minHPDHits = cms.int32(17),
   minRBXHits = cms.int32(999),
   minHPDNoOtherHits = cms.int32(10),
   minZeros = cms.int32(10),
   minHighEHitTime = cms.double(-9999.0),
   maxHighEHitTime = cms.double(9999.0),
   maxRBXEMF = cms.double(-999.0),
   minNumIsolatedNoiseChannels = cms.int32(10),
   minIsolatedNoiseSumE = cms.double(50.0),
   minIsolatedNoiseSumEt = cms.double(25.0),
   useTS4TS5 = cms.bool(False),
   useRBXRechitR45Loose = cms.bool(False),
   useRBXRechitR45Tight = cms.bool(False),
   IgnoreTS4TS5ifJetInLowBVRegion = cms.bool(True),
   jetlabel = cms.InputTag('ak5PFJets'),
   maxjetindex = cms.int32(0),
   maxNHF = cms.double(0.9)
   )
process.TFileService = cms.Service("TFileService",
   #fileName = cms.string("NoiseTree_Commissionig2014_HcalHPDNoise_v3_227489.root")
   fileName = cms.string("/tmp/fahim/NoiseTree_Commissionig2014_Cosmics_v3_231228_02.root")
   )
process.ExportTree = cms.EDAnalyzer("HcalNoiseAnalyzer",
  HBHERecHitCollection = cms.untracked.string('hbhereco'),
  IsCosmic             = cms.untracked.bool(True)
)
process.hcalNoiseAna = cms.EDAnalyzer('HcalNoiseHistogrammer',
    HBHERecHitCollection = cms.InputTag("hbhereco"),
    HBHEDigiCollection   = cms.InputTag("hcalDigis")
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionCosmics * process.HBHENoiseFilterResultProducer * process.ExportTree)
#process.reconstruction_step = cms.Path(process.reconstructionCosmics * process.HBHENoiseFilterResultProducer * process.hcalNoiseAna)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.endjob_step)


#===================== Message Logger =============================                                                                                         
##process.load("FWCore.MessageLogger.MessageLogger_cfi")
##process.MessageLogger.categories.append('PATSummaryTables')
##process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
##    limit = cms.untracked.int32(10),
##    reportEvery = cms.untracked.int32(1)
##    )
##process.options = cms.untracked.PSet(
##    wantSummary = cms.untracked.bool(True)
##    )
##process.MessageLogger.cerr.FwkReport.reportEvery = 1000

