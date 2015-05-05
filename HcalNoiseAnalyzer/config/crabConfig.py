from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Cosmics_Tree_2014-v3_RAW_HcalHPDNoise_227489'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HcalNoiseTree_RECO_RAW2DIGI_RECO_Cosmics.py'

config.section_("Data")
config.Data.inputDataset = '/HcalHPDNoise/Commissioning2014-v3/RAW'
config.Data.dbsUrl = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20 # 200
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.lumiMask = 'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt' # if you downloaded the file in the working directory
config.Data.runRange = '227489' # '193093-194075'
config.Data.publication = False
#config.Data.publishDbsUrl = 'phys03'
#config.Data.publishDataName = 'CRAB3_tutorial_Data_analysis_test5'

config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'
