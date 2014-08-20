#crab config file for generic MiniSelector job
#to use: change requestName, psetName, outputFiles, and inputDataset to desired values

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'razorJetAnalysis_TT40_logSquared_fixed'
config.General.workArea = 'crabMiniAOD'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/razorJetConfig.py'
config.JobType.outputFiles = ['razorJetAnalysis.root']

config.section_("Data")
#25ns scenario
#config.Data.inputDataset = '/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'

#50ns scenario
config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/StoreResults-Spring14dr_PU_S14_POSTLS170_V6AN1_miniAOD706p1_814812ec83fce2f620905d2bb30e9100-v2/USER'
#config.Data.inputDataset = '/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/duanders-Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1-814812ec83fce2f620905d2bb30e9100/USER'

#config.Data.dbsUrl = 'global'
config.Data.dbsUrl = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.publication = False
#config.Data.publishDbsUrl = 'phys03'
#config.Data.publishDataName = 'miniSelector_generic'
config.Data.ignoreLocality = True #enable AAA

config.section_("Site")
config.Site.storageSite = 'T2_US_Caltech'
