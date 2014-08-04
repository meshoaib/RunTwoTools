#crab config file for generic MiniSelector job
#to use: change requestName, psetName, outputFiles, and inputDataset to desired values

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'miniAOD_example'
config.General.workArea = 'crabMiniAOD'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/miniExampleConfig.py'
config.JobType.outputFiles = ['example.root']

config.section_("Data")
config.Data.inputDataset = '/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'

config.Data.dbsUrl = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
#config.Data.publishDbsUrl = 'phys03'
#config.Data.publishDataName = 'miniSelector_generic'
config.Data.ignoreLocality = True #enable AAA

config.section_("Site")
config.Site.storageSite = 'T2_US_Caltech'
