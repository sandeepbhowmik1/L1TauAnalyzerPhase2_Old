import FWCore.ParameterSet.Config as cms
process = cms.Process('Analyze')
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
	)
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3832.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3833.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3834.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3835.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3836.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3837.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3838.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3839.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3840.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3841.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3842.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3843.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3844.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3845.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3846.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3847.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3848.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3849.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3851.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3852.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3853.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3854.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3855.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3856.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3857.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3858.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3859.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3860.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3861.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3862.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3863.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3864.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3865.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3866.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3867.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3868.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3871.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3872.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3873.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3874.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3875.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3877.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3879.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3880.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3881.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3882.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3883.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3884.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3885.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3886.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3887.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3888.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3889.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3890.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3891.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3892.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3893.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3895.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3896.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3897.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3899.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3900.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3901.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3902.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3903.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3904.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3905.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3906.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3907.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3909.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3911.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3912.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3913.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3914.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3915.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3916.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3917.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3918.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3919.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3920.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3921.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3922.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3923.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3924.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3925.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3926.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3927.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3928.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3929.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3930.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3931.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3932.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3933.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3934.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3935.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3937.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3939.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3940.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3941.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3942.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3943.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3944.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3945.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3946.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3947.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3948.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3949.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3950.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3951.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3952.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3953.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3954.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3955.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3956.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3957.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3958.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3959.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3967.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3968.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3974.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3975.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3982.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3995.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4000.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4009.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4011.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4037.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4038.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4055.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4057.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4058.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4061.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4065.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4099.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4100.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4114.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4115.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4133.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4154.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4155.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4156.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4157.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4158.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4159.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4160.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4161.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4162.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4163.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4164.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4165.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4166.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4167.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4168.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4169.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4170.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4171.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4172.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4173.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4174.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4175.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4176.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4177.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4178.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4179.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4180.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4181.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4182.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4183.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4184.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4185.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4186.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4187.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4188.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4189.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4190.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4191.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4192.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4193.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4194.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4195.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4196.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4197.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4198.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4199.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4200.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4201.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4202.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4203.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4204.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4205.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4206.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4207.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4208.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4209.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4210.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4211.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4212.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4213.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4214.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4215.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4216.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4217.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4218.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4219.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4220.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4221.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4222.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4223.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4224.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4225.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4226.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4227.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4228.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4229.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4230.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4231.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4232.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4233.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4234.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4235.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4236.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4237.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4238.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4239.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4240.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4241.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4242.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4243.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4244.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4245.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4246.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4247.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4248.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4249.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4250.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4251.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4252.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4253.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4254.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4255.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4256.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4257.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4258.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4259.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4260.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4261.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4262.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4263.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4264.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0004/NTuple_L1TkMuonHPSPFTauProducer_4265.root',
	)
)
process.load("L1TauAnalyzerPhase2.L1TauAnalyzerPhase2.L1TkMuonHPSPFTauAnalyzer_cff")
process.L1TkMuonHPSPFTauAnalyzer.histRootFileName = cms.string("hist_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_15.root")
process.L1TkMuonHPSPFTauAnalyzer.bdtRootFileName = cms.string("bdt_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_15.root")
process.p = cms.Path(
	process.AnalyzerSeq
)
process.schedule = cms.Schedule(process.p)
process.TFileService=cms.Service('TFileService',fileName=cms.string("rootTree_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_15.root"))