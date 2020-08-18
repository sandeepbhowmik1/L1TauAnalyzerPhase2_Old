import FWCore.ParameterSet.Config as cms
process = cms.Process('Analyze')
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
	)
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3293.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3294.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3295.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3296.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3298.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3299.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3300.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3301.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3302.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3303.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3304.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3305.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3306.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3307.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3308.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3309.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3310.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3311.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3313.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3314.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3315.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3316.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3317.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3318.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3319.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3320.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3321.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3322.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3323.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3324.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3325.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3326.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3327.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3328.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3329.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3330.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3331.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3332.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3333.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3334.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3335.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3336.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3337.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3338.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3339.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3340.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3341.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3342.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3343.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3344.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3345.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3346.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3347.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3348.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3349.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3350.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3352.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3353.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3354.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3356.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3357.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3358.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3359.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3360.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3361.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3362.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3363.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3364.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3365.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3366.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3367.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3368.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3369.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3370.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3371.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3372.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3373.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3374.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3375.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3376.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3378.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3379.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3380.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3382.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3383.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3384.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3385.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3386.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3387.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3389.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3390.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3391.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3392.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3393.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3394.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3395.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3396.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3399.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3400.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3401.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3402.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3403.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3404.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3405.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3406.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3407.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3408.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3409.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3410.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3411.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3412.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3413.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3414.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3415.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3417.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3418.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3419.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3420.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3421.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3422.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3423.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3424.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3425.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3426.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3427.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3428.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3429.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3430.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3431.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3432.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3433.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3434.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3435.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3436.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3437.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3439.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3441.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3442.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3443.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3444.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3445.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3446.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3448.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3449.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3450.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3452.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3453.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3454.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3455.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3456.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3457.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3458.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3459.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3460.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3461.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3462.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3463.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3464.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3465.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3466.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3467.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3468.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3469.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3470.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3471.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3472.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3473.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3476.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3477.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3478.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3479.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3480.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3481.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3482.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3483.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3484.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3485.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3486.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3487.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3488.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3489.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3490.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3491.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3493.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3494.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3495.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3496.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3497.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3498.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3499.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3500.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3501.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3502.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3503.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3504.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3505.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3506.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3507.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3508.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3509.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3510.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3511.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3512.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3513.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3514.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3516.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3518.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3519.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3521.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3522.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3523.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3524.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3525.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3526.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3527.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3528.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3529.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3530.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3531.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3532.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3533.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3534.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3535.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3536.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3537.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3538.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3539.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3540.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3541.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3542.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3543.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3544.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3545.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3546.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3547.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3548.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3549.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3550.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3551.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3552.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3553.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3554.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3555.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3556.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3557.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3558.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3559.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3560.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3561.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3562.root',
	)
)
process.load("L1TauAnalyzerPhase2.L1TauAnalyzerPhase2.L1TkElectronHPSPFTauAnalyzer_cff")
process.L1TkElectronHPSPFTauAnalyzer.histRootFileName = cms.string("hist_test_L1TkElectronHPSPFTauAnalyzer_Background_20200803_part_13.root")
process.L1TkElectronHPSPFTauAnalyzer.bdtRootFileName = cms.string("bdt_test_L1TkElectronHPSPFTauAnalyzer_Background_20200803_part_13.root")
process.p = cms.Path(
	process.AnalyzerSeq
)
process.schedule = cms.Schedule(process.p)
process.TFileService=cms.Service('TFileService',fileName=cms.string("rootTree_test_L1TkElectronHPSPFTauAnalyzer_Background_20200803_part_13.root"))