import FWCore.ParameterSet.Config as cms
process = cms.Process('Analyze')
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
	)
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3024.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3026.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3027.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3028.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3029.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3030.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3031.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3032.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3033.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3034.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3035.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3036.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3037.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3038.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3039.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3040.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3041.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3042.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3043.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3044.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3045.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3046.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3047.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3049.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3050.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3051.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3052.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3053.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3055.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3056.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3057.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3058.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3059.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3060.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3061.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3062.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3063.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3064.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3065.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3066.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3067.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3068.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3069.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3070.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3071.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3072.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3073.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3074.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3075.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3076.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3077.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3078.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3079.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3080.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3081.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3082.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3083.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3084.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3085.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3086.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3087.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3088.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3089.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3090.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3091.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3092.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3093.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3094.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3095.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3096.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3097.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3098.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3099.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3101.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3102.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3103.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3104.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3105.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3106.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3107.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3109.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3110.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3111.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3112.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3114.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3115.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3116.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3117.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3118.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3119.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3120.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3121.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3122.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3123.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3124.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3125.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3126.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3127.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3128.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3129.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3130.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3131.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3132.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3133.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3134.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3135.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3138.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3139.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3140.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3142.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3143.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3144.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3146.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3147.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3148.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3149.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3150.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3151.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3152.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3153.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3154.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3155.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3156.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3158.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3159.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3160.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3161.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3162.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3163.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3164.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3165.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3166.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3168.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3169.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3170.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3171.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3172.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3173.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3174.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3176.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3177.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3178.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3179.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3180.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3181.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3182.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3183.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3184.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3185.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3186.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3187.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3188.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3189.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3190.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3192.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3193.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3194.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3195.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3196.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3197.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3198.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3200.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3201.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3202.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3203.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3204.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3205.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3206.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3207.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3208.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3209.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3210.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3211.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3212.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3213.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3214.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3215.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3216.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3217.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3218.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3219.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3220.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3221.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3222.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3223.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3224.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3225.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3226.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3227.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3228.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3229.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3231.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3232.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3233.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3234.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3235.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3236.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3237.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3238.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3239.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3240.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3241.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3242.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3243.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3244.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3245.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3247.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3248.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3249.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3250.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3251.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3252.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3253.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3254.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3255.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3256.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3257.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3258.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3259.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3260.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3261.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3262.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3263.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3264.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3265.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3266.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3267.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3268.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3269.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3270.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3271.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3272.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3274.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3275.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3276.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3277.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3278.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3279.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3280.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3281.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3282.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3284.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3285.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3286.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3287.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3288.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3289.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3290.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3291.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0003/NTuple_L1TkMuonHPSPFTauProducer_3292.root',
	)
)
process.load("L1TauAnalyzerPhase2.L1TauAnalyzerPhase2.L1TkMuonHPSPFTauAnalyzer_cff")
process.L1TkMuonHPSPFTauAnalyzer.histRootFileName = cms.string("hist_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_12.root")
process.L1TkMuonHPSPFTauAnalyzer.bdtRootFileName = cms.string("bdt_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_12.root")
process.p = cms.Path(
	process.AnalyzerSeq
)
process.schedule = cms.Schedule(process.p)
process.TFileService=cms.Service('TFileService',fileName=cms.string("rootTree_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_12.root"))
