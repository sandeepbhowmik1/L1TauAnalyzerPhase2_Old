import FWCore.ParameterSet.Config as cms
process = cms.Process('Analyze')
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
	)
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1379.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1380.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1381.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1382.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1383.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1384.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1385.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1386.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1387.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1388.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1390.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1391.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1392.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1393.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1394.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1395.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1396.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1397.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1399.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1400.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1401.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1402.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1403.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1404.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1405.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1406.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1407.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1408.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1409.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1411.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1412.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1414.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1416.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1417.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1418.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1419.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1420.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1421.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1422.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1423.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1424.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1425.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1426.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1427.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1428.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1429.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1430.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1431.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1432.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1433.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1434.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1436.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1437.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1438.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1439.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1440.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1441.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1442.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1443.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1444.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1445.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1446.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1447.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1448.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1449.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1450.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1451.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1452.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1453.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1454.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1455.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1456.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1457.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1458.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1459.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1460.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1461.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1462.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1463.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1464.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1465.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1466.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1467.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1468.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1469.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1470.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1471.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1472.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1473.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1474.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1475.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1476.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1477.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1478.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1480.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1481.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1482.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1483.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1484.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1487.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1489.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1490.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1491.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1492.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1493.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1494.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1495.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1496.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1497.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1498.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1499.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1500.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1501.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1502.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1503.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1504.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1505.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1506.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1507.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1508.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1509.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1510.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1512.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1514.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1515.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1516.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1517.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1518.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1519.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1521.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1522.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1523.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1525.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1526.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1527.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1528.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1529.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1530.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1531.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1532.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1533.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1535.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1536.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1537.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1538.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1539.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1540.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1541.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1542.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1543.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1544.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1545.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1546.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1548.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1549.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1550.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1551.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1553.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1554.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1555.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1556.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1557.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1558.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1560.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1561.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1562.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1563.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1564.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1565.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1567.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1568.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1569.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1570.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1571.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1573.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1574.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1575.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1576.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1577.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1578.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1579.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1580.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1581.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1582.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1583.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1584.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1585.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1586.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1587.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1588.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1589.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1590.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1591.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1592.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1593.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1594.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1595.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1596.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1597.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1598.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1599.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1600.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1601.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1602.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1603.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1604.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1605.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1606.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1607.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1608.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1609.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1610.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1611.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1613.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1615.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1617.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1619.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1620.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1621.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1622.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1623.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1624.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1625.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1626.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1627.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1629.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1630.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1631.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1632.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1633.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1634.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1635.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1636.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1637.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1638.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1640.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1641.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1642.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1643.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1644.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1645.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1646.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1647.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1648.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1649.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1650.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1651.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1652.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1654.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1655.root',
	)
)
process.load("L1TauAnalyzerPhase2.L1TauAnalyzerPhase2.L1TkElectronHPSPFTauAnalyzer_cff")
process.L1TkElectronHPSPFTauAnalyzer.histRootFileName = cms.string("hist_test_L1TkElectronHPSPFTauAnalyzer_Background_20200803_part_6.root")
process.L1TkElectronHPSPFTauAnalyzer.bdtRootFileName = cms.string("bdt_test_L1TkElectronHPSPFTauAnalyzer_Background_20200803_part_6.root")
process.p = cms.Path(
	process.AnalyzerSeq
)
process.schedule = cms.Schedule(process.p)
process.TFileService=cms.Service('TFileService',fileName=cms.string("rootTree_test_L1TkElectronHPSPFTauAnalyzer_Background_20200803_part_6.root"))
