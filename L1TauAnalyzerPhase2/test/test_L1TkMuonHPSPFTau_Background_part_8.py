import FWCore.ParameterSet.Config as cms
process = cms.Process('Analyze')
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
	)
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1915.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1916.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1917.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1920.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1921.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1922.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1923.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1924.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1925.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1926.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1927.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1929.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1930.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1931.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1932.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1933.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1934.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1935.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1936.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1937.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1938.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1940.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1941.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1942.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1943.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1944.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1945.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1946.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1947.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1948.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1949.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1950.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1952.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1953.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1954.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1955.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1956.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1958.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1959.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1960.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1961.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1962.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1963.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1964.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1966.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1967.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1968.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1969.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1970.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1971.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1972.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1973.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1974.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1975.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1977.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1978.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1979.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1980.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1981.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1982.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1983.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1984.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1985.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1986.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1987.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1988.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1990.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1992.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1993.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1994.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1995.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1996.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1997.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0001/NTuple_L1TkMuonHPSPFTauProducer_1999.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2001.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2002.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2004.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2005.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2006.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2007.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2008.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2009.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2010.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2011.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2012.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2013.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2014.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2015.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2016.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2017.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2018.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2019.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2020.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2021.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2022.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2023.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2024.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2025.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2026.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2027.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2028.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2030.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2031.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2033.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2034.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2036.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2037.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2038.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2039.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2040.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2041.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2042.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2043.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2044.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2045.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2046.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2047.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2048.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2049.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2051.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2052.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2054.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2055.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2056.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2057.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2059.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2060.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2061.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2062.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2063.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2064.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2065.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2067.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2068.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2069.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2070.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2071.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2072.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2073.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2074.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2075.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2076.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2077.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2078.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2079.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2080.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2081.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2083.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2084.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2085.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2086.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2087.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2088.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2089.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2091.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2092.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2093.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2094.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2096.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2097.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2098.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2099.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2100.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2102.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2103.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2104.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2105.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2107.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2108.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2110.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2111.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2113.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2114.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2115.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2116.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2117.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2118.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2119.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2120.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2121.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2122.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2123.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2124.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2125.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2127.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2128.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2129.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2130.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2131.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2132.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2133.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2134.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2135.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2136.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2137.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2139.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2140.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2141.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2142.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2143.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2144.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2145.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2146.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2147.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2148.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2149.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2150.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2151.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2152.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2153.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2154.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2155.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2156.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2157.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2158.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2159.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2160.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2161.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2162.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2163.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2164.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2165.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2166.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2167.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2168.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2169.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2170.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2171.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2172.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2174.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2175.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2176.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2177.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2178.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2180.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2181.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2182.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2183.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2184.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2186.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2187.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2188.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2189.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2190.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2191.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2192.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2193.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2194.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2195.root',
	'file:/hdfs/cms/store/user/sbhowmik/Nu_E10-pythia8-gun/PhaseIITDRSpring19MiniAOD_20200422/200423_140455/0002/NTuple_L1TkMuonHPSPFTauProducer_2196.root',
	)
)
process.load("L1TauAnalyzerPhase2.L1TauAnalyzerPhase2.L1TkMuonHPSPFTauAnalyzer_cff")
process.L1TkMuonHPSPFTauAnalyzer.histRootFileName = cms.string("hist_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_8.root")
process.L1TkMuonHPSPFTauAnalyzer.bdtRootFileName = cms.string("bdt_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_8.root")
process.p = cms.Path(
	process.AnalyzerSeq
)
process.schedule = cms.Schedule(process.p)
process.TFileService=cms.Service('TFileService',fileName=cms.string("rootTree_test_L1TkMuonHPSPFTauAnalyzer_Background_20200803_part_8.root"))
