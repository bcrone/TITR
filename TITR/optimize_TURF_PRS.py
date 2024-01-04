import numpy as np 
import pandas as pd
import glob
from subprocess import Popen, PIPE
import argparse
import itertools
from itertools import islice
from sklearn.linear_model import LinearRegression
from scipy.stats import t
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.diagnostic import compare_j
from statsmodels.stats.diagnostic import compare_cox
import pickle
import time
import multiprocessing
import scipy.stats as ss

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def calc_pval(model, x, y):
	y_len = len(y)
	beta_hat = [model.intercept_] + model.coef_.tolist()
	x1 = np.column_stack((np.ones(y_len), x))
	st_dev_noise = np.sqrt(np.sum(np.square(y - x1@beta_hat))/(y_len - x1.shape[1]))
	beta_cov = np.linalg.inv(x1.T@x1)
	t_val = beta_hat/(st_dev_noise*np.sqrt(np.diagonal(beta_cov)))
	p_vals = t.sf(np.abs(t_val), y_len-x1.shape[1])*2
	p = p_vals[-1]
	return p

def calc_R2(model, x, y, adj_R2):
	y_len = len(y)
	ss_residual = sum((y - model.predict(x))**2)
	ss_total = sum((y - np.mean(y))**2)
	R2 = 1 - float(ss_residual)/ss_total
	new_adj_R2 = 1 - (1 - R2)*(y_len-1)/(y_len-x.shape[1]-1)
	d_R2 = new_adj_R2 - adj_R2
	return new_adj_R2, d_R2

def rank_INT(series, c=3.0/8):
	np.random.seed(123)
	orig_idx = series.index
	series = series.loc[~pd.isnull(series)]
	rank = ss.rankdata(series, method="average")
	rank = pd.Series(rank, index=series.index)
	transformed = rank.apply(rank_to_normal, c=c, n=len(rank))
	return transformed[orig_idx]

def rank_to_normal(rank, c, n):
	x = (rank - c) / (n - 2*c + 1)
	return ss.norm.ppf(x)

def initializeNullModel(covariates_path, samples_path, trait, isQuant):
	null_pheno = pd.read_table(covariates_path, sep='\t', index_col='FID')
	samples = pd.read_table(samples_path, sep=' ', header=None)
	samples.columns = ['FID','IID']
	null_pheno = null_pheno[null_pheno['IID'].isin(samples['IID'])]
	null_pheno = null_pheno[null_pheno[trait].notna()]
	x = (null_pheno[[
		"Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
		]]).to_numpy()
	if isQuant:
		null_pheno["norm"] = rank_INT(null_pheno[trait])
		y = (null_pheno["norm"]).to_numpy()
	else:
		y = (null_pheno[f"{trait}"]).to_numpy()
	null_model = LinearRegression(n_jobs=-1).fit(x, y)
	null_adj_R2, null_d_R2 = calc_R2(null_model, x, y, 0)
	return [null_pheno, null_model, null_adj_R2]

def calculateCodingModel(coding_prefix, null_pheno, null_adj_R2, isQuant):
	results = []
	coding_results = pd.DataFrame(columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","SNP Count","Null Adj R2"])
	coding_files = glob.glob(f"{coding_prefix}.*.sscore")
	for coding_file in coding_files:
		try:
			coding_table = pd.read_table(coding_file, delim_whitespace=True, index_col="IID")
		except OSError as e:
			return
		p = f"0.{coding_file.split('.')[4]}"
		if p == "0.sscore":
			p = 1
		coding_count = max(coding_table["ALLELE_CT"])/2
		coding_pheno = pd.merge(null_pheno, coding_table["SCORE1_AVG"], on="IID")
		x = (coding_pheno[[
			"Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "SCORE1_AVG"
			]]).to_numpy()
		if isQuant:
			y = (coding_pheno["norm"]).to_numpy()
		else:
			y = (coding_pheno[f"{trait}"]).to_numpy()
		coding_model = LinearRegression(n_jobs=-1).fit(x, y)
		coding_adj_R2, coding_d_R2 = calc_R2(coding_model, x, y, null_adj_R2)
		coding_p = calc_pval(coding_model, x, y)
		results.append([p, "Coding", "NA", coding_adj_R2, coding_d_R2, coding_p, coding_count, null_adj_R2])
	return results

def addCodingModel(threshold, coding_prefix, null_pheno, null_adj_R2, isQuant):
	coding_file = f"{coding_prefix}.{threshold}.sscore"
	try:
		coding_table = pd.read_table(coding_file, delim_whitespace=True, index_col="IID")
	except OSError as e:
		return
	coding_count = max(coding_table['ALLELE_CT'])/2
	coding_pheno = pd.merge(null_pheno, coding_table["SCORE1_AVG"], on="IID")
	coding_pheno.rename(columns={'SCORE1_AVG':'CODING'},inplace=True)	
	if isQuant:
		coding_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING", data=coding_pheno).fit()
	else:
		coding_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING", data=coding_pheno).fit()
	coding_adj_R2 = coding_model.rsquared_adj
	coding_d_R2 = coding_adj_R2 - null_adj_R2
	coding_p = coding_model.pvalues['CODING']
	return [coding_pheno, coding_model, coding_adj_R2, coding_d_R2, coding_p, coding_count]

def updateNullModel(trait, iteration, traits_path, pheno, null_adj_R2, isQuant):
	null_table = pd.read_table(f"{traits_path}/{trait}/{trait}.TURF.ITERATION_{iteration}.sscore", delim_whitespace=True, index_col="IID")
	null_table.rename(columns={'SCORE1_AVG':'MASTER'},inplace=True)
	if iteration == 1: # 1st iteration logic
		master_pheno = pd.merge(pheno, null_table['MASTER'], on="IID")
	else:
		master_pheno = pd.merge(pheno.drop(columns='MASTER'), null_table['MASTER'], on="IID")
	if isQuant:
		master_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + MASTER", data=master_pheno).fit()
	else:
		master_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + MASTER", data=master_pheno).fit()
	master_adj_R2 = master_model.rsquared_adj
	master_adj_d_R2 = master_adj_R2 - null_adj_R2
	return [master_pheno, master_model, master_adj_R2, master_adj_d_R2]

def runComparison(trait, tissues_partitions, iteration, model, pheno, adj_R2, added_bin_paths, results_path, isQuant, lead_tissues):  # 1:514.4733331203461
	def expand_tissues(tissues, partitions):
		return list(itertools.chain.from_iterable(itertools.repeat(x, len(partitions)) for x in tissues))
	def get_TURF_results(file, iteration, isQuant, tissue, partition, adj_R2):
		p = f"0.{file.split('.')[5]}"
		if p == "0.sscore":
			p = 1
		TURF_table = pd.read_table(file, delim_whitespace=True, index_col='IID')
		if TURF_table['SCORE1_AVG'].isnull().values.any():
			return {'Threshold':p, 'Tissue':tissue, 'Model':partition, "Adj R2":adj_R2, 
				"Delta Adj R2":np.nan, "P":np.nan, "Count":np.nan}
		TURF_count = max(TURF_table["ALLELE_CT"])/2
		TURF_pheno = pd.merge(pheno, TURF_table["SCORE1_AVG"], on="IID")
		if iteration == 1: # 1st iteration logic
			x = (TURF_pheno[[
				"Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "CODING", "SCORE1_AVG"]]
				).to_numpy()
		else:
			x = (TURF_pheno[[
				"Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "CODING", "SCORE1_AVG"]]
				).to_numpy()
		if isQuant:
			y = (TURF_pheno["norm"]).to_numpy()
		else:
			y = (TURF_pheno[f"{trait}"]).to_numpy()
		TURF_model = LinearRegression(n_jobs=-1).fit(x, y)
		TURF_adj_R2, TURF_d_R2 = calc_R2(TURF_model, x, y, adj_R2)
		TURF_p = calc_pval(TURF_model, x, y)
		row = {'Threshold':p, 'Tissue':tissue, 'Model':partition, "Adj R2":TURF_adj_R2, 
			"Delta Adj R2":TURF_d_R2, "P":TURF_p, "Count":TURF_count}
		return row
	def get_max_TURF(tissue, partition, results_path, iteration, isQuant, adj_R2, model, pheno):
		tissue_path = f"{results_path}/{tissue}"
		tissue_files = glob.glob(f"{tissue_path}/{trait}.{tissue}.{partition}.ITERATION_{iteration}.*.sscore")
		TURF_result = pd.DataFrame(
			[get_TURF_results(file, iteration, isQuant, tissue, partition, adj_R2) for file in tissue_files],
			columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","Count"])
		TURF_result_pass = pd.DataFrame(columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","Count"])
		TURF_result_fail = pd.DataFrame(columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","Count","fwd","bwd"])
		TURF_result.dropna(inplace=True)
		for index, row in TURF_result.iterrows():
			if iteration==1:
				if runANOVA(row, model, pheno, trait, iteration, results_path, isQuant):
					TURF_result_pass = TURF_result_pass.append(row, ignore_index=True)
			else:
				[fwd, bwd] = runCoxtest(row, model, pheno, trait, iteration, results_path, isQuant)
				if fwd:
					TURF_result_pass = TURF_result_pass.append(row, ignore_index=True)
				else:
					row['fwd'] = fwd
					row['bwd'] = bwd
					TURF_result_fail = TURF_result_fail.append(row, ignore_index=True)
		if len(TURF_result_pass) == 0:
			TURF_result_fail.to_csv(f"{results_path}/traits/{trait}/{trait}.{tissue}.ITERATION_{iteration}.failures.tsv",sep="\t",index=False)
		TURF_bin_max = TURF_result_pass[TURF_result_pass['Adj R2']==TURF_result_pass['Adj R2'].max()]
		return TURF_result_pass
	TURF_head = pd.DataFrame(columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","Count"])
	TURF_max = [get_max_TURF(tissue, partition, results_path, iteration, isQuant, adj_R2, model, pheno) for tissue, partition in tissues_partitions.items() if f"{tissue}.{partition}" not in added_bin_paths and tissue in lead_tissues]
	return TURF_head.append(TURF_max, ignore_index=True)

def runCoxtest(max_model, model, pheno, trait, iteration, results_path, isQuant):
	path=f"{results_path}/{max_model.Tissue}"
	file = f"{path}/{trait}.{max_model.Tissue}.{max_model.Model}.ITERATION_{iteration}.{max_model.Threshold}.sscore"
	model_table = pd.read_table(file, delim_whitespace=True, index_col="IID")
	model_count = max(model_table["ALLELE_CT"])
	model_pheno = pd.merge(pheno, model_table["SCORE1_AVG"], on="IID")
	if isQuant:
		model_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + SCORE1_AVG", data=model_pheno).fit()
	else:
		model_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + SCORE1_AVG", data=model_pheno).fit()
	return [compare_cox(model_model, model)[1] <= 0.05, compare_cox(model, model_model)[1] > 0.05]

def runANOVA(max_model, model, pheno, trait, iteration, results_path, isQuant):
	path = f"{results_path}/{max_model.Tissue}"
	file = f"{path}/{trait}.{max_model.Tissue}.{max_model.Model}.ITERATION_{iteration}.{max_model.Threshold}.sscore"
	model_table = pd.read_table(file, delim_whitespace=True, index_col="IID")
	model_count = max(model_table["ALLELE_CT"])
	model_pheno = pd.merge(pheno, model_table["SCORE1_AVG"], on="IID")
	if iteration == 1: # 1st iteration logic
		if isQuant:
			model_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + SCORE1_AVG", data=model_pheno).fit()
		else:
			model_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + SCORE1_AVG", data=model_pheno).fit()
	else:
		if isQuant:
			model_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + MASTER + SCORE1_AVG", data=model_pheno).fit()
		else:
			model_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + MASTER + SCORE1_AVG", data=model_pheno).fit()
	return anova_lm(model,model_model)['Pr(>F)'][1] < 0.05

def addPartitionSNPs(master, trait, tissue, partition, threshold, iteration, GWAS_path, results_path):
	snp_pvalue = pd.read_table(f"{GWAS_path}/{trait}/{trait}.PLINK.TITR.SNP.pvalue", delim_whitespace=True)
	snp_pvalue.rename(columns={'SNP':'ID'}, inplace=True)
	partition_snp = pd.read_table(f"{results_path}/{tissue}/{trait}.{tissue}.{partition}.ITERATION_{iteration}.snplist",names=['ID'])
	partition_snp_pvalue = pd.merge(partition_snp, snp_pvalue, on="ID")
	partition_snp_threshold = partition_snp_pvalue[partition_snp_pvalue["P"] <= float(threshold)]
	[master.append(snp) for snp in partition_snp_threshold.ID.values.tolist() if snp not in master]

def tagClumpedSNPs(trait, tissue, partition, threshold, iteration, GWAS_path, results_path):
	clumped_snp = pd.read_table(f"{results_path}/{tissue}/{trait}.{tissue}.{partition}.ITERATION_{iteration}.clumped", delim_whitespace=True)
	clumped_snp_threshold = clumped_snp[clumped_snp["P"] < float(threshold)]
	clumped_snp_threshold['SP2'] = clumped_snp_threshold['SP2'].str.split(',').apply(lambda x: [e.strip('(1)') for e in x])
	return list(set(clumped_snp_threshold['SP2'].sum()))

def writeClumpedList(master, trait, ancestry, iteration, results_path):
	with open(f"{results_path}/traits/{trait}/{ancestry}.{trait}.ITERATION_{iteration}.SNPs", "w") as f:
		for snp in master:
			f.write("%s\n" % snp)

def writeExcludedSNPs(excluded_list, trait, ancestry, iteration, results_path):
	with open(f"{results_path}/traits/{trait}/{trait}.TURF.ITERATION_{iteration}.exclude", "w") as f:
		for snp in excluded_list:
			f.write("%s\n" % snp)

def writeModelTerminate(trait, results_path):
	open(f"{results_path}/{trait}/{trait}.terminate","a").close()

def writeMasterResults(master_results, trait, results_path):
	master_results.to_csv(f"{results_path}/{trait}/{trait}.TURF.master_results.tsv", sep="\t", index=False)

def readMasterResults(trait, results_path):
	master_results = pd.read_table(f"{results_path}/{trait}/{trait}.TURF.master_results.tsv")
	return master_results

def writeMasterPheno(master_pheno, trait, iteration, results_path):
	master_pheno.to_csv(f"{results_path}/{trait}/{trait}.TURF.master_pheno.ITERATION_{iteration}.tsv", sep="\t", index=False)

def writeTURFMax(TURF_max, trait, iteration, results_path):
	TURF_max.to_csv(f"{results_path}/{trait}/{trait}.TURF_max.ITERATION_{iteration}.tsv", sep="\t", index=False)

def readCoefTable(trait, results_path):
	coef_table = pd.read_table(f"{results_path}/{trait}/{trait}.TURF.master_coefs.tsv")
	return coef_table

def writeCoefTable(coef_table, trait, results_path):
	coef_table.to_csv(f"{results_path}/{trait}/{trait}.TURF.master_coefs.tsv", sep="\t", index=False)

def readMasterPheno(trait, iteration, results_path, isQuant):
	master_pheno = pd.read_table(f"{results_path}/{trait}/{trait}.TURF.master_pheno.ITERATION_{iteration-1}.tsv")
	if isQuant:
		master_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + MASTER", data=master_pheno).fit()
	else:
		master_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + MASTER", data=master_pheno).fit()
	return [master_pheno, master_model, master_model.rsquared_adj]

def readMasterSNPs(trait, iteration, results_path):
	with open(f"{results_path}/{trait}/{trait}.TURF.ITERATION_{iteration-1}.snplist") as f:
		master = f.readlines()
	master = list(map(lambda s: s.strip(), master))
	return master

def readExcludedSNPs(trait, iteration, results_path):
	with open(f"{results_path}/traits/{trait}/{trait}.TURF.ITERATION_{iteration-1}.exclude") as f:
		excluded_snps = f.readlines()
	excluded_snps = list(map(lambda s: s.strip(), excluded_snps))
	return excluded_snps

def readSNPList(trait, iteration, results_path):
	with open(f"{results_path}/{trait}/{trait}.TURF.ITERATION_{iteration}.snplist") as f:
		master = f.readlines()
	master = list(map(lambda s: s.strip(), master))
	return master

def readTissuesPartitions(path):
	tissues_partitions = {}
	with open(path,'r') as f:
		for line in f.readlines():
			tissues_partitions[line.split(",")[0]] = line.split(",")[1].strip()
	return tissues_partitions

def writeTissuesPartitions(tissues_partitions, path):
	with open(path,'w') as f:
		for tissue,partition in tissues_partitions.items():
			f.write(f'{tissue},{partition}\n')

def buildTissuePartitionList(master_results):
	tissues = master_results['Tissue'].tolist()
	partitions = master_results['Model'].tolist()
	return [str(a) + "." + str(b) for a,b in zip(tissues, partitions)]

def readLeadTissues(root,trait):
	with open(f"{trait}.lead_tissue","r") as f:
		lead_tissues = f.readlines()
	lead_tissues = list(map(lambda s: s.strip(), lead_tissues))
	return lead_tissues

def runPlinkScore(ancestry, trait, iteration, root):
	plink_path = "/nfs/turbo/boylelab/plink2/plink2"
	ukb_path = f'{root}/data/UKB'
	GWAS_path = f'{root}/data/GWAS'
	results_path = f'{root}/results/traits'

	err = {}
	returncode = {}

	cmd = f"{plink_path} \
	--bed {ukb_path}/{ancestry}/ukb_imp_chrALL_v3.bed \
	--bim {ukb_path}/{ancestry}/ukb_imp_chrALL_v3.bim \
	--fam {ukb_path}/phenos/{ancestry}/{ancestry}.{trait}.fam \
	--score {GWAS_path}/{trait}/{trait}.PLINK.TITR 1 2 3 header no-mean-imputation \
	--extract {results_path}/{trait}/{ancestry}.{trait}.ITERATION_{iteration}.SNPs \
	--keep {ukb_path}/phenos/{ancestry}/{ancestry}.{trait}.sample.IDs \
	--out {results_path}/{trait}/{trait}.TURF.ITERATION_{iteration} --allow-no-sex --write-snplist"

	p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)

	err['score'] = p.communicate()[1]
	returncode['score'] = p.returncode
	return [returncode,err]

def chunks(data, SIZE=10000):
	it = iter(data)
	for i in range(0, len(data), SIZE):
		yield {k:data[k] for k in islice(it, SIZE)}

def build_clumped_snps(chunk, snps_file):
	clumped_snps = set()
	for snps in chunk.values():
		clumped_snps = clumped_snps.union(set(snps).intersection(snps_file))
	return clumped_snps

def remove_clumps(chunk, candidate_snps):
	for index, snps in chunk.items():
		if index in candidate_snps:
			chunk = removekey(chunk, index)
		else:
			if set(candidate_snps).intersection(snps):
				chunk = removekey(chunk, index)
	return chunk

def write_clumps(trait, tissue, iteration, partition, GWAS_path, results_path):
	start = time.time()
	with open(f"{GWAS_path}/{trait}/1KG-clump/{trait}.TURF.ITERATION_{iteration-1}.clumps.pkl","rb") as f:
		clumps = pickle.load(f)

	with open(f"{results_path}/traits/{trait}/{trait}.TURF.ITERATION_{iteration}.snplist","r") as f:
		snps_file = set(line.strip() for line in f)

	index_snps = set(snps_file).intersection(clumps.keys())
	pool = multiprocessing.Pool(10)
	processes = [pool.apply_async(build_clumped_snps, args=(chunk, snps_file)) for chunk in chunks(clumps)]
	result = [p.get() for p in processes]
	clumped_snps = {i for l in result for i in l}
	current_snps = index_snps.union(clumped_snps)

	with open(f"{results_path}/traits/{trait}/{trait}.TURF.ITERATION_{iteration-1}.snplist","r") as f:
		master_snps = [s.strip() for s in f.readlines()]

	candidate_snps = list(current_snps.union(master_snps))

	processes = [pool.apply_async(remove_clumps, args=(chunk, candidate_snps)) for chunk in chunks(clumps)]
	result = [p.get() for p in processes]
	clumps_out = {}
	for d in result:
		clumps_out.update(d)

	with open(f"{GWAS_path}/{trait}/1KG-clump/{trait}.TURF.ITERATION_{iteration}.clumps.pkl","wb") as f:
		pickle.dump(clumps_out,f)

	end = time.time()

	with open(f"{results_path}/traits/{trait}/{trait}.write_blocks.time","a") as f:
		f.write(f"Iteration {iteration}: {(end-start)}\n")

def evaluateFirstIteration(ancestry, trait, iteration, isQuant, root, lead_tissues):
	# Intialize structures
	master = []
	results = []
	added_bin_paths = []
	master_results = pd.DataFrame(columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","SNP Count","Null Adj R2"])

	# Initialize paths
	results_path = f'{root}/results'
	traits_path = f'{results_path}/traits'
	GWAS_path = f'{root}/data/GWAS'
	covariates_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.covariates'
	samples_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.sample.IDs'
	coding_prefix = f'{root}/results/coding/{trait}/{trait}.{ancestry}.coding'

	# Build null model
	null_pheno, null_model, null_adj_R2 = initializeNullModel(covariates_path, samples_path, trait, isQuant)

	# Calculate coding models
	coding_results = pd.DataFrame(calculateCodingModel(coding_prefix, null_pheno, null_adj_R2, isQuant),
		columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","SNP Count","Null Adj R2"])

	# Select best coding model (greatest dR2 over null)
	best_coding_model = coding_results.iloc[coding_results['Delta Adj R2'].idxmax()]
	best_coding_threshold = float(best_coding_model['Threshold'])

	# Add coding model as covariate in model
	coding_pheno, coding_model, coding_adj_R2, coding_d_R2, coding_p, coding_count = addCodingModel(best_coding_threshold, coding_prefix, null_pheno, null_adj_R2, isQuant)
	row = {"Threshold":best_coding_threshold,"Tissue":"NULL","Model":"Coding","Adj R2":coding_adj_R2,"Delta Adj R2":coding_d_R2, 
		"P":coding_p, "SNP Count":coding_count, "Null Adj R2":null_adj_R2}
	master_results = master_results.append(row, ignore_index=True)
	
	# Read tissue/partition pairs
	tissues_partitions = readTissuesPartitions(f"{root}/repo/RegDB-tissue-heritability/PRS/{trait}.partitions")
	# Run comparisions of partitions
	TURF_max = runComparison(trait, tissues_partitions, iteration, coding_model, coding_pheno, coding_adj_R2, added_bin_paths, results_path, isQuant, lead_tissues)
	if len(TURF_max) == 0:
		print('Failed ANOVA!')
		writeMasterResult(master_results, trait, traits_path)
		writeModelTerminate(trait, traits_path)
		quit()
	TURF_max = TURF_max.sort_values(by="Delta Adj R2",ascending=False)
	TURF_max['Threshold'] = pd.to_numeric(TURF_max['Threshold'])
	writeTURFMax(TURF_max, trait, iteration, traits_path)
	TURF_max = TURF_max.sort_values(by="Threshold",ascending=True)
	max_model = TURF_max[TURF_max["Delta Adj R2"] == TURF_max["Delta Adj R2"].max()].iloc[0] #select lowest p-value threshold if tie
	print(max_model)
	addPartitionSNPs(master, trait, max_model.Tissue, max_model.Model, max_model.Threshold, iteration, GWAS_path, results_path)
	writeClumpedList(master, trait, ancestry, iteration, results_path)

	# Run master score
	isError, msg = runPlinkScore(ancestry, trait, iteration, root)
	if isError:
		print(msg)

	# Update null model/R2 and write coefficient table
	master_pheno, master_model, master_adj_R2, master_adj_d_R2 = updateNullModel(trait, iteration, traits_path, coding_pheno, coding_adj_R2, isQuant)
	coef_table = master_model.params.to_frame().transpose()
	writeCoefTable(coef_table, trait, traits_path)
	tissues_partitions[max_model.Tissue] = int(tissues_partitions[max_model.Tissue]) + 1
	
	# Write out SNPs to exclude from subsequent iterations
	write_clumps(trait, max_model.Tissue, iteration, max_model.Model, GWAS_path, results_path)
	snplist = readSNPList(trait, iteration, traits_path)

	# Write results report
	row = {"Threshold":max_model.Threshold,"Tissue":max_model.Tissue,"Model":max_model.Model,
	"Adj R2":master_adj_R2,"Delta Adj R2":master_adj_d_R2,"P":max_model['P'],"SNP Count":len(snplist),"Null Adj R2":null_adj_R2}

	master_results = master_results.append(row, ignore_index=True)
	print(f'Master results after {iteration} iterations for {trait}')
	print(master_results)
	writeMasterResults(master_results, trait, traits_path)
	writeMasterPheno(master_pheno, trait, iteration, traits_path)
	writeTissuesPartitions(tissues_partitions,f"{root}/repo/RegDB-tissue-heritability/PRS/{trait}.partitions")

def evaluateIteration(ancestry, trait, iteration, isQuant, root, lead_tissues):
	results_path = f'{root}/results'
	traits_path = f'{results_path}/traits'
	GWAS_path = f'{root}/data/GWAS'
	covariates_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.covariates'
	samples_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.sample.IDs'
	coding_prefix = f'{root}/results/coding/{trait}/{trait}.coding'

	master_results = readMasterResults(trait, traits_path)
	print(f'Master results after {iteration-1} iterations for {trait}')
	print(master_results)
	master = readMasterSNPs(trait, iteration, traits_path)

	null_pheno, null_model, null_adj_R2 = initializeNullModel(covariates_path, samples_path, trait, isQuant)
	master_pheno, master_model, master_adj_R2 = readMasterPheno(trait, iteration, traits_path, isQuant)
	coef_table = readCoefTable(trait, traits_path)
	added_bin_paths = buildTissuePartitionList(master_results)

	tissues_partitions = readTissuesPartitions(f"{root}/repo/RegDB-tissue-heritability/PRS/{trait}.partitions")

	TURF_max = runComparison(trait, tissues_partitions, iteration, master_model, master_pheno, master_adj_R2, added_bin_paths, results_path, isQuant, lead_tissues)
	if len(TURF_max) == 0:
		print('Exhausted partition list - model terminating')
		writeMasterResults(master_results, trait, traits_path)
		writeModelTerminate(trait, traits_path)
		quit()
	TURF_max = TURF_max.sort_values(by="Delta Adj R2",ascending=False)
	writeTURFMax(TURF_max, trait, iteration, traits_path)
	while TURF_max.shape[0] > 0:
		max_model = TURF_max[TURF_max["Delta Adj R2"] == TURF_max["Delta Adj R2"].max()].iloc[0]
		tmp = master.copy()
		addPartitionSNPs(tmp, trait, max_model.Tissue, max_model.Model, max_model.Threshold, iteration, GWAS_path, results_path)
		if len(tmp) == len(master):
			print(f"Model ({max_model.Tissue} {max_model.Model} {max_model.Threshold}) adds no unique SNPs - continuing")
			TURF_max = TURF_max[~((TURF_max['Tissue'] == max_model.Tissue) & (TURF_max['Threshold'] == max_model.Threshold))]
		else:
			break

	if len(TURF_max) == 0:
		print('Exhausted partition list - model terminating')
		writeMasterResults(master_results, trait, traits_path)
		writeModelTerminate(trait, traits_path)
		quit()
	
	print(f'Comparision complete - Max Model')
	print(max_model)
	addPartitionSNPs(master, trait, max_model.Tissue, max_model.Model, max_model.Threshold, iteration, GWAS_path, results_path)
	writeClumpedList(master, trait, ancestry, iteration, results_path)

	isError, msg = runPlinkScore(ancestry, trait, iteration, root)
	if isError:
		print(msg)

	master_pheno, master_model, master_adj_R2, master_adj_d_R2 = updateNullModel(trait, iteration, traits_path, master_pheno, master_adj_R2, isQuant)
	coef_table = pd.concat([coef_table,master_model.params.to_frame().transpose()])
	writeCoefTable(coef_table, trait, traits_path)
	tissues_partitions[max_model.Tissue] = int(tissues_partitions[max_model.Tissue]) + 1

	write_clumps(trait, max_model.Tissue, iteration, max_model.Model, GWAS_path, results_path)
	snplist = readSNPList(trait, iteration, traits_path)
	added_bin_paths.append(f"{max_model.Tissue}.{max_model.Model}")

	row = {"Threshold":max_model.Threshold,"Tissue":max_model.Tissue,"Model":max_model.Model,
	"Adj R2":master_adj_R2,"Delta Adj R2":master_adj_d_R2,"P":max_model['P'],"SNP Count":len(snplist),"Null Adj R2":null_adj_R2}

	master_results = master_results.append(row, ignore_index=True)
	print(f'Master results after {iteration} iterations for {trait}')
	print(master_results)
	writeMasterResults(master_results, trait, traits_path)
	writeMasterPheno(master_pheno, trait, iteration, traits_path)
	writeTissuesPartitions(tissues_partitions,f"{root}/repo/RegDB-tissue-heritability/PRS/{trait}.partitions")

def main(ancestry, trait, iteration, isQuant, root):

	print(f'Evaluating iteration {iteration} for {trait}')

	lead_tissues = readLeadTissues(root,trait)

	# First iteration, follow first iteration logic
	if iteration == 1:
		evaluateFirstIteration(ancestry, trait, iteration, isQuant, root, lead_tissues)
	# Subsequent iterations
	else:
		evaluateIteration(ancestry, trait, iteration, isQuant, root, lead_tissues)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Optimize TURF PRS')
	parser.add_argument('--ancestry', type=str)
	parser.add_argument('--trait', type=str)
	parser.add_argument('--iteration', type=int)
	parser.add_argument('--isQuant', type=bool)
	parser.add_argument('--root', type=str)

	args = parser.parse_args()
	main(ancestry=args.ancestry, trait=args.trait, iteration=args.iteration, isQuant=args.isQuant, root=args.root)
