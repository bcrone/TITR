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
import pickle
import time
import multiprocessing

def calc_pval(model, x, y):
	y_len = len(y)
	beta_hat = [model.intercept_] + model.coef_.tolist()
	x1 = np.column_stack((np.ones(y_len), x))
	st_dev_noise = np.sqrt(np.sum(np.square(y - x1@beta_hat))/(y_len - x1.shape[1]))
	if np.linalg.det(x1.T@x1) == 0:
		return 0
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

def initializeNullModel(covariates_path, samples_path, trait, isQuant):
	null_pheno = pd.read_table(covariates_path, sep='\t', index_col='FID')
	samples = pd.read_table(samples_path, sep=' ', header=None)
	samples.columns = ['FID','IID']
	null_pheno = null_pheno[null_pheno['IID'].isin(samples['IID'])]
	null_pheno = null_pheno[null_pheno[trait].notna()]
	if isQuant:
		null_pheno["norm"] = (null_pheno[trait] - null_pheno[trait].mean()) / null_pheno[trait].std()
		null_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", data=null_pheno).fit()
	else:
		null_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", data=null_pheno).fit()
	null_adj_R2 = null_model.rsquared_adj
	return [null_pheno, null_model, null_adj_R2]

def calculateCodingModel(coding_prefix, null_pheno, null_adj_R2, isQuant):
	results = []
	coding_results = pd.DataFrame(columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","SNP Count","Null Adj R2"])
	coding_files = glob.glob(f"{coding_prefix}.*.sscore")
	for coding_file in coding_files:
		try:
			coding_table = pd.read_table(coding_file, delim_whitespace=True, index_col="IID")
		except OSError as e:
			continue
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

def scoreStandardModel(trait, ancestry, results_path, adj_R2, isQuant, pheno, iteration, threshold=None):
	def get_score_results(file, iteration, adj_R2, isQuant, pheno):
		p = f"0.{file.split('.')[4]}"
		if p == "0.sscore":
			p = 1
		score_table = pd.read_table(file, delim_whitespace=True, index_col='IID')
		if score_table['SCORE1_AVG'].isnull().values.any():
			return {'Model':'Standard', 'Threshold':p, 'Partition':"NA", 'Adj R2':adj_R2,
				'Delta Adj R2':np.nan, 'P':np.nan, 'SNP Count':np.nan, 'Iteration':iteration}
		score_count = max(score_table['ALLELE_CT']) / 2
		score_pheno = pd.merge(pheno, score_table['SCORE1_AVG'], on='IID')
		x = (score_pheno[[
			"Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "CODING", "SCORE1_AVG"]]
			).to_numpy()
		if isQuant:
			y = (score_pheno["norm"]).to_numpy()
		else:
			y = (score_pheno[f"{trait}"]).to_numpy()
		score_model = LinearRegression(n_jobs=-1).fit(x, y)
		score_adj_R2, score_d_R2 = calc_R2(score_model, x, y, adj_R2)
		score_p = calc_pval(score_model, x, y)
		row = {'Model':'Standard', 'Threshold':p, 'Partition':"NA", 'Adj R2':score_adj_R2,
			'Delta Adj R2':score_d_R2, 'P':score_p, 'SNP Count':score_count, 'Iteration':iteration}
		return row
	def get_max_score(trait, ancestry, results_path, adj_R2, isQuant, iteration, threshold=None):
		score_path = f"{results_path}/standard"
		if threshold:
			score_files = glob.glob(f"{score_path}/{trait}.{ancestry}.standard.{threshold}.sscore")
		else:
			score_files = glob.glob(f"{score_path}/{trait}.{ancestry}.standard.*.sscore")
		score_results = pd.DataFrame(
			[get_score_results(file, iteration, adj_R2, isQuant, pheno) for file in score_files],
			columns=["Model","Threshold","Partition","Adj R2","Delta Adj R2","P","SNP Count","Iteration"])
		score_results.dropna(inplace=True)
		score_max = score_results[score_results["Delta Adj R2"]==score_results["Delta Adj R2"].max()]
		return score_max
	score_head = pd.DataFrame(columns=["Model","Threshold","Partition","Adj R2","Delta Adj R2","P","SNP Count","Iteration"])
	score_max = [get_max_score(trait, ancestry, results_path, adj_R2, isQuant, iteration, threshold)]
	return score_head.append(score_max, ignore_index=True)

def scoreModel(model, trait, ancestry, results_path, adj_R2, isQuant, pheno, iteration,threshold=None):  # 1:514.4733331203461
	def get_score_results(file, model, partition, adj_R2, isQuant, pheno, iteration):
		p = f"0.{file.split('.')[5]}"
		if p == "0.sscore":
			p = 1
		score_table = pd.read_table(file, delim_whitespace=True, index_col='IID')
		if score_table['SCORE1_AVG'].isnull().values.any():
			return {'Model':model, 'Thereshold':p, 'Partition':partition, 'Adj R2':adj_R2, 
				'Delta Adj R2':np.nan, 'P':np.nan, 'SNP Count':np.nan, 'Iteration':iteration}
		score_count = max(score_table["ALLELE_CT"])/2
		score_pheno = pd.merge(pheno, score_table["SCORE1_AVG"], on="IID")
		x = (score_pheno[[
			"Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "CODING", "SCORE1_AVG"]]
			).to_numpy()
		if isQuant:
			y = (score_pheno["norm"]).to_numpy()
		else:
			y = (score_pheno[f"{trait}"]).to_numpy()
		score_model = LinearRegression(n_jobs=-1).fit(x, y)
		score_adj_R2, score_d_R2 = calc_R2(score_model, x, y, adj_R2)
		score_p = calc_pval(score_model, x, y)
		row = {'Model':model, 'Threshold':p, 'Partition':partition, 'Adj R2':score_adj_R2, 
			'Delta Adj R2':score_d_R2, 'P':score_p, 'SNP Count':score_count, 'Iteration':iteration}
		return row
	def get_max_score(model, partition, ancestry, results_path, isQuant, adj_R2, iteration, threshold=None):
		score_path = f"{results_path}/{model}"
		if threshold:
			score_files = glob.glob(f"{score_path}/{trait}.{ancestry}.{model}.{partition}.{threshold}.sscore")
		else:
			score_files = glob.glob(f"{score_path}/{trait}.{ancestry}.{model}.{partition}.*.sscore")
		score_results = pd.DataFrame(
			[get_score_results(file, model, partition, adj_R2, isQuant, pheno, iteration) for file in score_files],
			columns=["Model","Threshold","Partition","Adj R2","Delta Adj R2","P","SNP Count","Iteration"])
		score_results.dropna(inplace=True)
		score_max = score_results[score_results['Delta Adj R2']==score_results['Delta Adj R2'].max()]
		return score_max
	score_head = pd.DataFrame(columns=["Model","Threshold","Partition","Adj R2","Delta Adj R2","P","SNP Count","Iteration"])
	if threshold:
		score_max = [get_max_score(model, partition, ancestry, results_path, isQuant, adj_R2, iteration, threshold[partition]) for partition in [10,50,100,200,500]]
	else:
		score_max = [get_max_score(model, partition, ancestry, results_path, isQuant, adj_R2, iteration, None) for partition in [10,50,100,200,500]]
	return score_head.append(score_max, ignore_index=True)

def main(ancestry, trait, isQuant, root):
	results_path = f'{root}/results/single_tissue/{trait}'
	standard_path = f'{results_path}/standard'
	IMPACT_path = f'{results_path}/IMPACT'
	SURF_path = f'{results_path}/SURF'
	TURF_path = f'{results_path}/TURF'

	covariates_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.covariates'
	samples_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.sample.IDs'
	coding_prefix = f'{root}/results/coding/{trait}/{trait}.{ancestry}.coding'

	master_results = pd.DataFrame(columns=["Model","Threshold","Partition","Adj R2","Delta Adj R2","P","SNP Count","Iteration"])

	null_pheno, null_model, null_adj_R2 = initializeNullModel(covariates_path, samples_path, trait, isQuant)
	results = calculateCodingModel(coding_prefix, null_pheno, null_adj_R2, isQuant)
	coding_results = pd.DataFrame(results,columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","SNP Count","Null Adj R2"])
	best_coding_model = coding_results.iloc[coding_results['Delta Adj R2'].idxmax()]
	best_coding_threshold = best_coding_model['Threshold']
	coding_pheno, coding_model, coding_adj_R2, coding_d_R2, coding_p, coding_count = addCodingModel(best_coding_threshold, coding_prefix, null_pheno, null_adj_R2, isQuant)
	row = {"Model":"Coding", "Threshold":best_coding_threshold, "Partition":np.nan, "Adj R2":coding_adj_R2,
		   "Delta Adj R2":coding_d_R2, "P":coding_p, "SNP Count":coding_count, "Iteration":np.nan}
	master_results = master_results.append(row, ignore_index=True)

	standard_max = scoreStandardModel(trait, ancestry, results_path, coding_adj_R2, isQuant, coding_pheno, 0, None)
	IMPACT_max = scoreModel("IMPACT", trait, ancestry, results_path, coding_adj_R2, isQuant, coding_pheno, 0, None)
	SURF_max = scoreModel("SURF", trait, ancestry, results_path, coding_adj_R2, isQuant, coding_pheno, 0, None)
	TURF_max = scoreModel("TURF", trait, ancestry, results_path, coding_adj_R2, isQuant, coding_pheno, 0, None)

	standard_threshold = standard_max["Threshold"].values[0]
	IMPACT_threshold = {partition:IMPACT_max.loc[IMPACT_max["Partition"]==partition]["Threshold"].values[0] for partition in [10,50,100,200,500]}
	SURF_threshold = {partition:SURF_max.loc[SURF_max["Partition"]==partition]["Threshold"].values[0] for partition in [10,50,100,200,500]}
	TURF_threshold = {partition:TURF_max.loc[TURF_max["Partition"]==partition]["Threshold"].values[0] for partition in [10,50,100,200,500]}

	for i in range(1,1001):
		print(trait, i)
		sample_pheno = coding_pheno.sample(coding_pheno.shape[0], replace=True, random_state=i)
		standard_max = standard_max.append(scoreStandardModel(trait, ancestry, results_path, coding_adj_R2, isQuant, sample_pheno, i, standard_threshold), ignore_index=True)
		IMPACT_max = IMPACT_max.append(scoreModel("IMPACT", trait, ancestry, results_path, coding_adj_R2, isQuant, sample_pheno, i, IMPACT_threshold), ignore_index=True)
		SURF_max = SURF_max.append(scoreModel("SURF", trait, ancestry, results_path, coding_adj_R2, isQuant, sample_pheno, i, SURF_threshold), ignore_index=True)
		TURF_max = TURF_max.append(scoreModel("TURF", trait, ancestry, results_path, coding_adj_R2, isQuant, sample_pheno, i, TURF_threshold), ignore_index=True)

	standard_max.to_csv(f"{standard_path}/{trait}.{ancestry}.standard.bootstrap",sep='\t',header=True,index=False)
	IMPACT_max.to_csv(f"{IMPACT_path}/{trait}.{ancestry}.IMPACT.bootstrap",sep='\t',header=True,index=False)
	SURF_max.to_csv(f"{SURF_path}/{trait}.{ancestry}.SURF.bootstrap",sep='\t',header=True,index=False)
	TURF_max.to_csv(f"{TURF_path}/{trait}.{ancestry}.TURF.bootstrap",sep='\t',header=True,index=False)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Assess single tissue model - standard/IMPACT/SURF/TURF')
	parser.add_argument('--ancestry', type=str)
	parser.add_argument('--trait', type=str)
	parser.add_argument('--isQuant', type=bool)
	parser.add_argument('--root', type=str)

	args = parser.parse_args()
	main(ancestry=args.ancestry, trait=args.trait, isQuant=args.isQuant, root=args.root)