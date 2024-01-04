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
	null_pheno = null_pheno[~null_pheno['IID'].isin(samples['IID'])]
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

def addCodingModel(coding_file, null_pheno, null_adj_R2, isQuant):
	try:
		coding_table = pd.read_table(coding_file, delim_whitespace=True, index_col="IID")
	except OSError as e:
		return
	coding_count = max(coding_table['ALLELE_CT'])
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

def evaluateIteration(ancestry, trait, iteration, clump, coding_threshold, isQuant, root):
	covariates_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.covariates'
	samples_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.NA'
	coding_file = f'{root}/results/coding/{trait}/{trait}.{ancestry}.coding.{coding_threshold}.sscore'
	iteration_file = f'{root}/results/traits/{trait}/{trait}.{ancestry}.ITERATION_{iteration}.sscore'
	null_pheno, null_model, null_adj_R2 = initializeNullModel(covariates_path, samples_path, trait, isQuant)
	coding_pheno, coding_model, coding_adj_R2, coding_d_R2, coding_p, coding_count = addCodingModel(coding_file, null_pheno, null_adj_R2, isQuant)
	iteration_table = pd.read_table(iteration_file, delim_whitespace=True, index_col='IID')
	iteration_count = int(max(iteration_table["ALLELE_CT"])/2)
	iteration_pheno = pd.merge(coding_pheno, iteration_table["SCORE1_AVG"], on="IID")
	iteration_pheno.rename(columns={'SCORE1_AVG':'SCORE'}, inplace=True)
	if isQuant:
		iteration_model = ols("norm ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + SCORE", data=iteration_pheno).fit()
	else:
		iteration_model = ols(f"{trait} ~ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + CODING + SCORE", data=iteration_pheno).fit()
	iteration_adj_R2 = iteration_model.rsquared_adj
	iteration_adj_d_R2 = iteration_adj_R2 - null_adj_R2
	return [iteration_pheno, iteration_model, iteration_adj_R2, iteration_adj_d_R2, iteration_count, null_adj_R2]

def main(ancestry, trait, iteration, clump, isQuant, root):
	covariates_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.covariates'
	samples_path = f'{root}/data/UKB/phenos/{ancestry}/{ancestry}.{trait}.NA'
	coding_prefix = f'{root}/results/coding/{trait}/{trait}.{ancestry}.coding'
	null_pheno, null_model, null_adj_R2 = initializeNullModel(covariates_path, samples_path, trait, isQuant)

	optimized_results = pd.read_table(f'{root}/results/traits/{trait}/{trait}.TURF.master_results.tsv')
	results = calculateCodingModel(coding_prefix, null_pheno, null_adj_R2, isQuant)
	coding_results = pd.DataFrame(results,columns=["Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","SNP Count","Null Adj R2"])
	best_coding_model = coding_results.iloc[coding_results['Delta Adj R2'].idxmax()]
	coding_threshold = best_coding_model['Threshold']
	#coding_threshold = optimized_results.loc[optimized_results['Model']=='Coding']['Threshold'].values[0]

	master_results = pd.DataFrame(columns=["Iteration","Threshold","Tissue","Model","Adj R2","Delta Adj R2","P","SNP Count","Null Adj R2"])

	row = {'Iteration':0, 'Threshold':'Null', 'Tissue':'Null', 'Model':0, 'Adj R2':null_adj_R2, 'Delta Adj R2':0, "P":0, "SNP Count":0, "Null Adj R2":null_adj_R2}
	master_results = master_results.append(row, ignore_index=True)

	iterations = range(1,iteration+1)

	for i in iterations:
		tissue = optimized_results.iloc[i]['Tissue']
		model = optimized_results.iloc[i]['Model']
		threshold = optimized_results.iloc[i]['Threshold']
		print(tissue)
		iteration_pheno, iteration_model, iteration_adj_R2, iteration_adj_d_R2, iteration_count, null_adj_R2 = evaluateIteration(ancestry, trait, i, clump, coding_threshold, isQuant, root)
		row = {'Iteration':i, 'Threshold':threshold, 'Tissue':tissue, 'Model':model, 'Adj R2':iteration_adj_R2, 'Delta Adj R2': iteration_adj_d_R2, 'P':iteration_model.pvalues['SCORE'], 'SNP Count':iteration_count, 'Null Adj R2':null_adj_R2}
		master_results = master_results.append(row, ignore_index=True)

	master_results.to_csv(f"{root}/results/traits/{trait}/{trait}.{ancestry}.master_results.tsv",sep="\t",index=False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Optimize TURF PRS')
	parser.add_argument('--ancestry', type=str)
	parser.add_argument('--trait', type=str)
	parser.add_argument('--iteration', type=int)
	parser.add_argument('--clump', type=float)
	parser.add_argument('--isQuant', type=bool)
	parser.add_argument('--root', type=str)

	args = parser.parse_args()
	main(ancestry=args.ancestry, trait=args.trait, iteration=args.iteration, clump=args.clump, isQuant=args.isQuant, root=args.root)
