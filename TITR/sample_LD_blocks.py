import pickle
import argparse
import pandas as pd
import time
from itertools import islice
import multiprocessing

def chunks(data, SIZE=10000):
	it = iter(data)
	for i in range(0, len(data), SIZE):
		yield {k:data[k] for k in islice(it, SIZE)}

def build_clumped_snps(chunk, index_snps, snps_file, snps_P):
	clumped_snps = []
	for index, snps in chunk.items():
		if index in index_snps:
			continue
		snp_intersect = snps_file.intersection(set(snps))
		if len(snp_intersect) == 1:
			clumped_snps.extend(snp_intersect)
		elif len(snp_intersect) > 1:
			clumped_snps.extend([select_min_P(snp_intersect, snps_P)])
	return clumped_snps

def select_min_P(snps, snps_P):
	a = snps_P.loc[snps_P["SNP"].isin(snps)]
	return(a.loc[a["P"] == min(a["P"])].reset_index(drop=True)['SNP'][0])

def main(trait, tissue, iteration, partition, data_path, results_path):
	with open(f"{data_path}/GWAS/{trait}/1KG-clump/{trait}.TURF.ITERATION_{iteration-1}.clumps.pkl","rb") as f:
		clumps = pickle.load(f)

	with open(f"{data_path}/RegDB/{tissue}/1000G_phase3_master_scores_1KG-pruned_quantile_normalized_chrALL.{tissue}.{partition}.SNPs","r") as f:
		snps_file = set(line.strip() for line in f)

	snps_P = pd.read_table(f"{data_path}/GWAS/{trait}/{trait}.PLINK.TITR.SNP.pvalue",sep=" ")

	index_snps = set(snps_file).intersection(clumps.keys())

	pool = multiprocessing.Pool(10)
	start_time = time.perf_counter()
	processes = [pool.apply_async(build_clumped_snps, args=(chunk, index_snps, snps_file, snps_P)) for chunk in chunks(clumps)]
	result = [p.get() for p in processes]
	clumped_snps = {i for l in result for i in l}
	finish_time = time.perf_counter()
	multi_time = finish_time - start_time

	current_snps = index_snps.union(clumped_snps)

	with open(f"{results_path}/traits/{trait}/{trait}.TURF.ITERATION_{iteration-1}.snplist","r") as f:
		master_snps = [s.strip() for s in f.readlines()]

	candidate_snps = list(current_snps.union(master_snps))

	with open(f"{results_path}/{tissue}/{trait}.{tissue}.{partition}.ITERATION_{iteration}.SNPs","w") as f:
		for snp in candidate_snps:
			f.write(snp + "\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--trait', type=str)
	parser.add_argument('--tissue', type=str)
	parser.add_argument('--iteration', type=int)
	parser.add_argument('--partition', type=int)

	args = parser.parse_args()

	data_path = "/path/to/data"
	results_path = "/path/to/results"
	start = time.time()
	main(trait=args.trait, tissue=args.tissue, iteration=args.iteration, partition=args.partition, data_path=data_path, results_path=results_path)
	end = time.time()
	with open(f"{results_path}/{args.tissue}/{args.trait}.{args.tissue}.sample_LD.time","a") as f:
		f.write(f"Iteration {args.iteration}: {(end-start)}\n")

