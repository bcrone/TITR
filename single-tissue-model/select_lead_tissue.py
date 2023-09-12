import argparse
import pandas as pd
import glob

def merge_tau_star_results(trait, path):
	summary = pd.DataFrame(columns=["Category","Prop._SNPs","Prop._h2","Prop._h2_std_error",
									"Enrichment","Enrichment_std_error","Enrichment_p","Coefficient",
									"Coefficient_std_error","Coefficient_z-score","Coefficient_p",
									"Tissue_M","Tissue_h2g","tau_star","tau_star_se","tau_star_p"])
	results = glob.glob(f"{path}/{trait}/{trait}.*.tau_star.results")
	for file in results:
		result = pd.read_table(file)
		summary = summary.append(result, ignore_index=True)
	return summary

def select_lead_tau_star(t):
	return (t.loc[t["tau_star_p"] < (1-(1-0.05)**(1/t.shape[0]))].sort_values(by="tau_star",ascending=False)).iloc[0]

def main(trait, path):
	summary = merge_tau_star_results(trait, path)
	summary = summary.sort_values(by="Category")
	summary.to_csv(f'{path}/{trait}/{trait}.tau_star.summary', sep='\t', index=False)

	lead_tau_star = select_lead_tau_star(summary)
	lead_tau_star["Category"].to_csv(f"{path}/{trait}/{trait}.lead_tissue", index=False, header=False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--trait', type=str, required=True)
	parser.add_argument('--path', type=str, required=True)
	
	args = parser.parse_args()

	main(trait=args.trait, path=args.path)