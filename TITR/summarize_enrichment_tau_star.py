import pandas as pd
import argparse
import scipy.stats

def parse_log_file(log_path):
	with open(log_path,'r') as f:
		log = f.readlines()
		for line in log:
			if "chi^2" in line:
				M = int(line.split("(")[1].split(" ")[0])
			if "Total Observed scale" in line:
				h2g = float(line.split(":")[1].strip().split(" ")[0])
	return (M, h2g)

def parse_results_file(results_path, M, h2g):
	results_table = pd.read_table(results_path)
	result = results_table.iloc[-1]
	tau_star = calculate_tau_star(M, h2g, result["Category"], result["Coefficient"])
	tau_star_se = calculate_tau_star_se(M, h2g, result["Category"], result["Coefficient_std_error"])
	tau_star_p = calculate_tau_star_p(tau_star, tau_star_se)
	row = {"Category":result["Category"].split("L2_0")[0],"Prop._SNPs":result["Prop._SNPs"],"Prop._h2":result["Prop._h2"],"Prop._h2_std_error":result["Prop._h2_std_error"],
		   "Enrichment":result["Enrichment"],"Enrichment_std_error":result["Enrichment_std_error"],"Enrichment_p":result["Enrichment_p"],"Coefficient":result["Coefficient"],
		   "Coefficient_std_error":result["Coefficient_std_error"],"Coefficient_z-score":result["Coefficient_z-score"],"Coefficient_p":scipy.stats.norm.sf(abs(result["Coefficient_z-score"])),
		   "Tissue_M":M,"Tissue_h2g":h2g,"tau_star":tau_star,"tau_star_se":tau_star_se,"tau_star_p":tau_star_p}
	return row

def get_tissue_sd(tissue):
	tissue_sd_table = pd.read_table("../data/RegDB.tissue.sd")
	tissue_sd = tissue_sd_table.loc[tissue_sd_table['TISSUE'] == tissue]
	return tissue_sd.iloc[0]["SD"]

def calculate_tau_star(M, h2g, tissue, tau):
	sd = get_tissue_sd(tissue)
	return ((M * sd)/h2g) * tau

def calculate_tau_star_se(M, h2g, tissue, tau_se):
	sd = get_tissue_sd(tissue)
	return ((M * sd)/h2g) * tau_se

def calculate_tau_star_p(tau_star, tau_star_se):
	return (1 - scipy.stats.norm.cdf(abs(tau_star), loc = 0, scale = tau_star_se))

def select_lead_tau_star(t):
	return (t.loc[t["tau_star_p"] < (1-(1-0.05)**(1/t.shape[0]))].sort_values(by="tau_star",ascending=False))

def main(path, prefix, trait, iteration):

	TISSUES = ["ADIPOSE_TISSUE","ADRENAL_GLAND","ARTERIAL_BLOOD_VESSEL","BLOOD","BLOOD_VESSEL","BONE_ELEMENT",
		   "BONE_MARROW","BRAIN","BREAST","COLON","CONNECTIVE_TISSUE","EAR","EMBRYO","ENDOCRINE_GLAND",
		   "EPITHELIUM","ESOPHAGUS","EXOCRINE_GLAND","EXTRAEMBRYONIC_COMPONENT","EYE","GONAD",
		   "HEART","IMMUNE_ORGAN","INTESTINE","KIDNEY","LARGE_INTESTINE","LIMB","LIVER","LUNG","LYMPHOID_TISSUE",
		   "LYMPH_NODE","MAMMARY_GLAND","MOUTH","MUSCULATURE_OF_BODY","NERVE","OVARY","PANCREAS","PENIS","PLACENTA",
		   "PROSTATE_GLAND","SKIN_OF_BODY","SKIN_OF_PREPUCE_OF_PENIS","SMALL_INTESTINE","SPINAL_CORD","SPLEEN",
		   "STOMACH","TESTIS","THYMUS","THYROID_GLAND","UTERUS","VAGINA","VASCULATURE"]
	print(f"Starting summary for {prefix}")

	summary = pd.DataFrame(columns=["Category","Prop._SNPs","Prop._h2","Prop._h2_std_error",
									"Enrichment","Enrichment_std_error","Enrichment_p","Coefficient",
									"Coefficient_std_error","Coefficient_z-score","Coefficient_p",
									"Tissue_M","Tissue_h2g","tau_star","tau_star_se","tau_star_p"])

	for tissue in TISSUES:
		print(f"Summarizing tissue {tissue}")
		log_path = f"{path}/{tissue}/{prefix}.{tissue}.log"
		results_path = f"{path}/{tissue}/{prefix}.{tissue}.results"
		M, h2g = parse_log_file(log_path)
		row = parse_results_file(results_path, M, h2g)
		summary = summary.append(row,ignore_index=True)

	out_path = f"{path}/traits/{trait}/{prefix}.ITERATION_{iteration}.summary"
	summary.to_csv(out_path, sep= '\t', index=False)

	lead_tau_star = select_lead_tau_star(summary)
	lead_tau_star["Category"].to_csv(f"{trait}.lead_tissue", index=False, header=False)

	print(f"Summary for {prefix} complete")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--path', type=str, required=True)
	parser.add_argument('--prefix', type=str, required=True)
	parser.add_argument('--trait', type=str, required=True)
	parser.add_argument('--iteration', type=str, required=True)

	args = parser.parse_args()

	main(path=args.path, prefix=args.prefix, trait=args.trait, iteration=args.iteration)
