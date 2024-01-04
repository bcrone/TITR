import argparse
import pandas
import pickle

def build_1KG_clumps(root, trait):
	GWAS_path = f'{root}/data/GWAS'

	with open(f'{GWAS_path}/{trait}/1KG-clump/{trait}.r_0.2.clumped','r') as f:
		r02_file = f.readlines()

	r02_file = r02_file[1:]

	r02_clumps = {}
	for i in range(len(r02_file)):
		try:
			index = r02_file[i].split()[2]
			if index not in r02_clumps.keys():
				r02_clumps[index] = [s.replace('(1)','') for s in r02_file[i].split()[11].split(',')]
		except:
			continue

	with open(f'{GWAS_path}/{trait}/1KG-clump/{trait}.master_clumps.r_0.2.pkl','wb') as f:
		pickle.dump(r02_clumps,f)

	with open(f'{GWAS_path}/{trait}/1KG-clump/{trait}.r_0.8.clumped','r') as f:
		r08_file = f.readlines()

	r08_file = r08_file[1:]

	r08_clumps = {}

	for i in range(len(r08_file)):
		try:
			index = r08_file[i].split()[2]
			if index not in r08_clumps.keys():
				r08_clumps[index] = [s.replace('(1)','') for s in r08_file[i].split()[11].split(',')]
		except:
			continue

	with open(f'{GWAS_path}/{trait}/1KG-clump/{trait}.master_clumps.r_0.8.pkl','wb') as f:
		pickle.dump(r08_clumps,f)

	with open(f'{GWAS_path}/{trait}/1KG-clump/{trait}.r_1.clumped','r') as f:
		r1_file = f.readlines()

	r1_file = r1_file[1:]

	r1_clumps = {}

	for i in range(len(r1_file)):
		try:
			index = r1_file[i].split()[2]
			if index not in r1_clumps.keys():
				r1_clumps[index] = [s.replace('(1)','') for s in r1_file[i].split()[11].split(',')]
		except:
			continue

	with open(f'{GWAS_path}/{trait}/1KG-clump/{trait}.master_clumps.r_1.pkl','wb') as f:
		pickle.dump(r1_clumps,f)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Build 1KG clumps for r=0.2 and r=0.8')
	parser.add_argument('--root', type=str)
	parser.add_argument('--trait', type=str)
	args = parser.parse_args()
	build_1KG_clumps(args.root, args.trait)