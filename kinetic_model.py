# Euler's integration technique
# Built for kinetic model
# Generic kinetic model 
# Edit, writing skip

import sys
import math

in_file=open(sys.argv[1], 'r')
out_file=open("kinetic_model.out", 'w')
temp_file=open("temp_km.out", 'w')
log_file=open("log_kinetic.log",'w')
log_file.writelines(f"input file {sys.argv[1]} \n")
n_comp=int(in_file.readline().split()[0])
log_file.writelines(f"Total number of componets {n_comp} \n")
log_file.writelines(f"initial concentration \n")
conc=[]
for i in range(n_comp):
	conc.append(float(in_file.readline().split()[0]))
	log_file.writelines(f"Comp {i+1} -> {conc[i]} \n")
in_file.readline()
n_rxn =int(in_file.readline().split()[0])
log_file.writelines(f"Total number of reactions {n_rxn} \n")
log_file.writelines(f"Rate constants \n")
kconst=[]
for i in range(n_rxn):
	line=in_file.readline().split()
	kconst.append(float(line[0]))
	log_file.writelines(f"Rxn {i+1} -> {kconst[i]} \n")
in_file.readline()
co=[[0]*n_comp for i in range(n_rxn)]
log_file.writelines(f"Role of components in reaction (reactant -1, product 1, none 0)\n")
for i in range(n_rxn):
	line=in_file.readline().split()
	log_file.writelines(f"Reaction {i+1} : \n")
	for j in range(n_comp):
		co[i][j]=float(line[j])
		log_file.writelines(f"Comp {j+1}: {co[i][j]} \t")
	log_file.writelines(f"\n")
in_file.readline()
power=[[0]*n_comp for i in range(n_rxn)]

log_file.writelines(f"\n Coefficient of components in reaction \n")
for i in range(n_rxn):
	line=in_file.readline().split()
	log_file.writelines(f"Reaction {i+1} : \n")
	for j in range(n_comp):
		power[i][j]=float(line[j])
		log_file.writelines(f"Comp{j+1}: {power[i][j]} \t")
	log_file.writelines(f"\n")
in_file.readline()
dt=float(in_file.readline().split()[0])
tot=float(in_file.readline().split()[0])
skip_level=int(in_file.readline().split()[0])
del_rxn=[[0.0]*n_rxn]
log_file.writelines(f" Time step  : {dt} \n")
log_file.writelines(f" Total time : {tot} in steps ={int(tot/dt)} \n")
log_file.writelines(f" Skip level : Output after {skip_level} iterations. \n")



out_file.writelines(f"# kinetic model for {n_rxn} reaction with {n_comp} \n")
#out_file.writelines(f"{t} {conc[:]} \n")
t=0
conct=[0.0]*n_comp
#out_file.writelines(f"{t}\t\t") # initial concentration
#for j in range(n_comp):
#	out_file.writelines(f"{conc[j]}\t")
#	out_file.writelines(f"\n")

for i in range(int(tot/dt)):
	temp_file.writelines(f"{t}\t \t")
	for j in range(n_comp):
		temp_file.writelines(f"{conc[j]}\t")
	temp_file.writelines(f"\n")
	del_rxn=[0.0]*n_rxn
	for k in range(n_rxn):
		del_rxn[k] =kconst[k]
		for j in range(n_comp):
			if(co[k][j]<0.0) :
				del_rxn[k] *= pow(conc[j],power[k][j])
				#print(k, j, conc[j], power[k][j], del_rxn[k])
		if(del_rxn[k] < 0.0):
			print("Reduce the time step")
			exit()
	for j in range(n_comp):
		conct[j]=conc[j]
		for k in range(n_rxn):
			if(power[k][j] != 0.0):
				conct[j] += dt* co[k][j]*del_rxn[k]/power[k][j]

	if(i%skip_level == 0) :
		out_file.writelines(f"{t}\t\t")
		for j in range(n_comp):
			out_file.writelines(f"{conc[j]}\t")
		out_file.writelines(f"\n")
		temp_file.close()
		temp_file=open("temp_km.out", 'w')
	conc[:] =conct[:]
	t +=dt
out_file.writelines(f"# END")
out_file.close()
