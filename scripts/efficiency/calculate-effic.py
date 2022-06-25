




import csv
import argparse
import os

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--atp')
parser.add_argument('--carbon')
parser.add_argument('--tempid', help='04C or 15C or 20C')
parser.add_argument('--datadir', help='')
args = parser.parse_args()


atp_dict = {}
with open(args.atp) as f:
    for row in csv.reader(f, delimiter='\t'):
        atp_dict[row[0]] = float(row[2])

carbon_dict = {}
with open(args.carbon) as f:
    for row in csv.reader(f, delimiter='\t'):
        carbon_dict[row[0]] = float(row[1])

def replace_dir(x):
    x = x.replace('_forward', '')
    x = x.replace('_reverse', '')
    return x

CUE_l = []
prod_l = []
effic_l = []
for f in os.listdir('{}'.format(args.datadir)):
    if args.tempid in f:
        atp_fluxes = {}
        carbon_fluxes = {}
        biomass_flux = None
        carbon_uptake = None
        for row in csv.reader(open('{}'.format(args.datadir)+f), delimiter='\t'):
            if row[0] == 'Flux':
                if row[1] == 'ATPM':
                    atpm_flux = float(row[2])
                if row[1] == 'Biomass_WP2':
                    biomass_flux = float(row[2])
                if replace_dir(row[1]) in carbon_dict.keys():
                    if row[1] == 'EX_cpd_acgam[e]':
                        carbon_uptake = float(row[2])
                    carbon_fluxes[row[1]] = float(row[2])*carbon_dict[row[1]]
                if replace_dir(row[1]) in atp_dict.keys():
                    atp_fluxes[row[1]] = float(row[2])

        #calculate CUE
        carbon_prod = 0
        carbon_cons = 0
        for key, value in carbon_fluxes.items():
            # print(key, value)
            if value < 0:
                carbon_cons += float(value)*-1
            elif value > 0:
                carbon_prod += float(value)
        #print(carbon_cons, carbon_prod)
        CUE = (carbon_cons - carbon_prod) / carbon_cons
        # print('CUE\t{}'.format(str(CUE)))
        CUE_l.append(CUE)
        # print(atp_fluxes)
        # Fix issue with resolving directionality and ATP value. clean up this part of code.
        atp_prod = 0
        atp_con = 0
        for key, value in atp_fluxes.items():
            atp_flux = (value / carbon_uptake*-1) * atp_dict[replace_dir(key)]
            atp_fluxes[key] = atp_flux

            #print(key, value, atp_dict[key], carbon_uptake, atp_flux)
            if atp_flux > 0:
                atp_prod += atp_flux
            if atp_flux < 0:
                atp_con += value
        # print(atp_con, biomass_flux)
        atp_con_normalized = atp_con/biomass_flux

        #Calculate ATP production
        # print(atp_prod, atp_con_normalized)
        prod_l.append(atp_prod)
        effic_l.append(atp_con_normalized)

# def get_mean_std(l):
#     mean = sum(l) / len(l)
#     variance = sum([((x - mean) ** 2) for x in l]) / len(l)
#     res = variance ** 0.5
#     return mean, res

# print('CUE',get_mean_std(CUE_l))
# print('Prod',get_mean_std(prod_l))
# print('cons',get_mean_std(effic_l))

for i in CUE_l:
    print('{}\tCUE\t{}'.format(args.tempid, i))
for i in prod_l:
    print('{}\tATPP\t{}'.format(args.tempid, i))
for i in effic_l:
    print('{}\tATPC\t{}'.format(args.tempid, i))
