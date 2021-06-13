import csv
from decimal import *
from collections import defaultdict, Counter
import argparse




parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--four', type=str,
					help='4C simulation file')
parser.add_argument('--fifteen', type=str,
					help='15C simulation file')
parser.add_argument('--twenty', type=str,
					help='20C simulation file')
parser.add_argument('--out', type=str,
					help='Output files prefix')
args = parser.parse_args()

def is_number(x):
	try:
		float(x)
		return True
	except:
		return False

def compare_differences_absMin(con1, con2, con1_name, con2_name):
	if is_number(con1) and is_number(con2):
		con1 = abs(round(Decimal(con1), 6))
		con2 = abs(round(Decimal(con2), 6))
		if con1 == 0 and con2 == 0:
			comp = 'fixed zero'
		elif con1 == 0 or con2 == 0:
			comp = 'one zero'
		elif con1 == con2:
			comp = 'fixed non-zero'
		else:
			comp = 'diff values'

		if con1 < con2:
			return '{} < {}, {}'.format(con1_name, con2_name, comp)
		elif con2 < con1:
			return '{} < {}, {}'.format(con2_name, con1_name, comp)
		elif con1 == con2:
			return '{} == {}, {}'.format(con1_name, con2_name, comp)
	else:
		return 'No Comparison'


def compare_differences(con1l, con1u, con1r, con2l, con2u, con2r, con1_name, con2_name):
	if is_number(con1l) and is_number(con2l):
		# lower comparison
		if con1l == con2l:
			l = '{}L == {}L'.format(con1_name, con2_name)
		if con1l < con2l:
			l = '{}L < {}L'.format(con1_name, con2_name)
		elif con2l < con1l:
			l = '{}L < {}L'.format(con2_name, con1_name)

		# upper comparison
		if con1u == con2u:
			u = '{}U == {}U'.format(con1_name, con2_name)
		if con1u < con2u:
			u = '{}U < {}U'.format(con1_name, con2_name)
		elif con2u < con1u:
			u = '{}U < {}U'.format(con2_name, con1_name)

		# range comparison
		if con1r == con2r:
			r = '{}R == {}R'.format(con1_name, con2_name)
		if con1r < con2r:
			r = '{}R < {}R'.format(con1_name, con2_name)
		elif con2r < con1r:
			r = '{}R < {}R'.format(con2_name, con1_name)
		comparison = ' | '.join([l, u, r])
		return comparison
	else:
		return 'No Comparison'

def absolute_minimum(l, u):
	if l < 0 < u:
		return 0
	elif l < 0 and u < 0:
		return max([l, u])
	elif l > 0 and u > 0:
		return min([l, u])
	elif l == 0 or u == 0:
		return 0

def combine_fluxes(f_l, f_u, r_l, r_u):
	# print(f_l, f_u, r_l, r_u)
	if f_l > 0:
		if r_l and r_u != 0:
			print(f_l, f_u, r_l, r_u)
			raise ValueError
			quit()
		else:
			return f_l, f_u
	elif r_l > 0:
		if f_l and f_u != 0:
			print(f_l, f_u, r_l, r_u)
			raise ValueError
			quit()
		else:
			return -r_u, -r_l
	elif f_l == 0.0:
		if r_l == 0.0:
			return -r_u, f_u


def combine_drg(f_l, f_u, r_l, r_u):
	# print(f_l, f_u, r_l, r_u)
	if f_u != -r_l:
		raise ValueError
		quit()
	elif f_l != -r_u:
		raise ValueError
		quit()
	else:
		return f_l, f_u

file_list = [(args.four, 4),
			 (args.fifteen, 15),
			 (args.twenty, 20)]

var_list_flux = set()
var_list_drg = set()
var_list_cpd = set()

flux_entry_dict = defaultdict(dict)
drg_entry_dict = defaultdict(dict)
cpd_entry_dict = defaultdict(dict)

reversed_reactions = set()
reverse_orig_reactions = set()
for (f, t) in file_list:
	for row in csv.reader(open(f, mode='rU'), delimiter='\t'):
		if row[0] == 'Flux':
			if '_forward' in row[1]:
				reversed_reactions.add(row[1])
				reverse_orig_reactions.add(row[1][:-8])
			if '_reverse' in row[1]:
				reversed_reactions.add(row[1])
				reverse_orig_reactions.add(row[1][:-8])

for (f, t) in file_list:
	reversed_flux_dict = {}
	reversed_drg_dict = {}
	for row in csv.reader(open(f, mode='rU'), delimiter='\t'):
		if row[0] == 'Flux':
			if row[1] not in reversed_reactions:
				var_list_flux.add(row[1])
				flux_entry_dict[t][row[1]] = (round(float(row[2]), 8), round(float(row[3]), 8))
			else:
				reversed_flux_dict[row[1]] = (round(float(row[2]), 8), round(float(row[3]), 8))

		if row[0] == 'DGR':
			if row[1] not in reversed_reactions:
				var_list_drg.add(row[1])
				drg_entry_dict[t][row[1]] = (row[2], row[3])
			else:
				reversed_drg_dict[row[1]] = (row[2], row[3])

		if row[0] == 'CONC':
			var_list_cpd.add(row[1])
			cpd_entry_dict[t][row[1]] = (row[2], row[3])

	for reaction in reversed_reactions:
		if reaction in reversed_flux_dict.keys():
			new_reaction = reaction[:-8]
			(for_flux_l, for_flux_u) = reversed_flux_dict['{}_forward'.format(new_reaction)]
			(rev_flux_l, rev_flux_u) = reversed_flux_dict['{}_reverse'.format(new_reaction)]
			# print(reaction, for_flux_l, for_flux_u, rev_flux_l, rev_flux_u)
			new_l, new_u = combine_fluxes(for_flux_l, for_flux_u, rev_flux_l, rev_flux_u)
			flux_entry_dict[t][new_reaction] = (new_l, new_u)
			var_list_flux.add(new_reaction)
			var_list_drg.add(new_reaction)

			(for_drg_l, for_drg_u) = reversed_drg_dict['{}_forward'.format(new_reaction)]
			(rev_drg_l, rev_drg_u) = reversed_drg_dict['{}_reverse'.format(new_reaction)]

			new_drg_l, new_drg_u = combine_drg(round(float(for_drg_l), 4), round(float(for_drg_u), 4), round(float(rev_drg_l), 4), round(float(rev_drg_u), 4))
			drg_entry_dict[t][new_reaction] = (new_drg_l, new_drg_u)


cpd_15_20 = []
cpd_04_20 = []
cpd_04_15 = []


compound_file = open('{}_Compund_Variability.tsv'.format(args.out), mode='w')
compound_file.write('Compound_Variability\tcompound\tlower_4\tupper_4\trange_4\tlower_15\tupper_15\trange_15\tlower_20\tupper_20\trange_20\tfour_fifteen\tfour_twenty\tfifteen_twenty\n')
for compound in sorted(var_list_cpd):
	try:
		(lower_4, upper_4) = cpd_entry_dict[4][compound]
		range_4 = float(upper_4) - float(lower_4)
	except KeyError:
		lower_4, upper_4 = 'NA', 'NA'
		range_4 =  'NA'

	try:
		(lower_15, upper_15) = cpd_entry_dict[15][compound]
		range_15 = float(upper_15) - float(lower_15)
	except KeyError:
		lower_15, upper_15 = 'NA', 'NA'
		range_15 = 'NA'

	try:
		(lower_20, upper_20) = cpd_entry_dict[20][compound]
		range_20 = float(upper_20) - float(lower_20)
	except KeyError:
		lower_20, upper_20 = 'NA', 'NA'
		range_20 =  'NA'

	fifteen_twenty = compare_differences(lower_15, upper_15, range_15, lower_20, upper_20, range_20, '15C', '20C')
	four_twenty = compare_differences(lower_4, upper_4, range_4, lower_20, upper_20, range_20, '04C', '20C')
	four_fifteen = compare_differences(lower_4, upper_4, range_4, lower_15, upper_15, range_15, '04C', '15C')
	cpd_15_20.append(fifteen_twenty)
	cpd_04_20.append(four_twenty)
	cpd_04_15.append(four_fifteen)

	compound_file.write('Compound_Variability\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(compound, lower_4, upper_4, range_4, lower_15, upper_15, range_15,
																			lower_20, upper_20, range_20, four_fifteen, four_twenty, fifteen_twenty))

compound_file.close()

flux_file = open('{}_Flux_Variability.tsv'.format(args.out), mode='w')

flux_15_20 = []
flux_04_20 = []
flux_04_15 = []
flux_file.write('Flux_Variability\tID\tlower_4\tupper_4\trange_4\tlower_15\tupper_15\trange_15\tlower_20\tupper_20\trange_20\tfour_fifteen\tfour_twenty\tfifteen_twenty\n')
for reaction in sorted(var_list_flux):
	try:
		(lower_4, upper_4) = flux_entry_dict[4][reaction]
		range_4 = float(upper_4) - float(lower_4)
	except (KeyError, ValueError):
		lower_4, upper_4 = 'NA', 'NA'
		range_4 =  'NA'

	try:
		(lower_15, upper_15) = flux_entry_dict[15][reaction]
		range_15 = float(upper_15) - float(lower_15)
	except (KeyError, ValueError):
		lower_15, upper_15 = 'NA', 'NA'
		range_15 = 'NA'

	try:
		(lower_20, upper_20) = flux_entry_dict[20][reaction]
		range_20 = float(upper_20) - float(lower_20)
	except (KeyError, ValueError):
		lower_20, upper_20 = 'NA', 'NA'
		range_20 =  'NA'

	fifteen_twenty = compare_differences(lower_15, upper_15, range_15, lower_20, upper_20, range_20, '15C', '20C')
	four_twenty = compare_differences(lower_4, upper_4, range_4, lower_20, upper_20, range_20, '04C', '20C')
	four_fifteen = compare_differences(lower_4, upper_4, range_4, lower_15, upper_15, range_15, '04C', '15C')

	flux_15_20.append(fifteen_twenty)
	flux_04_20.append(four_twenty)
	flux_04_15.append(four_fifteen)

	flux_file.write('Flux_Variability\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(reaction, lower_4, upper_4, range_4, lower_15, upper_15, range_15,
																			lower_20, upper_20, range_20, four_fifteen, four_twenty, fifteen_twenty))
flux_file.close()

deltaG_15_20 = []
deltaG_04_20 = []
deltaG_04_15 = []

deltaG_file = open('{}_deltaG_Variability.tsv'.format(args.out), mode='w')

deltaG_file.write('DeltaG_Variability\tID\tlower_4\tupper_4\trange_4\tlower_15\tupper_15\trange_15\tlower_20\tupper_20\trange_20\tfour_fifteen\tfour_twenty\tfifteen_twenty\n')
for reaction in sorted(var_list_drg):
	try:
		(lower_4, upper_4) = drg_entry_dict[4][reaction]
		range_4 = float(upper_4) - float(lower_4)
	except (KeyError, ValueError):
		lower_4, upper_4 = 'NA', 'NA'
		range_4 =  'NA'

	try:
		(lower_15, upper_15) = drg_entry_dict[15][reaction]
		range_15 = float(upper_15) - float(lower_15)
	except (KeyError, ValueError):
		lower_15, upper_15 = 'NA', 'NA'
		range_15 = 'NA'

	try:
		(lower_20, upper_20) = drg_entry_dict[20][reaction]
		range_20 = float(upper_20) - float(lower_20)
	except (KeyError, ValueError):
		lower_20, upper_20 = 'NA', 'NA'
		range_20 =  'NA'

	fifteen_twenty = compare_differences(lower_15, upper_15, range_15, lower_20, upper_20, range_20, '15C', '20C')
	four_twenty = compare_differences(lower_4, upper_4, range_4, lower_20, upper_20, range_20, '04C', '20C')
	four_fifteen = compare_differences(lower_4, upper_4, range_4, lower_15, upper_15, range_15, '04C', '15C')

	deltaG_15_20.append(fifteen_twenty)
	deltaG_04_20.append(four_twenty)
	deltaG_04_15.append(four_fifteen)

	deltaG_file.write('DeltaG_Variability\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(reaction, lower_4, upper_4, range_4, lower_15, upper_15, range_15,
																			lower_20, upper_20, range_20, four_fifteen, four_twenty, fifteen_twenty))

deltaG_file.close()
(biomass_4, _) = flux_entry_dict[4]['Biomass_WP2']
(biomass_15, _) = flux_entry_dict[15]['Biomass_WP2']
(biomass_20, _) = flux_entry_dict[20]['Biomass_WP2']


absMin_15_20 = []
absMin_04_20 = []
absMin_04_15 = []

absMin_file = open('{}_absMin_Variability.tsv'.format(args.out), mode='w')

absMin_file.write('AbsMin_Comparisons\tID\tlower_4\tupper_4\tabsMin_4\tlower_15\tupper_15\tabsMin_15\tlower_20\tupper_20\tabsMin_20\tfour_fifteen\tfour_twenty\tfifteen_twenty\n')
for reaction in sorted(var_list_flux):
		try:
			(lower_4, upper_4) = flux_entry_dict[4][reaction]
			normalized_4l, normalized_4u = round(Decimal(lower_4 / biomass_4), 6), round(Decimal(upper_4 / biomass_4), 6)
			absMin4 = absolute_minimum(normalized_4l, normalized_4u)
		except (KeyError, ValueError):
			lower_4, upper_4 = 'NA', 'NA'
			normalized_4l, normalized_4u =  'NA', 'NA'
			absMin4 = 'NA'

		try:
			(lower_15, upper_15) = flux_entry_dict[15][reaction]
			normalized_15l, normalized_15u = round(Decimal(lower_15 / biomass_15), 6), round(Decimal(upper_15 / biomass_15), 6)
			absMin15 = absolute_minimum(normalized_15l, normalized_15u)
		except (KeyError, ValueError):
			lower_15, upper_15 = 'NA', 'NA'
			normalized_15l, normalized_15u =  'NA', 'NA'
			absMin15 = 'NA'

		try:
			(lower_20, upper_20) = flux_entry_dict[20][reaction]
			normalized_20l, normalized_20u = round(Decimal(lower_20 / biomass_20), 6), round(Decimal(upper_20 / biomass_20), 6)
			absMin20 = absolute_minimum(normalized_20l, normalized_20u)
		except (KeyError, ValueError):
			lower_20, upper_20 = 'NA', 'NA'
			normalized_20l, normalized_20u =  'NA', 'NA'
			absMin20 = 'NA'

		fifteen_twenty = compare_differences_absMin(absMin15, absMin20, '15C', '20C')
		four_twenty = compare_differences_absMin(absMin4, absMin20, '4C', '20C')
		four_fifteen = compare_differences_absMin(absMin4, absMin15, '4C', '15C')

		absMin_15_20.append(fifteen_twenty)
		absMin_04_20.append(four_twenty)
		absMin_04_15.append(four_fifteen)

		absMin_file.write('AbsMin_Comparisons\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(reaction, normalized_4l,
													normalized_4u, absMin4, normalized_15l, normalized_15u,
													absMin15, normalized_20l, normalized_20u, absMin20,
													four_fifteen, four_twenty, fifteen_twenty))

absMin_file.close()


carbNorm_file = open('{}_carbNorm_Variability.tsv'.format(args.out), mode='w')

carbNorm_15_20 = []
carbNorm_04_20 = []
carbNorm_04_15 = []

(_, carb_uptake_4) = flux_entry_dict[4]['EX_cpd_acgam[e]']
(_, carb_uptake_15) = flux_entry_dict[15]['EX_cpd_acgam[e]']
(_, carb_uptake_20) = flux_entry_dict[20]['EX_cpd_acgam[e]']
carb_uptake_4 = abs(carb_uptake_4)
carb_uptake_15 = abs(carb_uptake_15)
carb_uptake_20 = abs(carb_uptake_20)
carbNorm_file.write('carbNorm_Comparisons\tID\tlower_4\tupper_4\tcarbNorm_4\tlower_15\tupper_15\tcarbNorm_15\tlower_20\tupper_20\tcarbNorm_20\tfour_fifteen\tfour_twenty\tfifteen_twenty\n')
for reaction in sorted(var_list_flux):
		try:
			(lower_4, upper_4) = flux_entry_dict[4][reaction]
			normalized_4l, normalized_4u = round(Decimal(lower_4 / carb_uptake_4), 6), round(Decimal(upper_4 / carb_uptake_4), 6)
			carbNorm4 = absolute_minimum(normalized_4l, normalized_4u)
		except (KeyError, ValueError):
			lower_4, upper_4 = 'NA', 'NA'
			normalized_4l, normalized_4u =  'NA', 'NA'
			carbNorm4 = 'NA'

		try:
			(lower_15, upper_15) = flux_entry_dict[15][reaction]
			normalized_15l, normalized_15u = round(Decimal(lower_15 / carb_uptake_15), 6), round(Decimal(upper_15 / carb_uptake_15), 6)
			carbNorm15 = absolute_minimum(normalized_15l, normalized_15u)
		except (KeyError, ValueError):
			lower_15, upper_15 = 'NA', 'NA'
			normalized_15l, normalized_15u =  'NA', 'NA'
			carbNorm15 = 'NA'

		try:
			(lower_20, upper_20) = flux_entry_dict[20][reaction]
			normalized_20l, normalized_20u = round(Decimal(lower_20 / carb_uptake_20), 6), round(Decimal(upper_20 / carb_uptake_20), 6)
			carbNorm20 = absolute_minimum(normalized_20l, normalized_20u)
		except (KeyError, ValueError):
			lower_20, upper_20 = 'NA', 'NA'
			normalized_20l, normalized_20u =  'NA', 'NA'
			carbNorm20 = 'NA'

		fifteen_twenty = compare_differences_absMin(carbNorm15, carbNorm20, '15C', '20C')
		four_twenty = compare_differences_absMin(carbNorm4, carbNorm20, '4C', '20C')
		four_fifteen = compare_differences_absMin(carbNorm4, carbNorm15, '4C', '15C')

		carbNorm_15_20.append(fifteen_twenty)
		carbNorm_04_20.append(four_twenty)
		carbNorm_04_15.append(four_fifteen)

		carbNorm_file.write('carbNorm_Comparisons\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(reaction, normalized_4l,
													normalized_4u, carbNorm4, normalized_15l, normalized_15u,
													carbNorm15, normalized_20l, normalized_20u, carbNorm20,
													four_fifteen, four_twenty, fifteen_twenty))

carbNorm_file.close()

summary_file = open('{}_SummaryCounts.tsv'.format(args.out), mode='w')
flux_04_15_c, flux_04_20_c, flux_15_20_c = Counter(flux_04_15), Counter(flux_04_20), Counter(flux_15_20)
flux_keys = set(list(flux_04_15_c.keys())+list(flux_04_20_c.keys())+list(flux_15_20_c.keys()))
summary_file.write('Flux_Category\t04_15\t04_20\t15_20\n'.format())
for key in sorted(flux_keys):
	summary_file.write('{}\t{}\t{}\t{}\n'.format(key, flux_04_15_c.get(key, 0), flux_04_20_c.get(key, 0), flux_15_20_c.get(key, 0)))

summary_file.write('\n')

deltaG_04_15_c, deltaG_04_20_c, deltaG_15_20_c = Counter(deltaG_04_15), Counter(deltaG_04_20), Counter(deltaG_15_20)
deltaG_keys = set(list(deltaG_04_15_c.keys())+list(deltaG_04_20_c.keys())+list(deltaG_15_20_c.keys()))
summary_file.write('DeltaG_Category\t04_15\t04_20\t15_20\n'.format())
for key in sorted(deltaG_keys):
	summary_file.write('{}\t{}\t{}\t{}\n'.format(key, deltaG_04_15_c.get(key, 0), deltaG_04_20_c.get(key, 0), deltaG_15_20_c.get(key, 0)))

summary_file.write('\n')

cpd_04_15_c, cpd_04_20_c, cpd_15_20_c = Counter(cpd_04_15), Counter(cpd_04_20), Counter(cpd_15_20)
cpd_keys = set(list(cpd_04_15_c.keys())+list(cpd_04_20_c.keys())+list(cpd_15_20_c.keys()))
summary_file.write('Cpd_Category\t04_15\t04_20\t15_20\n'.format())
for key in sorted(cpd_keys):
	summary_file.write('{}\t{}\t{}\t{}\n'.format(key, cpd_04_15_c.get(key, 0), cpd_04_20_c.get(key, 0), cpd_15_20_c.get(key, 0)))

summary_file.write('\n')

absMin_04_15_c, absMin_04_20_c, absMin_15_20_c = Counter(absMin_04_15), Counter(absMin_04_20), Counter(absMin_15_20)
absMin_keys = set(list(absMin_04_15_c.keys())+list(absMin_04_20_c.keys())+list(absMin_15_20_c.keys()))
summary_file.write('absMin_Category\t04_15\t04_20\t15_20\n'.format())
for key in sorted(absMin_keys):
	summary_file.write('{}\t{}\t{}\t{}\n'.format(key, absMin_04_15_c.get(key, 0), absMin_04_20_c.get(key, 0), absMin_15_20_c.get(key, 0)))
