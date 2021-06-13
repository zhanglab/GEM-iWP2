from psamm.commands.tmfa import *
import yaml
import argparse
from psamm.datasource import native
from psamm.lpsolver import generic
import numpy as np

solver = generic.Solver()

parser = argparse.ArgumentParser()
parser.add_argument('--config', type=argparse.FileType('r'),
                    help='Config file for TMFA settings')
parser.add_argument('--model')
parser.add_argument('--biomass')
parser.add_argument('--varying')
parser.add_argument('--temp', help='Temperature in Celsius',
                    default=25)
parser.add_argument('--phin', default='4,11',
                            help='Allowable range for internal compartment pH '
                            'format: Lower,Upper')
parser.add_argument('--phout', default='4,11',
                            help='Allowable range for external compartment pH '
                            'format: Lower,Upper')
args = parser.parse_args()


model_reader = native.ModelReader.reader_from_path(args.model)
model = model_reader.create_model()
mm = model.create_metabolic_model()

phin_s = args.phin
phout_s = args.phout
phinsplit = phin_s.split(',')
phoutsplit = phout_s.split(',')
phin = (float(phinsplit[0]), float(phinsplit[1]))
phout = (float(phoutsplit[0]), float(phoutsplit[1]))
ph_bounds = [phin, phout]

config_dict = yaml.safe_load(args.config)
# Parse exclude file provided through command line.
base_exclude_list = []
if config_dict.get('exclude') is not None:
    for line in open(config_dict.get('exclude')).readlines():
        base_exclude_list.append(line.rstrip())

# Parse set of reactions that will have their deltaG
# constrained only by deltaG of the transport component.
ph_difference_rxn = []
if config_dict.get('rxn_conc_only') is not None:
    for line in open(config_dict.get('rxn_conc_only')).readlines():
        rxn = line.rstrip()
        ph_difference_rxn.append(rxn)
        ph_difference_rxn.append(u'{}_forward'.format(rxn))
        ph_difference_rxn.append(u'{}_reverse'.format(rxn))

# Parse out the file of lumped reactions or set the lumps
# dictionaries to be empty
if config_dict.get('lumped_reactions') is not None:
    lumpid_to_rxn, rxn_lumpid, lump_to_rxnids, \
        lump_to_rxnids_dir = lump_parser(
            open(config_dict.get('lumped_reactions')))
else:
    logger.info('No lumped reactions were provided')
    lumpid_to_rxn = {}
    lump_to_rxnids = {}
    lump_to_rxnids_dir = {}

# Add exchange reactions to the exclude list
for reaction in mm.reactions:
    if mm.is_exchange(reaction):
        base_exclude_list.append(reaction)

# these lines will add ph_difference list to the base exclude list.
base_exclude_list = base_exclude_list + ph_difference_rxn

# exclude lump_list is a list of all excluded reactions and
# lump reactions. exclude_lump_unkown list is a list of all
# excluded reactions, lump reactions and lumped reactions.
exclude_lump_list = set(base_exclude_list)
exclude_unknown_list = set(base_exclude_list)
for lump_id, sub_rxn_list in iteritems(lump_to_rxnids):
    exclude_lump_list.add(lump_id)
    for sub_rxn in sub_rxn_list:
        exclude_unknown_list.add(sub_rxn)
# Make list of all excluded, lumped, and unknown dgr.
exclude_lump_unkown = exclude_unknown_list.union(exclude_lump_list)

cpd_conc_dict = {}
# Parse file containing set concentrations and ranges.
if config_dict.get('concentrations') is not None:
    for row in csv.reader(
        open(config_dict.get('concentrations')),
            delimiter=str('\t')):
        cpd, lower, upper = row
        cpd_conc_dict[convert_to_unicode(cpd)] = [lower, upper]

# Add lump reactions to the mm_irreversible model
for lump_id, rxn in iteritems(lumpid_to_rxn):
    reaction = parse_reaction(rxn)
    mm.database.set_reaction(lump_id, reaction)
    mm.add_reaction(lump_id)
    mm.limits[lump_id].lower = 0
    mm.limits[lump_id].upper = 0

# Remove limits on varying reaction
mm.limits[args.varying].lower = 0
mm.limits[args.varying].upper = 100

# make an irreversible version of the metabolic model:
mm_irreversible, _, split_reversible, \
    reversible_lump_to_rxn_dict = mm.make_irreversible(
        gene_dict={}, exclude_list=exclude_lump_unkown,
        lumped_rxns=lump_to_rxnids_dir,
        all_reversible=False)

split_reactions = []
for_rev_reactions = []
for (i, j) in split_reversible:
    for_rev_reactions.append(i)
    for_rev_reactions.append(j)
    split_reactions.append(i[:-8])

for rx in [i for i in model.reactions]:
    if rx.id in split_reactions:
        f = ReactionEntry(
            dict(id=u'{}_forward'.format(rx.id), genes=rx.genes))
        r = ReactionEntry(
            dict(id=u'{}_reverse'.format(rx.id), genes=rx.genes))
        model.reactions.add_entry(f)
        model.reactions.add_entry(r)

# update these lists for the new reversible lumps and constituent
# reactions. exclude lump_list is a list of all excluded reactions
# and lump reactions. exclude_lump_unkown list is a list of
# all excluded reactions, lump reactions and lumped reactions.
for lump_id, sub_rxn_list in iteritems(reversible_lump_to_rxn_dict):
    exclude_lump_list.add(lump_id)
    for sub_rxn in sub_rxn_list:
        exclude_unknown_list.add(sub_rxn)
# Make list of all excluded, lumped, and unknown dgr.
exclude_lump_unkown = exclude_unknown_list.union(exclude_lump_list)

# Read DeltaG information from file.
dgr_dict = {}
if config_dict.get('deltaG') is not None:
    dgr_dict = parse_dgr_file(open(
        config_dict.get('deltaG')), mm_irreversible)
elif config_dict.get('deltaGf') is not None:
    # Parse the deltaGf values for all metabolites
    # from a supplied file.
    dgf_dict = parse_dgf(mm_irreversible, open(
        config_dict.get('deltaGf')))
    # Calculate the deltaGr values for all reactions
    dgr_dict = calculate_dgr(
        mm_irreversible, dgf_dict, exclude_unknown_list,
        open(config_dict.get('transporters')),
        ph_difference_rxn, open(config_dict.get('scaled_compounds')))
    print('using dgf file')

if config_dict.get('transporters') is not None:
    transport_parameters = parse_tparam_file(
        open(config_dict.get('transporters')))

prob, v, zi, dgri, xij, cp_list = make_tmfa_problem(
    mm_irreversible, solver)

prob, cpd_xij_dict = add_conc_constraints(
    xij, prob, cpd_conc_dict, cp_list,
    config_dict.get('water'), config_dict.get('proton-in'),
    config_dict.get('proton-out'), config_dict.get('proton-other'))
testing_list_tmp = list(mm_irreversible.reactions)
prob, excluded_compounds = add_reaction_constraints(
    prob, v, zi, dgri, xij, mm_irreversible, exclude_lump_list,
    exclude_unknown_list, exclude_lump_unkown, dgr_dict,
    reversible_lump_to_rxn_dict, split_reversible,
    transport_parameters, testing_list_tmp,
    config_dict.get('scaled_compounds'), config_dict.get('water'),
    config_dict.get('proton-in'), config_dict.get('proton-out'),
    config_dict.get('proton-other'), ph_bounds,
    args.temp, False)


max = get_var_bound(prob, v(args.varying),
                                lp.ObjectiveSense.Maximize)
max_vary = np.round(max, 1)
try:
    for i in np.arange(0, max_vary+0.1, 0.05):
        c, = prob.add_linear_constraints(v(args.varying) == i)
        bio_flux = get_var_bound(prob, v(args.biomass),
                                        lp.ObjectiveSense.Maximize)
        print('{}\t{}'.format(i, bio_flux))
        c.delete()
except:
    pass
