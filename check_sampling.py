#!/usr/bin/env python3
"""
Compare random vs. stratified term sampling across knowledge graphs.

For each KG, this script:
  1. Loads the dataset via PyKEEN
  2. Assigns entity types (same logic as prepare.py)
  3. Shows the full type distribution
  4. Compares random vs. stratified sampling at DOWNSAMPLE_SIZE
  5. Reports representation ratios and lost/gained types

Usage:
    ./check_sampling.py [dataset ...]
    ./check_sampling.py                   # all four
    ./check_sampling.py hetionet primekg  # specific ones

No training triples or S3 metadata needed for type assignment:
  - Hetionet: type is embedded in entity name prefix
  - BioKG: type is the PyKEEN entity prefix (before the colon)
  - OpenBioLink: type is the PyKEEN entity prefix (before the colon)
  - PrimeKG: type is the PyKEEN entity prefix (before the colon)

If S3 metadata files are available, the script will use them for richer
type labels (matching prepare.py exactly). Set USE_S3_METADATA=True below.
"""

import sys
import random
import numpy as np
import pandas as pd
from collections import Counter

RANDOM_SEED = 42
DOWNSAMPLE_SIZE = 10_000
USE_S3_METADATA = True
S3_BASE_PATH = 's3://alethiotx-artemis/data/kgs/raw'

SUPPORTED = ['hetionet', 'biokg', 'openbiolink', 'primekg']


def set_seeds(seed=RANDOM_SEED):
    random.seed(seed)
    np.random.seed(seed)


def get_all_entities(kg):
    """Collect all entities across train/test/val splits."""
    entities = set()
    for split in [kg.training, kg.testing, kg.validation]:
        entities.update(split.entity_to_id.keys())
    return sorted(entities)


# ── Type assignment ──────────────────────────────────────────────────────────

def assign_types_hetionet(entities):
    """Hetionet entities are 'Type::Name', e.g. 'Gene::BRCA1'."""
    prefixes = [
        'Biological Process', 'Pathway', 'Cellular Component',
        'Molecular Function', 'Anatomy', 'Compound', 'Disease', 'Gene'
    ]
    types = {}
    for e in entities:
        matched = False
        for p in prefixes:
            if e.startswith(f'{p}::'):
                types[e] = p
                matched = True
                break
        if not matched:
            types[e] = 'Other'
    return types


def assign_types_from_prefix(entities):
    """
    For BioKG, OpenBioLink, and PrimeKG the PyKEEN entity IDs often carry a
    type indicator.  This function uses heuristics on entity name patterns.

    BioKG entities look like: 'Q9Y6K9' (protein), 'DOID:...' (disease), etc.
    OpenBioLink: 'GENE:1234', 'GO:...', 'DIS:...', etc.
    PrimeKG: 'aspirin' (drug), 'BRCA1' (gene/protein), etc.

    Because the exact mapping depends on the KG and its metadata, we fall
    back to the S3 metadata logic from prepare.py when USE_S3_METADATA is set.
    Otherwise, we use the PyKEEN entity-to-id mapping plus relation endpoints
    to infer types.
    """
    # Generic fallback: classify by first component before ':'
    types = {}
    for e in entities:
        if '::' in e:
            types[e] = e.split('::')[0]
        elif ':' in e:
            types[e] = e.split(':')[0]
        else:
            types[e] = 'Unknown'
    return types


def assign_types_biokg_from_metadata(entities):
    """Use S3 metadata files (matches prepare.py logic)."""
    biokg_path = f'{S3_BASE_PATH}/biokg'
    entity_set = set(entities)
    types = {}

    protein = pd.read_csv(f'{biokg_path}/biokg.metadata.protein.tsv', sep="\t", header=None)
    human_ids = set(protein[(protein[1] == 'SPECIES') & (protein[2] == 'HUMAN')][0])
    protein_ids = human_ids & entity_set
    for e in protein_ids:
        types[e] = 'Protein'

    disease = pd.read_csv(f'{biokg_path}/biokg.metadata.disease.tsv', sep="\t", header=None)
    for e in set(disease[disease[1] == 'NAME'][0]) & entity_set:
        types.setdefault(e, 'Disease')

    pathway = pd.read_csv(f'{biokg_path}/biokg.metadata.pathway.tsv', sep="\t", header=None)
    for e in set(pathway[pathway[1] == 'NAME'][0]) & entity_set:
        types.setdefault(e, 'Pathway')

    drug = pd.read_csv(f'{biokg_path}/biokg.metadata.drug.tsv', sep="\t", header=None)
    for e in set(drug[drug[1] == 'NAME'][0]) & entity_set:
        types.setdefault(e, 'Drug')

    links = pd.read_csv(f'{biokg_path}/biokg.links.tsv', sep="\t", header=None)
    for e in set(links[links[1] == 'MEMBER_OF_COMPLEX'][2]) & entity_set:
        types.setdefault(e, 'Complex')
    for e in set(links[links[1] == 'RELATED_GENETIC_DISORDER'][2]) & entity_set:
        types.setdefault(e, 'Genetic Disorder')

    for e in entities:
        types.setdefault(e, 'Other')
    return types


def assign_types_openbiolink_from_metadata(entities):
    """Use S3 nodes.csv (matches prepare.py logic)."""
    obl_path = f'{S3_BASE_PATH}/openbiolink'
    nodes = pd.read_csv(f'{obl_path}/graph_files/nodes.csv', sep="\t", header=None)
    type_map = dict(zip(nodes[0], nodes[1]))
    return {e: type_map.get(e, 'Unknown') for e in entities}


def assign_types_primekg_from_metadata(entities):
    """Use S3 nodes.csv (matches prepare.py logic)."""
    primekg_path = f'{S3_BASE_PATH}/primekg'
    nodes = pd.read_csv(f'{primekg_path}/nodes.csv')
    type_map = dict(zip(nodes['node_name'], nodes['node_type']))
    return {e: type_map.get(e, 'Unknown') for e in entities}


def assign_types_via_relations(kg, entities):
    """
    Infer entity types from the relation endpoints in the KG triples.
    For each entity, find all relations where it appears as head or tail,
    and use the most frequent relation as a proxy for type.
    """
    entity_set = set(entities)
    triples = kg.training.triples
    head_rels = {}
    tail_rels = {}

    for h, r, t in triples:
        if h in entity_set:
            head_rels.setdefault(h, []).append(r)
        if t in entity_set:
            tail_rels.setdefault(t, []).append(r)

    types = {}
    for e in entities:
        rels = head_rels.get(e, []) + tail_rels.get(e, [])
        if rels:
            most_common = Counter(rels).most_common(1)[0][0]
            types[e] = most_common
        else:
            types[e] = 'Unknown'
    return types


def get_entity_types(dataset_name, kg, entities):
    """Return {entity: type} dict for the given KG."""
    if dataset_name == 'hetionet':
        return assign_types_hetionet(entities)

    if USE_S3_METADATA:
        if dataset_name == 'biokg':
            return assign_types_biokg_from_metadata(entities)
        elif dataset_name == 'openbiolink':
            return assign_types_openbiolink_from_metadata(entities)
        elif dataset_name == 'primekg':
            return assign_types_primekg_from_metadata(entities)

    # Try prefix-based assignment first
    types = assign_types_from_prefix(entities)
    unique_types = set(types.values())

    # If prefix-based gives only 'Unknown', fall back to relation-based
    if unique_types == {'Unknown'}:
        print("  Prefix-based typing found no types, falling back to relation-based inference...")
        types = assign_types_via_relations(kg, entities)

    return types


# ── Sampling ─────────────────────────────────────────────────────────────────

def random_sample(entities, n, seed=RANDOM_SEED):
    rng = random.Random(seed)
    return rng.sample(entities, min(n, len(entities)))


def stratified_sample(entities, types, n, seed=RANDOM_SEED):
    """
    Proportional stratified sample: each type gets floor(proportion * n)
    slots, remainder allocated round-robin to under-represented types.
    """
    rng = random.Random(seed)
    by_type = {}
    for e in entities:
        by_type.setdefault(types[e], []).append(e)

    total = len(entities)
    n = min(n, total)
    sampled = []
    quotas = {}

    # Proportional allocation
    for t, members in by_type.items():
        quota = int(len(members) / total * n)
        quotas[t] = quota

    # Distribute remainder to largest groups first
    remainder = n - sum(quotas.values())
    sorted_types = sorted(by_type.keys(), key=lambda t: len(by_type[t]), reverse=True)
    for i in range(remainder):
        quotas[sorted_types[i % len(sorted_types)]] += 1

    # Sample within each type
    for t, members in by_type.items():
        q = min(quotas[t], len(members))
        sampled.extend(rng.sample(members, q))

    return sampled


# ── Reporting ────────────────────────────────────────────────────────────────

def type_distribution(entities, types):
    counts = Counter(types[e] for e in entities)
    return counts


def print_comparison(dataset_name, entities, types):
    total = len(entities)
    full_dist = type_distribution(entities, types)

    print(f"\n{'='*70}")
    print(f"  {dataset_name.upper()}  —  {total:,} total entities, {len(full_dist)} types")
    print(f"{'='*70}")

    if total <= DOWNSAMPLE_SIZE:
        print(f"  ⓘ  Only {total:,} entities — below {DOWNSAMPLE_SIZE:,} threshold, no downsampling needed.")
        print_dist_table(full_dist, full_dist, full_dist, total)
        return

    rand = random_sample(entities, DOWNSAMPLE_SIZE)
    strat = stratified_sample(entities, types, DOWNSAMPLE_SIZE)

    rand_dist = type_distribution(rand, types)
    strat_dist = type_distribution(strat, types)

    print_dist_table(full_dist, rand_dist, strat_dist, total)

    # Summary statistics
    print(f"\n  Representation error (avg |sample% - full%| across types):")
    rand_err = mean_pct_error(full_dist, rand_dist, total, DOWNSAMPLE_SIZE)
    strat_err = mean_pct_error(full_dist, strat_dist, total, DOWNSAMPLE_SIZE)
    print(f"    Random:     {rand_err:.3f} pp")
    print(f"    Stratified: {strat_err:.3f} pp")

    # Types lost
    rand_lost = set(full_dist) - set(rand_dist)
    strat_lost = set(full_dist) - set(strat_dist)
    if rand_lost:
        print(f"\n  ⚠  Types lost in random sample: {', '.join(sorted(rand_lost))}")
    if strat_lost:
        print(f"  ⚠  Types lost in stratified sample: {', '.join(sorted(strat_lost))}")
    if not rand_lost and not strat_lost:
        print(f"\n  ✓  No types lost in either sampling method.")


def print_dist_table(full, rand, strat, total_full):
    total_rand = sum(rand.values())
    total_strat = sum(strat.values())

    all_types = sorted(set(full) | set(rand) | set(strat), key=lambda t: -full.get(t, 0))

    print(f"\n  {'Type':<30} {'Full':>8} {'%':>7}  {'Random':>8} {'%':>7}  {'Strat.':>8} {'%':>7}")
    print(f"  {'─'*30} {'─'*8} {'─'*7}  {'─'*8} {'─'*7}  {'─'*8} {'─'*7}")

    for t in all_types:
        fc = full.get(t, 0)
        rc = rand.get(t, 0)
        sc = strat.get(t, 0)
        fp = fc / total_full * 100 if total_full else 0
        rp = rc / total_rand * 100 if total_rand else 0
        sp = sc / total_strat * 100 if total_strat else 0
        print(f"  {t:<30} {fc:>8,} {fp:>6.2f}%  {rc:>8,} {rp:>6.2f}%  {sc:>8,} {sp:>6.2f}%")

    print(f"  {'─'*30} {'─'*8} {'─'*7}  {'─'*8} {'─'*7}  {'─'*8} {'─'*7}")
    print(f"  {'TOTAL':<30} {total_full:>8,}         {total_rand:>8,}         {total_strat:>8,}")


def mean_pct_error(full_dist, sample_dist, total_full, total_sample):
    all_types = set(full_dist) | set(sample_dist)
    errors = []
    for t in all_types:
        fp = full_dist.get(t, 0) / total_full * 100
        sp = sample_dist.get(t, 0) / total_sample * 100
        errors.append(abs(fp - sp))
    return np.mean(errors) if errors else 0.0


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    set_seeds()

    from pykeen.datasets import Hetionet, BioKG, OpenBioLink, PrimeKG

    dataset_map = {
        'hetionet': Hetionet,
        'biokg': BioKG,
        'openbiolink': OpenBioLink,
        'primekg': PrimeKG,
    }

    datasets = sys.argv[1:] if len(sys.argv) > 1 else SUPPORTED
    for name in datasets:
        name = name.lower()
        if name not in dataset_map:
            print(f"Unknown dataset: {name}, skipping")
            continue

        print(f"\nLoading {name}...")
        kg = dataset_map[name]()
        entities = get_all_entities(kg)
        print(f"  {len(entities):,} entities extracted")

        types = get_entity_types(name, kg, entities)
        print_comparison(name, entities, types)


if __name__ == '__main__':
    main()
