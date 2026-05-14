from .resolve_perm import resolve_permutation
from .filter_synteny_blocks import get_perm_from_psl
import os
import pandas as pd


def inv_resolver(input_file, output_dir, query, target, seed,
                 overlap_length, overlap_percent, remove_overlaps,
                 chain_step, use_chaining, scaled_v, compact_v):

    current_dir = os.getcwd()
    new_dir_path = os.path.join(current_dir, output_dir)
    os.makedirs(new_dir_path, exist_ok=True)
    dir = new_dir_path.replace("\\", "/")

    perms = get_perm_from_psl(input_file, dir, query, target,
                              overlap_length, overlap_percent, remove_overlaps,
                              chain_step, use_chaining
                              )
    chr_chains = pd.read_csv(f"{dir}/{query}.to.{target}.chains.csv")

    for perm in perms:
        resolve_permutation(perm, ".".join(perm.split(".")[0:5]), seed, query, target, scaled_v, compact_v, chr_chains)
