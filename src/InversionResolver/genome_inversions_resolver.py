from .resolve_perm import resolve_permutation
from .filter_synteny_blocks import get_perm_from_psl
import os


def inv_resolver(input_file, output_dir, species, seed,
                 overlap_length, overlap_percent, remove_overlaps,
                 chain_step, use_chaining):

    current_dir = os.getcwd()
    new_dir_path = os.path.join(current_dir, output_dir)
    os.makedirs(new_dir_path, exist_ok=True)

    perms = get_perm_from_psl(input_file, new_dir_path.replace("\\", "/"), species)

    for perm in perms:
        print(perm)
        # resolve_permutation(perm, perm.split(".")[0:5], seed)
