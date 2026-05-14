import random
import sys
import pandas as pd

from .perm_sorting import get_reversal_sequences
from .visual import draw_bezier_curves
from .reversals_order import get_perm_for_image


def resolve_permutation(input_file: str, output: str, seed_val, query, target, scaled_v, compact_v, chr_chains=None):
    random.seed(seed_val)

    try:
        with open(input_file, "r") as f_in:
            perm = [int(x) for x in f_in.readline().split()]

        print(input_file)
        perm_global_list, blocks_global_list, id_global_list = get_reversal_sequences(perm, 10)

        if not output:
            output = ".".join(input_file.split(".")[0:5])

        output_file_id = "{}.id".format(output)
        output_file_perm = "{}.permseq".format(output)
        scaled_output_png = "{}_scaled.png".format(output)
        compact_output_png = "{}_compact.png".format(output)
        scaled_compact_output_png = "{}_scaled_compact.png".format(output)
        output_png = "{}.png".format(output)

        with open(output_file_perm, 'w') as f_out_perm:
            f_out_perm.write('\n'.join(' '.join(str(x) for x in row) for row in perm_global_list))

        with open(output_file_id, 'w') as f_out_perm:
            f_out_perm.write('\n'.join(' '.join(str(x) for x in row) for row in id_global_list))

        print(f"Successfully processed '{input_file}'")
        if not blocks_global_list:
            print("Number of reversals: 0")
            return

        print(f"Number of reversals: {str(len(blocks_global_list[0]))}")
        print(f"Number of retained reversal sequences: {len(perm_global_list)}")

        try:
            synteny_block_names_f = "{}.synteny_blocks_metadata.csv".format(input_file.split(".filtered_synteny_blocks")[0])
            synteny_block = pd.read_csv(synteny_block_names_f)
            synteny_block_names = synteny_block['synteny_block_id'].tolist()
            synteny_block_len = synteny_block['tHitLen'].tolist()

            synteny_block_names_dict = {abs(perm[i+1]): f"{synteny_block_names[i]}" for i in range(len(synteny_block_names))}
            synteny_block_names_dict[0] = "out"
            synteny_block_names_dict[len(synteny_block_names)+1] = "out"

            blocks_len = {abs(perm[i + 1]): synteny_block_len[i] for i in range(len(synteny_block_names))}
            blocks_len[0] = 50000
            blocks_len[len(synteny_block_names)+1] = 50000

            with open(output_file_id, 'w') as f_out_id:
                for option in id_global_list:
                    blocks = []
                    for row in option:
                        blocks.append([synteny_block_names_dict[i] for i in row])
                    f_out_id.write(' '.join(str(block) for block in blocks) + "\n")

            synteny_order = {abs(perm[i+1]): f"{abs(perm[i+1])}" for i in range(len(synteny_block_names))}
            synteny_order[0] = "0"
            synteny_order[len(synteny_block_names)+1] = f"{len(synteny_block_names) + 1}"

            chr_list_q = []
            chr_list_t = []
            if chr_chains is not None:
                target_chr = input_file.split(".")[4]
                query_chr = input_file.split(".")[1]
                chr_list_q = chr_chains[chr_chains["qName"] == query_chr].sort_values(by="qStart")["synteny_block_id"].unique().tolist()
                chr_list_t = chr_chains[chr_chains["tName"] == target_chr].sort_values(by="tStart")["synteny_block_id"].unique().tolist()

            print("Drawing")
            perms_for_image = perm_global_list[0]

            draw_bezier_curves(perms_for_image, output_png, synteny_block_names_dict,
                               synteny_order, query, target, None, chr_list_q, chr_list_t, synteny_block_names)

            if scaled_v:
                draw_bezier_curves(perms_for_image, scaled_output_png, synteny_block_names_dict,
                               synteny_order, query, target, blocks_len, chr_list_q, chr_list_t, synteny_block_names)
            if compact_v:
                perms_for_image = get_perm_for_image(perms_for_image, id_global_list[0])
                draw_bezier_curves(perms_for_image, compact_output_png, synteny_block_names_dict,
                                   synteny_order, query, target, None, chr_list_q, chr_list_t, synteny_block_names)
            if scaled_v and compact_v:
                perms_for_image = get_perm_for_image(perms_for_image, id_global_list[0])
                draw_bezier_curves(perms_for_image, scaled_compact_output_png, synteny_block_names_dict,
                                   synteny_order, query, target, blocks_len, chr_list_q, chr_list_t, synteny_block_names)


            # i = 0
            # for perm_i in perm_global_list:
            #     if i > 10:
            #         break
            #     i += 1
            #     # draw_bezier_curves(perm_i, output_png, output)
            #     draw_bezier_curves(perm_i, output + "_" + str(i) + ".png", output, synteny_block_names_dict, synteny_order)

        except FileNotFoundError:
            n = len(perm_global_list[0][0])
            synteny_block_names_dict = {i: f"SB_{i}" for i in range(n - 1)}
            synteny_block_names_dict[0] = "out"
            synteny_block_names_dict[n - 1] = "out"

            synteny_order = {abs(perm[i + 1]): f"{abs(perm[i + 1])}" for i in range(n - 1)}
            synteny_order[0] = "0"
            synteny_order[n - 1] = f"{n - 1}"

            print("Drawing")
            draw_bezier_curves(perm_global_list[0], output_png, synteny_block_names_dict, synteny_order)

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied when accessing files.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

