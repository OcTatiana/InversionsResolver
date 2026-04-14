import random
import argparse
import sys
import pandas as pd

from perm_sorting import get_reversal_sequences
from visual import draw_bezier_curves
# from reversals_order import canonicalize


def resolve_permutation(input_file: str, output: str, seed_val):
    random.seed(seed_val)

    try:
        with open(input_file, "r") as f_in:
            perm = [int(x) for x in f_in.readline().split()]

        # print(perm)
        perm_global_list, id_global_list = get_reversal_sequences(perm, 100)

        if not output:
            output = ".".join(input_file.split(".")[0:5])

        output_file_id = "{}.id".format(output)
        output_file_perm = "{}.permseq".format(output)
        scaled_output_png = "{}_scaled.png".format(output)
        output_png = "{}.png".format(output)

        with open(output_file_perm, 'w') as f_out_perm:
            f_out_perm.write('\n'.join(' '.join(str(x) for x in row) for row in perm_global_list))

        print(f"Successfully processed '{input_file}' and saved to '{output}'")
        print(f"Number of reversals: {str(len(id_global_list[0]))}")
        print(f"Number of retained reversal sequences: {len(perm_global_list)}")

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
        synteny_order[0] = 0
        synteny_order[len(synteny_block_names)+1] = len(synteny_block_names) + 1

        print("Drawing")
        draw_bezier_curves(perm_global_list[0], scaled_output_png, output.split("\\")[-1], synteny_block_names_dict,
                           synteny_order, blocks_len)
        draw_bezier_curves(perm_global_list[0], output_png, output.split("\\")[-1], synteny_block_names_dict,
                           synteny_order)

        # i = 0
        # for perm_i in perm_global_list:
        #     if i > 10:
        #         break
        #     i += 1
        #     # draw_bezier_curves(perm_i, output_png, output)
        #     draw_bezier_curves(perm_i, output + "_" + str(i) + ".png", output, synteny_block_names_dict, synteny_order)


    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied when accessing files.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Read a permutation from input file and resolve it',
        epilog='Example: python script.py -i input.txt [-o output_name]'
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to the input file to read'
    )

    parser.add_argument(
        '-o', '--output',
        required=False,
        help='Name the output file to write'
    )

    args = parser.parse_args()
    resolve_permutation(args.input, args.output)


if __name__ == "__main__":
    main()

