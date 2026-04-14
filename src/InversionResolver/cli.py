import argparse
from .resolve_perm import resolve_permutation
from .genome_inversions_resolver import inv_resolver


def main():
    parser = argparse.ArgumentParser(
        prog="InversionResolver",
        description="Tool for inversion resolving"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # =========================================================
    # Command 1: SORT (algorithm only)
    # =========================================================
    sort_parser = subparsers.add_parser(
        "sort",
        help="Resolve single permutation"
    )

    sort_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file with space-separated signed permutation"
    )

    sort_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file name"
    )

    # sort_parser.add_argument(
    #     "-v", "--visualize",
    #     action="store_true",
    #     help="Enable visualization"
    # )

    sort_parser.add_argument(
        "-s", "--seed",
        type=int,
        default=30,
        help="Random seed (default 30)"
    )

    # =========================================================
    # Command 2: GENOME (pipeline)
    # =========================================================
    genome_parser = subparsers.add_parser(
        "genome",
        help="Run full pipeline: preprocessing + algorithm"
    )

    genome_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input raw data file (psl)"
    )

    genome_parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory"
    )

    genome_parser.add_argument(
        "-sp", "--species",
        required=True,
        help="Species name (used in file naming)"
    )

    # -------- Optional parameters --------

    # genome_parser.add_argument(
    #     "-v", "--visualize",
    #     action="store_true",
    #     help="Enable visualization"
    # )

    genome_parser.add_argument(
        "-s", "--seed",
        type=int,
        default=30,
        help="Random seed (default: 30)"
    )

    genome_parser.add_argument(
        "-ol", "--overlap-length",
        type=int,
        default=100000,
        help="Overlap length (default: 10^5)"
    )

    genome_parser.add_argument(
        "-op", "--overlap-percent",
        type=float,
        default=0.1,
        help="Overlap percent (default: 0.1)"
    )

    genome_parser.add_argument(
        "-ro", "--remove-overlaps",
        action="store_true",
        help="Remove all overlaps (default: False)"
    )

    genome_parser.add_argument(
        "-cs", "--chain-step",
        type=int,
        default=1000000,
        help="Chaining step (default: 10^6)"
    )

    genome_parser.add_argument(
        "-nc", "--no-chaining",
        action="store_true",
        help="Disable chaining (enabled by default)"
    )

    # =========================================================
    # Parse and dispatch
    # =========================================================
    args = parser.parse_args()

    if args.command == "sort":
        resolve_permutation(
            input_file=args.input,
            output=args.output,
            seed_val=args.seed
        )

    elif args.command == "genome":
        inv_resolver(
            input_file=args.input,
            output_dir=args.output_dir,
            species=args.species,
            seed=args.seed,
            overlap_length=args.overlap_length,
            overlap_percent=args.overlap_percent,
            remove_overlaps=args.remove_overlaps,
            chain_step=args.chain_step,
            use_chaining=not args.no_chaining
        )


if __name__ == "__main__":
    main()