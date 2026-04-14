import pandas as pd
import os


OVERLAP_PERCENT = 0.1
OVERLAP_LEN = 10 ** 5
REMOVE_ALL = False
CHAINING_STEP = 10 ** 6
CHAINING_COUNT = 1


def is_nested(start1, end1, start2, end2):
    if start1 > end1 or start2 > end2:
        # raise ValueError("Start coordinate is bigger than end one")
        print("Start coordinate is bigger than end one")
        return None
    if start1 < start2:
        return end2 <= end1
    elif start1 > start2:
        return end1 <= end2
    else:
        return True


def is_overlapped(start1, end1, start2, end2):
    if start1 > end1 or start2 > end2:
        # raise ValueError("Start coordinate is bigger than end one")
        print("Start coordinate is bigger than end one")
        return None
    if start1 < start2:
        return start2 < end1 < end2
    elif start1 > start2:
        return start1 < end2 < end1
    else:
        return False


def chaining(chrom, species, dir_path):
    global CHAINING_COUNT
    chaining_flag = True
    chaining_blocks = []

    blocks = []
    t_coord = chrom.loc[0, "tStart"]

    for i in range(chrom.shape[0]):
        if not chaining_flag:
            if blocks:
                chaining_blocks.append((blocks, t_coord))
                blocks = []
                t_coord = chrom.loc[i, "tStart"]
            chaining_flag = True

        if not blocks:
            blocks.append(chrom.loc[i, "synteny_block_id"])
        end_t = chrom.loc[i, "tEnd"]
        if i + 1 == chrom.shape[0]:
            break
        start_t = chrom.loc[i + 1, "tStart"]
        strand_cur = chrom.loc[i, "strand"]
        strand_next = chrom.loc[i + 1, "strand"]

        if strand_cur != strand_next:
            chaining_flag = False

        elif start_t - end_t < CHAINING_STEP:
            if strand_cur == "+":
                end_q = chrom.loc[i, "qEnd"]
                start_q = chrom.loc[i + 1, "qStart"]
                end_tmp = chrom.loc[i + 1, "qEnd"]
                if start_q - end_q < CHAINING_STEP:
                    blocks.append(chrom.loc[i + 1, "synteny_block_id"])
                else:
                    chaining_flag = False

            if strand_cur == "-":
                end_q = chrom.loc[i + 1, "qEnd"]
                start_q = chrom.loc[i, "qStart"]
                end_tmp = chrom.loc[i, "qEnd"]
                if start_q - end_q < CHAINING_STEP:
                    blocks.append(chrom.loc[i + 1, "synteny_block_id"])
                else:
                    chaining_flag = False

        else:
            chaining_flag = False

    chaining_blocks.append((blocks, t_coord))
    chaining_blocks = sorted(chaining_blocks, key=lambda x: x[1])
    # print(chaining_blocks)

    with open(f"{dir_path}/{species}.to.MFOI.chainID", "a") as f:
        for chain in chaining_blocks:
            chain = chain[0]
            if len(chain) > 1:
                n = len(chrom)
                chrom = pd.concat([chrom, chrom[chrom["synteny_block_id"] == chain[0]]], ignore_index=True)
                if chrom.loc[n, "strand"] == "+":
                    chrom.loc[n, "qEnd"] = chrom[chrom["synteny_block_id"] == chain[-1]]["qEnd"].reset_index(drop=True)[0]
                    chrom.loc[n, "tEnd"] = chrom[chrom["synteny_block_id"] == chain[-1]]["tEnd"].reset_index(drop=True)[0]

                if chrom.loc[n, "strand"] == "-":
                    chrom.loc[n, "qStart"] = chrom[chrom["synteny_block_id"] == chain[-1]]["qStart"].reset_index(drop=True)[0]
                    chrom.loc[n, "tEnd"] = chrom[chrom["synteny_block_id"] == chain[-1]]["tEnd"].reset_index(drop=True)[0]

                chrom.loc[n, "tHitLen"] = sum(chrom[chrom["synteny_block_id"].isin(chain)]["tHitLen"])
                chrom.loc[n, "qHitLen"] = sum(chrom[chrom["synteny_block_id"].isin(chain)]["qHitLen"])
                chrom.loc[n, "synteny_block_id"] = f"SC_{CHAINING_COUNT}"
                f.write(f"SC_{CHAINING_COUNT}," + ",".join(chain) + "\n")
                CHAINING_COUNT += 1
                # chr = chr[~chr["synteny_block_id"].isin(chain)].reset_index(drop=True)

            else:
                idx = chrom[chrom["synteny_block_id"] == chain[0]].index[0]
                chrom.loc[idx, "synteny_block_id"] = f"SC_{CHAINING_COUNT}"
                f.write(f"SC_{CHAINING_COUNT}," + ",".join(chain) + "\n")
                CHAINING_COUNT += 1

    for chain in chaining_blocks:
        chain = chain[0]
        if len(chain) > 1:
            chrom = chrom[~chrom["synteny_block_id"].isin(chain)]

    chrom = chrom.reset_index(drop=True)

    return chrom


def filter_parsed_psl(chrom, species, chr_name, target_chr_name, dir_path):
    log_file = open(f"{dir_path}/{species}.{chr_name}.to.MFOI.{target_chr_name}.filtered_synteny_blocks.log", "w")
    log_file.write("\t".join(
        ["strand", "qName", "qStart", "qEnd", "tName", "tStart", "tEnd", "tHitLen", "qHitLen", "synteny_block_id",
         "translocationChr", "translocationLen",
         "qNested", "tNested",
         "qOverlapLen", "qOverlapID", "tOverlapLen", "tOverlapID"]) + "\n")

    to_remove = []

    for i in range(chrom.shape[0]):
        start1 = chrom.loc[i, "qStart"]
        end1 = chrom.loc[i, "qEnd"]

        chr2 = chrom.loc[i, "tName"]
        if chr2 != target_chr_name:
            to_remove.append(chrom.loc[i, "synteny_block_id"])
            report = chrom.iloc[i, ].tolist()[1:11]
            report.append(chr2)
            report.append(end1 - start1 + 1)
            for _ in range(6):
                report.append(".")
            log_file.write("\t".join([str(r) for r in report]) + "\n")

        for j in range(i + 1, chrom.shape[0]):
            start2 = chrom.loc[j, "qStart"]
            end2 = chrom.loc[j, "qEnd"]

            if is_nested(start1, end1, start2, end2):
                # print("Nested", start1, end1, start2, end2)
                len1 = end1 - start1 + 1
                len2 = end2 - start2 + 1
                if len1 > len2:
                    report = chrom.iloc[j, ].tolist()[1:11]
                    for _ in range(2):
                        report.append(".")
                    report.append("+")
                    for _ in range(5):
                        report.append(".")
                    to_remove.append(chrom.loc[j, "synteny_block_id"])
                    log_file.write("\t".join([str(r) for r in report]) + "\n")

                else:
                    report = chrom.iloc[i, ].tolist()[1:11]
                    for _ in range(2):
                        report.append(".")
                    report.append("+")
                    for _ in range(5):
                        report.append(".")
                    to_remove.append(chrom.loc[i, "synteny_block_id"])
                    log_file.write("\t".join([str(r) for r in report]) + "\n")
                    to_remove.append(i)

            elif is_overlapped(start1, end1, start2, end2):
                # print("Overlapping", start1, end1, start2, end2)
                overlap_len = end1 - start2 + 1
                len1 = end1 - start1 + 1
                len2 = end2 - start2 + 1

                if overlap_len > OVERLAP_LEN:
                    if len1 > len2:
                        d, s = j, i
                    else:
                        d, s = i, j

                    to_remove.append(chrom.loc[d, "synteny_block_id"])

                    report = chrom.iloc[d, ].tolist()[1:11]
                    for _ in range(4):
                        report.append(".")
                    report.append(overlap_len)
                    report.append(chrom.loc[s, "synteny_block_id"])
                    for _ in range(2):
                        report.append(".")
                    log_file.write("\t".join([str(r) for r in report]) + "\n")

                    if REMOVE_ALL:
                        to_remove.append(chrom.loc[s, "synteny_block_id"])

                        report = chrom.iloc[s, ].tolist()[1:11]
                        for _ in range(4):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[d, "synteny_block_id"])
                        for _ in range(2):
                            report.append(".")
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                else:
                    if len1 * OVERLAP_PERCENT < overlap_len and len2 * OVERLAP_PERCENT < overlap_len:
                        if len1 > len2:
                            min_len, max_len, d, s = len2, len1, j, i
                        else:
                            min_len, max_len, d, s = len1, len2, i, j

                        to_remove.append(chrom.loc[d, "synteny_block_id"])
                        report = chrom.iloc[d, ].tolist()[1:11]
                        for _ in range(4):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[s, "synteny_block_id"])
                        for _ in range(2):
                            report.append(".")
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                    if REMOVE_ALL:
                        to_remove.append(chrom.loc[s, "synteny_block_id"])
                        report = chrom.iloc[s, ].tolist()[1:11]
                        for _ in range(4):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[d, "synteny_block_id"])
                        for _ in range(2):
                            report.append(".")
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                    elif len1 * OVERLAP_PERCENT < overlap_len:
                        to_remove.append(chrom.loc[i, "synteny_block_id"])
                        report = chrom.iloc[i, ].tolist()[1:11]
                        for _ in range(4):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[j, "synteny_block_id"])
                        for _ in range(2):
                            report.append(".")
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                    elif len2 * OVERLAP_PERCENT < overlap_len:
                        to_remove.append(chrom.loc[j, "synteny_block_id"])
                        report = chrom.iloc[j, ].tolist()[1:11]
                        for _ in range(4):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[i, "synteny_block_id"])
                        for _ in range(2):
                            report.append(".")
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

            else:
                break

    chrom = chrom.sort_values(by="tStart").reset_index(drop=True)

    for i in range(chrom.shape[0]):
        start1 = chrom.loc[i, "tStart"]
        end1 = chrom.loc[i, "tEnd"]

        chr2 = chrom.loc[i, "tName"]
        if chr2 != target_chr_name:
            continue

        for j in range(i + 1, chrom.shape[0]):
            start2 = chrom.loc[j, "tStart"]
            end2 = chrom.loc[j, "tEnd"]

            if is_nested(start1, end1, start2, end2):
                # print("Nested", start1, end1, start2, end2)
                len1 = end1 - start1 + 1
                len2 = end2 - start2 + 1
                if len1 > len2:
                    report = chrom.iloc[j, ].tolist()[1:11]
                    for _ in range(3):
                        report.append(".")
                    report.append("+")
                    for _ in range(4):
                        report.append(".")
                    to_remove.append(chrom.loc[j, "synteny_block_id"])
                    log_file.write("\t".join([str(r) for r in report]) + "\n")

                else:
                    report = chrom.iloc[i, ].tolist()[1:11]
                    for _ in range(3):
                        report.append(".")
                    report.append("+")
                    for _ in range(4):
                        report.append(".")
                    to_remove.append(chrom.loc[i, "synteny_block_id"])
                    log_file.write("\t".join([str(r) for r in report]) + "\n")

            elif is_overlapped(start1, end1, start2, end2):
                overlap_len = end1 - start2 + 1
                len1 = end1 - start1 + 1
                len2 = end2 - start2 + 1

                if overlap_len > OVERLAP_LEN:
                    if len1 > len2:
                        d, s = j, i
                    else:
                        d, s = i, j

                    to_remove.append(chrom.loc[d, "synteny_block_id"])

                    report = chrom.iloc[d, ].tolist()[1:11]
                    for _ in range(6):
                        report.append(".")
                    report.append(overlap_len)
                    report.append(chrom.loc[s, "synteny_block_id"])
                    log_file.write("\t".join([str(r) for r in report]) + "\n")

                    if REMOVE_ALL:
                        to_remove.append(chrom.loc[s, "synteny_block_id"])

                        report = chrom.iloc[s, ].tolist()[1:11]
                        for _ in range(6):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[d, "synteny_block_id"])
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                else:
                    if len1 * OVERLAP_PERCENT < overlap_len and len2 * OVERLAP_PERCENT < overlap_len:
                        if len1 > len2:
                            min_len, max_len, d, s = len2, len1, j, i
                        else:
                            min_len, max_len, d, s = len1, len2, i, j

                        to_remove.append(chrom.loc[d, "synteny_block_id"])
                        report = chrom.iloc[d, ].tolist()[1:11]
                        for _ in range(6):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[s, "synteny_block_id"])
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                    if REMOVE_ALL:
                        to_remove.append(chrom.loc[s, "synteny_block_id"])
                        report = chrom.iloc[s, ].tolist()[1:11]
                        for _ in range(6):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[d, "synteny_block_id"])
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                    elif len1 * OVERLAP_PERCENT < overlap_len:
                        to_remove.append(chrom.loc[i, "synteny_block_id"])
                        report = chrom.iloc[i, ].tolist()[1:11]
                        for _ in range(6):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[j, "synteny_block_id"])
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

                    elif len2 * OVERLAP_PERCENT < overlap_len:
                        to_remove.append(chrom.loc[j, "synteny_block_id"])
                        report = chrom.iloc[j, ].tolist()[1:11]
                        for _ in range(6):
                            report.append(".")
                        report.append(overlap_len)
                        report.append(chrom.loc[i, "synteny_block_id"])
                        log_file.write("\t".join([str(r) for r in report]) + "\n")

            else:
                break

    log_file.close()

    to_remove = list(set(to_remove))
    chrom = chrom[~chrom["synteny_block_id"].isin(to_remove)].reset_index(drop=True)
    chrom = chrom.sort_values(by="qStart").reset_index(drop=True)
    chrom = chrom.sort_values(by="tStart")

    return chrom


def get_perm_from_psl(input_file, dir_path, species):
    table = pd.read_csv(input_file, sep="\t")

    query_chr = table['qName'].unique()
    target_chr = table['tName'].unique()
    file_perm_names = []

    output_csv = f"{dir_path}/{species}.to.MFOI.chains.csv"

    with open(f"{dir_path}/{species}.to.MFOI.general_stats.stats", "w") as f:
        f.write("\t".join(["qName", "tName", "N", "qSumLen", "tSumLen"]) + "\n")
        for q in query_chr:
            for t in target_chr:
                chrom = table[(table['qName'] == q) & (table['tName'] == t)]
                f.write("\t".join(
                    [q, t, str(chrom.shape[0]), str(sum(chrom["qHitLen"])), str(sum(chrom["tHitLen"]))]) + "\n")

    with open(f"{dir_path}/{species}.to.MFOI.chainID", "w") as f:
        pass

    for chr_name in query_chr:
        chrom = table[table['qName'] == chr_name].reset_index(drop=True)
        chrom = chrom.sort_values(by="qStart")
        target_chr_name = chrom['tName'].value_counts().sort_values(ascending=False).reset_index().loc[0, "tName"]

        chrom = filter_parsed_psl(chrom, species, chr_name, target_chr_name, dir_path)
        chrom = chaining(chrom, species, dir_path)
        chrom = chrom.sort_values(by="tStart").reset_index(drop=True)

        chrom.to_csv(output_csv, mode='a', header=not os.path.exists(output_csv), index=False)
        chrom = chrom.sort_values(by="qStart")

        perm = chrom.index.tolist()
        signs = chrom['strand'].tolist()
        for i, sign in enumerate(signs):
            perm[i] += 1
            if sign == "-":
                perm[i] *= -1
        perm.insert(0, 0)
        perm.append(len(perm))

        chrom[['synteny_block_id', 'tHitLen']].to_csv(
            f"{dir_path}/{species}.{chr_name}.to.MFOI.{target_chr_name}.synteny_blocks_metadata.csv",
            index=False)

        # synteny_block_names = chrom['synteny_block_id'].tolist()
        # with open(f"{dir_path}/{species}.{chr_name}.to.MFOI.{target_chr_name}.synteny_blocks_metadata.csv", "w") as f:
        #     f.write(",".join(synteny_block_names))

        perm_file = open(f"{dir_path}/{species}.{chr_name}.to.MFOI.{target_chr_name}.filtered_synteny_blocks.perm",
                         "w")
        file_perm_names.append(f"{dir_path}/{species}.{chr_name}.to.MFOI.{target_chr_name}.filtered_synteny_blocks.perm")
        perm_file.write(" ".join([str(i) for i in perm]))
        perm_file.close()

    all_chains = pd.read_csv(f"{dir_path}/{species}.to.MFOI.chains.csv")
    all_chains.to_excel(f"{dir_path}/{species}.to.MFOI.chains.xlsx", index=False)

    return file_perm_names
