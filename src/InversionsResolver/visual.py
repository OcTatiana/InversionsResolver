import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, PathPatch
from matplotlib.path import Path


def find_reversal_boundaries(p1, p2):
    boundaries = []
    if p1 == p2:
        return boundaries

    n = len(p1)
    i = 0

    while i < n:
        if p1[i] != p2[i]:
            start = i
            j = i
            while j < n and p1[j] != p2[j]:
                j += 1
            k = start
            while k < j:
                for end in range(k, j):
                    length = end - k + 1
                    is_rev = True
                    for offset in range(length):
                        if p1[k + offset] != -p2[end - offset]:
                            is_rev = False
                            break
                    if is_rev:
                        boundaries.append((abs(p1[k]), abs(p1[end])))
                        k = end + 1
                        break
                else:
                    k += 1
            i = j
        else:
            i += 1

    return boundaries


def draw_bezier_curves(perm_list, filepath, synteny_names, synteny_order,
                       query, target,
                       blocks_len=None, chr_list_q=None, chr_list_t=None, synteny_block_names=None):

    if blocks_len is None:
        blocks_len = {i: 1 for i in range(len(perm_list[0]))}

    n_steps = len(perm_list)
    n_genes = len(perm_list[0])
    y_positions = []
    shift = 0

    if chr_list_t:
        fig, ax = plt.subplots(figsize=(n_steps*2+6, n_genes*2))
        shift = 2
    else:
        fig, ax = plt.subplots(figsize=(n_steps*2, n_genes*2))

    norm_coef = sum(blocks_len[i] for i in blocks_len) / n_genes

    if chr_list_q:
        n = len(chr_list_q)
        if n < 10:
            f = 16
        elif n < 20:
            f = 12
        else:
            f = 8
        width = 1
        height = n_genes*2
        x0, y0 = -0.5, -n_genes-2

        # Draw outer rounded rectangle (chromosome shape)
        chromosome = FancyBboxPatch(
            (x0, y0),
            width,
            height,
            boxstyle="round,pad=0.02,rounding_size=0.5",
            linewidth=2,
            edgecolor='black',
            facecolor='none'
        )
        ax.add_patch(chromosome)

        for i in range(n):
            block_name = chr_list_q[n - i - 1]
            y = y0 + i*n_genes*2/n

            if block_name in synteny_block_names:
                rect = Rectangle(
                    (x0, y),
                    width,
                    n_genes*2/n,
                    facecolor='violet',
                    alpha=0.5,
                    edgecolor='none'
                )
                rect.set_clip_path(chromosome)
                ax.add_patch(rect)

            if i != 0:
                ax.plot([x0, x0 + width], [y, y], color='black', linewidth=1)

            ax.text(
                x0 + width / 2,
                y + n_genes*2/(2*n),
                block_name,
                ha='center',
                va='center',
                fontsize=15
            )

        ax.text(
            x0 + 0.5,
            y0 + height + 1,
            f'{query}\n{filepath.split(".")[1]}',
            ha='center',
            va='center',
            fontsize=f
        )

    for i, perm in enumerate(perm_list):
        x = i * 2.0 + shift
        y = n_genes
        y_position = {}
        for j, gene in enumerate(perm):
            # color = colors[abs(gene)]
            color = "grey"
            direction = np.sign(gene)
            if gene == 0:
                direction = 1

            length = 2 * blocks_len[abs(gene)] / norm_coef

            # print(length)

            if direction > 0:
                y_start = y - length
            else:
                y_start = y

            ax.arrow(
                x,
                y_start,
                0,
                length * direction,
                head_width=0.2,
                head_length=0.2,
                fc=color,
                ec=color,
                length_includes_head=True
            )

            if length >= 0:
                if direction > 0:
                    y_text = y_start + length / 2
                else:
                    y_text = y_start - length / 2

                ax.text(x-0.6, y_text-0.2, synteny_names[abs(gene)],
                        ha='center', va='center', color='grey', fontsize=15, fontweight='bold')

                ax.text(x-0.6, y_text+0.2, synteny_order[abs(gene)],
                        ha='center', va='center', color='black', fontsize=15, fontweight='bold')

            y_position[f"left_{abs(gene)}"] = y
            y_position[f"right_{abs(gene)}"] = y - length
            y -= length + 0.5

        y_positions.append(y_position)

    for step in range(n_steps - 1):
        p1 = perm_list[step]
        p2 = perm_list[step + 1]
        x1 = step * 2.0 + shift
        x2 = (step + 1) * 2.0 + shift

        for left_gene, right_gene in find_reversal_boundaries(p1, p2):
            # print(left_gene, right_gene)
            y1_left = y_positions[step][f"left_{abs(left_gene)}"]
            y2_left = y_positions[step+1][f"right_{abs(left_gene)}"]

            y1_right = y_positions[step][f"right_{abs(right_gene)}"]
            y2_right = y_positions[step+1][f"left_{abs(right_gene)}"]

            for (y1, y2, gene) in [(y1_left, y2_left, left_gene), (y1_right, y2_right, right_gene)]:
                left_sign = np.sign(next(g for g in p1 if abs(g) == gene))
                right_sign = np.sign(next(g for g in p2 if abs(g) == gene))
                curve_color = "red" if left_sign != right_sign else "gray"
                verts = [
                    (x1 + 0.3, y1),
                    (x1 + 0.7, (y1 + y2)/2),
                    (x2 - 0.7, (y1 + y2)/2),
                    (x2 - 0.3, y2)
                ]
                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
                path = Path(verts, codes)
                patch = PathPatch(path, facecolor='none', lw=1.5, edgecolor=curve_color, alpha=1)
                ax.add_patch(patch)

    # ax.set_xlim(-1, 2.0 * (n_steps - 1) + 1)
    # ax.set_ylim(-1, n_genes)
    ax.axis('off')

    # Add frame
    vertices = [(-0.2 + shift, y_position[f"left_{abs(1)}"] + 0.1), (-0.2 + shift, y_position[f"right_{abs(n_genes-2)}"] - 0.1),
                (n_steps*2 - 1.8 + shift, y_position[f"right_{abs(n_genes-2)}"] - 0.1), (n_steps*2 - 1.8 + shift, y_position[f"left_{abs(1)}"] + 0.1)]
    x = [v[0] for v in vertices]
    y = [v[1] for v in vertices]
    x.append(vertices[0][0])
    y.append(vertices[0][1])
    plt.plot(x, y, linestyle=':', color='green', linewidth=2)

    # Add chromosome diagrams
    if chr_list_t:
        n = len(chr_list_t)
        if n < 10:
            f = 16
        elif n < 20:
            f = 12
        else:
            f = 8
        width = 1
        height = n_genes * 2
        x0, y0 = 2 * n_steps + 1, -n_genes - 2

        # Draw outer rounded rectangle (chromosome shape)
        chromosome = FancyBboxPatch(
            (x0, y0),
            width,
            height,
            boxstyle="round,pad=0.02,rounding_size=0.5",
            linewidth=2,
            edgecolor='black',
            facecolor='none'
        )
        ax.add_patch(chromosome)

        for i in range(n):
            block_name = chr_list_t[n - i - 1]
            y = y0 + i*n_genes*2/n

            if block_name in synteny_block_names:
                rect = Rectangle(
                    (x0, y),
                    width,
                    n_genes*2/n,
                    facecolor='violet',
                    alpha=0.5,
                    edgecolor='none'
                )
                rect.set_clip_path(chromosome)
                ax.add_patch(rect)

            if i != 0:
                ax.plot([x0, x0 + width], [y, y], color='black', linewidth=1)

            ax.text(
                x0 + width / 2,
                y + n_genes*2/(2*n),
                block_name,
                ha='center',
                va='center',
                fontsize=f
            )

        ax.text(
            x0 + 0.5,
            y0 + height + 1,
            f'{target}\n{filepath.split(".")[4].split("_")[0]}',
            ha='center',
            va='center',
            fontsize=15
        )

    # ax.set_title(f"Chromosomal Reversals Visualization {title}",
    #             fontsize=14, fontweight='bold')

    # plt.tight_layout()
    # plt.show()

    plt.savefig(filepath, bbox_inches='tight')