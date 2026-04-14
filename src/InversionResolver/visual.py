import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, PathPatch
from matplotlib.path import Path


def find_reversal_boundaries(p1, p2):
    boundaries = []
    if p1 == p2:
        return boundaries

    diff = [i for i in range(len(p1)) if p1[i] != p2[i]]
    if not diff:
        return boundaries

    left, right = diff[0], diff[-1]
    boundaries.append((abs(p1[left]), abs(p1[right])))
    return boundaries


def draw_bezier_curves(perm_list, filepath, title, synteny_names, synteny_order, blocks_len = None):
    if blocks_len is None:
        blocks_len = {i: 1 for i in range(len(perm_list[0]))}

    n_steps = len(perm_list)
    n_genes = len(perm_list[0])
    y_positions = []

    fig, ax = plt.subplots(figsize=(n_steps*2, n_genes*2))

    norm_coef = sum(blocks_len[i] for i in blocks_len) / n_genes

    for i, perm in enumerate(perm_list):
        x = i * 2.0
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

                ax.text(x - 0.6, y_text+0.2, synteny_order[abs(gene)],
                        ha='center', va='center', color='black', fontsize=15, fontweight='bold')

            y_position[f"left_{abs(gene)}"] = y
            y_position[f"right_{abs(gene)}"] = y - length
            y -= length + 0.5

        y_positions.append(y_position)

    for step in range(n_steps - 1):
        p1 = perm_list[step]
        p2 = perm_list[step + 1]
        x1 = step * 2.0
        x2 = (step + 1) * 2.0

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

    ax.set_xlim(-1, 2.0 * (n_steps - 1) + 1)
    # ax.set_ylim(-1, n_genes)
    ax.axis('off')
    ax.set_title(f"Chromosomal Reversals Visualization {title}",
                fontsize=14, fontweight='bold')

    # plt.tight_layout()
    # plt.show()

    plt.savefig(filepath, bbox_inches='tight')