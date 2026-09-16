"""Assemble the supplied DFT images without changing their numerical contents."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


def assemble(only=None):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                         'axes.labelsize': 8, 'axes.titlesize': 9,
                         'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'legend.fontsize': 7, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'axes.linewidth': 0.6})
    def panels(items, name, ncols=3, width=6.8):
        if only is not None and name != only:
            return
        nrows = int(np.ceil(len(items) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(width, 2.38 * nrows), squeeze=False)
        for index, (filename, title) in enumerate(items):
            ax = axes.flat[index]
            ax.imshow(Image.open(ROOT / 'data/original_panels' / filename), interpolation='none')
            ax.axis('off')
            ax.set_title(f'({chr(97 + index)}) {title}', fontsize=9, pad=2)
        for ax in list(axes.flat)[len(items):]:
            ax.axis('off')
        fig.subplots_adjust(left=.005, right=.995, bottom=.005, top=.94, hspace=.18, wspace=.015)
        fig.savefig(ROOT / 'figures' / name, dpi=450)
        plt.close(fig)

    # The legacy 2-85 filename contains the coordinate-verified 2.84-A series.
    szv = [('DOS_3-30_theta30-60.jpg', '3.30'), ('DOS_3-0_theta30-60.jpg', '3.00'),
           ('DOS_2-85_theta30-60.jpg', '2.84'), ('DOS_2-70_theta30-60.jpg', '2.70'),
           ('DOS_2-55_theta30-60.jpg', '2.55'), ('DOS_2-30_theta30-60_l.jpg', '2.30'),
           ('DOS_2-25_theta30-60_l.jpg', '2.25'), ('DOS_2-20_theta30-60_l.jpg', '2.20'),
           ('DOS_2-15_theta30-60_l.jpg', '2.15'), ('DOS_2-0_theta30-60_l.jpg', '2.00')]
    panels([(szv[i][0], f'$d={szv[i][1]}$ Å') for i in [0, 4, 8]], 'dft_main.pdf')
    panels([(filename, f'$d={distance}$ Å') for filename, distance in szv[:6]], 'dft_survey_1.pdf')
    panels([(filename, f'$d={distance}$ Å') for filename, distance in szv[6:]], 'dft_survey_2.pdf', ncols=2)
    items = []
    for basis in ['SZV', 'DZVP']:
        for distance in ['25', '20', '15']:
            filename = f'DOS_2-{distance}_theta30-60_l' + ('_DZVP' if basis == 'DZVP' else '') + '.jpg'
            items.append((filename, f'{basis}, $d=2.{distance}$ Å'))
    panels(items, 'basis_comparison.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', choices=['dft_main.pdf', 'dft_survey_1.pdf',
                                         'dft_survey_2.pdf', 'basis_comparison.pdf'])
    assemble(parser.parse_args().only)
