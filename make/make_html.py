#!/usr/bin/env python

'Make an HTML page of all the given figures.'

import argparse
import sys


_MODULE_DOC = __doc__
TEMPLATE = """<p>
<h3>{fig}</h3>
<a href="{fig}"><img src="{fig}" style="width: 500px;"></a>
</p>
"""


def make_figure_html(figures):
    html = ''

    for fig in figures:
        html += TEMPLATE.format(fig=fig)

    return html


def save_html(html, outfile):
    with open('template.html', 'r') as t:
        template = t.read()

    with open(outfile, 'w') as out:
        out.write(template.format(html=html))


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=_MODULE_DOC)
    parser.add_argument('outfile', type=str, help='Name of output HTML file.')
    parser.add_argument('figures', type=str, nargs='+',
                        help='Figures to include in the webpage.')
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    html = make_figure_html(args.figures)
    save_html(html, args.outfile)


if __name__ == '__main__':
    sys.exit(main())
