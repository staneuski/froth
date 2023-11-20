#!/usr/bin/env pvbatch
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from sys import version_info
from typing import Callable

import paraview.simple as pv

LOGGING_FORMAT = '[foamio:%(levelname)s] (called at %(asctime)s) %(message)s'


def _load(args: argparse.Namespace) -> None:
    """https://kitware.github.io/paraview-docs/latest/python/_modules/paraview/simple.html#LoadState
    """

    logging.info(f'loading {args.state=} at {args.case=}')
    pv.LoadState(
        str(args.state),
        data_directory=str(args.case) if not args.case is None else None,
    )


def add_args(parser: argparse.ArgumentParser, pv_func: Callable) -> None:
    parser.add_argument('out', type=Path)
    parser.add_argument(
        '--dict',
        type=Path,
        metavar='JSON',
        help=f'dictionary file with any `{pv_func.__name__}` keyword arguments',
    )


def timesteps(args: argparse.Namespace) -> list:
    _load(args)
    return pv.GetTimeKeeper().TimestepValues


def visualise(args: argparse.Namespace, pv_func: Callable) -> None:
    _load(args)

    # Ensure that write directory exists
    args.out.parent.mkdir(parents=True, exist_ok=True)

    if args.dict is None:
        pv_func(str(args.out))
        return

    with open(args.dict) as f:
        kwargs = json.load(f)
        logging.info(f'{kwargs=}')
        pv_func(str(args.out), **kwargs)


def main() -> argparse.Namespace:
    parent_parser = argparse.ArgumentParser(
        description='Missing ParaView CLI.',
        epilog=f'foamio pvbatch'
        f' [Python {version_info.major}.{version_info.minor}.{version_info.micro}]'
        '\nCopyright (c) 2021-2023 Stanislau Stasheuski',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parent_parser.add_argument(
        '-v',
        '--verbose',
        help='output information',
        action='store_const',
        dest='loglevel',
        const=logging.INFO,
    )

    parent_parser.add_argument(
        '-c',
        '--case',
        type=Path,
        help='directory from where to load files'
        ' https://kitware.github.io/paraview-docs/latest/python/_modules/paraview/simple.html#LoadState',
    )

    parent_parser.add_argument(
        'state',
        type=Path,
        metavar='PVSM',
        help='path to ParaView state',
    )

    subparsers = parent_parser.add_subparsers(
        title='subcommands',
        dest='command',
        required=True,
    )

    #: TimestepValues {{{
    parser = subparsers.add_parser(
        'timesteps',
        help='output time-steps'
        'https://kitware.github.io/paraview-docs/latest/python/_modules/paraview/simple.html#GetTimeKeeper',
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help='suppress indices',
    )
    parser.set_defaults(func=lambda args: print(
        timesteps(args) if args.list else {
            i: time
            for i, time in enumerate(timesteps(args))
        }))
    #: }}}

    #: SaveScreenshot {{{
    parser = subparsers.add_parser(
        'screenshot',
        aliases=['s'],
        help='https://kitware.github.io/paraview-docs/latest/python/_modules/paraview/simple.html#SaveScreenshot',
    )
    add_args(parser, pv.SaveScreenshot)
    parser.set_defaults(func=lambda args: visualise(args, pv.SaveScreenshot))
    #: }}}

    #: SaveAnimation {{{
    parser = subparsers.add_parser(
        'animate',
        aliases=['a'],
        help='https://kitware.github.io/paraview-docs/latest/python/_modules/paraview/simple.html#SaveAnimation',
    )
    add_args(parser, pv.SaveAnimation)
    parser.set_defaults(func=lambda args: visualise(args, pv.SaveAnimation))
    #: }}}

    args = parent_parser.parse_args()
    logging.basicConfig(level=args.loglevel, format=LOGGING_FORMAT)

    args.func(args)


if __name__ == '__main__':
    main()
