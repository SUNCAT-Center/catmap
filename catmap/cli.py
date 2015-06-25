import os
import shutil

usage = {}

usage['import'] = """catmap import <mkm-file>
    Open a *.mkm project file and work with it interactively.
"""


def get_options(args=None, get_parser=False):
    import optparse
    import os
    from glob import glob
    import catmap

    parser = optparse.OptionParser(
        'Usage: %prog [help] ('
        + '|'.join(sorted(usage.keys()))
        + ') [options]',
        version=catmap.__version__)

    if args is not None:
        options, args = parser.parse_args(args.split())
    else:
        options, args = parser.parse_args()
    if len(args) < 1:
        parser.error('Command expected')
    if get_parser:
        return options, args, parser
    else:
        return options, args


def match_keys(arg, usage, parser):
    """Try to match part of a command against
       the set of commands from usage. Throws
       an error if not successful.

    """
    possible_args = [key for key in usage if key.startswith(arg)]
    if len(possible_args) == 0:
        parser.error('Command "%s" not understood.' % arg)
    elif len(possible_args) > 1:
        parser.error(('Command "%s" ambiguous.\n'
                      'Could be one of %s\n\n') % (arg, possible_args))
    else:
        return possible_args[0]


def main(args=None):
    """The CLI main entry point function.

    The optional argument args, can be used to
    directly supply command line argument like

    $ catmap <args>

    otherwise args will be taken from STDIN.

    """
    from glob import glob

    options, args, parser = get_options(args, get_parser=True)

    if not args[0] in usage.keys():
        args[0] = match_keys(args[0], usage, parser)

    elif args[0] == 'import':
        if len(args) < 2:
            parser.error('mkm filename expected.')

        from catmap import ReactionModel
        mkm_file = args[1]
        global model
        model = ReactionModel(setup_file=mkm_file)
        sh(banner='Note: model = catmap.ReactionModel(setup_file=\'%s\')\n# do model.run()\nfor a fully initialized model.' %
           args[1])


def sh(banner):
    """Wrapper around interactive ipython shell
    that factors out ipython version depencies.

    """

    from distutils.version import LooseVersion
    import IPython
    if hasattr(IPython, 'release'):
        try:
            from IPython.terminal.embed import InteractiveShellEmbed
            InteractiveShellEmbed(banner1=banner)()

        except ImportError:
            try:
                from IPython.frontend.terminal.embed \
                    import InteractiveShellEmbed
                InteractiveShellEmbed(banner1=banner)()

            except ImportError:
                from IPython.Shell import IPShellEmbed
                IPShellEmbed(banner=banner)()
    else:
        from IPython.Shell import IPShellEmbed
        IPShellEmbed(banner=banner)()
