import click
from click.exceptions import UsageError
from click._compat import get_text_stderr
import sys

import run
from version import __version__ as version

def _show_usage_error(self, file=None):
    if file is None:
        file = get_text_stderr()
    color = None
    if self.ctx is not None:
        color = self.ctx.color
        click.echo(self.ctx.get_help() + '\n', file=file, color=color)
    click.echo('ERROR: %s' % self.format_message(), file=file, color=color)


UsageError.show = _show_usage_error
CONTEXT_SETTINGS = dict(
    help_option_names=['-h', '--help'],
    
)

@click.group(no_args_is_help=True, context_settings=CONTEXT_SETTINGS, invoke_without_command=True)
@click.version_option(version, message='%(prog)s %(version)s')
# @click.pass_context
def cli():
    pass


cli.add_command(run.run)




if __name__ == '__main__':
    cli()
