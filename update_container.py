'''
Given a new amrfinderplus version or DB, update the container
'''

import pathlib
import click
import jinja2
import toml
import pendulum


def load_template(name):
    '''
    Return the singularity recipe template as unicode text
    '''
    template = pathlib.Path(name).read_text()
    return template


@click.command()
@click.option("--amrfinderplus_version", default=None)
@click.option("--db_version", is_flag=True)
@click.option("--author", default=None)
@click.option("-c", "--config", default="config.toml")
def update_amrfinderplus_singularity(amrfinderplus_version, db_version, author, config):
    '''
    Use the config.toml, or override any of the options via the command line
    '''
    # load the params
    config = toml.load(config)
    if amrfinderplus_version is not None:
        config['amrfinderplus_version'] = mlst_version
    if author is not None:
        config['author'] = author
    if db_version:
        config['db_version'] = False
    # load the template
    loader = jinja2.FunctionLoader(load_template)
    env = jinja2.Environment(loader=loader)
    SINGULARITY_RECIPE = env.get_template("amrfinderplus.singularity").render(config)
    # prepare the folders
    version_path = pathlib.Path(f'v{config["amrfinderplus_version"]}')
    if not version_path.exists():
        version_path.mkdir()
    today = pendulum.today().format('YYYYMMDD')
    subfolder_path = version_path / today
    if not subfolder_path.exists():
        subfolder_path.mkdir()
    # create local version
    local_recipe = subfolder_path / \
        f'Singularity.v{config["amrfinderplus_version"]}_{today}'
    local_recipe.write_text(SINGULARITY_RECIPE)
    # create global version
    global_recipe = pathlib.Path("Singularity")
    global_recipe.write_text(SINGULARITY_RECIPE)


if __name__ == "__main__":
    update_amrfinderplus_singularity()