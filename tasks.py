"""
Automate deployment to PyPi
"""

import invoke


@invoke.task
def deploy(ctx):
    """
    Automate deployment
    rm -rf build/* dist/*
    bumpversion patch --verbose
    python3 setup.py sdist bdist_wheel
    twine upload dist/*
    git push --tags
    """
    ctx.run("rm -rf build/* dist/*")
    # ctx.run("bumpversion {bump} --verbose")
    ctx.run("python3 setup.py sdist bdist_wheel")
    ctx.run("python3 -m twine check dist/*")
    ctx.run("python3 -m twine upload dist/*")
    ctx.run("git push origin --tags")
    # ctx.run("git push kristy --tags")

@invoke.task
def gitpush(ctx, message):
    """
    Automate push for minor changes including typos and reversions
    """
    message = ' '.join(message.split('_'))
    ctx.run("git add -A")
    ctx.run(f"git commit -m '{message}'")
    ctx.run("git push origin")
    # ctx.run("git push kristy")

@invoke.task
def gittag(ctx):
    """
    Automate push and tagging
    for any change that will change behaviour
    """
    # message = ' '.join(message.split('_'))
    # ctx.run("git add -A")
    # ctx.run(f"git commit -m '{message}'")
    # ctx.run(f"git tag -a {tag} -m {message}")
    ctx.run("git push origin --tags")
    # ctx.run("git push kristy --tags")

import subprocess, datetime, yaml, toml

# from mdu_writer.verifications.meningotype_write import WriteMenigotypeVerify

@task
def update_singularity(c):
    date = datetime.datetime.today().strftime("%d_%m_%y")
    c.run("python3 update_container.py")
    c.run("git add *")
    c.run(f"git commit -m updated singulairty {date}")
    c.run("git push")

@task
def build_container(c):
    config = toml.load("config.toml")
    date = datetime.datetime.today().strftime("%d_%m_%y")
    archive_dir = f"{config['archive_dir']}/{config['version']}_DB{date}"
    print("Building container")
    c.sudo(f'singularity build salmonella_typing.simg Singularity')
    c.sudo(f"mkdir {archive_dir}") # make archive directory for this image
    c.sudo(f"cp {config['container_dir']}/salmonella_typing.simg {archive_dir}") #copy to the archive directory
    c.sudo(f"mv salmonella_typing.simg {config['container_dir']}") #move to the execution directory
