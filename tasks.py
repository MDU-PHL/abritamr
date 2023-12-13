"""
Automate deployment to PyPi
"""

import invoke,json,pathlib


def get_token():

    with open(f"{pathlib.Path(__file__).parent /'token.json' }", "r") as cfg:

        return cfg['token']

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
    ctx.run(f"python3 -m twine upload --config-file {pathlib.Path(__file__).parent / '.pypirc' } dist/*")
    # ctx.run("git push origin --tags")
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

