import click
from nanotext.cli.embed import train
from nanotext.cli.similarity import search, compare, taxonomy, lookup
from nanotext.cli.predict import predict

# TODO: bash completion
# http://click.palletsprojects.com/en/7.x/bashcomplete/


@click.group()
def cli():
    pass


# embed
cli.add_command(train)

# similarity
cli.add_command(search)
cli.add_command(compare)
cli.add_command(taxonomy)
cli.add_command(lookup)

# predict
cli.add_command(predict)
