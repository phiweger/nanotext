from setuptools import setup, find_packages
# https://packaging.python.org/tutorials/packaging-projects/
# https://docs.pytest.org/en/latest/goodpractices.html


with open('README.md', 'r') as fh:
    long_description = fh.read()


setup(
    name='nanotext',
    version='0.1',
    author='Adrian Viehweger',
    author_email='foo@bar.com',
    description='Domains as words, genomes as documents.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    # url='',
    license='BSD 3-clause',
    
    # next 2 lines for pytest
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

    install_requires=[
        'click',
        'numpy',
        'pandas',
        'tqdm',
        # 'pysam==0.15.1',
        'scikit-learn',
        # 'pybedtools==0.7.10',
        # relies on pysam
        # https://github.com/pysam-developers/pysam/issues/697
        'gensim',
        # current version ==3.7.1
        # 'faiss-cpu==1.5.0',
    ],

    # https://click.palletsprojects.com/en/7.x/setuptools/#testing-the-script
    # https://click.palletsprojects.com/en/7.x/setuptools/#scripts-in-packages
    packages=find_packages(),
    entry_points='''
        [console_scripts]
        nanotext=nanotext.__main__:cli
    ''',
)