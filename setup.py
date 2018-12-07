import setuptools
# https://packaging.python.org/tutorials/packaging-projects/
# https://docs.pytest.org/en/latest/goodpractices.html

with open('README.md', 'r') as fh:
    long_description = fh.read()


setuptools.setup(
    name='nanotext',
    version='0.0.1',
    author='Adrian Viehweger',
    author_email='adrian.viehweger@googlemail.com',
    description='Domains as words, genomes as documents.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    # next 2 lines for pytest
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)