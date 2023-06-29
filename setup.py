from setuptools import setup, find_packages

setup(
    name='megs',
    version='0.1.0.dev0',
    author='Ufuk Cakir',
    author_email='ufuk.cakir@stud.uni-heidelberg.de',
    description='MEGS: Morphological Evaluation of Galactic Structure',
    url='https://github.com/ufuk-cakir/MEGS',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
)
