from setuptools import setup, find_packages
import pathlib

import pkg_resources
import setuptools
import codecs
import os

# here = os.path.abspath(os.path.dirname(__file__))

# with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
#     long_description = "\n" + fh.read()

VERSION = '0.0.2'
DESCRIPTION = 'Visualize seismicity'
LONG_DESCRIPTION = 'Visualize seismicity through high-resolution maps and profiles.'

req_path = os.path.join(os.path.dirname(__file__),"requirements.txt")
with pathlib.Path('requirements.txt').open() as requirements_txt:
    install_requires = [
        str(requirement)
        for requirement
        in pkg_resources.parse_requirements(requirements_txt)
    ]

# Setting up
setup(
    name="EQViewer",
    version=VERSION,
    author="ecastillot (Emmanuel Castillo)",
    author_email="<ecastillot@unal.edu.co>",
    url="https://github.com/ecastillot/EQViewer",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=install_requires,
    keywords=['python', "eqviewer","earthquakes","seismology"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ],
    python_requires='>=3.8'
)

# python setup.py sdist bdist_wheel
# twine upload dist/*