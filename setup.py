import os
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), "co1_fish", "version.py")
version = open(version_py).read().strip().split("=")[-1].replace('"', "").strip()
long_description = """
``co1_fish``
"""
# NTS -> Fix all this when you actually have a working piece of software

HERE = os.path.dirname(__file__)

with open(os.path.join(HERE, "requirements.txt"), "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
    name="co1_fish",
    version=version,
    install_requires=install_requires,
    requires=["python (>=3.9)"],
    packages=["co1_fish"],
    author="Sam AJ Wilkinson",
    description="",
    long_description=long_description,
    url="",  # Add later
    package_dir={"co1_fish": "co1_fish"},
    package_data={"co1_fish": []},
    zip_safe=False,
    include_package_data=True,
    entry_points={"console_scripts": ["co1_fish=co1_fish.co1_fish_cli:main"],},
    author_email="s.a.j.wilkinson@bham.ac.uk",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
