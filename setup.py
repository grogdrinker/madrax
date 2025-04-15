#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   Copyright 2021 Gabriele Orlando
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


from setuptools import setup

#with open("README.md", "r") as fh:

#    long_description = fh.read()
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
setup(
     name='madrax',
     zip_safe=False,
     include_package_data=True,

     version='0.0.6',

     author="Gabriele Orlando",

     author_email="gabriele.orlando@umontpellier.fr",

     description="A differentiable force field implemented in pytorch which can be integrated in neural networks in an end-to-end way",

     long_description=long_description,

     long_description_content_type="text/markdown",

     url="https://github.com/grogdrinker/madrax",
     
     packages=['madrax','madrax.mutate','madrax.sources','madrax.sources.kde'],
     package_dir={'madrax': 'madrax/','madrax.sources':'madrax/sources','madrax.mutate':'madrax/mutate','madrax.sources.kde':'madrax/sources/kde'},
     package_data={'madrax': ['parameters/**', '*.so',"parameters/weightsKDE/*"]},

     python_requires=">=3.7,<=3.12",

     install_requires=["torch","numpy","matplotlib"],

    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: POSIX :: Linux",
    ],

 )
