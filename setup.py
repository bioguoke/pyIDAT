#coding=utf-8
#!/usr/bin/python
# Created on 
# Author: Zhihua Pei
# Organization: www.bioguoke.com
# website: www.zilhua.com
# github: github.com/zilhua

from setuptools import setup,find_packages

setup(
    name="pyIDAT",
    version="1.0",
    description="analysis for illumina beadarray .idat format file",
    author="zilhua",
    author_email="andyzilhua@qq.com",
    url="https://github.com/bioguoke/pyIDAT",
    license="GPLv3",
    keywords=["pyIDAT", "idat"],
    packages=["."],
    install_requires=["pandas"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ]
)