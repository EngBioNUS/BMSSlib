# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 14:18:45 2018

@author: JingWui
"""

from distutils.core import setup
setup(
  name = 'BMSSlib',         # How you named your package folder (BMSSlib)
  packages = ['BMSSlib'],   # Chose the same as "name"
  version = '1.0',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Bio-Model Selection System',   
  author = 'Jing Wui Yeoh', 'Chueh Loo Poh',   
  author_email = 'yeohjingwui@gmail.com', 'poh.chuehloo@nus.edu.sg',      
  url = 'https://github.com/EngBioNUS/BMSSlib',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/EngBioNUS/BMSSlib/archive/v_01.tar.gz',   
  keywords = ['Genetic Circuit Modelling', 'Model Selection', 'Model Fitting', 'Automation', 'Characterization Data'],   # Keywords that define your package best
  install_requires=[            
          'tabulate',
          'tesbml',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3.6', # 3.6.5
      ],
)

# https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56
